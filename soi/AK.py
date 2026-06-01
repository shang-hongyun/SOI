import sys
import os
import time
import itertools
import random
import math
import re
from collections import Counter, defaultdict, OrderedDict
from typing import Optional, Set, Dict, List, Tuple
import networkx as nx

from .mcscan import ColinearGroups, Gff, GffGraph, SyntenyGraph, GffLine
from .hog import HOG, HOGrecord
from .RunCmdsMP import logger
from .ak_dotplot import draw_akr_dotplots, draw_dotplot, aag_to_karyo

try:
    import pulp
except ImportError:
    pulp = None


# =====================
# 重排事件记录
# =====================

class RearrangementEvent:
    """记录一个重排事件的详细信息"""
    EVENT_TYPES = {
        'tandem_dup', 'proximal_dup', 'dispersed_dup',
        'indel', 'inversion', 'internal_inversion', 'telomere_inversion',
        'unidir_trans',
        'reciprocal_translocation',
        'unbalanced_reciprocal_translocation',
        'eej',
        'fission',
        'ncf',
    }

    def __init__(self, event_type, node,
                 genes_involved=None, parent_chroms=None,
                 desc='', support=None, child_source=None, og_supported=None):
        if event_type not in self.EVENT_TYPES:
            logger.warning('Unknown event type: {}'.format(event_type))
        self.event_type = event_type
        self.node = node
        self.genes_involved = genes_involved or []
        self.parent_chroms = parent_chroms or []
        self.desc = desc
        self.support = support
        self.child_source = child_source  # which child branch this event belongs to
        self.og_supported = og_supported  # whether supported by outgroup evidence

    def __repr__(self):
        return '<{} on {}: {}>'.format(self.event_type, self.node, self.desc)

    def to_dict(self):
        d = {
            'event_type': self.event_type,
            'node': self.node,
            'genes_involved': [str(g) for g in self.genes_involved],
            'parent_chroms': self.parent_chroms,
            'desc': self.desc,
            'support': self.support,
            'child_source': self.child_source,
        }
        if self.og_supported is not None:
            d['og_supported'] = self.og_supported
        return d


# =====================
# 邻接图 + 断点图融合：AncestralAdjacencyGraph
# =====================

class AncestralAdjacencyGraph:
    """
    融合邻接图与断点图思想的祖先邻接图。
    以带端粒的 GffGraph 为骨架，为每个节点（基因）和端粒连接
    维护跨物种的观测信息，便于外类群极化。
    """

    def __init__(self, node_id, species_set=None):
        self.node_id = node_id
        self.species_set = species_set or set()
        self.graph = nx.DiGraph()
        self.telomeres = set()
        self.gene_nodes = set()
        self.hog_map = {}
        self.chrom_map = {}
        self.events = []
        # 染色体来源追踪：记录每个节点/边来自哪个子图（source_id）
        self.node_sources = {}   # node -> source_id (e.g., child node name)
        self.edge_sources = {}   # (n1, n2) -> source_id
        # 染色体物种来源：chrom_index -> set of species
        self.chrom_species = {}  # chrom_index -> {species}

    @classmethod
    def from_gffgraph(cls, node_id, gffG, add_telomere=True):
        aag = cls(node_id=node_id)
        for n in gffG.nodes:
            aag.graph.add_node(n, **dict(gffG.nodes[n]))
            if isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R'):
                aag.telomeres.add(n)
            else:
                aag.gene_nodes.add(n)
                if hasattr(n, 'chrom'):
                    aag.chrom_map[n] = n.chrom
        for n1, n2 in gffG.edges:
            aag.graph.add_edge(n1, n2, **dict(gffG[n1][n2]))
        if add_telomere and not aag.telomeres:
            aag._add_telomeres()
        return aag

    def _add_telomeres(self):
        seen_chroms = set()
        # 先收集所有起始点，避免迭代时修改 graph
        starts = list(self.starts)
        for idx, start in enumerate(starts):
            chrom_nodes = list(self.iter_chrom(start))
            if not chrom_nodes:
                continue
            first = chrom_nodes[0]
            last = chrom_nodes[-1]
            # 确定染色体标识：基因对象用 chrom 属性，HOG ID 字符串用虚拟名
            if hasattr(first, 'chrom'):
                chrom = first.chrom
            elif hasattr(last, 'chrom'):
                chrom = last.chrom
            else:
                chrom = "{}_chrom_{}".format(self.node_id, idx)
            if chrom in seen_chroms:
                continue
            seen_chroms.add(chrom)
            left_tel = (chrom, 'L')
            right_tel = (chrom, 'R')
            self.graph.add_node(left_tel, telomere=True)
            self.graph.add_node(right_tel, telomere=True)
            self.graph.add_edge(left_tel, first)
            self.graph.add_edge(last, right_tel)
            self.telomeres.add(left_tel)
            self.telomeres.add(right_tel)

    def rebuild_edges_from_chrom_hogs(self):
        """从 chrom_hogs 重建图边，移除自环和重复边。"""
        ch = getattr(self, 'chrom_hogs', None)
        if not ch:
            return
        chrom_ends = {}
        self.graph.clear_edges()
        for ci, hogs in ch.items():
            gene_hogs = [h for h in hogs if h not in self.telomeres]
            for i in range(len(gene_hogs) - 1):
                self.graph.add_edge(gene_hogs[i], gene_hogs[i + 1])
            if gene_hogs:
                chrom_ends[ci] = (gene_hogs[0], gene_hogs[-1])
        # 重连端粒
        for tel in self.telomeres:
            if not isinstance(tel, tuple) or len(tel) != 2:
                continue
            chrom_name, end = tel
            for ci, (first, last) in chrom_ends.items():
                if f'chrom_{ci}' in chrom_name:
                    if end == 'L':
                        self.graph.add_edge(tel, first)
                    elif end == 'R':
                        self.graph.add_edge(last, tel)
                    break

    @property
    def starts(self):
        for node, pred in self.graph.pred.items():
            if not pred and (node in self.gene_nodes or node in self.telomeres):
                yield node

    def iter_chrom(self, node):
        current = node
        yield current
        visited = {current}
        while True:
            succs = list(self.graph.successors(current))
            if not succs:
                break
            current = succs[0]
            if current in visited:
                break
            visited.add(current)
            yield current
            if current in self.telomeres and current != node:
                break

    @property
    def chromosomes(self):
        seen = set()
        for start in self.starts:
            chrom = []
            for n in self.iter_chrom(start):
                chrom.append(n)
                if n in seen and n not in self.telomeres:
                    break
            if chrom:
                real_nodes = [n for n in chrom if n not in self.telomeres]
                if any(n in seen for n in real_nodes):
                    continue
                seen.update(real_nodes)
                yield chrom

    def get_adjacencies(self, include_telomere=False):
        adjs = set()
        for n1, n2 in self.graph.edges():
            if not include_telomere:
                if n1 in self.telomeres or n2 in self.telomeres:
                    continue
            adjs.add((n1, n2))
        return adjs

    def remove_nodes(self, nodes):
        for n in list(nodes):
            if n in self.graph:
                preds = list(self.graph.predecessors(n))
                succs = list(self.graph.successors(n))
                for p in preds:
                    for s in succs:
                        self.graph.add_edge(p, s)
                self.graph.remove_node(n)
            self.gene_nodes.discard(n)
            self.telomeres.discard(n)

    def add_path(self, path):
        for i in range(len(path) - 1):
            self.graph.add_edge(path[i], path[i + 1])
            self.gene_nodes.add(path[i])
        if path:
            self.gene_nodes.add(path[-1])

    def to_gfa(self, fout):
        """直接输出 GFA 格式，去重双向边，并为端粒/普通节点着色"""
        # Segment 行：带颜色标签 CL:Z:#RRGGBB
        for node in self.graph.nodes():
            data = dict(self.graph.nodes[node])
            if data.get('telomere'):
                color = '#FF0000'  # 端粒红色
            else:
                color = '#808080'  # 普通节点灰色
            line = ['S', node, '*', 'CL:Z:{}'.format(color)]
            fout.write('\t'.join(map(str, line)) + '\n')
        # Link 行：去重，只输出一次
        seen = set()
        for n1, n2 in self.graph.edges():
            if (n2, n1) in seen:
                continue
            seen.add((n1, n2))
            line = ['L', n1, '+', n2, '+', '0M']
            fout.write('\t'.join(map(str, line)) + '\n')


# =====================
# AKR 主类
# =====================

class AKR:
    """
    Ancestral Karyotype Reconstruction (AKR)
    基于 HOG 层级同源群 + 邻接/断点图融合 + 外类群极化 的祖先核型重建
    核心关注 telomere-centric 重排：EEJ、Fission、NCF、相互易位
    """

    def __init__(self,
                 ogfile=None,
                 orthfiles=None,
                 sptreefile=None,
                 outpre="AKR",
                 paralog=False,
                 gfffile=None,
                 spsd=None,
                 rounds=3,
                 chrom_list=None,
                 min_genes=0,
                 timeout=600,
                 node_timeout=300,
                 use_ilp_sa=True,
                 sa_iterations=5000,
                 reconstruction_algorithm='v4_colored',
                 **kargs):

        self.ogfile = ogfile
        self.orthfiles = orthfiles
        self.sptreefile = sptreefile
        self.outpre = outpre
        self.paralog = paralog

        self.gfffile = gfffile
        self.spsd = spsd
        self.rounds = rounds
        self.chrom_list = chrom_list
        self.min_genes = min_genes
        self.timeout = timeout
        self.node_timeout = node_timeout
        self.use_ilp_sa = use_ilp_sa and pulp is not None
        self.sa_iterations = sa_iterations
        self._start_time = None

        self.hog = None
        self.tree = None
        self.leaf_graphs = {}
        self.anc_graphs = {}
        self.pre_wgd_graphs = {}  # WGD节点的pre-WGD图
        self.events = defaultdict(list)

        self.hogs_by_node = defaultdict(list)
        self.gene_to_hog = {}
        self.node_gene_to_hog = defaultdict(dict)
        self.gene_to_all_hogs = defaultdict(list)  # gene -> [hog_id, ...]（所有层级）
        self.hog_parent = {}           # child_hog_id -> parent_hog_id
        self.hog_children = defaultdict(list)  # parent_hog_id -> [child_hog_id]
        self.hog_species = defaultdict(set)    # hog_id -> {species}
        self._hog_node_cache = {}      # hog_id -> node_id（预解析，加速查找）
        self.ploidy_map = {}           # node_name -> ploidy_factor（从树文件解析）

    def _parse_ploidy_map(self):
        """从物种树文件中解析 [p=N] 多倍体标注"""
        if not self.sptreefile:
            return
        try:
            with open(self.sptreefile) as f:
                nw = f.read()
            for m in re.finditer(r'([\w. -]+)\[p=(\d+)\]', nw):
                name = m.group(1).strip()
                ploidy = int(m.group(2))
                if ploidy > 1:
                    self.ploidy_map[name] = ploidy
            if self.ploidy_map:
                logger.info("Ploidy annotations: {}".format(self.ploidy_map))
        except Exception as e:
            logger.warning("Failed to parse ploidy map: {}".format(e))

    def run(self):
        """Main entry point"""
        logger.info("=== Ancestral Karyotype Reconstruction (AKR) started ===")
        self._start_time = time.time()

        self._build_hogs()
        self._parse_ploidy_map()
        self._build_leaf_graphs()

        # v4_colored event-driven reconstruction
        from soi.takr_event_driven import reconstruct_event_driven_v2
        anc_graphs = reconstruct_event_driven_v2(self)
        self.anc_graphs = anc_graphs

        self._export_results()

        logger.info("=== Ancestral Karyotype Reconstruction (AKR) completed ===")
        return self.anc_graphs

    def _build_hogs(self):
        """Run HOG pipeline and build indexes"""
        self.hog = HOG(
            ogfile=self.ogfile,
            orthfiles=self.orthfiles,
            sptreefile=self.sptreefile,
            outpre=self.outpre,
            paralog=self.paralog
        )
        self.hog.pipe(write_tsv=True)
        self.tree = self.hog.tree

        for hog_id, rec in self.hog.all_hogs.items():
            self.hogs_by_node[rec.node_id].append(rec)
            parts = hog_id.split('.')
            self._hog_node_cache[hog_id] = parts[1] if len(parts) >= 3 else ''
            for gene in rec.genes:
                self.gene_to_hog[gene] = hog_id
                self.node_gene_to_hog[rec.node_id][gene] = hog_id
                self.gene_to_all_hogs[gene].append(hog_id)
                if '|' in gene:
                    self.hog_species[hog_id].add(gene.split('|')[0])
                else:
                    self.hog_species[hog_id].add('unknown')
            if rec.parent:
                self.hog_parent[hog_id] = rec.parent
                self.hog_children[rec.parent].append(hog_id)

        logger.info("Indexed {} genes into HOGs".format(len(self.gene_to_hog)))
        logger.info("HOGs per node: {}".format(
            {k: len(v) for k, v in self.hogs_by_node.items()}))

    def _build_leaf_graphs(self):
        """Build AncestralAdjacencyGraph for each leaf species"""
        if self.gfffile is None:
            raise ValueError("gfffile is required for AKR")

        # 读取目标染色体列表（如果提供）
        allowed_chroms = None
        if self.chrom_list:
            with open(self.chrom_list) as f:
                allowed_chroms = {line.strip().split()[0] for line in f}

        gff = Gff(self.gfffile)
        d_chrom = OrderedDict()
        gff_species = set()
        for line in gff:
            sp = line.species
            gff_species.add(sp)
            if sp not in self.hog.species:
                continue
            try:
                d_chrom[(sp, line.chrom)] += [line]
            except KeyError:
                d_chrom[(sp, line.chrom)] = [line]

        skipped_sps = gff_species - set(self.hog.species)
        if skipped_sps:
            logger.info("GFF contains {} species not in tree, skipped: {}".format(
                len(skipped_sps), ', '.join(sorted(skipped_sps))))

        for (sp, chrom), lines in d_chrom.items():
            lines = sorted(lines, key=lambda x: x.start)
            for i, line in enumerate(lines):
                line.index = i

        # 诊断：检查 gene_to_hog 与 GFF 基因 ID 的格式匹配（样本交集比例）
        sample_hog_keys = list(self.gene_to_hog.keys())[:5]
        sample_gff_ids = []
        for line in Gff(self.gfffile):
            if line.species in self.hog.species:
                sample_gff_ids.append(line.id)
            if len(sample_gff_ids) >= 100:
                break
        has_pipe_hog = any('|' in k for k in sample_hog_keys)
        has_pipe_gff = any('|' in k for k in sample_gff_ids)
        matched = sum(1 for gid in sample_gff_ids if gid in self.gene_to_hog)
        ratio = matched / len(sample_gff_ids) * 100 if sample_gff_ids else 0
        logger.info("Format check: sample intersection ratio {}/{} = {:.1f}% (pipe_hog={}, pipe_gff={})".format(
            matched, len(sample_gff_ids), ratio, has_pipe_hog, has_pipe_gff))

        # 自动检测并修复格式不匹配
        if has_pipe_hog and not has_pipe_gff:
            logger.info("Format mismatch detected: HOG has SPECIES|GENEID, GFF has plain GENEID. Auto-fixing...")
            self._gene_to_hog_fixed = {}
            for gid, hid in self.gene_to_hog.items():
                plain = gid.split('|', 1)[1] if '|' in gid else gid
                self._gene_to_hog_fixed[plain] = hid
            self.gene_to_hog = self._gene_to_hog_fixed
            for node_id, mapping in self.node_gene_to_hog.items():
                fixed = {}
                for gid, hid in mapping.items():
                    plain = gid.split('|', 1)[1] if '|' in gid else gid
                    fixed[plain] = hid
                self.node_gene_to_hog[node_id] = fixed
        elif not has_pipe_hog and has_pipe_gff:
            logger.info("Format mismatch detected: HOG has plain GENEID, GFF has SPECIES|GENEID. Auto-fixing...")
            self._gene_to_hog_fixed = {}
            for gid, hid in self.gene_to_hog.items():
                self._gene_to_hog_fixed[gid] = hid
            self.gene_to_hog = self._gene_to_hog_fixed

        for sp in self.hog.species:
            gg = GffGraph()
            sp_chroms = {k: v for k, v in d_chrom.items() if k[0] == sp}
            if not sp_chroms:
                logger.warning("Species {} not found in GFF, skipping".format(sp))
                continue
            for (sp_name, chrom), lines in sp_chroms.items():
                if allowed_chroms is not None and chrom not in allowed_chroms:
                    continue
                if len(lines) < self.min_genes:
                    continue
                gg.add_path(lines)

            # 删除不在任何 HOG 中的孤儿基因（orphan genes），使用 GffGraph.remove_internals 两侧重连
            # 注意：用 gene_to_all_hogs 判断（包含所有层级），避免被高层级 HOG 覆盖导致误判
            orphan_genes = [line for line in gg.nodes() if line.id not in self.gene_to_all_hogs]
            if orphan_genes:
                gg.remove_internals(orphan_genes)
                logger.info("Species {}: removed {} orphan genes from GffGraph".format(
                    sp, len(orphan_genes)))

            # 构建带端粒的叶子邻接图
            aag = AncestralAdjacencyGraph.from_gffgraph(sp, gg, add_telomere=True)
            aag.species_set = {sp}

            # 建立基因到 HOG 的映射（叶子层级）
            mapped = 0
            for gene in aag.gene_nodes:
                if gene.id in self.gene_to_all_hogs:
                    # 取该基因的第一个 HOG（通常是叶子层级）作为起点映射
                    hog_id = self.gene_to_all_hogs[gene.id][0]
                    aag.hog_map[gene] = hog_id
                    gene.hog_id = hog_id
                    mapped += 1
            self.leaf_graphs[sp] = aag
            self.anc_graphs[sp] = aag
            logger.info("Built leaf graph for {}: {} genes, {} chromosomes".format(
                sp, len(aag.gene_nodes), len(list(aag.chromosomes))))

    def _draw_step_dotplots(self, parent_aag, mapped_children, child_graphs,
                            child_source_ids, target_node_id):
        """每重建一步输出父节点与每个子节点的dotplot。

        Parameters:
            parent_aag: 重建后的祖先图
            mapped_children: 映射后的子图列表（含hog_map）
            child_graphs: 原始子图列表
            child_source_ids: 子节点ID列表
            target_node_id: 当前重建的节点ID
        """
        try:
            import matplotlib
            matplotlib.use('Agg')
        except ImportError:
            return

        outdir = os.path.join(os.path.dirname(self.outpre) or '.', 'dotplots')
        os.makedirs(outdir, exist_ok=True)

        # 祖先核型：使用HOG ID
        parent_karyo = aag_to_karyo(parent_aag)
        if not parent_karyo:
            return

        for i, (mc, cg, cid) in enumerate(zip(mapped_children, child_graphs, child_source_ids)):
            # 子节点核型：叶子用gene ID，祖先用hog_id
            child_karyo = aag_to_karyo(cg)
            if not child_karyo:
                continue

            # 构建子→父映射：子节点的gene/hog ID → 父节点HOG ID
            gene_mapping = {}
            for gene, hog_rec in mc.hog_map.items():
                if hasattr(gene, 'hog_id'):
                    key = gene.hog_id
                else:
                    key = str(gene)
                gene_mapping[key] = hog_rec.hog_id

            outpath = os.path.join(outdir,
                                   "dotplot_{}_vs_{}.png".format(target_node_id, cid))
            draw_dotplot(parent_karyo, child_karyo, outpath,
                         ref_name=target_node_id,
                         query_name=cid,
                         gene_mapping=gene_mapping,
                         dpi=150)
            logger.info("  Dotplot: {} vs {} -> {}".format(target_node_id, cid, outpath))

        for path in paths:
            # 添加节点（保留strand等属性）
            for h in path:
                if h not in seen_nodes:
                    new_graph.add_node(h, **dict(merged.graph.nodes[h]))
                    result.gene_nodes.add(h)
                    result.hog_map[h] = h
                    seen_nodes.add(h)
                    # 保留node_sources
                    if h in merged.node_sources:
                        result.node_sources[h] = set(merged.node_sources[h])

            # 添加有向边
            for i in range(len(path) - 1):
                new_graph.add_edge(path[i], path[i + 1])
                # 保留edge_sources
                key_fwd = (path[i], path[i + 1])
                key_rev = (path[i + 1], path[i])
                if key_fwd in merged.edge_sources:
                    result.edge_sources[key_fwd] = set(merged.edge_sources[key_fwd])
                if key_rev in merged.edge_sources:
                    result.edge_sources[key_rev] = set(merged.edge_sources[key_rev])

        result.graph = new_graph
        return result

    def _map_child_node_to_hog(self, node_id, child_graph):
        """
        将子节点图的节点映射到当前节点(node_id)的 HOGrecord 对象。
        对于祖先子节点：沿 hog_parent 链向上爬，直到找到属于当前节点的 HOG。
        对于叶子子节点：通过 node_gene_to_hog[node_id] 直接查找。
        返回 dict{原节点: HOGrecord}。
        """
        node_hog_map = self.node_gene_to_hog.get(node_id, {})
        node_hog_records = {rec.hog_id: rec for rec in self.hogs_by_node.get(node_id, [])}
        mapping = {}
        for n in child_graph.gene_nodes:
            if isinstance(n, HOGrecord):
                # 子节点是祖先图，节点是 HOGrecord：沿 parent 链向上爬
                current = n.hog_id
                while current:
                    if current in node_hog_records:
                        mapping[n] = node_hog_records[current]
                        break
                    current = self.hog_parent.get(current)
            elif isinstance(n, str):
                # 子节点是祖先图，节点是字符串 HOG ID：沿 parent 链向上爬
                current = n
                while current:
                    if current in node_hog_records:
                        mapping[n] = node_hog_records[current]
                        break
                    current = self.hog_parent.get(current)
            else:
                # 子节点是叶子图，节点是基因对象
                gid = getattr(n, 'id', str(n))
                if gid in node_hog_map:
                    hid = node_hog_map[gid]
                    if hid in node_hog_records:
                        mapping[n] = node_hog_records[hid]
        return mapping

    def _map_to_parent_hogs(self, node_id, graph, source_id=None):
        """
        将子图（叶子图或祖先图）的节点映射到指定 node_id 的 HOGrecord，
        返回一个新的 AncestralAdjacencyGraph，节点全部为 HOGrecord 对象。
        边按照原图的邻接关系复制，去除自环和端粒边。
        source_id: 标记映射后图中节点/边的来源（如子节点名称），
                   用于合并后追踪染色体来源。
        """
        mapping = self._map_child_node_to_hog(node_id, graph)
        mapped = AncestralAdjacencyGraph(node_id="{}_mapped".format(graph.node_id))
        mapped.species_set = set(graph.species_set)
        src = source_id or graph.node_id

        # 添加节点：只添加成功映射的 HOGrecord
        for rec in set(mapping.values()):
            mapped.graph.add_node(rec, hog=rec)
            mapped.gene_nodes.add(rec)
            mapped.node_sources[rec] = {src}

        # 保留原始映射（基因->HOG），用于后续 duplication 检测
        for gene, rec in mapping.items():
            mapped.hog_map[gene] = rec

        # 记录每个 HOG 的链向（来自子节点原始基因或祖先图节点属性）
        # 对映射到同一 HOG 的多个基因进行投票，保留多数链向
        strand_votes = defaultdict(lambda: {'+': 0, '-': 0})
        for gene, rec in mapping.items():
            # 优先从图节点属性读取（祖先图已设置 strand），
            # 否则从基因对象本身读取（GffLine 有 strand 属性）
            strand = graph.graph.nodes[gene].get('strand')
            if strand is None:
                strand = getattr(gene, 'strand', '+')
            strand_votes[rec][strand] += 1
        for rec, votes in strand_votes.items():
            consensus = '+' if votes['+'] >= votes['-'] else '-'
            mapped.graph.nodes[rec]['strand'] = consensus

        # 记录每个HOG在子节点中的端点信息（telomere-centric追踪）
        # hog_endpoints: hog_rec -> {'left': [(chrom_id, tel_left), ...], 'right': [...]}
        mapped.hog_endpoints = defaultdict(lambda: {'left': [], 'right': []})
        # 记录每个HOG在子节点中的染色体位置信息（从原始图计算，不依赖mc.chromosomes遍历）
        # hog_to_chrom_idx: hog_rec -> chrom_idx (在该子节点中的染色体编号)
        # hog_chrom_ends: hog_rec -> 'left'/'right'/'both'/'internal' (在该染色体上的位置)
        # hog_to_pos: hog_rec -> int (在该染色体映射后HOG中的相对位置)
        # chrom_hogs: chrom_idx -> [hog_rec, ...] (按原始基因顺序排列的映射后HOG列表)
        mapped.hog_to_chrom_idx = {}
        mapped.hog_chrom_ends = {}
        mapped.hog_to_pos = {}
        mapped.chrom_hogs = {}
        for cci, chrom in enumerate(graph.chromosomes):
            if not chrom:
                continue
            genes = [n for n in chrom if n not in graph.telomeres]
            if not genes:
                continue
            left_tel = chrom[0] if chrom[0] in graph.telomeres else None
            right_tel = chrom[-1] if chrom[-1] in graph.telomeres else None
            # 左端点：从最左端开始找第一个映射成功的基因
            left_hog = None
            for g in genes:
                if g in mapping:
                    left_hog = mapping[g]
                    break
            if left_hog is not None:
                mapped.hog_endpoints[left_hog]['left'].append((graph.node_id, left_tel))
            # 右端点：从最右端开始找第一个映射成功的基因
            right_hog = None
            for g in reversed(genes):
                if g in mapping:
                    right_hog = mapping[g]
                    break
            if right_hog is not None:
                mapped.hog_endpoints[right_hog]['right'].append((graph.node_id, right_tel))
            # 记录每个映射成功的HOG在该染色体上的位置
            mapped_genes = []
            for j, g in enumerate(genes):
                if g in mapping:
                    hog = mapping[g]
                    mapped_genes.append(hog)
                    mapped.hog_to_chrom_idx[hog] = cci
            # 计算端点类型和位置：基于映射后的HOG在该染色体上的相对位置
            if mapped_genes:
                n_mapped = len(mapped_genes)
                mapped.chrom_hogs[cci] = mapped_genes
                for k, hog in enumerate(mapped_genes):
                    mapped.hog_to_pos[hog] = k
                    if n_mapped == 1:
                        mapped.hog_chrom_ends[hog] = 'both'
                    elif k == 0:
                        mapped.hog_chrom_ends[hog] = 'left'
                    elif k == n_mapped - 1:
                        mapped.hog_chrom_ends[hog] = 'right'
                    else:
                        mapped.hog_chrom_ends[hog] = 'internal'

        # 复制边：按原始染色体方向添加有向边
        # 统计每个有向 HOG 对的支持数（保留方向信息）
        edge_support = Counter()
        tandem_dup_counts = Counter()
        for n1, n2 in graph.get_adjacencies(include_telomere=False):
            if n1 in mapping and n2 in mapping:
                h1 = mapping[n1]
                h2 = mapping[n2]
                if h1 == h2:
                    tandem_dup_counts[h1] += 1
                    # 保留自环边，让 dedup 阶段处理
                edge_support[(h1, h2)] += 1

        mapped.tandem_dup_counts = tandem_dup_counts

        for (h1, h2), count in edge_support.items():
            mapped.graph.add_edge(h1, h2, support=count)
            mapped.edge_sources[(h1, h2)] = {src}
            mapped.edge_sources[(h2, h1)] = {src}

        # 记录原始染色体数（映射前后不变）
        mapped.n_chromosomes = len(list(graph.chromosomes))

        # 重建端粒，使 chromosomes property 正确返回染色体列表
        mapped._add_telomeres()

        return mapped

    def _merge_child_graphs(self, node_id, mapped_children, child_source_ids=None):
        """
        合并多个已映射到同一层 HOG 的邻接图。
        节点取并集，边按支持度投票：
        - support=子图数：所有子图都支持，高置信度保留
        - support<子图数：仅部分子图支持，标记为冲突边
        同时标记 conflict 边是否涉及子节点染色体端点（端点邻接更可能是重排断点），
        以及冲突边在其他子节点中是否位于同一染色体（同一染色体更可能是祖先邻接断裂）。
        返回合并后的 AncestralAdjacencyGraph。
        """
        n_children = len(mapped_children)
        if child_source_ids is None:
            child_source_ids = [mc.node_id for mc in mapped_children]

        merged = AncestralAdjacencyGraph(node_id=node_id)
        merged.species_set = set()
        for g in mapped_children:
            merged.species_set |= g.species_set

        # 节点并集
        all_hogs = set()
        for g in mapped_children:
            all_hogs |= set(g.gene_nodes)
        for rec in all_hogs:
            merged.graph.add_node(rec, hog=rec)
            merged.gene_nodes.add(rec)
            merged.hog_map[rec] = rec

        # 记录节点来源（可能来自多个子图）
        for i, g in enumerate(mapped_children):
            src = child_source_ids[i]
            for rec in g.gene_nodes:
                if rec not in merged.node_sources:
                    merged.node_sources[rec] = set()
                merged.node_sources[rec].add(src)

        # 预计算每个子图的 HOG -> 染色体索引 / 位置
        child_chrom_info = []
        for g in mapped_children:
            h2c = {}
            h2p = {}
            for ci, chrom in enumerate(g.chromosomes):
                genes = [n for n in chrom if n not in g.telomeres]
                for pos, h in enumerate(genes):
                    h2c[h] = ci
                    h2p[h] = pos
            child_chrom_info.append((h2c, h2p))

        # 边支持度统计（无向，按 hog_id 字符串排序作为 key）
        edge_support = Counter()
        edge_endpoint = {}  # key -> bool (是否涉及子节点端点)
        edge_children = defaultdict(set)  # key -> {child_index}
        for i, g in enumerate(mapped_children):
            hog_endpoints = getattr(g, 'hog_endpoints', {})
            for n1, n2 in g.get_adjacencies(include_telomere=False):
                key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                edge_support[key] += 1
                edge_children[key].add(i)
                # 若任一端点是子节点染色体端点，标记为端点边
                for h in (n1, n2):
                    ep = hog_endpoints.get(h, {})
                    if ep.get('left') or ep.get('right'):
                        edge_endpoint[key] = True

        # 记录边来源
        for key, child_indices in edge_children.items():
            h1, h2 = key
            srcs = {child_source_ids[i] for i in child_indices}
            merged.edge_sources[(h1, h2)] = srcs
            merged.edge_sources[(h2, h1)] = srcs

        # 为 conflict 边计算：在其他子节点中是否位于同一染色体及距离
        same_chrom_other = {}
        dist_other = {}
        for key, count in edge_support.items():
            if count >= n_children:
                continue
            child_idx = list(edge_children[key])[0]
            h1, h2 = key
            # 在所有其他子图中检查是否同染色体，取最小距离
            min_dist = 999999
            any_same = False
            for other_idx in range(n_children):
                if other_idx == child_idx:
                    continue
                other_h2c = child_chrom_info[other_idx][0]
                other_h2p = child_chrom_info[other_idx][1]
                c1 = other_h2c.get(h1)
                c2 = other_h2c.get(h2)
                if c1 is not None and c2 is not None and c1 == c2:
                    any_same = True
                    p1 = other_h2p.get(h1, 0)
                    p2 = other_h2p.get(h2, 0)
                    d = abs(p1 - p2)
                    if d < min_dist:
                        min_dist = d
            same_chrom_other[key] = any_same
            dist_other[key] = min_dist if any_same else 999999

        for (h1, h2), count in edge_support.items():
            is_endpoint = edge_endpoint.get((h1, h2), False)
            merged.graph.add_edge(h1, h2, support=count,
                                  conflict=(count < n_children), endpoint=is_endpoint)
            merged.graph.add_edge(h2, h1, support=count,
                                  conflict=(count < n_children), endpoint=is_endpoint)

        # 合并 hog_endpoints，保留端粒位置信息供后续端粒约束使用
        merged.hog_endpoints = defaultdict(lambda: {'left': [], 'right': []})
        for g in mapped_children:
            if hasattr(g, 'hog_endpoints'):
                for h, ends in g.hog_endpoints.items():
                    if h in merged.gene_nodes:
                        merged.hog_endpoints[h]['left'].extend(ends['left'])
                        merged.hog_endpoints[h]['right'].extend(ends['right'])

        # 将 conflict 边排序辅助信息附加到 merged 图对象上，供 _linearize_graph 使用
        merged.same_chrom_other = same_chrom_other
        merged.dist_other = dist_other

        return merged

    def _find_descendant_hog(self, source_hog_id, target_node_id):
        """
        从 source_hog_id 向下 BFS 查找属于 target_node_id 的后代 HOG。
        用于外类群反向映射：外类群基因属于高层 HOG，需找到当前节点层级的对应 HOG。
        """
        if self._hog_node_cache.get(source_hog_id) == target_node_id:
            return source_hog_id
        queue = [source_hog_id]
        visited = {source_hog_id}
        while queue:
            current = queue.pop(0)
            if self._hog_node_cache.get(current) == target_node_id:
                return current
            for child in self.hog_children.get(current, []):
                if child not in visited:
                    visited.add(child)
                    queue.append(child)
        return None

    def _map_outgroup_to_current_hogs(self, node_id, og_graph):
        """
        将外类群图映射到当前节点的 HOGrecord（反向映射），返回新的 AncestralAdjacencyGraph。
        外类群不是当前节点的后代，其基因不在 node_gene_to_hog[node_id] 中。
        策略：
        1. 若外类群是祖先图（节点为 HOGrecord），沿 hog_parent 向上爬直到当前节点。
        2. 若外类群是叶子图（节点为基因），通过 gene_to_all_hogs 找到该基因的所有层级 HOG，
           再向下 BFS 查找属于当前节点的后代 HOG。
        映射后，按原图的邻接关系重建边（去除自环和端粒）。
        """
        node_hog_records = {rec.hog_id: rec for rec in self.hogs_by_node.get(node_id, [])}
        mapping = {}
        for n in og_graph.gene_nodes:
            if isinstance(n, HOGrecord):
                # 外类群是祖先图：沿 parent 链向上爬
                current = n.hog_id
                while current:
                    if current in node_hog_records:
                        mapping[n] = node_hog_records[current]
                        break
                    current = self.hog_parent.get(current)
            else:
                # 外类群是叶子图：自上而下反向查找
                gid = getattr(n, 'id', str(n))
                source_hogs = self.gene_to_all_hogs.get(gid, [])
                target_hog = None
                for shog in source_hogs:
                    target = self._find_descendant_hog(shog, node_id)
                    if target:
                        target_hog = target
                        break
                if target_hog and target_hog in node_hog_records:
                    mapping[n] = node_hog_records[target_hog]

        # 构建映射后的图
        mapped = AncestralAdjacencyGraph(node_id="{}_og_mapped".format(og_graph.node_id))
        mapped.species_set = set(og_graph.species_set)
        for rec in set(mapping.values()):
            mapped.graph.add_node(rec, hog=rec)
            mapped.gene_nodes.add(rec)
            mapped.hog_map[rec] = rec
        seen_edges = set()
        for n1, n2 in og_graph.get_adjacencies(include_telomere=False):
            if n1 in mapping and n2 in mapping:
                h1 = mapping[n1]
                h2 = mapping[n2]
                if h1 == h2:
                    continue
                key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
                if key in seen_edges:
                    continue
                seen_edges.add(key)
                mapped.graph.add_edge(h1, h2)
        return mapped

    def _filter_fragments(self, result, aag, target_chromosomes):
        """
        端粒约束后处理：丢弃过短的碎片路径。
        核心原则：只有NCF/EEJ/fission应改变染色体数量。
        保留主染色体（最长的 target 条），并额外保留达到一定长度的路径，
        避免过度丢弃基因。短碎片不判定为染色体。
        """
        if target_chromosomes is None:
            return result

        paths = list(result.chromosomes)
        n_paths = len(paths)
        if n_paths <= target_chromosomes:
            return result

        path_info = []
        total_genes = len(aag.gene_nodes)
        for p in paths:
            n_genes = len([n for n in p if n not in result.telomeres])
            path_info.append((n_genes, p))
        path_info.sort(key=lambda x: x[0], reverse=True)

        # 动态保留策略：优先满足目标染色体数，同时保证基本覆盖率
        # 若路径数远超目标，采用更激进的过滤；若接近目标，放宽保留
        excess = n_paths - target_chromosomes
        max_paths = max(target_chromosomes * 2, target_chromosomes + 5)
        if excess > target_chromosomes * 3:
            # 严重碎片化：严格限制额外路径，阈值提高
            max_paths = min(max_paths, target_chromosomes + 3)
            threshold = max(5, total_genes // (target_chromosomes * 4))
        elif excess > target_chromosomes:
            max_paths = min(max_paths, target_chromosomes + 4)
            threshold = max(5, total_genes // (target_chromosomes * 6))
        else:
            threshold = max(5, total_genes // (target_chromosomes * 8))

        keep_paths = []
        kept_genes = 0
        for n_genes, p in path_info:
            if len(keep_paths) < target_chromosomes:
                keep_paths.append(p)
                kept_genes += n_genes
            elif len(keep_paths) >= max_paths:
                break
            elif n_genes >= threshold:
                keep_paths.append(p)
                kept_genes += n_genes
            else:
                break

        # 若主染色体覆盖率过低(<55%)，按长度顺序继续补充，最多到max_paths
        if kept_genes < total_genes * 0.55 and len(keep_paths) < len(path_info):
            for n_genes, p in path_info[len(keep_paths):]:
                if len(keep_paths) >= max_paths:
                    break
                keep_paths.append(p)
                kept_genes += n_genes
                if kept_genes >= total_genes * 0.55:
                    break

        keep_nodes = set()
        for p in keep_paths:
            keep_nodes.update(p)

        # 重建有向图
        filtered = nx.DiGraph()
        for n in keep_nodes:
            filtered.add_node(n, **dict(result.graph.nodes[n]))
        for u, v in result.graph.edges():
            if u in keep_nodes and v in keep_nodes:
                filtered.add_edge(u, v)

        result.graph = filtered
        result.gene_nodes = keep_nodes - result.telomeres
        result.events = list(aag.events)
        result.chrom_map = {}
        for n in result.gene_nodes:
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom

        discarded = total_genes - len(result.gene_nodes)
        logger.info("  Telomere filter: {} paths -> {} paths (discarded {} orphan genes, coverage {:.1%})".format(
            n_paths, len(keep_paths), discarded, len(result.gene_nodes) / max(total_genes, 1)))
        return result

    def _export_results(self):
        """导出祖先核型和重排事件"""
        for node_id, aag in self.anc_graphs.items():
            if node_id in self.leaf_graphs:
                continue
            prefix = "{}.{}".format(self.outpre, node_id)
            # 直接为祖先图生成 wgdi 格式的 gff 和 lens 文件
            try:
                self._export_aag_wgdi(aag, prefix)
                logger.info("Exported {} to {}".format(node_id, prefix))
            except Exception as e:
                logger.warning("Failed to export {}: {}".format(node_id, e))

        # Build parent_of mapping for branch-level output
        parent_of = {}
        if self.tree is not None:
            for node in self.tree.traverse("preorder"):
                if not node.is_root():
                    parent_of[node.name] = node.up.name

        event_file = "{}.events.tsv".format(self.outpre)
        with open(event_file, 'w') as f:
            header = ["branch", "event_type", "genes", "chroms", "desc", "support", "child_source"]
            f.write('\t'.join(header) + '\n')

            def _write_events_from_graph(node_id, aag, parent):
                """Helper: write events from a single graph, compute branch."""
                for ev in aag.events:
                    # v4 events use TAKREvent with ev.branch, v3 uses node_id
                    if hasattr(ev, 'branch') and ev.branch:
                        branch = ev.branch
                    else:
                        branch = "{}-{}".format(parent, node_id)
                    genes = ','.join(str(g) for g in ev.genes_involved)
                    # parent_chroms only exists on v3 RearrangementEvent
                    if hasattr(ev, 'parent_chroms') and ev.parent_chroms:
                        chroms = ','.join(str(c) for c in ev.parent_chroms)
                    else:
                        chroms = ''
                    cs = ev.child_source if hasattr(ev, 'child_source') and ev.child_source else ''
                    line = [branch, ev.event_type, genes, chroms, ev.desc, str(ev.support), cs]
                    f.write('\t'.join(line) + '\n')

            # Ancestral nodes (internal)
            for node_id, aag in self.anc_graphs.items():
                if node_id in self.leaf_graphs:
                    continue
                parent = parent_of.get(node_id, '?')
                _write_events_from_graph(node_id, aag, parent)

            # Pre-WGD pseudo-nodes — events detected between pre and post
            for pre_name, aag in self.pre_wgd_graphs.items():
                if aag.events:
                    parent = parent_of.get(pre_name, '?')
                    _write_events_from_graph(pre_name, aag, parent)
            for leaf_name, aag in self.leaf_graphs.items():
                if aag.events:
                    parent = parent_of.get(leaf_name, '?')
                    _write_events_from_graph(leaf_name, aag, parent)
        logger.info("Exported rearrangement events to {}".format(event_file))

        # Draw dotplots
        outdir = os.path.dirname(self.outpre) or '.'
        try:
            draw_akr_dotplots(self, outdir)
            logger.info("Exported dotplots to {}".format(outdir))
        except Exception as e:
            logger.warning("Failed to export dotplots: {}".format(e))

    def _export_aag_wgdi(self, aag, prefix):
        """将 AncestralAdjacencyGraph 导出为 wgdi 格式的 gff 和 lens"""
        fgff = open(prefix + '.gff', 'w')
        flen = open(prefix + '.lens', 'w')
        chrom_idx = 0
        for chrom in aag.chromosomes:
            real_nodes = [n for n in chrom if n not in aag.telomeres]
            if not real_nodes:
                continue
            chrom_idx += 1
            chrom_name = "{}_chrom_{}".format(aag.node_id, chrom_idx)
            chrom_end = 0
            for i, node in enumerate(real_nodes):
                idx = i + 1
                start = idx * 100 - 99
                end = idx * 100
                chrom_end = end
                if isinstance(node, HOGrecord):
                    gene_id = node.hog_id
                    strand = getattr(node, 'strand', '+')
                else:
                    gene_id = str(node)
                    strand = '+'
                line = [chrom_name, gene_id, start, end, strand, idx, gene_id]
                fgff.write('\t'.join(map(str, line)) + '\n')
            flen.write('\t'.join(map(str, [chrom_name, chrom_end, len(real_nodes)])) + '\n')
        fgff.close()
        flen.close()


