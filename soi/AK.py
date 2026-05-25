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
from .chromosome_path_cover import ChromosomePathCover, GreedyPathCover, ORTOOLS_AVAILABLE
from .telomere_constraint import TelomereConstraint
from .validator import TAKRValidator

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

    def get_chromosome_ends(self):
        ends = []
        for chrom in self.chromosomes:
            if not chrom:
                continue
            left = chrom[0] if chrom[0] in self.telomeres else None
            right = chrom[-1] if chrom[-1] in self.telomeres else None
            genes = [n for n in chrom if n not in self.telomeres]
            first_gene = genes[0] if genes else None
            last_gene = genes[-1] if genes else None
            ends.append((left, first_gene, last_gene, right))
        return ends

    def get_adjacencies(self, include_telomere=False):
        adjs = set()
        for n1, n2 in self.graph.edges():
            if not include_telomere:
                if n1 in self.telomeres or n2 in self.telomeres:
                    continue
            adjs.add((n1, n2))
        return adjs

    def get_telomere_adjacencies(self):
        adjs = []
        for n1, n2 in self.graph.edges():
            if n1 in self.telomeres or n2 in self.telomeres:
                adjs.append((n1, n2))
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

    def to_gffgraph(self):
        """转换为 GffGraph，自动去重双向边（只保留单向）"""
        gg = GffGraph()
        for n in self.graph.nodes():
            gg.add_node(n, **dict(self.graph.nodes[n]))
        seen = set()
        for n1, n2 in self.graph.edges():
            if (n2, n1) not in seen:
                gg.add_edge(n1, n2, **dict(self.graph[n1][n2]))
                seen.add((n1, n2))
        return gg

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
                 use_v3=True,
                 use_v4=False,
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
        self.timeout = timeout  # 全局总超时秒数
        self.node_timeout = node_timeout  # 单节点重建超时秒数
        self.use_ilp_sa = use_ilp_sa and pulp is not None  # 是否使用ILP+SA混合线性化
        self.sa_iterations = sa_iterations
        # Unified algorithm selector
        # reconstruction_algorithm: 'v3' (CP-SAT), 'v4' (event-driven), 'v4_colored' (ColoredGraph)
        # Falls back to use_v3/use_v4 for backward compatibility
        self.use_v4_colored = False
        if reconstruction_algorithm is not None:
            if reconstruction_algorithm == 'v3':
                self.use_v3 = True and ORTOOLS_AVAILABLE
                self.use_v4 = False
            elif reconstruction_algorithm == 'v4':
                self.use_v3 = False
                self.use_v4 = True
            elif reconstruction_algorithm == 'v4_colored':
                self.use_v3 = False
                self.use_v4 = False
                self.use_v4_colored = True
            else:
                raise ValueError("Unknown reconstruction_algorithm: %s (use 'v3', 'v4', or 'v4_colored')" % reconstruction_algorithm)
        else:
            self.use_v3 = use_v3 and ORTOOLS_AVAILABLE
            self.use_v4 = use_v4
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

        # Handle polyploid leaves (v4_colored handles this internally)
        if not self.use_v4_colored:
            for leaf_name, ploidy in self.ploidy_map.items():
                if leaf_name in self.leaf_graphs and ploidy > 1:
                    self._collapse_polyploid_leaf(leaf_name)

        # Post-order traversal: reconstruct speciation nodes first, then pre-WGD pseudo-nodes
        for node in self.tree.traverse(strategy="postorder"):
            if node.is_leaf():
                continue
            elapsed = time.time() - self._start_time
            if self.timeout > 0 and elapsed > self.timeout:
                logger.warning("Global timeout ({}s) reached at {:.1f}s, skipping remaining nodes".format(
                    self.timeout, elapsed))
                break
            # Speciation node reconstruction
            if self.use_v4_colored:
                # ColoredGraph event-driven handles all nodes in one pass
                from soi.takr_event_driven import reconstruct_event_driven_v2
                anc_graphs = reconstruct_event_driven_v2(self)
                self.anc_graphs = anc_graphs
                break
            elif self.use_v4:
                # Event-driven reconstruction handles all nodes in one pass
                from soi.takr_event_driven import reconstruct_event_driven
                anc_graphs = reconstruct_event_driven(self)
                self.anc_graphs = anc_graphs
                break  # event-driven handles all nodes at once
            elif self.use_v3:
                self._reconstruct_node_v3(node)
            else:
                self._reconstruct_node(node)

            # WGD nodes: additional pre-WGD pseudo-node reconstruction (unified log format)
            if node.name in self.ploidy_map and self.ploidy_map[node.name] > 1:
                if self.use_v3:
                    self._collapse_wgd_v3(node)
                else:
                    self._collapse_wgd(node)

        # Skip optimization + top-down detection for v4/v4_colored (handled internally)
        if not self.use_v4 and not self.use_v4_colored:
            for i in range(self.rounds):
                logger.info('Optimization round {}'.format(i))
                self._optimize_round()

            self._detect_events_topdown()
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

    def _reconstruct(self, target_node_id, child_graphs, child_source_ids=None,
                     hog_level=None, outgroup_graphs=None,
                     target_chromosomes=None, is_pre_wgd=False):
        """
        统一重建方法：物种分化节点和 pre-WGD 节点共用同一流程。
        核心流程：
        1. 将每个子节点图映射到目标 HOG 层级
        2. 合并映射后的邻接图（产生冲突边）
        3. 外类群反向映射到目标 HOG（如有）
        4. 基于外类群投票解决合并冲突
        5. 线性化并重建端粒结构
        6. 逐个解决小规模重排
        7. 解决 telomere-centric 大规模重排

        Parameters:
            target_node_id: 重建目标节点名称（如 "N4" 或 "pre-WGD N4"）
            child_graphs: 子节点图列表（AncestralAdjacencyGraph）
            child_source_ids: 子图来源标识列表（用于追踪染色体来源）
            hog_level: 映射目标 HOG 层级的节点名（None 则用 target_node_id）
            outgroup_graphs: 外类群图列表 [(graph, weight), ...]（None 则不使用外类群）
            target_chromosomes: 目标染色体数目（None 则自动估计）
            is_pre_wgd: 是否为 pre-WGD 重建（影响部分后处理）
        """
        if child_source_ids is None:
            child_source_ids = [cg.node_id for cg in child_graphs]
        if hog_level is None:
            hog_level = target_node_id

        prefix = "pre-WGD " if is_pre_wgd else ""
        logger.info("Reconstructing node {}{}".format(prefix, target_node_id))
        t0 = time.time()

        # 打印子节点信息
        for i, cg in enumerate(child_graphs):
            n_chrom = len(list(cg.chromosomes))
            n_genes = len(cg.gene_nodes)
            n_edges = len(list(cg.get_adjacencies(include_telomere=False)))
            logger.info("  Child {}: {} genes, {} chromosomes, {} adjacencies".format(
                child_source_ids[i], n_genes, n_chrom, n_edges))

        # ============================================================
        # Step 1: 将每个子节点映射到目标 HOGrecord
        # ============================================================
        t1 = time.time()
        mapped_children = [self._map_to_parent_hogs(hog_level, cg, source_id=sid)
                           for cg, sid in zip(child_graphs, child_source_ids)]
        for i, mc in enumerate(mapped_children):
            n_edges = len(list(mc.get_adjacencies(include_telomere=False)))
            n_chrom = mc.n_chromosomes
            logger.info("  Mapped child {} -> {} HOGs, {} adjacencies, {} chromosomes ({:.2f}s)".format(
                child_source_ids[i], len(mc.gene_nodes), n_edges, n_chrom, time.time()-t1))

        # ============================================================
        # Step 2: 合并映射后的邻接图
        # ============================================================
        t1 = time.time()
        merged = self._merge_child_graphs(target_node_id, mapped_children, child_source_ids)
        n_conflict = sum(1 for u, v, d in merged.graph.edges(data=True) if d.get('conflict'))
        logger.info("  After merge: {} HOGs, {} conflict edges ({:.2f}s)".format(
            len(merged.gene_nodes), n_conflict, time.time()-t1))

        # 根据子节点投票确定每个 HOG 的共识链向
        strand_votes = defaultdict(lambda: {'+': 0, '-': 0})
        for mc in mapped_children:
            for rec in mc.gene_nodes:
                strand = mc.graph.nodes[rec].get('strand')
                if strand is None:
                    strand = getattr(rec, 'strand', '+')
                strand_votes[rec][strand] += 1
        for rec in merged.gene_nodes:
            votes = strand_votes.get(rec, {'+': 1, '-': 0})
            consensus = '+' if votes['+'] >= votes['-'] else '-'
            merged.graph.nodes[rec]['strand'] = consensus

        # 导出重建前 GFA
        t1 = time.time()
        merge_gfa = "{}.{}.merge.gfa".format(self.outpre, target_node_id)
        with open(merge_gfa, 'w') as fout:
            merged.to_gfa(fout)
        logger.info("  Exported pre-reconstruction GFA: {} ({:.2f}s)".format(merge_gfa, time.time()-t1))

        # ============================================================
        # Step 3: 外类群反向映射到目标 HOG（带系统发育距离权重）
        # ============================================================
        t1 = time.time()
        mapped_outgroups = []
        if outgroup_graphs:
            mapped_outgroups = [(self._map_outgroup_to_current_hogs(hog_level, og), weight)
                                for og, weight in outgroup_graphs]
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        logger.info("  Mapped {} outgroup graphs, total_weight={:.2f} ({:.2f}s)".format(
            len(mapped_outgroups), total_weight, time.time()-t1))

        # ============================================================
        # Step 4: 基于外类群投票解决合并冲突 + 全局线性化
        # ============================================================
        t1 = time.time()
        merged = self._resolve_merge_conflicts(merged, mapped_children, mapped_outgroups)
        t2 = time.time()
        if target_chromosomes is None:
            target_chromosomes = max(len(list(child_graphs[0].chromosomes)), 1) if child_graphs else 1
        if self.use_ilp_sa:
            merged = self._linearize_graph_ilp_sa(merged, target_chromosomes=target_chromosomes)
        else:
            merged = self._linearize_graph(merged, target_chromosomes=target_chromosomes)
        t3 = time.time()
        merged._add_telomeres()
        logger.info("  After conflict resolution: {} HOGs, {} chromosomes (conflict {:.2f}s, linearize {:.2f}s, telomere {:.2f}s)".format(
            len(merged.gene_nodes), len(list(merged.chromosomes)), t2-t1, t3-t2, time.time()-t3))

        # ============================================================
        # Step 5: 逐个解决小规模重排
        # ============================================================
        t1 = time.time()
        merged = self._resolve_indels(merged, mapped_children, mapped_outgroups)
        t2 = time.time()
        merged = self._resolve_duplications(merged, mapped_children)
        logger.info("  Small rearrangements: indel {:.2f}s, dup {:.2f}s".format(
            t2-t1, time.time()-t2))

        # ============================================================
        # Step 6: 检测重排事件（输入 vs 输出差异）
        # ============================================================
        t1 = time.time()
        self._detect_bottomup_events(merged, mapped_children, child_source_ids)
        logger.info("  Bottom-up event detection: {:.2f}s".format(time.time()-t1))

        # 记录染色体物种来源
        for ci, chrom in enumerate(merged.chromosomes):
            genes = [n for n in chrom if n not in merged.telomeres]
            species = set()
            for g in genes:
                srcs = merged.node_sources.get(g, set())
                for src in srcs:
                    # 从 child_graphs 中查找 source_id 对应的物种集
                    for j, sid in enumerate(child_source_ids):
                        if sid == src and j < len(child_graphs):
                            species |= child_graphs[j].species_set
            merged.chrom_species[ci] = species

        # 事件分类统计
        n_chrom = len(list(merged.chromosomes))
        all_event_types = ['indel', 'tandem_dup', 'proximal_dup', 'dispersed_dup',
                          'inversion', 'internal_inversion', 'telomere_inversion',
                          'fission', 'ncf', 'eej', 'unidir_trans',
                          'translocation', 'reciprocal_translocation']
        evt_counts = Counter(e.event_type for e in merged.events)
        evt_summary = ', '.join('{}:{}'.format(k, evt_counts.get(k, 0)) for k in all_event_types)

        n_fission = evt_counts.get('fission', 0)
        n_ncf = evt_counts.get('ncf', 0)
        n_eej = evt_counts.get('eej', 0)
        chrom_delta = n_fission - n_ncf - n_eej
        child_chroms = [len(list(cg.chromosomes)) for cg in child_graphs]
        avg_child_chrom = sum(child_chroms) / len(child_chroms) if child_chroms else n_chrom
        expected_delta = avg_child_chrom - n_chrom
        logger.info("Node {} final: {} chromosomes, events: {}".format(
            target_node_id, n_chrom, evt_summary))
        logger.info("  Chromosome validation: fission({}) - ncf({}) - eej({}) = delta {}".format(
            n_fission, n_ncf, n_eej, chrom_delta))
        logger.info("  Expected delta (child {:.1f} - anc {}) = {:.1f}, match={}".format(
            avg_child_chrom, n_chrom, expected_delta,
            "YES" if abs(chrom_delta - expected_delta) < 1.5 else "NO"))

        node_elapsed = time.time() - t0
        if self.node_timeout > 0 and node_elapsed > self.node_timeout:
            logger.warning("Node {} reconstruction took {:.1f}s, exceeding node_timeout {}s".format(
                target_node_id, node_elapsed, self.node_timeout))

        # 导出重建后 GFA
        anc_gfa = "{}.{}.anc.gfa".format(self.outpre, target_node_id)
        with open(anc_gfa, 'w') as fout:
            merged.to_gfa(fout)
        logger.info("  Exported post-reconstruction GFA: {}".format(anc_gfa))

        # 输出父节点与子节点的dotplot
        self._draw_step_dotplots(merged, mapped_children, child_graphs,
                                 child_source_ids, target_node_id)

        return merged

    def _reconstruct_node(self, node):
        """
        重建一个物种分化内部节点的祖先核型。
        委托给统一的 _reconstruct 方法。
        """
        node_id = node.name
        children = node.children
        if len(children) != 2:
            logger.warning("Node {} has {} children; expecting binary tree".format(
                node_id, len(children)))

        # 收集子节点图（叶子图或已重建的祖先图）
        # WGD节点优先使用pre-WGD图参与父节点重建
        child_graphs = []
        for child in children:
            child_id = child.name
            if child_id in self.pre_wgd_graphs:
                child_graphs.append(self.pre_wgd_graphs[child_id])
            elif child_id in self.anc_graphs:
                child_graphs.append(self.anc_graphs[child_id])
            elif child.is_leaf() and child_id in self.leaf_graphs:
                child_graphs.append(self.leaf_graphs[child_id])

        if not child_graphs:
            logger.warning("No child graphs for {}".format(node_id))
            return

        # 获取外类群图
        outgroup_graphs = self._get_outgroup_graphs(node)

        # 估计目标染色体数（使用外类群信号）
        target_chromosomes = self._estimate_target_chromosomes(node, outgroup_graphs)

        merged = self._reconstruct(
            target_node_id=node_id,
            child_graphs=child_graphs,
            hog_level=node_id,
            outgroup_graphs=outgroup_graphs,
            target_chromosomes=target_chromosomes,
            is_pre_wgd=False
        )

        self.anc_graphs[node_id] = merged

    def _reconstruct_node_v3(self, node):
        """
        使用v3 CP-SAT路径覆盖重建物种分化内部节点。
        映射、合并、外类群同现有流程，但线性化步骤替换为CP-SAT。
        """
        node_id = node.name
        children = node.children
        if len(children) != 2:
            logger.warning("Node {} has {} children; expecting binary tree".format(
                node_id, len(children)))

        child_graphs = []
        child_source_ids = []
        for child in children:
            child_id = child.name
            if child_id in self.pre_wgd_graphs:
                child_graphs.append(self.pre_wgd_graphs[child_id])
            elif child_id in self.anc_graphs:
                child_graphs.append(self.anc_graphs[child_id])
            elif child.is_leaf() and child_id in self.leaf_graphs:
                child_graphs.append(self.leaf_graphs[child_id])
            child_source_ids.append(child_id)

        if not child_graphs:
            logger.warning("No child graphs for {}".format(node_id))
            return

        outgroup_graphs = self._get_outgroup_graphs(node)
        target_chromosomes = self._estimate_target_chromosomes(node, outgroup_graphs)

        merged = self._reconstruct_v3(
            target_node_id=node_id,
            child_graphs=child_graphs,
            child_source_ids=child_source_ids,
            hog_level=node_id,
            outgroup_graphs=outgroup_graphs,
            target_chromosomes=target_chromosomes,
            is_pre_wgd=False
        )

        self.anc_graphs[node_id] = merged

    def _collapse_wgd_v3(self, node):
        """使用v3 CP-SAT路径覆盖重建pre-WGD基因组。"""
        node_id = node.name
        post_graph = self.anc_graphs.get(node_id)
        if post_graph is None:
            return
        parent = node.up
        if not parent:
            return

        parent_id = parent.name
        ploidy = self.ploidy_map.get(node_id, 2)
        target = max(len(list(post_graph.chromosomes)) // ploidy, 1)

        pre_id = "pre-WGD {}".format(node_id)
        pre_graph = self._reconstruct_v3(
            target_node_id=pre_id,
            child_graphs=[post_graph],
            child_source_ids=[node_id],
            hog_level=parent_id,
            outgroup_graphs=None,
            target_chromosomes=target,
            is_pre_wgd=True
        )

        pre_graph.node_id = "{}_pre".format(node_id)
        self.pre_wgd_graphs[node_id] = pre_graph
        n_pre = len(list(pre_graph.chromosomes))
        logger.info("  WGD pre-WGD v3 for {}: {} post chroms -> {} pre chroms (target={})".format(
            node_id, len(list(post_graph.chromosomes)), n_pre, target))

    def _reconstruct_v3(self, target_node_id, child_graphs, child_source_ids=None,
                        hog_level=None, outgroup_graphs=None,
                        target_chromosomes=None, is_pre_wgd=False):
        """
        统一重建方法v3：使用CP-SAT路径覆盖替代贪心线性化。

        Steps 1-3（映射、合并、外类群）同 _reconstruct，
        Step 4 替换为 CP-SAT 路径覆盖 + 端粒约束。
        """
        if child_source_ids is None:
            child_source_ids = [cg.node_id for cg in child_graphs]
        if hog_level is None:
            hog_level = target_node_id

        prefix = "pre-WGD " if is_pre_wgd else ""
        logger.info("Reconstructing node {}{} [v3]".format(prefix, target_node_id))
        t0 = time.time()

        # Step 1: 映射子节点到目标 HOGrecord
        t1 = time.time()
        mapped_children = [self._map_to_parent_hogs(hog_level, cg, source_id=sid)
                           for cg, sid in zip(child_graphs, child_source_ids)]
        for i, mc in enumerate(mapped_children):
            n_edges = len(list(mc.get_adjacencies(include_telomere=False)))
            n_chrom = mc.n_chromosomes
            logger.info("  Mapped child {} -> {} HOGs, {} adjacencies, {} chromosomes ({:.2f}s)".format(
                child_source_ids[i], len(mc.gene_nodes), n_edges, n_chrom, time.time()-t1))

        # Step 2: 合并映射后的邻接图
        t1 = time.time()
        merged = self._merge_child_graphs(target_node_id, mapped_children, child_source_ids)
        n_conflict = sum(1 for u, v, d in merged.graph.edges(data=True) if d.get('conflict'))
        logger.info("  After merge: {} HOGs, {} conflict edges ({:.2f}s)".format(
            len(merged.gene_nodes), n_conflict, time.time()-t1))

        # 共识链向
        strand_votes = defaultdict(lambda: {'+': 0, '-': 0})
        for mc in mapped_children:
            for rec in mc.gene_nodes:
                strand = mc.graph.nodes[rec].get('strand')
                if strand is None:
                    strand = getattr(rec, 'strand', '+')
                strand_votes[rec][strand] += 1
        for rec in merged.gene_nodes:
            votes = strand_votes.get(rec, {'+': 1, '-': 0})
            consensus = '+' if votes['+'] >= votes['-'] else '-'
            merged.graph.nodes[rec]['strand'] = consensus

        # Step 3: 外类群反向映射
        t1 = time.time()
        mapped_outgroups = []
        if outgroup_graphs:
            mapped_outgroups = [(self._map_outgroup_to_current_hogs(hog_level, og), weight)
                                for og, weight in outgroup_graphs]
        logger.info("  Mapped {} outgroup graphs ({:.2f}s)".format(
            len(mapped_outgroups), time.time()-t1))

        # Step 4: 外类群投票解决冲突
        t1 = time.time()
        merged = self._resolve_merge_conflicts(merged, mapped_children, mapped_outgroups)
        logger.info("  Conflict resolution: {:.2f}s".format(time.time()-t1))

        # Step 5: 构建端粒约束
        n_children = len(child_graphs)
        ploidy_ctx = self.ploidy_map.get(target_node_id) if is_pre_wgd else None
        tc = TelomereConstraint(merged.hog_endpoints, n_children=n_children,
                                ploidy_context=ploidy_ctx)
        logger.info("  {}".format(tc.summary()))

        # Step 6: 构建边权重和支持度
        edge_weights = {}
        support_counts = {}
        for u, v, d in merged.graph.edges(data=True):
            if u in merged.telomeres or v in merged.telomeres:
                continue
            if u == v:
                continue
            key = (u, v) if u.hog_id < v.hog_id else (v, u)
            if key not in edge_weights:
                support = d.get('support', 1)
                support_counts[key] = support
                og_weight = d.get('og_weight', 0)
                if support >= n_children:
                    weight = 200
                else:
                    weight = 100 + min(og_weight * 50, 100)
                edge_weights[key] = weight

        # Step 7: CP-SAT求解
        if target_chromosomes is None:
            target_chromosomes = max(
                len(list(child_graphs[0].chromosomes)), 1) if child_graphs else 1

        telomere_hogs = tc.get_telomere_hogs()
        telomere_weights = tc.get_telomere_weights()
        cpc = ChromosomePathCover(
            edge_weights=edge_weights,
            support_counts=support_counts,
            telomere_hogs=telomere_hogs,
            target_chromosomes=target_chromosomes,
            n_children=n_children,
            cpsat_timeout=120 if is_pre_wgd else 60,
            telomere_weights=telomere_weights,
        )
        paths, stats = cpc.solve()

        # Step 8: Fallback
        if paths is None:
            logger.warning("CP-SAT failed for {}, fallback to greedy".format(
                target_node_id))
            greedy = GreedyPathCover(edge_weights, telomere_hogs, target_chromosomes,
                                     telomere_weights=telomere_weights)
            paths, stats = greedy.solve()

        if paths is None:
            logger.error("Both CP-SAT and greedy failed for {}".format(
                target_node_id))
            return merged

        logger.info("  Path cover: {} paths, status={}".format(
            len(paths), stats.get('status', 'UNKNOWN')))

        # Step 9: 构建结果图
        result = self._paths_to_aag_v3(paths, merged, target_chromosomes)

        # Step 10: Small-scale rearrangement resolution (indel/duplication)
        # NOTE: v3 uses _resolve_indels_v3 to preserve CP-SAT path structure
        # Do NOT call _resolve_indels (it calls _linearize_graph which destroys paths)
        paths, indel_events, n_indel_paths_removed = self._resolve_indels_v3(
            paths, merged, mapped_children, mapped_outgroups)

        # Rebuild result from updated paths (paths may have changed after indel removal)
        result = self._paths_to_aag_v3(paths, merged, target_chromosomes)
        result.events.extend(indel_events)
        result = self._resolve_duplications(result, mapped_children)

        # Step 11: 事件检测
        result._add_telomeres()
        self._detect_bottomup_events(result, mapped_children, child_source_ids, mapped_outgroups)

        # Step 12: 染色体物种来源
        for ci, chrom in enumerate(result.chromosomes):
            genes = [n for n in chrom if n not in result.telomeres]
            species = set()
            for g in genes:
                srcs = result.node_sources.get(g, set())
                for src in srcs:
                    for j, sid in enumerate(child_source_ids):
                        if sid == src and j < len(child_graphs):
                            species |= child_graphs[j].species_set
            result.chrom_species[ci] = species

        # Step 13: 校验
        validator = TAKRValidator()
        child_chrom_counts = [mc.n_chromosomes for mc in mapped_children]
        is_ok, exp, act = validator.validate_chrom_event_balance(
            len(list(result.chromosomes)), result.events, child_chrom_counts,
            child_source_ids=child_source_ids,
            ploidy_context=ploidy_ctx, is_pre_wgd=is_pre_wgd,
            n_indel_paths_removed=n_indel_paths_removed)
        if not is_ok:
            logger.warning("  Chrom-event imbalance: expected={}, actual={}".format(
                exp, act))

        # 事件统计
        n_chrom = len(list(result.chromosomes))
        all_event_types = ['indel', 'tandem_dup', 'proximal_dup', 'dispersed_dup',
                          'inversion', 'internal_inversion', 'telomere_inversion',
                          'fission', 'ncf', 'eej', 'unidir_trans',
                          'reciprocal_translocation', 'unbalanced_reciprocal_translocation']
        evt_counts = Counter(e.event_type for e in result.events)
        evt_summary = ', '.join('{}:{}'.format(k, evt_counts.get(k, 0))
                                 for k in all_event_types)
        logger.info("  Node {}{} final [v3]: {} chromosomes, events: {}".format(
            prefix, target_node_id, n_chrom, evt_summary))

        # 导出GFA
        anc_gfa = "{}.{}.anc.gfa".format(self.outpre, target_node_id)
        with open(anc_gfa, 'w') as fout:
            result.to_gfa(fout)

        # Step 14: 输出父节点与子节点的dotplot
        self._draw_step_dotplots(result, mapped_children, child_graphs,
                                 child_source_ids, target_node_id)

        return result

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

    def _paths_to_aag_v3(self, paths, merged, target_chromosomes):
        """将CP-SAT路径列表转为AncestralAdjacencyGraph。"""
        result = AncestralAdjacencyGraph(node_id=merged.node_id)
        result.species_set = merged.species_set
        result.hog_endpoints = getattr(merged, 'hog_endpoints', {})
        result.same_chrom_other = getattr(merged, 'same_chrom_other', {})
        result.dist_other = getattr(merged, 'dist_other', {})
        result.events = list(merged.events)

        new_graph = nx.DiGraph()
        seen_nodes = set()

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

    def _estimate_target_chromosomes(self, node, outgroup_graphs=None):
        """
        估计目标染色体数目。
        使用子节点染色体数的最大值作为基础（祖先通常 ≥ 每个子节点），
        因为 NCF/EEJ 总是减少染色体数，fission 较少见。
        使用外类群染色体数作为第三信号，辅助推断是否发生了多次融合。
        """
        child_counts = []
        for child in node.children:
            child_name = child.name
            if child_name in self.anc_graphs:
                child_counts.append(len(list(self.anc_graphs[child_name].chromosomes)))
            elif child_name in self.leaf_graphs:
                child_counts.append(len(list(self.leaf_graphs[child_name].chromosomes)))

        clade_leaves = set(node.get_leaf_names())
        leaf_counts = []
        for leaf in clade_leaves:
            if leaf in self.leaf_graphs:
                leaf_counts.append(len(list(self.leaf_graphs[leaf].chromosomes)))

        all_counts = child_counts + leaf_counts
        if not all_counts:
            return 1

        # 使用子节点最大值作为基础（祖先 ≥ 每个子节点）
        base = max(child_counts) if child_counts else max(all_counts)

        # 如果两个子节点染色体数不同，较少的分支发生了更多融合，
        # 祖先可能比 max(child) 稍多
        if len(child_counts) >= 2:
            c_min, c_max = min(child_counts), max(child_counts)
            if c_max - c_min >= 2:
                # 差异较大时加1
                base += 1

        # 使用后代叶子最大值：祖先不应少于叶子最大值
        # （除非叶子发生了 WGD，但那是特殊情况）
        if leaf_counts:
            leaf_max = max(leaf_counts)
            if base < leaf_max:
                base = leaf_max

        # 使用外类群辅助推断
        if outgroup_graphs and len(outgroup_graphs) > 0:
            og_counts = [len(list(og.chromosomes)) for og, w in outgroup_graphs if og]
            if og_counts and len(og_counts) >= 2:
                # 外类群中位数（鲁棒估计，避免极端值影响）
                sorted_og = sorted(og_counts)
                n = len(sorted_og)
                if n % 2 == 0:
                    og_median = (sorted_og[n//2 - 1] + sorted_og[n//2]) / 2.0
                else:
                    og_median = float(sorted_og[n//2])
                og_median = int(og_median + 0.5)  # round

                # 如果子节点最大值 << 外类群中位数，暗示双方都发生了融合
                # 需要提升祖先染色体数估计
                child_max = max(child_counts) if child_counts else 0
                deficit = og_median - child_max
                if deficit >= 2:
                    # 双方都发生了大量融合：提升至外类群中位数
                    # 但受限于叶子最大值（避免过度推断）
                    if base < og_median:
                        logger.info(
                            "  Outgroup signal: child_max={} vs og_median={}, "
                            "raising ancestor estimate from {} to {}".format(
                                child_max, og_median, base, og_median))
                        base = og_median
                elif deficit >= 1 and len(child_counts) >= 2:
                    # 轻微差异：如果两个子节点都略低于外类群，加1
                    c_min2 = min(child_counts)
                    if c_min2 < og_median:
                        base += 1

        return max(base, 1)

    def _linearize_graph(self, aag, target_chromosomes=None):
        """
        基于目标染色体数的 Kruskal 风格线性化。
        1. 优先保留 support=2 边（高置信度祖先邻接）。
        2. 按外类群权重排序 conflict 边，仅当连接不同连通分量且度<2时才加入。
        3. 当连通分量数达到 target_chromosomes 时停止，避免过度连接。
        4. 对每个无向路径分量进行一致定向，解决方向反转导致的碎片化。
        """
        undirected = aag.graph.to_undirected()
        undirected.remove_nodes_from(aag.telomeres)
        new_graph = nx.DiGraph()

        if not undirected.nodes():
            result = AncestralAdjacencyGraph(node_id=aag.node_id)
            result.species_set = aag.species_set
            result.hog_map = dict(aag.hog_map)
            result.graph = new_graph
            result.events = list(aag.events)
            return result

        # 计算边权重，去重无向边
        edge_weights = {}
        og_weights = {}
        for u, v, d in aag.graph.edges(data=True):
            if u == v:
                continue
            if u in aag.telomeres or v in aag.telomeres:
                continue
            support = d.get('support', 1)
            og_weight = d.get('og_weight', 0)
            key = (u, v) if u.hog_id < v.hog_id else (v, u)
            if support >= 2:
                weight = 200
            else:
                # 外类群支持越强，权重越高（100 ~ 200）
                weight = 100 + min(og_weight * 50, 100)
            edge_weights[key] = weight
            og_weights[key] = og_weight

        # 分离 support=2 和 conflict 边
        support2_edges = []
        conflict_edges = []
        for key, weight in edge_weights.items():
            if weight >= 200:
                support2_edges.append((key, weight))
            else:
                conflict_edges.append((key, weight))

        # 按权重降序排序 conflict 边；涉及子节点端点的边降级（更可能是重排断点）
        # 优先选择在另一子节点中位于同一染色体且距离近的冲突边，这些边更可能是
        # 祖先邻接在另一子节点中被局部重排打断，而非衍生融合/易位。
        # 距离阈值设为 50：同一染色体但距离 > 50 的边被视为可疑，与不同染色体边同级。
        endpoint_penalty = {}
        for u, v, d in aag.graph.edges(data=True):
            if d.get('conflict') and d.get('endpoint'):
                key = (u, v) if u.hog_id < v.hog_id else (v, u)
                endpoint_penalty[key] = True
        same_chrom_other = getattr(aag, 'same_chrom_other', {})
        dist_other = getattr(aag, 'dist_other', {})

        def _conflict_edge_priority(key):
            same = same_chrom_other.get(key, False)
            dist = dist_other.get(key, 999999)
            endpoint = endpoint_penalty.get(key, False)
            # Tier 0: same chrom, close (< 50), non-endpoint
            # Tier 1: same chrom, close (< 50), endpoint
            # Tier 2: same chrom, far (>= 50), non-endpoint
            # Tier 3: same chrom, far (>= 50), endpoint  OR  diff chrom, non-endpoint
            # Tier 4: diff chrom, endpoint
            if same and dist < 50:
                tier = 0 if not endpoint else 1
            elif same:
                tier = 2 if not endpoint else 3
            else:
                tier = 3 if not endpoint else 4
            return tier

        conflict_edges.sort(key=lambda x: (
            -x[1],                          # 权重高优先
            _conflict_edge_priority(x[0]),  # 质量 tier
            dist_other.get(x[0], 999999),   # 距离近优先（仅对同 tier 有效）
            x[0][0].hog_id, x[0][1].hog_id  # 确定性兜底
        ))

        # 将 conflict 边分为高质量和低质量两组，先尽量用高质量边，若仍无法达到目标再用低质量边
        high_quality = [e for e in conflict_edges if _conflict_edge_priority(e[0]) <= 1]
        low_quality = [e for e in conflict_edges if _conflict_edge_priority(e[0]) > 1]

        # Union-Find
        parent = {n: n for n in undirected.nodes()}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        def count_components():
            return len({find(n) for n in undirected.nodes()})

        degree = {n: 0 for n in undirected.nodes()}
        selected = set()

        # Phase 1: 加入所有 support=2 边
        for (u, v), weight in support2_edges:
            if degree[u] < 2 and degree[v] < 2:
                selected.add((u, v))
                degree[u] += 1
                degree[v] += 1
                union(u, v)

        initial_components = count_components()

        # Phase 2: Kruskal 加入 conflict 边，直到达到目标染色体数
        # 2a: 先使用高质量冲突边（同染色体、距离近）
        if target_chromosomes is None or initial_components > target_chromosomes:
            for (u, v), weight in high_quality:
                if degree[u] >= 2 or degree[v] >= 2:
                    continue
                if find(u) == find(v):
                    continue
                selected.add((u, v))
                degree[u] += 1
                degree[v] += 1
                union(u, v)
                if target_chromosomes is not None:
                    current = count_components()
                    if current <= target_chromosomes:
                        break

        # 2b: 若仍未达到目标，回退到低质量冲突边
        # 当无目标染色体数时，也加入低质量边以尽可能补全图
        need_low_quality = target_chromosomes is None
        if target_chromosomes is not None and count_components() > target_chromosomes:
            need_low_quality = True
        if need_low_quality:
            for (u, v), weight in low_quality:
                if degree[u] >= 2 or degree[v] >= 2:
                    continue
                if find(u) == find(v):
                    continue
                selected.add((u, v))
                degree[u] += 1
                degree[v] += 1
                union(u, v)
                if target_chromosomes is not None:
                    current = count_components()
                    if current <= target_chromosomes:
                        break

        final_components = count_components()
        logger.info("  Linearization stats: initial={}, target={}, final={}".format(
            initial_components, target_chromosomes, final_components))
        if target_chromosomes is not None and final_components > target_chromosomes:
            logger.info("  Linearization: target {} chromosomes, reached {} (insufficient edges)".format(
                target_chromosomes, final_components))

        # Phase 3: 强制连接孤立节点（度为0）到最佳邻居，避免单基因染色体 artifact
        isolated = [n for n in undirected.nodes() if degree.get(n, 0) == 0]
        for h in isolated:
            best_n = None
            best_w = -1
            for n in undirected.neighbors(h):
                if degree.get(n, 0) >= 2:
                    continue
                key = (h, n) if h.hog_id < n.hog_id else (n, h)
                w = edge_weights.get(key, 0)
                if w > best_w:
                    best_w = w
                    best_n = n
            if best_n is not None:
                selected.add((h, best_n))
                degree[h] = 1
                degree[best_n] += 1

        # 构建无向子图
        subg = nx.Graph()
        for u, v in selected:
            subg.add_edge(u, v)

        # 对每个连通分量进行一致定向
        seen_nodes = set()
        for comp_nodes in nx.connected_components(subg):
            if len(comp_nodes) == 1:
                n = list(comp_nodes)[0]
                new_graph.add_node(n)
                seen_nodes.add(n)
                continue

            comp = subg.subgraph(comp_nodes).copy()
            # 检测并打破环（理论上 Kruskal 不会产生环，但 support=2 边可能已有环）
            if comp.number_of_edges() == len(comp_nodes):
                # 找到权重最低的边删除
                weakest = None
                min_w = float('inf')
                for u, v in comp.edges():
                    key = (u, v) if u.hog_id < v.hog_id else (v, u)
                    w = edge_weights.get(key, 0)
                    if w < min_w:
                        min_w = w
                        weakest = (u, v)
                if weakest:
                    comp.remove_edge(*weakest)

            # 找到路径端点（度为1的节点）
            endpoints = [n for n in comp_nodes if comp.degree(n) == 1]
            if endpoints:
                start = endpoints[0]
            else:
                start = list(comp_nodes)[0]

            # 沿路径行走并定向
            curr = start
            prev = None
            while True:
                nbrs = [n for n in comp.neighbors(curr) if n != prev]
                if not nbrs:
                    break
                nxt = nbrs[0]
                new_graph.add_edge(curr, nxt)
                prev = curr
                curr = nxt

            seen_nodes.update(comp_nodes)

        # 补充孤立节点
        for n in undirected.nodes():
            if n not in seen_nodes:
                new_graph.add_node(n)

        # 复制原图的节点属性（如 strand）到新图
        for n in new_graph.nodes():
            if n in aag.graph:
                for k, v in aag.graph.nodes[n].items():
                    if k not in new_graph.nodes[n]:
                        new_graph.nodes[n][k] = v

        result = AncestralAdjacencyGraph(node_id=aag.node_id)
        result.species_set = aag.species_set
        result.hog_map = dict(aag.hog_map)
        result.graph = new_graph
        result.gene_nodes = set(new_graph.nodes())
        result.events = list(aag.events)
        for n in new_graph.nodes():
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom

        # 端粒约束后处理：丢弃碎片路径
        result = self._filter_fragments(result, aag, target_chromosomes)
        return result

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

    def _resolve_merge_conflicts(self, merged, mapped_children, mapped_outgroups):
        """
        解决合并冲突。
        核心思想：外类群投票用于为 conflict 边加权，而不是直接删除边。
        所有边保留，由 _linearize_graph 基于全局权重选择最优边集。
        这样可以避免 short-cut 检测过度删除导致的碎片化。
        """
        # 收集外类群映射图中的所有边（无向），按系统发育距离加权
        outgroup_edge_weights = defaultdict(float)
        for og, weight in mapped_outgroups:
            for n1, n2 in og.get_adjacencies(include_telomere=False):
                key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                outgroup_edge_weights[key] += weight
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        weight_threshold = total_weight / 3.0 if total_weight > 0 else 0

        # 为每条 conflict 边标记外类群支持权重，不删除任何边
        n_conflict = 0
        n_outgroup_supported = 0
        for n1, n2 in list(merged.graph.edges()):
            if n1 in merged.telomeres or n2 in merged.telomeres:
                continue
            data = merged.graph[n1][n2]
            if data.get('conflict'):
                n_conflict += 1
                key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                edge_weight = outgroup_edge_weights.get(key, 0)
                merged.graph[n1][n2]['og_weight'] = edge_weight
                merged.graph[n1][n2]['og_supported'] = edge_weight >= weight_threshold
                if edge_weight >= weight_threshold:
                    n_outgroup_supported += 1

        if n_conflict > 0:
            logger.info("  Conflict stats: {} conflict edges, {} outgroup-supported".format(
                n_conflict, n_outgroup_supported))
        return merged

    def _prune_branches(self, merged, mapped_outgroups):
        """
        分支修剪：将度>2的节点修剪到度<=2。
        优先保留 support=2 的边，其次是外类群支持的 support=1 边，
        最后按路径长度保留最长的两条分支（保留主干，丢弃短侧枝）。
        这能显著减少 _linearize_graph 产生的染色体碎片。
        """
        # 收集外类群边，按系统发育距离加权
        outgroup_edge_weights = defaultdict(float)
        for og, weight in mapped_outgroups:
            for n1, n2 in og.get_adjacencies(include_telomere=False):
                key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                outgroup_edge_weights[key] += weight
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        weight_threshold = total_weight / 3.0 if total_weight > 0 else 0

        def edge_quality(n1, n2):
            """边质量：2=双支持, 1=单支持+外类群加权支持, 0=单支持无外类群"""
            d1 = merged.graph.get_edge_data(n1, n2, default={})
            d2 = merged.graph.get_edge_data(n2, n1, default={})
            support = max(d1.get('support', 1), d2.get('support', 1))
            if support >= 2:
                return 2
            key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
            return 1 if outgroup_edge_weights.get(key, 0) >= weight_threshold else 0

        def path_length_from(start, avoid):
            """从 start 出发，不经过 avoid，沿度<=2的链走多远"""
            visited = {avoid}
            curr = start
            length = 1
            visited.add(curr)
            while True:
                nbrs = set()
                for n in merged.graph.successors(curr):
                    if n not in merged.telomeres and n not in visited:
                        nbrs.add(n)
                for n in merged.graph.predecessors(curr):
                    if n not in merged.telomeres and n not in visited:
                        nbrs.add(n)
                if len(nbrs) != 1:
                    break
                curr = nbrs.pop()
                length += 1
                visited.add(curr)
            return length

        edges_removed = 0
        while True:
            changed = False
            for node in list(merged.gene_nodes):
                if node in merged.telomeres:
                    continue
                # 收集所有非端粒邻居
                neighbors = set()
                for succ in merged.graph.successors(node):
                    if succ not in merged.telomeres:
                        neighbors.add(succ)
                for pred in merged.graph.predecessors(node):
                    if pred not in merged.telomeres:
                        neighbors.add(pred)

                if len(neighbors) <= 2:
                    continue

                # 按（质量, 路径长度）排序，保留最优的两条边
                scored = []
                for nb in neighbors:
                    qual = edge_quality(node, nb)
                    plen = path_length_from(nb, node)
                    scored.append((qual, plen, nb))
                scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
                keep = {nb for _, _, nb in scored[:2]}

                for nb in neighbors:
                    if nb in keep:
                        continue
                    if merged.graph.has_edge(node, nb):
                        merged.graph.remove_edge(node, nb)
                        edges_removed += 1
                    if merged.graph.has_edge(nb, node):
                        merged.graph.remove_edge(nb, node)
                        edges_removed += 1
                    changed = True

            if not changed:
                break

        if edges_removed > 0:
            logger.info("  Pruned branches: removed {} edges, all nodes now degree <= 2".format(edges_removed))
        return merged

    def _resolve_indels_v3(self, paths, merged, mapped_children, mapped_outgroups):
        """
        Resolve indels for v3 without destroying CP-SAT path structure.

        "挖洞"策略：从路径中移除 indel HOG，但不断开路径。
        相邻 HOG 直接相连，保持染色体数与 CP-SAT 输出一致。
        只有整条路径全为 indel 时才丢弃该路径。

        Returns:
            (updated_paths, indel_events)
        """
        node_id = merged.node_id

        # Same presence/absence logic as v2
        hog_to_children = defaultdict(set)
        for mc in mapped_children:
            for rec in mc.gene_nodes:
                hog_to_children[rec].add(mc.node_id)

        outgroup_hog_weights = defaultdict(float)
        for og, weight in mapped_outgroups:
            for rec in og.gene_nodes:
                outgroup_hog_weights[rec] += weight
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        weight_threshold = total_weight / 3.0 if total_weight > 0 else 0
        outgroup_hogs = {hog for hog, w in outgroup_hog_weights.items()
                         if w >= weight_threshold}

        children_with_data = [mc for mc in mapped_children if len(mc.gene_nodes) > 0]
        n_children_data = len(children_with_data)

        genes_to_remove = set()
        for hog_rec in list(merged.gene_nodes):
            child_present = len(hog_to_children.get(hog_rec, set()))
            outgroup_present = hog_rec in outgroup_hogs
            if n_children_data > 0 and child_present < n_children_data and not outgroup_present:
                genes_to_remove.add(hog_rec)

        if not genes_to_remove:
            return paths, [], 0

        # 挖洞：从路径中移除 indel HOG，但保持路径连续
        new_paths = []
        events = []
        n_removed_total = 0
        n_paths_removed = 0
        for path in paths:
            remaining = [h for h in path if h not in genes_to_remove]
            removed_in_path = [h for h in path if h in genes_to_remove]
            n_removed_total += len(removed_in_path)

            if not remaining:
                # 整条路径都是 indel → 丢弃，记录事件
                n_paths_removed += 1
                events.append(RearrangementEvent(
                    'indel', node_id,
                    genes_involved=path,
                    desc="indel: entire path ({} HOGs) removed".format(len(path)),
                    support=len(path)
                ))
            else:
                # 保留路径，挖洞连接
                new_paths.append(remaining)
                if removed_in_path:
                    events.append(RearrangementEvent(
                        'indel', node_id,
                        genes_involved=removed_in_path,
                        desc="indel: {} HOGs removed from path (punched through)".format(
                            len(removed_in_path)),
                        support=len(removed_in_path)
                    ))

        logger.info("  [v3] After indel resolution: removed {} HOGs, {} paths -> {} paths (punch-through, {} paths dropped), {} indel events".format(
            n_removed_total, len(paths), len(new_paths), n_paths_removed, len(events)))

        return new_paths, events, n_paths_removed
        
    def _resolve_indels(self, merged, mapped_children, mapped_outgroups):
        """
        Resolve indels (gene loss/insertion) for v2 path.
        If a HOG appears in only some children and is absent from outgroups,
        it is inferred as an insertion in those children and removed from ancestor.
        Note: children without data (e.g., outgroup leaves without genes) should not be counted.
        """
        node_id = merged.node_id
        # Count which children each HOG appears in
        hog_to_children = defaultdict(set)
        for mc in mapped_children:
            for rec in mc.gene_nodes:
                hog_to_children[rec].add(mc.node_id)

        # Collect outgroup HOGs (for inferring ancestral presence), weighted by phylogenetic distance
        # HOG weight >= 1/3 of total weight is considered outgroup-supported
        outgroup_hog_weights = defaultdict(float)
        for og, weight in mapped_outgroups:
            for rec in og.gene_nodes:
                outgroup_hog_weights[rec] += weight
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        weight_threshold = total_weight / 3.0 if total_weight > 0 else 0
        outgroup_hogs = {hog for hog, w in outgroup_hog_weights.items()
                         if w >= weight_threshold}

        # Only count children with data
        children_with_data = [mc for mc in mapped_children if len(mc.gene_nodes) > 0]
        n_children_data = len(children_with_data)

        genes_to_remove = set()
        for hog_rec in list(merged.gene_nodes):
            child_present = len(hog_to_children.get(hog_rec, set()))
            outgroup_present = hog_rec in outgroup_hogs
            # If HOG is not present in all children with data and absent from outgroups, mark as indel
            if n_children_data > 0 and child_present < n_children_data and not outgroup_present:
                genes_to_remove.add(hog_rec)
                merged.events.append(RearrangementEvent(
                    'indel', node_id,
                    genes_involved=[hog_rec],
                    desc="HOG {} present in only {}/{} children, absent in outgroup".format(
                        hog_rec.hog_id, child_present, n_children_data),
                    support="{}/{} children".format(child_present, n_children_data)
                ))

        if genes_to_remove:
            merged.remove_nodes(genes_to_remove)
            merged = self._linearize_graph(merged)
            merged._add_telomeres()
            logger.info("  After indel resolution: removed {} HOGs".format(len(genes_to_remove)))
        return merged

    def _resolve_duplications(self, merged, mapped_children):
        """
        Resolve duplications. If a child maps multiple original genes to one HOG,
        record duplication events. Ancestor graph keeps one representative per HOG.
        """
        node_id = merged.node_id
        for mc in mapped_children:
            # 统计该子节点中每个 HOG 对应多少个原始基因
            # 需要反向建立 HOG -> [原始基因] 的映射
            hog_to_genes = defaultdict(list)
            for gene, rec in mc.hog_map.items():
                if isinstance(gene, HOGrecord):
                    continue  # 祖先图节点不统计
                hog_to_genes[rec].append(gene)

            for hog_rec, copies in hog_to_genes.items():
                if len(copies) <= 1:
                    continue
                copies_sorted = sorted(copies, key=lambda x: getattr(x, 'index', 0))
                min_dist = None
                for c1, c2 in itertools.combinations(copies_sorted, 2):
                    if hasattr(c1, 'index') and hasattr(c2, 'index') \
                            and getattr(c1, 'chrom', None) == getattr(c2, 'chrom', None):
                        d = abs(c1.index - c2.index)
                        if min_dist is None or d < min_dist:
                            min_dist = d
                if min_dist == 1:
                    etype = 'tandem_dup'
                elif min_dist is not None and min_dist <= 5:
                    etype = 'proximal_dup'
                else:
                    etype = 'dispersed_dup'
                gene_names = ','.join(str(getattr(c, 'id', str(c))) for c in copies[:3])
                if len(copies) > 3:
                    gene_names += '...'
                merged.events.append(RearrangementEvent(
                    etype, node_id,
                    genes_involved=copies,
                    desc="{}: HOG {} x{} in {} (genes: {})".format(
                        etype, hog_rec.hog_id, len(copies), mc.node_id, gene_names),
                    support=min_dist
                ))
        return merged

    def _detect_bottomup_events(self, merged, mapped_children, child_source_ids, mapped_outgroups=None):
        """
        Detect bottom-up rearrangement events: compare each child mapped graph (input)
        vs the linearized ancestor graph (output).
        Each conflict type corresponds to a rearrangement type, detected in order:

        1. indel: HOG missing in some children
        2. duplication: HOG has multiple copies in a child
        3. inversion: strand mismatch (child vs ancestor consensus)
        4. unidir_trans: segment moved to another chromosome (non-telomere exchange)
        5. reciprocal_translocation: two adjacencies break and cross-connect (A-B,C-D → A-D,C-B)
        6. fission: ancestor adjacency breaks in child, ends become different chromosome endpoints
        7. NCF: chromosome fusion (one-end or non-endpoint adjacency change)
        8. EEJ: new child adjacency connects two ancestor chromosome endpoints

        Outgroup adjacencies are used for:
        - Validating anc_only edges: if outgroup also has the adjacency, the child truly lost it
        - Validating child_only edges: if outgroup also has the adjacency, ancestor may have missed it
        """
        # Build outgroup adjacency set (weighted)
        og_adj_weights = defaultdict(float)
        if mapped_outgroups:
            for og, weight in mapped_outgroups:
                for n1, n2 in og.get_adjacencies(include_telomere=False):
                    key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                    og_adj_weights[key] += weight
            total_og_weight = sum(weight for _, weight in mapped_outgroups)
            og_threshold = total_og_weight / 3.0 if total_og_weight > 0 else 0
            og_supported = {k for k, v in og_adj_weights.items() if v >= og_threshold}
        else:
            og_supported = set()
        node_id = merged.node_id

        # 构建祖先图位置索引
        merged_hog_to_chrom = {}
        merged_hog_to_pos = {}
        merged_chrom_ends = {}  # hog -> 'left'/'right'/'both'/'internal'
        merged_adjs_norm = set()
        for ci, chrom in enumerate(merged.chromosomes):
            genes = [n for n in chrom if n not in merged.telomeres]
            for j, h in enumerate(genes):
                merged_hog_to_chrom[h] = ci
                merged_hog_to_pos[h] = j
                if j == 0 and j == len(genes) - 1:
                    merged_chrom_ends[h] = 'both'
                elif j == 0:
                    merged_chrom_ends[h] = 'left'
                elif j == len(genes) - 1:
                    merged_chrom_ends[h] = 'right'
                else:
                    merged_chrom_ends[h] = 'internal'

        # Build support count map for merged adjacencies
        merged_adj_support = {}
        for n1, n2, d in merged.graph.edges(data=True):
            if n1 in merged.telomeres or n2 in merged.telomeres:
                continue
            if n1 == n2:
                continue
            key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
            merged_adj_support[key] = d.get('support', 1)

        for n1, n2 in merged.get_adjacencies(include_telomere=False):
            key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
            merged_adjs_norm.add(key)

        def get_merged_strand(h):
            s = merged.graph.nodes[h].get('strand')
            return s if s else getattr(h, 'strand', '+')

        for ci_idx, mc in enumerate(mapped_children):
            src = child_source_ids[ci_idx] if ci_idx < len(child_source_ids) else mc.node_id

            # 子节点位置索引：优先使用 _map_to_parent_hogs 中从原始图计算的属性
            # 回退到 mc.chromosomes 遍历（兼容旧代码路径）
            if hasattr(mc, 'hog_to_chrom_idx') and mc.hog_to_chrom_idx:
                mc_hog_to_chrom = mc.hog_to_chrom_idx
                mc_chrom_ends = mc.hog_chrom_ends
            else:
                mc_hog_to_chrom = {}
                mc_hog_to_pos = {}
                mc_chrom_ends = {}
                for cci, chrom in enumerate(mc.chromosomes):
                    genes = [n for n in chrom if n not in mc.telomeres]
                    for j, h in enumerate(genes):
                        mc_hog_to_chrom[h] = cci
                        mc_hog_to_pos[h] = j
                        if j == 0 and j == len(genes) - 1:
                            mc_chrom_ends[h] = 'both'
                        elif j == 0:
                            mc_chrom_ends[h] = 'left'
                        elif j == len(genes) - 1:
                            mc_chrom_ends[h] = 'right'
                        else:
                            mc_chrom_ends[h] = 'internal'
            mc_adjs_norm = set()
            for n1, n2 in mc.get_adjacencies(include_telomere=False):
                key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
                mc_adjs_norm.add(key)

            def get_mc_strand(h):
                s = mc.graph.nodes[h].get('strand')
                return s if s else getattr(h, 'strand', '+')

            # 邻接差异集合
            # 祖先有但子节点无 → 子节点分支上断裂
            anc_only = merged_adjs_norm - mc_adjs_norm
            # 子节点有但祖先无 → 子节点分支上创造新连接
            child_only = mc_adjs_norm - merged_adjs_norm

            # 过滤重建伪影：anc_only中support=1的边是重建伪影
            # 当祖先只采用了一个子节点的邻接（support=1），另一个子节点会
            # 显示该邻接为"丢失"，但这并非真实的演化事件
            # 只有support>=2的anc_only边才是真实的断裂事件
            anc_only = {key for key in anc_only
                        if merged_adj_support.get(key, 1) >= 2}

            # 过滤重建伪影：child_only中support=1的边也是重建伪影
            # 当祖先拒绝了一个子节点独有的邻接（support=1），这个邻接在
            # 该子节点中显示为"新产生"，但祖先可能是有意拒绝的（识别为
            # 子节点特有重排），不应计为演化事件
            child_only = {key for key in child_only
                          if merged_adj_support.get(key, 1) >= 2}

            # ============================================================
            # 1. Inversion 检测（小重排）
            # 子节点链向与祖先共识不一致的连续HOG块
            # 分类: internal_inversion / telomere_inversion
            # ============================================================
            # 使用 chrom_hogs（从原始图计算的映射后HOG列表），回退到 mc.chromosomes
            if hasattr(mc, 'chrom_hogs') and mc.chrom_hogs:
                chrom_iter = [(cci, hogs) for cci, hogs in sorted(mc.chrom_hogs.items())]
            else:
                chrom_iter = [(cci, [n for n in chrom if n not in mc.telomeres])
                              for cci, chrom in enumerate(mc.chromosomes)]
            for cci, genes in chrom_iter:
                genes = [h for h in genes if h in merged_hog_to_chrom]
                if len(genes) < 3:
                    continue
                blocks = []
                start = None
                for i, h in enumerate(genes):
                    mc_s = get_mc_strand(h)
                    anc_s = get_merged_strand(h)
                    if mc_s != anc_s:
                        if start is None:
                            start = i
                    else:
                        if start is not None:
                            blocks.append((start, i))
                            start = None
                if start is not None:
                    blocks.append((start, len(genes)))

                min_block_size = max(3, int(len(genes) * 0.02) + 1)
                for bstart, bend in blocks:
                    block = genes[bstart:bend]
                    if len(block) < min_block_size:
                        continue
                    is_left = (bstart == 0)
                    is_right = (bend == len(genes))

                    if is_left or is_right:
                        etype = 'telomere_inversion'
                    else:
                        etype = 'internal_inversion'
                    merged.events.append(RearrangementEvent(
                        etype, node_id,
                        genes_involved=block,
                        desc="{}: {} HOGs [{}..{}] in {}".format(
                            etype, len(block),
                            block[0].hog_id, block[-1].hog_id, src),
                        support=len(block),
                        child_source=src
                    ))

                    # 移除倒位边界邻接，避免被误判为 NCF/unidir_trans
                    # 祖先顺序: left_nb - block[-1] - ... - block[0] - right_nb
                    # 子节点顺序: left_nb - block[0] - ... - block[-1] - right_nb
                    if bstart > 0:
                        left_nb = genes[bstart - 1]
                        # 祖先邻接: left_nb - block[-1] (在 anc_only 中)
                        inv_anc = (left_nb, block[-1]) if left_nb.hog_id < block[-1].hog_id \
                            else (block[-1], left_nb)
                        anc_only.discard(inv_anc)
                        # 子节点邻接: left_nb - block[0] (在 child_only 中)
                        inv_child = (left_nb, block[0]) if left_nb.hog_id < block[0].hog_id \
                            else (block[0], left_nb)
                        child_only.discard(inv_child)
                    if bend < len(genes):
                        right_nb = genes[bend]
                        # 祖先邻接: block[0] - right_nb (在 anc_only 中)
                        inv_anc = (block[0], right_nb) if block[0].hog_id < right_nb.hog_id \
                            else (right_nb, block[0])
                        anc_only.discard(inv_anc)
                        # 子节点邻接: block[-1] - right_nb (在 child_only 中)
                        inv_child = (block[-1], right_nb) if block[-1].hog_id < right_nb.hog_id \
                            else (right_nb, block[-1])
                        child_only.discard(inv_child)

            # ============================================================
            # 2. Reciprocal Translocation 检测（先于 unidir_trans）
            # 特征：祖先邻接 A-B 和 C-D 在子节点中同时断裂，
            #        子节点产生 A-D 和 C-B（交叉连接）
            # 子类型区分:
            #   - RT: 两个断点都在染色体内部（平衡交换）
            #   - URT: 一个断点在端粒，另一个在内部（不平衡交换）
            # 必须在 unidir_trans 之前检测，否则交叉连接边会被
            # unidir_trans 移除，导致 RT 无法识别。
            # ============================================================
            rt_matched_anc = set()
            rt_matched_child = set()
            for h1, h2 in list(anc_only):
                if (h1, h2) in rt_matched_anc:
                    continue
                a_partners_in_child = {e for e in child_only
                                       if (e[0] == h1 or e[1] == h1)
                                       and e not in rt_matched_child}
                for ae in a_partners_in_child:
                    x = ae[1] if ae[0] == h1 else ae[0]
                    for h3, h4 in list(anc_only):
                        if (h3, h4) in rt_matched_anc:
                            continue
                        cross_in_child = None
                        other_break = None
                        if x == h3:
                            cross_in_child = (h4, h2) if h4.hog_id < h2.hog_id else (h2, h4)
                            other_break = (h3, h4)
                        elif x == h4:
                            cross_in_child = (h3, h2) if h3.hog_id < h2.hog_id else (h2, h3)
                            other_break = (h3, h4)
                        if cross_in_child and cross_in_child in child_only \
                                and cross_in_child not in rt_matched_child:
                            c_ab1 = merged_hog_to_chrom.get(h1)
                            c_ab2 = merged_hog_to_chrom.get(h2)
                            c_cd1 = merged_hog_to_chrom.get(h3)
                            c_cd2 = merged_hog_to_chrom.get(h4)
                            if c_ab1 is not None and c_ab2 is not None \
                                    and c_cd1 is not None and c_cd2 is not None \
                                    and c_ab1 != c_ab2 and c_cd1 != c_cd2:
                                # 判断 RT vs URT: 检查是否有断点涉及端粒
                                # 在祖先图中，端点HOG标记为 left/right/both
                                h1_end = merged_chrom_ends.get(h1, 'internal')
                                h2_end = merged_chrom_ends.get(h2, 'internal')
                                h3_end = merged_chrom_ends.get(h3, 'internal')
                                h4_end = merged_chrom_ends.get(h4, 'internal')
                                # 如果任一断点涉及端粒（非internal），则为URT
                                is_unbalanced = any(
                                    e != 'internal' for e in [h1_end, h2_end, h3_end, h4_end]
                                )
                                if is_unbalanced:
                                    etype = 'unbalanced_reciprocal_translocation'
                                    desc = "URT: {}-{} and {}-{} crossed in {} (unbalanced)".format(
                                        h1.hog_id, h2.hog_id, h3.hog_id, h4.hog_id, src)
                                else:
                                    etype = 'reciprocal_translocation'
                                    desc = "RT: {}-{} and {}-{} crossed in {}".format(
                                        h1.hog_id, h2.hog_id, h3.hog_id, h4.hog_id, src)
                                merged.events.append(RearrangementEvent(
                                    etype, node_id,
                                    genes_involved=[h1, h2, h3, h4],
                                    desc=desc,
                                    support=2,
                                    child_source=src
                                ))
                                rt_matched_anc.add((h1, h2))
                                rt_matched_anc.add(other_break)
                                rt_matched_child.add(ae)
                                rt_matched_child.add(cross_in_child)
                                break
                    else:
                        continue
                    break

            # 去掉已匹配为 RT 的边
            anc_only_rem = anc_only - rt_matched_anc
            child_only_rem = child_only - rt_matched_child

            # ============================================================
            # 3. unidir_trans 检测（中等重排）
            # 片段从一条染色体移到另一条，两端不涉及染色体端点
            # 在 RT 检测之后执行，剩余边中不会再有交叉连接
            # ============================================================
            unidir_anc = set()   # 已匹配为 unidir_trans 的 anc_only 边
            unidir_child = set()  # 已匹配为 unidir_trans 的 child_only 边

            # anc_only_rem 中：两个HOG在子节点分属不同染色体，且都是内部节点
            for h1, h2 in list(anc_only_rem):
                e1 = mc_chrom_ends.get(h1, 'internal')
                e2 = mc_chrom_ends.get(h2, 'internal')
                c1 = mc_hog_to_chrom.get(h1)
                c2 = mc_hog_to_chrom.get(h2)
                if c1 is not None and c2 is not None and c1 != c2 \
                        and e1 == 'internal' and e2 == 'internal':
                    merged.events.append(RearrangementEvent(
                        'unidir_trans', node_id,
                        genes_involved=[h1, h2],
                        desc="unidir_trans: {}-{} adjacent in ancestor, "
                             "on diff chroms in {} (internal)".format(
                            h1.hog_id, h2.hog_id, src),
                        support=1,
                        child_source=src
                    ))
                    unidir_anc.add((h1, h2))

            # child_only_rem 中：都不是祖先染色体端点
            for h1, h2 in list(child_only_rem):
                e1 = merged_chrom_ends.get(h1, 'internal')
                e2 = merged_chrom_ends.get(h2, 'internal')
                if e1 == 'internal' and e2 == 'internal':
                    merged.events.append(RearrangementEvent(
                        'unidir_trans', node_id,
                        genes_involved=[h1, h2],
                        desc="unidir_trans: {}-{} new adjacency in {} "
                             "(internal in ancestor)".format(
                            h1.hog_id, h2.hog_id, src),
                        support=1,
                        child_source=src
                    ))
                    unidir_child.add((h1, h2))

            # 从剩余边中去掉已匹配为 unidir_trans 的边
            anc_only_final = anc_only_rem - unidir_anc
            child_only_final = child_only_rem - unidir_child

            # ============================================================
            # 4-5. Fission / NCF / EEJ 检测 + 合并
            # 先收集所有断裂和连接事件，再匹配配对（同一重排的断裂+连接），
            # 合并为单个事件，避免重复计数。
            # ============================================================

            # 收集断裂事件（祖先邻接在子节点中消失）
            break_events = []  # [(h1, h2, event_type)]
            for h1, h2 in anc_only_final:
                e1 = mc_chrom_ends.get(h1, 'internal')
                e2 = mc_chrom_ends.get(h2, 'internal')
                c1 = mc_hog_to_chrom.get(h1)
                c2 = mc_hog_to_chrom.get(h2)

                if c1 is not None and c2 is not None and c1 != c2:
                    # 染色体间断裂：两个HOG在子节点分属不同染色体
                    if e1 != 'internal' and e2 != 'internal':
                        break_events.append((h1, h2, 'fission'))
                    elif e1 != 'internal' or e2 != 'internal':
                        # 一端是子节点端点 → 可能是 translocation 的一部分
                        break_events.append((h1, h2, 'endpoint_break'))
                elif c1 is not None and c2 is not None and c1 == c2:
                    # 染色体内断裂：两个HOG在子节点同一染色体上但不再相邻
                    # 这是 NCF 插入位点的标志（recipient 内部邻接被 donor 打断）
                    if e1 == 'internal' and e2 == 'internal':
                        break_events.append((h1, h2, 'ncf_insertion'))

            # 收集连接事件（子节点中有但祖先中没有的邻接）
            join_events = []  # [(h1, h2, event_type)]
            for h1, h2 in child_only_final:
                e1 = merged_chrom_ends.get(h1, 'internal')
                e2 = merged_chrom_ends.get(h2, 'internal')

                if e1 != 'internal' and e2 != 'internal':
                    join_events.append((h1, h2, 'eej'))
                elif e1 != 'internal' or e2 != 'internal':
                    join_events.append((h1, h2, 'ncf'))

            # 配对 NCF 连接事件：
            # NCF（嵌套融合）在 child_only 中产生两个连接事件：
            #   (recipient_P, donor_endpoint_A) 和 (donor_endpoint_C, recipient_Q)
            # 两个 donor 端点来自同一祖先染色体，两个 recipient HOG 在祖先中相邻（插入位点）。
            # 这两个事件应合并为一个 NCF（净 -1 染色体），而非两个（净 -2）。
            ncf_paired_joins = set()  # 已配对的 NCF join_event 索引
            for ji1 in range(len(join_events)):
                if ji1 in ncf_paired_joins:
                    continue
                jh1_1, jh2_1, jt1 = join_events[ji1]
                if jt1 != 'ncf':
                    continue
                # 识别 donor 端点 HOG（祖先染色体端点）和 recipient 内部 HOG
                e1_1 = merged_chrom_ends.get(jh1_1, 'internal')
                if e1_1 != 'internal':
                    donor1, recip1 = jh1_1, jh2_1
                else:
                    donor1, recip1 = jh2_1, jh1_1
                for ji2 in range(ji1 + 1, len(join_events)):
                    if ji2 in ncf_paired_joins:
                        continue
                    jh1_2, jh2_2, jt2 = join_events[ji2]
                    if jt2 != 'ncf':
                        continue
                    e1_2 = merged_chrom_ends.get(jh1_2, 'internal')
                    if e1_2 != 'internal':
                        donor2, recip2 = jh1_2, jh2_2
                    else:
                        donor2, recip2 = jh2_2, jh1_2
                    # 条件1：两个 donor 端点来自同一祖先染色体
                    dc1 = merged_hog_to_chrom.get(donor1)
                    dc2 = merged_hog_to_chrom.get(donor2)
                    if dc1 is None or dc2 is None or dc1 != dc2:
                        continue
                    # 条件2：两个 recipient HOG 在祖先中相邻（插入位点断裂）
                    rk = (recip1, recip2) if recip1.hog_id < recip2.hog_id \
                        else (recip2, recip1)
                    if rk in merged_adjs_norm:
                        # 确认这是一次 NCF 事件
                        all_hogs = list(set([donor1, donor2, recip1, recip2]))
                        merged.events.append(RearrangementEvent(
                            'ncf', node_id,
                            genes_involved=all_hogs,
                            desc="NCF: donor({},{}) from chrom {} inserted "
                                 "between {}-{} in {}".format(
                                donor1.hog_id, donor2.hog_id, dc1,
                                recip1.hog_id, recip2.hog_id, src),
                            support=2,
                            child_source=src
                        ))
                        ncf_paired_joins.add(ji1)
                        ncf_paired_joins.add(ji2)
                        break

            # 合并断裂-连接配对：
            # 如果一个断裂的端点HOG出现在一个连接中，说明这是同一次重排的两个方面，
            # 应合并为一个事件（净染色体变化为0，而非-1-1=-2）。
            matched_breaks = set()
            matched_joins = set()
            for bi, (h1, h2, bt) in enumerate(break_events):
                if bi in matched_breaks:
                    continue
                # 找与这个断裂共享端点HOG的连接（跳过已配对的 NCF）
                for ji, (jh1, jh2, jt) in enumerate(join_events):
                    if ji in matched_joins or ji in ncf_paired_joins:
                        continue
                    # 端点匹配：断裂的端点HOG出现在连接中
                    if h1 == jh1 or h1 == jh2 or h2 == jh1 or h2 == jh2:
                        # 判断是 RT 还是 UT：检查是否有端点参与
                        h1_end = merged_chrom_ends.get(h1, 'internal')
                        h2_end = merged_chrom_ends.get(h2, 'internal')
                        jh1_end = merged_chrom_ends.get(jh1, 'internal')
                        jh2_end = merged_chrom_ends.get(jh2, 'internal')
                        is_unbalanced = any(
                            e != 'internal' for e in [h1_end, h2_end, jh1_end, jh2_end]
                        )
                        if is_unbalanced:
                            etype = 'unbalanced_reciprocal_translocation'
                            desc = "URT: {}-{} broken, {}-{} joined in {}".format(
                                h1.hog_id, h2.hog_id, jh1.hog_id, jh2.hog_id, src)
                        else:
                            etype = 'reciprocal_translocation'
                            desc = "RT: {}-{} broken, {}-{} joined in {}".format(
                                h1.hog_id, h2.hog_id, jh1.hog_id, jh2.hog_id, src)
                        all_hogs = list(set([h1, h2, jh1, jh2]))
                        merged.events.append(RearrangementEvent(
                            etype, node_id,
                            genes_involved=all_hogs,
                            desc=desc,
                            support=2,
                            child_source=src
                        ))
                        matched_breaks.add(bi)
                        matched_joins.add(ji)
                        break

            # 配对跨染色体 endpoint_break 事件：
            # 两个 endpoint_break 来自不同祖先染色体 → 相互易位的两端
            # 净染色体变化=0（不是两个 NCF = -2）
            ep_breaks = [(bi, h1, h2) for bi, (h1, h2, bt) in enumerate(break_events)
                         if bt == 'endpoint_break' and bi not in matched_breaks]
            for i in range(len(ep_breaks)):
                bi1, h1_1, h2_1 = ep_breaks[i]
                if bi1 in matched_breaks:
                    continue
                c1_1 = merged_hog_to_chrom.get(h1_1)
                c2_1 = merged_hog_to_chrom.get(h2_1)
                if c1_1 is None or c2_1 is None:
                    continue
                for j in range(i + 1, len(ep_breaks)):
                    bi2, h1_2, h2_2 = ep_breaks[j]
                    if bi2 in matched_breaks:
                        continue
                    c1_2 = merged_hog_to_chrom.get(h1_2)
                    c2_2 = merged_hog_to_chrom.get(h2_2)
                    if c1_2 is None or c2_2 is None:
                        continue
                    # 两个断裂来自不同的祖先染色体对
                    chroms_1 = {c1_1, c2_1}
                    chroms_2 = {c1_2, c2_2}
                    if chroms_1 == chroms_2 and c1_1 != c2_1:
                        # RT 模式：两个断裂都跨越相同的两条祖先染色体
                        # 判断 RT vs URT
                        h1_1_end = merged_chrom_ends.get(h1_1, 'internal')
                        h2_1_end = merged_chrom_ends.get(h2_1, 'internal')
                        h1_2_end = merged_chrom_ends.get(h1_2, 'internal')
                        h2_2_end = merged_chrom_ends.get(h2_2, 'internal')
                        is_unbalanced = any(
                            e != 'internal' for e in [h1_1_end, h2_1_end, h1_2_end, h2_2_end]
                        )
                        if is_unbalanced:
                            etype = 'unbalanced_reciprocal_translocation'
                            desc = "URT: {}-{} and {}-{} exchanged in {}".format(
                                h1_1.hog_id, h2_1.hog_id, h1_2.hog_id, h2_2.hog_id, src)
                        else:
                            etype = 'reciprocal_translocation'
                            desc = "RT: {}-{} and {}-{} exchanged in {}".format(
                                h1_1.hog_id, h2_1.hog_id, h1_2.hog_id, h2_2.hog_id, src)
                        all_hogs = list(set([h1_1, h2_1, h1_2, h2_2]))
                        merged.events.append(RearrangementEvent(
                            etype, node_id,
                            genes_involved=all_hogs,
                            desc=desc,
                            support=2,
                            child_source=src
                        ))
                        matched_breaks.add(bi1)
                        matched_breaks.add(bi2)
                        break

            # 配对跨染色体 ncf_insertion + endpoint_break 事件：
            # 当 ncf_insertion（染色体内断裂）和 endpoint_break（染色体间断裂）
            # 配对时，表示 RT 的一部分：染色体内断裂是 RT 的插入位点，
            # 染色体间断裂是 RT 的交换端点。净染色体变化=0。
            ncf_ins_breaks = [(bi, h1, h2) for bi, (h1, h2, bt) in enumerate(break_events)
                              if bt == 'ncf_insertion' and bi not in matched_breaks]
            ep_brks = [(bi, h1, h2) for bi, (h1, h2, bt) in enumerate(break_events)
                       if bt == 'endpoint_break' and bi not in matched_breaks]
            for ni, (bi_n, nh1, nh2) in enumerate(ncf_ins_breaks):
                if bi_n in matched_breaks:
                    continue
                nc_n = merged_hog_to_chrom.get(nh1)  # ncf_insertion: same ancestor chrom
                for ei, (bi_e, eh1, eh2) in enumerate(ep_brks):
                    if bi_e in matched_breaks:
                        continue
                    # endpoint_break spans two different ancestor chromosomes
                    ec1 = merged_hog_to_chrom.get(eh1)
                    ec2 = merged_hog_to_chrom.get(eh2)
                    if ec1 is None or ec2 is None:
                        continue
                    # 如果 endpoint_break 跨越的染色体之一包含 ncf_insertion，
                    # 说明这是 RT 的一部分
                    if nc_n is not None and (ec1 == nc_n or ec2 == nc_n):
                        # 判断 RT vs URT
                        nh1_end = merged_chrom_ends.get(nh1, 'internal')
                        nh2_end = merged_chrom_ends.get(nh2, 'internal')
                        eh1_end = merged_chrom_ends.get(eh1, 'internal')
                        eh2_end = merged_chrom_ends.get(eh2, 'internal')
                        is_unbalanced = any(
                            e != 'internal' for e in [nh1_end, nh2_end, eh1_end, eh2_end]
                        )
                        if is_unbalanced:
                            etype = 'unbalanced_reciprocal_translocation'
                            desc = "URT: {}-{} insertion site and {}-{} exchange in {}".format(
                                nh1.hog_id, nh2.hog_id, eh1.hog_id, eh2.hog_id, src)
                        else:
                            etype = 'reciprocal_translocation'
                            desc = "RT: {}-{} insertion site and {}-{} exchange in {}".format(
                                nh1.hog_id, nh2.hog_id, eh1.hog_id, eh2.hog_id, src)
                        all_hogs = list(set([nh1, nh2, eh1, eh2]))
                        merged.events.append(RearrangementEvent(
                            etype, node_id,
                            genes_involved=all_hogs,
                            desc=desc,
                            support=2,
                            child_source=src
                        ))
                        matched_breaks.add(bi_n)
                        matched_breaks.add(bi_e)
                        break

            # 配对两个 ncf_insertion 事件（来自不同祖先染色体）：
            # 两个 ncf_insertion 来自不同祖先染色体 → 交叉插入（RT 模式）
            # 净染色体变化=0
            for i in range(len(ncf_ins_breaks)):
                bi1, h1_1, h2_1 = ncf_ins_breaks[i]
                if bi1 in matched_breaks:
                    continue
                nc1 = merged_hog_to_chrom.get(h1_1)
                for j in range(i + 1, len(ncf_ins_breaks)):
                    bi2, h1_2, h2_2 = ncf_ins_breaks[j]
                    if bi2 in matched_breaks:
                        continue
                    nc2 = merged_hog_to_chrom.get(h1_2)
                    if nc1 is not None and nc2 is not None and nc1 != nc2:
                        # 两个 ncf_insertion 在不同祖先染色体上 → RT
                        # 判断 RT vs URT
                        h1_1_end = merged_chrom_ends.get(h1_1, 'internal')
                        h2_1_end = merged_chrom_ends.get(h2_1, 'internal')
                        h1_2_end = merged_chrom_ends.get(h1_2, 'internal')
                        h2_2_end = merged_chrom_ends.get(h2_2, 'internal')
                        is_unbalanced = any(
                            e != 'internal' for e in [h1_1_end, h2_1_end, h1_2_end, h2_2_end]
                        )
                        if is_unbalanced:
                            etype = 'unbalanced_reciprocal_translocation'
                            desc = "URT: {}-{} and {}-{} cross-insertion in {}".format(
                                h1_1.hog_id, h2_1.hog_id, h1_2.hog_id, h2_2.hog_id, src)
                        else:
                            etype = 'reciprocal_translocation'
                            desc = "RT: {}-{} and {}-{} cross-insertion in {}".format(
                                h1_1.hog_id, h2_1.hog_id, h1_2.hog_id, h2_2.hog_id, src)
                        all_hogs = list(set([h1_1, h2_1, h1_2, h2_2]))
                        merged.events.append(RearrangementEvent(
                            etype, node_id,
                            genes_involved=all_hogs,
                            desc=desc,
                            support=2,
                            child_source=src
                        ))
                        matched_breaks.add(bi1)
                        matched_breaks.add(bi2)
                        break

            # 发出未匹配的断裂事件
            for bi, (h1, h2, bt) in enumerate(break_events):
                if bi in matched_breaks:
                    continue
                if bt == 'fission':
                    merged.events.append(RearrangementEvent(
                        'fission', node_id,
                        genes_involved=[h1, h2],
                        desc="Fission: {}-{} adjacent in ancestor, "
                             "separated in {}".format(
                            h1.hog_id, h2.hog_id, src),
                        support=1,
                        child_source=src
                    ))
                elif bt == 'ncf_insertion':
                    # NCF 插入位点断裂：recipient 内部邻接被 donor 打断
                    # 只有当配对确认时才计为 NCF（影响染色体数），
                    # 未配对的不应影响染色体计数（可能是倒位残留或 RT 部分）
                    # 仍记录为事件用于注释，但不计入 NCF
                    merged.events.append(RearrangementEvent(
                        'unidir_trans', node_id,
                        genes_involved=[h1, h2],
                        desc="Insertion site: {}-{} adjacent in ancestor, "
                             "interrupted in {} (unpaired)".format(
                            h1.hog_id, h2.hog_id, src),
                        support=1,
                        child_source=src
                    ))
                else:
                    # endpoint_break：染色体间断裂，一端是子节点端点
                    # 未配对的可能是 EEJ 的一部分或孤立事件，按 EEJ 处理
                    merged.events.append(RearrangementEvent(
                        'eej', node_id,
                        genes_involved=[h1, h2],
                        desc="EEJ: {}-{} adjacent in ancestor, "
                             "broken at end in {}".format(
                            h1.hog_id, h2.hog_id, src),
                        support=1,
                        child_source=src
                    ))

            # 发出未匹配的连接事件（跳过已配对的 NCF）
            for ji, (jh1, jh2, jt) in enumerate(join_events):
                if ji in matched_joins or ji in ncf_paired_joins:
                    continue
                if jt == 'eej':
                    merged.events.append(RearrangementEvent(
                        'eej', node_id,
                        genes_involved=[jh1, jh2],
                        desc="EEJ: {}-{} joined in {} "
                             "(separate ends in ancestor)".format(
                            jh1.hog_id, jh2.hog_id, src),
                        support=1,
                        child_source=src
                    ))
                else:
                    merged.events.append(RearrangementEvent(
                        'ncf', node_id,
                        genes_involved=[jh1, jh2],
                        desc="NCF: {}-{} adjacent in {} "
                             "near end in ancestor".format(
                            jh1.hog_id, jh2.hog_id, src),
                        support=1,
                        child_source=src
                    ))

    def _detect_events_topdown(self):
        """
        自顶向下事件检测：将每个节点（内部节点和叶子）与其父节点比较，
        通过断点图（breakpoint graph）思想检测 EEJ、Fission、NCF、
        Translocation 和 Inversion 事件。
        """
        for node in self.tree.traverse(strategy="preorder"):
            if node.is_root():
                continue
            node_id = node.name
            parent_id = node.up.name
            # 对内部节点使用 anc_graphs，对叶子使用 leaf_graphs
            if node.is_leaf():
                if node_id not in self.leaf_graphs or parent_id not in self.anc_graphs:
                    continue
                anc = self.leaf_graphs[node_id]
                par = self.anc_graphs[parent_id]
            else:
                if node_id not in self.anc_graphs or parent_id not in self.anc_graphs:
                    continue
                anc = self.anc_graphs[node_id]
                par = self.anc_graphs[parent_id]
            self._compare_ancestor_to_parent(anc, par, node_id)

    def _compare_ancestor_to_parent(self, anc, par, node_id):
        """
        比较祖先图 anc 与其父节点图 par，检测重排事件并写入 anc.events。
        注意：anc 和 par 使用不同节点层级的 HOGrecord（hog_id 包含节点名），
        需要通过 hog_parent 映射建立对应关系后再比较。
        """
        # 0. 建立 anc HOG -> par HOG 的映射
        anc_to_par = {}
        par_hog_by_id = {h.hog_id: h for h in par.gene_nodes}
        for h in anc.gene_nodes:
            parent_hog_id = self.hog_parent.get(h.hog_id)
            if parent_hog_id and parent_hog_id in par_hog_by_id:
                anc_to_par[h] = par_hog_by_id[parent_hog_id]

        # 1. 邻接集合比较（映射到父节点 HOG 空间）
        par_adjs = set(par.get_adjacencies(include_telomere=False))
        anc_adjs_mapped = set()
        for h1, h2 in anc.get_adjacencies(include_telomere=False):
            p1 = anc_to_par.get(h1)
            p2 = anc_to_par.get(h2)
            if p1 and p2 and p1 != p2:
                key = (p1, p2) if p1.hog_id < p2.hog_id else (p2, p1)
                anc_adjs_mapped.add(key)
        par_adjs_norm = set()
        for h1, h2 in par_adjs:
            key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
            par_adjs_norm.add(key)

        lost = par_adjs_norm - anc_adjs_mapped   # 父节点有，anc 没有
        gained = anc_adjs_mapped - par_adjs_norm  # anc 有，父节点没有

        # 1.5 倒位检测与边界邻接移除
        # 遍历祖先染色体，映射到父节点 HOG 空间，检测倒位块，
        # 移除倒位产生的边界邻接，避免被误判为 NCF/unidir_trans
        par_hog_to_pos = {}
        for ci, chrom in enumerate(par.chromosomes):
            genes = [n for n in chrom if n not in par.telomeres]
            for j, h in enumerate(genes):
                par_hog_to_pos[h] = (ci, j)

        for chrom in anc.chromosomes:
            anc_genes = [n for n in chrom if n not in anc.telomeres]
            # 映射到父节点 HOG 并记录其在父节点中的位置
            par_mapped = []
            for h in anc_genes:
                ph = anc_to_par.get(h)
                if ph and ph in par_hog_to_pos:
                    par_mapped.append(ph)
                else:
                    par_mapped.append(None)

            # 检测倒位块：连续 HOG 在父节点中属于同一染色体但位置递减
            inv_blocks = []
            start = None
            for i in range(1, len(par_mapped)):
                prev, curr = par_mapped[i - 1], par_mapped[i]
                if prev is None or curr is None:
                    if start is not None:
                        inv_blocks.append((start, i))
                        start = None
                    continue
                prev_ci, prev_pos = par_hog_to_pos[prev]
                curr_ci, curr_pos = par_hog_to_pos[curr]
                if prev_ci == curr_ci and curr_pos < prev_pos:
                    # 位置递减 → 可能是倒位
                    if start is None:
                        start = i - 1
                else:
                    if start is not None:
                        inv_blocks.append((start, i))
                        start = None
            if start is not None:
                inv_blocks.append((start, len(par_mapped)))

            for bstart, bend in inv_blocks:
                block = par_mapped[bstart:bend]
                if len(block) < 3:
                    continue
                # 移除倒位边界邻接
                # 父节点顺序: left_nb - block[-1] - ... - block[0] - right_nb
                # 祖先顺序:   left_nb - block[0] - ... - block[-1] - right_nb
                if bstart > 0 and par_mapped[bstart - 1] is not None:
                    left_nb = par_mapped[bstart - 1]
                    # lost: left_nb - block[-1] (父节点有，祖先没有)
                    lk = (left_nb, block[-1]) if left_nb.hog_id < block[-1].hog_id \
                        else (block[-1], left_nb)
                    lost.discard(lk)
                    # gained: left_nb - block[0] (祖先有，父节点没有)
                    gk = (left_nb, block[0]) if left_nb.hog_id < block[0].hog_id \
                        else (block[0], left_nb)
                    gained.discard(gk)
                if bend < len(par_mapped) and par_mapped[bend] is not None:
                    right_nb = par_mapped[bend]
                    # lost: block[0] - right_nb (父节点有，祖先没有)
                    lk = (block[0], right_nb) if block[0].hog_id < right_nb.hog_id \
                        else (right_nb, block[0])
                    lost.discard(lk)
                    # gained: block[-1] - right_nb (祖先有，父节点没有)
                    gk = (block[-1], right_nb) if block[-1].hog_id < right_nb.hog_id \
                        else (right_nb, block[-1])
                    gained.discard(gk)

        # 将规范化后的邻接还原为 HOGrecord 对
        lost_pairs = list(lost)
        gained_pairs = list(gained)

        # 2. 父节点端点 HOG + 染色体索引
        par_ends = set()
        par_hog_to_chrom = {}
        par_chrom_ends = {}  # hog -> 'left'/'right'/'both'/'internal'
        for chrom in par.chromosomes:
            genes = [n for n in chrom if n not in par.telomeres]
            if genes:
                par_ends.add(genes[0])
                par_ends.add(genes[-1])
        for ci, chrom in enumerate(par.chromosomes):
            genes = [n for n in chrom if n not in par.telomeres]
            for j, h in enumerate(genes):
                par_hog_to_chrom[h] = ci
                if j == 0 and j == len(genes) - 1:
                    par_chrom_ends[h] = 'both'
                elif j == 0:
                    par_chrom_ends[h] = 'left'
                elif j == len(genes) - 1:
                    par_chrom_ends[h] = 'right'
                else:
                    par_chrom_ends[h] = 'internal'

        # 3. 分类 gained adjacencies（新连接）
        # 先收集 NCF 类型的 join 事件，再配对
        gained_eej = []      # [(p1, p2)] 两个父节点端点连接
        gained_ncf = []      # [(p1, p2)] 一个父节点端点 + 一个内部
        gained_other = []    # [(p1, p2)] 两个内部（不涉及端点）
        for p1, p2 in gained_pairs:
            e1 = par_chrom_ends.get(p1, 'internal')
            e2 = par_chrom_ends.get(p2, 'internal')
            if e1 != 'internal' and e2 != 'internal':
                gained_eej.append((p1, p2))
            elif e1 != 'internal' or e2 != 'internal':
                gained_ncf.append((p1, p2))
            else:
                gained_other.append((p1, p2))

        # 配对 NCF 连接事件：两个 NCF join 来自同一 donor 染色体
        # 且 recipient HOG 在父节点中相邻 → 合并为一个 NCF 事件
        par_adjs_norm_set = set()
        for n1, n2 in par.get_adjacencies(include_telomere=False):
            key = (n1, n2) if n1.hog_id < n2.hog_id else (n2, n1)
            par_adjs_norm_set.add(key)

        ncf_paired = set()  # indices of paired NCF gains
        for gi1 in range(len(gained_ncf)):
            if gi1 in ncf_paired:
                continue
            p1_1, p2_1 = gained_ncf[gi1]
            e1_1 = par_chrom_ends.get(p1_1, 'internal')
            if e1_1 != 'internal':
                donor1, recip1 = p1_1, p2_1
            else:
                donor1, recip1 = p2_1, p1_1
            for gi2 in range(gi1 + 1, len(gained_ncf)):
                if gi2 in ncf_paired:
                    continue
                p1_2, p2_2 = gained_ncf[gi2]
                e1_2 = par_chrom_ends.get(p1_2, 'internal')
                if e1_2 != 'internal':
                    donor2, recip2 = p1_2, p2_2
                else:
                    donor2, recip2 = p2_2, p1_2
                # 条件1：两个 donor 端点来自同一父节点染色体
                dc1 = par_hog_to_chrom.get(donor1)
                dc2 = par_hog_to_chrom.get(donor2)
                if dc1 is None or dc2 is None or dc1 != dc2:
                    continue
                # 条件2：两个 recipient HOG 在父节点中相邻
                rk = (recip1, recip2) if recip1.hog_id < recip2.hog_id \
                    else (recip2, recip1)
                if rk in par_adjs_norm_set:
                    all_hogs = list(set([donor1, donor2, recip1, recip2]))
                    anc.events.append(RearrangementEvent(
                        'ncf', node_id,
                        genes_involved=all_hogs,
                        desc="NCF: donor({},{}) from chrom {} inserted "
                             "between {}-{}".format(
                            donor1.hog_id, donor2.hog_id, dc1,
                            recip1.hog_id, recip2.hog_id),
                        support=2
                    ))
                    ncf_paired.add(gi1)
                    ncf_paired.add(gi2)
                    break

        # 发出 EEJ 事件
        for p1, p2 in gained_eej:
            anc.events.append(RearrangementEvent(
                'eej', node_id,
                genes_involved=[p1, p2],
                desc="EEJ: {} and {} joined (separate ends in parent)".format(
                    p1.hog_id, p2.hog_id),
                support=1
            ))

        # 发出未配对的 NCF 连接事件
        for gi, (p1, p2) in enumerate(gained_ncf):
            if gi in ncf_paired:
                continue
            anc.events.append(RearrangementEvent(
                'ncf', node_id,
                genes_involved=[p1, p2],
                desc="NCF: {}-{} near end in parent".format(
                    p1.hog_id, p2.hog_id),
                support=1
            ))

        # 4. 分类 lost adjacencies（断裂）-> fission
        # 将父节点的 lost adjacency 映射回 anc 的染色体位置
        par_to_anc = {p: a for a, p in anc_to_par.items()}
        anc_hog_to_chrom = {}
        anc_chrom_ends = set()  # anc 中端点 HOG 集合
        for i, chrom in enumerate(anc.chromosomes):
            genes = [n for n in chrom if n not in anc.telomeres]
            for j, n in enumerate(genes):
                anc_hog_to_chrom[n] = i
                if j == 0 or j == len(genes) - 1:
                    anc_chrom_ends.add(n)

        # 统计断裂情况，用于判断是否为广泛碎片化
        n_lost_cross_chrom = 0
        n_lost_endpoints = 0
        for p1, p2 in lost_pairs:
            a1 = par_to_anc.get(p1)
            a2 = par_to_anc.get(p2)
            c1 = anc_hog_to_chrom.get(a1) if a1 else None
            c2 = anc_hog_to_chrom.get(a2) if a2 else None
            if c1 is not None and c2 is not None and c1 != c2:
                n_lost_cross_chrom += 1
                if a1 in anc_chrom_ends and a2 in anc_chrom_ends:
                    n_lost_endpoints += 1

        # 如果祖先染色体数远大于父节点，且大多数断裂不涉及端点，
        # 说明是重建过度碎片化，放宽fission检测标准
        anc_chrom_count = len(list(anc.chromosomes))
        par_chrom_count = len(list(par.chromosomes))
        fragmented = anc_chrom_count > par_chrom_count * 2 and n_lost_cross_chrom > par_chrom_count * 2

        for p1, p2 in lost_pairs:
            a1 = par_to_anc.get(p1)
            a2 = par_to_anc.get(p2)
            c1 = anc_hog_to_chrom.get(a1) if a1 else None
            c2 = anc_hog_to_chrom.get(a2) if a2 else None
            if c1 is not None and c2 is not None and c1 != c2:
                # 染色体间断裂
                if a1 in anc_chrom_ends and a2 in anc_chrom_ends:
                    anc.events.append(RearrangementEvent(
                        'fission', node_id,
                        genes_involved=[p1, p2],
                        desc="Fission: {} and {} separated".format(p1.hog_id, p2.hog_id),
                        support=1
                    ))
                elif fragmented:
                    anc.events.append(RearrangementEvent(
                        'fission', node_id,
                        genes_involved=[p1, p2],
                        desc="Fission (fragmented): {} and {} separated".format(p1.hog_id, p2.hog_id),
                        support=1
                    ))
            elif c1 is not None and c2 is not None and c1 == c2:
                # 染色体内断裂：NCF 插入位点
                # 父节点中 P-Q 相邻，anc 中被 donor 打断
                anc.events.append(RearrangementEvent(
                    'ncf', node_id,
                    genes_involved=[p1, p2],
                    desc="NCF insertion site: {}-{} interrupted".format(
                        p1.hog_id, p2.hog_id),
                    support=1
                ))

        # 5. Inversion 检测：链向比较（使用映射后的 HOG）
        # 预计算父节点染色体长度
        par_chrom_lengths = {}
        par_hog_to_chrom = {}
        for pci, chrom in enumerate(par.chromosomes):
            pgenes = [n for n in chrom if n not in par.telomeres]
            par_chrom_lengths[pci] = len(pgenes)
            for g in pgenes:
                par_hog_to_chrom[g] = pci

        for chrom in anc.chromosomes:
            genes = [n for n in chrom if n not in anc.telomeres and n in anc_to_par]
            if len(genes) < 3:
                continue
            # 找出链向与父节点不一致的连续块
            blocks = []
            start = 0
            for i in range(len(genes)):
                anc_s = anc.graph.nodes[genes[i]].get('strand')
                if anc_s is None:
                    anc_s = getattr(genes[i], 'strand', '+')
                par_h = anc_to_par[genes[i]]
                par_s = par.graph.nodes[par_h].get('strand')
                if par_s is None:
                    par_s = getattr(par_h, 'strand', '+')
                if anc_s == par_s:
                    if i > start:
                        blocks.append(genes[start:i])
                    start = i + 1
            if start < len(genes):
                blocks.append(genes[start:])

            # 阈值：固定最小 3，且占染色体长度 >= 2%（避免大染色体漏掉小 inversion）
            min_block_size = max(3, int(len(genes) * 0.02) + 1)
            for block in blocks:
                if len(block) < min_block_size:
                    continue
                # 过滤：不能占父节点对应染色体的绝大部分（排除 NCF/EEJ/RT 导致的整臂翻转）
                par_chroms = {par_hog_to_chrom.get(anc_to_par[h]) for h in block}
                if len(par_chroms) != 1 or None in par_chroms:
                    continue
                par_ci = list(par_chroms)[0]
                if len(block) >= par_chrom_lengths.get(par_ci, 0) * 0.9:
                    continue
                is_left = block[0] == genes[0]
                is_right = block[-1] == genes[-1]

                if is_left or is_right:
                    etype = 'telomere_inversion'
                else:
                    etype = 'internal_inversion'
                block_par_hogs = [anc_to_par[h] for h in block]
                anc.events.append(RearrangementEvent(
                    etype, node_id,
                    genes_involved=block_par_hogs,
                    desc="{}: {} HOGs from {} to {}".format(
                        etype, len(block), block[0].hog_id, block[-1].hog_id),
                    support=len(block)
                ))

    def _get_outgroup_graphs(self, node):
        """
        获取外类群图及其系统发育距离权重。
        距离越近的外类群在投票中权重越高（反比距离加权）。
        返回 [(graph, weight), ...] 列表。
        """
        ingroup_leaves = set(node.get_leaf_names())
        all_leaves = set(self.tree.get_leaf_names())
        outgroup_leaves = all_leaves - ingroup_leaves
        result = []
        for leaf in outgroup_leaves:
            if leaf not in self.leaf_graphs:
                continue
            leaf_node = self.tree.search_nodes(name=leaf)
            if not leaf_node:
                continue
            leaf_node = leaf_node[0]
            dist = node.get_distance(leaf_node)
            # 反比距离加权：距离越近权重越高
            # 使用 1/(1+dist) 使权重范围在 (0,1] 之间，避免极端权重差异
            weight = 1.0 / (1.0 + dist)
            result.append((self.leaf_graphs[leaf], weight))
        return result

    def _optimize_round(self):
        """迭代优化：清理孤立节点"""
        for node_id, aag in self.anc_graphs.items():
            if aag.node_id in self.leaf_graphs:
                continue
            to_remove = {n for n in aag.gene_nodes if aag.graph.degree(n) == 0}
            if to_remove:
                aag.remove_nodes(to_remove)
                if self.use_ilp_sa:
                    new_aag = self._linearize_graph_ilp_sa(aag)
                else:
                    new_aag = self._linearize_graph(aag)
                new_aag._add_telomeres()
                self.anc_graphs[node_id] = new_aag

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


    def _collapse_polyploid_leaf(self, leaf_id):
        """
        将多倍体叶子基因组推断为pre-WGD基因组。
        直接将post-WGD图映射到父节点HOG，由support计数自然区分
        亚基因组内保留邻接(support=ploidy)与亚基因组间假邻接(support=1)，
        无需亚基因组拆分。
        """
        leaf_graph = self.leaf_graphs.get(leaf_id)
        if leaf_graph is None:
            return
        ploidy = self.ploidy_map.get(leaf_id, 2)

        leaf_node = self.tree.search_nodes(name=leaf_id)
        if not leaf_node:
            return
        parent = leaf_node[0].up
        if not parent:
            return
        parent_id = parent.name

        # 目标染色体数 = post-WGD染色体数 / 倍性
        target = max(len(list(leaf_graph.chromosomes)) // ploidy, 1)

        # 将整个post-WGD图作为单个子图，映射到父节点HOG
        # 支持计数自然区分：亚基因组内保守邻接 support=ploidy，
        # 亚基因组间假邻接 support=1
        pre_id = "pre-WGD {}".format(leaf_id)
        if self.use_v3:
            pre_graph = self._reconstruct_v3(
                target_node_id=pre_id,
                child_graphs=[leaf_graph],
                child_source_ids=[leaf_id],
                hog_level=parent_id,
                outgroup_graphs=None,
                target_chromosomes=target,
                is_pre_wgd=True
            )
        else:
            pre_graph = self._reconstruct(
                target_node_id=pre_id,
                child_graphs=[leaf_graph],
                child_source_ids=[leaf_id],
                hog_level=parent_id,
                outgroup_graphs=None,   # 多倍体叶子无外类群
                target_chromosomes=target,
                is_pre_wgd=True
            )

        pre_graph.node_id = "{}_pre".format(leaf_id)
        self.pre_wgd_graphs[leaf_id] = pre_graph
        n_pre = len(list(pre_graph.chromosomes))
        logger.info("  Leaf {} pre-WGD: {} post chroms -> {} pre chroms (target={})".format(
            leaf_id, len(list(leaf_graph.chromosomes)), n_pre, target))

    def _collapse_wgd(self, node):
        """
        将WGD节点的post-WGD基因组推断为pre-WGD基因组。
        直接将post-WGD图映射到父节点HOG，由support计数自然区分
        亚基因组内保留邻接与亚基因组间假邻接，
        无需亚基因组拆分。
        """
        node_id = node.name
        post_graph = self.anc_graphs.get(node_id)
        if post_graph is None:
            return
        parent = node.up
        if not parent:
            return

        parent_id = parent.name
        ploidy = self.ploidy_map.get(node_id, 2)

        # 目标染色体数 = post-WGD染色体数 / 倍性
        target = max(len(list(post_graph.chromosomes)) // ploidy, 1)

        # 将整个post-WGD图作为单个子图，映射到父节点HOG
        pre_id = "pre-WGD {}".format(node_id)
        pre_graph = self._reconstruct(
            target_node_id=pre_id,
            child_graphs=[post_graph],
            child_source_ids=[node_id],
            hog_level=parent_id,
            outgroup_graphs=None,   # pre-WGD无外类群
            target_chromosomes=target,
            is_pre_wgd=True
        )

        pre_graph.node_id = "{}_pre".format(node_id)
        self.pre_wgd_graphs[node_id] = pre_graph
        n_pre = len(list(pre_graph.chromosomes))
        logger.info("WGD pre-WGD for {}: {} post chroms -> {} pre chroms (target={})".format(
            node_id, len(list(post_graph.chromosomes)), n_pre, target))

    # =====================
    # ILP + SA 混合线性化
    # =====================

    def _linearize_graph_ilp_sa(self, aag, target_chromosomes=None):
        """
        ILP松弛 + 模拟退火混合线性化。
        Phase 1: ILP求解度<=2、边数=N-target_chrom的松弛问题（可能含环）
        Phase 2: 以ILP解为初始解，SA消除环并优化目标
        """
        if pulp is None:
            logger.warning("pulp not available, falling back to greedy linearization")
            return self._linearize_graph(aag, target_chromosomes=target_chromosomes)

        undirected = aag.graph.to_undirected()
        undirected.remove_nodes_from(aag.telomeres)
        nodes = list(undirected.nodes())

        if not nodes:
            result = AncestralAdjacencyGraph(node_id=aag.node_id)
            result.species_set = aag.species_set
            result.hog_map = dict(aag.hog_map)
            result.graph = nx.DiGraph()
            result.events = list(aag.events)
            return result

        # 计算边权重
        edge_weights = {}
        for u, v, d in aag.graph.edges(data=True):
            if u in aag.telomeres or v in aag.telomeres:
                continue
            key = (u, v) if u.hog_id < v.hog_id else (v, u)
            if key not in edge_weights:
                support = d.get('support', 1)
                og_weight = d.get('og_weight', 0)
                weight = 200 if support >= 2 else 100 + min(og_weight * 50, 100)
                edge_weights[key] = weight

        if not edge_weights:
            return self._linearize_graph(aag, target_chromosomes=target_chromosomes)

        # Phase 1: ILP松弛
        t0 = time.time()
        prob = pulp.LpProblem("Linearize", pulp.LpMaximize)

        x = {}
        for key in edge_weights:
            u, v = key
            x[key] = pulp.LpVariable("x_{}_{}".format(u.hog_id, v.hog_id), lowBound=0, upBound=1, cat='Binary')

        prob += pulp.lpSum(edge_weights[key] * x[key] for key in edge_weights)

        node_edges = defaultdict(list)
        for key in edge_weights:
            u, v = key
            node_edges[u].append(x[key])
            node_edges[v].append(x[key])

        for n in nodes:
            if n in node_edges:
                prob += pulp.lpSum(node_edges[n]) <= 2, "degree_{}".format(n.hog_id)

        n_nodes = len(nodes)
        # 孤立节点无法参与任何边，需要从边数目标中排除
        isolated_nodes = [n for n in nodes if n not in node_edges]
        n_isolated = len(isolated_nodes)
        effective_nodes = n_nodes - n_isolated

        if target_chromosomes is not None and target_chromosomes > 0:
            # 孤立节点必然各自成为单节点路径，调整目标
            adjusted_target = max(target_chromosomes - n_isolated, 0)
            target_edges = max(effective_nodes - adjusted_target, 0)
            if target_edges > 0:
                # 使用上界约束，允许更少的边（SA会优化 toward target）
                prob += pulp.lpSum(x[key] for key in edge_weights) <= target_edges, "edge_count_max"

        try:
            solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=20)
            prob.solve(solver)
        except Exception as e:
            logger.warning("ILP solve failed: {}, falling back to greedy".format(e))
            return self._linearize_graph(aag, target_chromosomes=target_chromosomes)

        if pulp.LpStatus[prob.status] not in ['Optimal', 'Feasible']:
            logger.warning("ILP status: {}, falling back to greedy".format(pulp.LpStatus[prob.status]))
            return self._linearize_graph(aag, target_chromosomes=target_chromosomes)

        selected_edges = set()
        for key, var in x.items():
            if var.value() is not None and var.value() > 0.5:
                selected_edges.add(key)

        logger.info("  ILP solved: {} edges in {:.2f}s".format(len(selected_edges), time.time() - t0))

        # Phase 2: 模拟退火
        t0 = time.time()
        current_edges = set(selected_edges)
        current_score = self._evaluate_edge_set(current_edges, edge_weights, nodes, target_chromosomes,
                                                 hog_endpoints=getattr(aag, 'hog_endpoints', None))
        best_edges = set(current_edges)
        best_score = current_score

        temperature = 100.0
        cooling_rate = 0.9995
        min_temp = 0.01
        all_edges = list(edge_weights.keys())

        iteration = 0
        max_iter = self.sa_iterations

        # 预计算节点到边的映射，加速remove_add
        node_to_edges = defaultdict(list)
        for e in all_edges:
            node_to_edges[e[0]].append(e)
            node_to_edges[e[1]].append(e)

        # 预计算current_edges中节点的度
        def compute_deg(edges):
            deg = defaultdict(int)
            for u, v in edges:
                deg[u] += 1
                deg[v] += 1
            return deg

        while temperature > min_temp and iteration < max_iter:
            iteration += 1
            op = random.choice(['remove_add', 'swap'])
            new_edges = None

            if op == 'remove_add':
                if not current_edges:
                    continue
                e_remove = random.choice(list(current_edges))
                u_r, v_r = e_remove
                base_deg = compute_deg(current_edges)
                base_deg[u_r] -= 1
                base_deg[v_r] -= 1
                # 只考虑与u_r或v_r相连的边，避免全局遍历
                candidate_set = set(node_to_edges.get(u_r, []) + node_to_edges.get(v_r, []))
                candidates = []
                for e in candidate_set:
                    if e in current_edges:
                        continue
                    u_a, v_a = e
                    if base_deg[u_a] + 1 <= 2 and base_deg[v_a] + 1 <= 2:
                        candidates.append(e)
                if not candidates:
                    continue
                e_add = random.choice(candidates)
                new_edges = set(current_edges)
                new_edges.remove(e_remove)
                new_edges.add(e_add)

            elif op == 'swap':
                if len(current_edges) < 2:
                    continue
                edges_list = list(current_edges)
                e1 = random.choice(edges_list)
                e2 = random.choice(edges_list)
                if e1 == e2:
                    continue
                u1, v1 = e1
                u2, v2 = e2
                swaps = []
                if u1 != u2 and v1 != v2:
                    ne1 = (u1, u2) if u1.hog_id < u2.hog_id else (u2, u1)
                    ne2 = (v1, v2) if v1.hog_id < v2.hog_id else (v2, v1)
                    if ne1 in edge_weights and ne2 in edge_weights:
                        swaps.append((ne1, ne2))
                if u1 != v2 and v1 != u2:
                    ne1 = (u1, v2) if u1.hog_id < v2.hog_id else (v2, u1)
                    ne2 = (v1, u2) if v1.hog_id < u2.hog_id else (u2, v1)
                    if ne1 in edge_weights and ne2 in edge_weights:
                        swaps.append((ne1, ne2))
                if not swaps:
                    continue
                ne1, ne2 = random.choice(swaps)
                new_edges = set(current_edges)
                new_edges.remove(e1)
                new_edges.remove(e2)
                new_edges.add(ne1)
                new_edges.add(ne2)

            if new_edges is None:
                continue

            new_score = self._evaluate_edge_set(new_edges, edge_weights, nodes, target_chromosomes,
                                                 hog_endpoints=getattr(aag, 'hog_endpoints', None))
            delta = new_score - current_score

            if delta > 0 or random.random() < math.exp(delta / max(temperature, 0.001)):
                current_edges = new_edges
                current_score = new_score
                if current_score > best_score:
                    best_edges = set(current_edges)
                    best_score = current_score

            temperature *= cooling_rate

        logger.info("  SA: {} iters, score {:.1f} -> {:.1f} in {:.2f}s".format(
            iteration, self._evaluate_edge_set(selected_edges, edge_weights, nodes, target_chromosomes,
                                               hog_endpoints=getattr(aag, 'hog_endpoints', None)),
            best_score, time.time() - t0))

        # Phase 3: 构建结果图
        return self._edges_to_aag(aag, best_edges, nodes, edge_weights, target_chromosomes=target_chromosomes)

    def _evaluate_edge_set(self, edges, edge_weights, nodes, target_chromosomes, hog_endpoints=None):
        """评估边集质量：权重和 - 惩罚项
        端粒约束：只有NCF/EEJ/fission应改变染色体数量。
        对过多的连通分量、小片段、孤立节点施加强惩罚。
        新增：端粒HOG应优先作为路径端点（度<=1），内部使用（度==2）暗示重排断点。
        """
        score = sum(edge_weights.get(e, 0) for e in edges)

        degree = {n: 0 for n in nodes}
        for u, v in edges:
            degree[u] = degree.get(u, 0) + 1
            degree[v] = degree.get(v, 0) + 1

        penalty = 0
        endpoint_set = set(hog_endpoints.keys()) if hog_endpoints else set()
        for n, d in degree.items():
            if d > 2:
                penalty += (d - 2) * 500

        subg = nx.Graph()
        for n in nodes:
            subg.add_node(n)
        for u, v in edges:
            subg.add_edge(u, v)

        # 组件数惩罚：强制接近目标染色体数
        if target_chromosomes is not None:
            n_components = nx.number_connected_components(subg)
            penalty += abs(n_components - target_chromosomes) * 800

            # 碎片惩罚：小片段（<3节点）强烈惩罚，不应被判定为染色体
            for comp in nx.connected_components(subg):
                comp_size = len(comp)
                if comp_size == 1:
                    penalty += 300  # 孤立节点
                elif comp_size == 2:
                    penalty += 150  # 2节点碎片

            # 端粒端点支持度评分
            # 路径端点（度==1）应尽可能对应子染色体端点（hog_endpoints），
            # 否则意味着发生了裂变（生物学罕见）
            supported_ends = 0
            unsupported_ends = 0
            for comp in nx.connected_components(subg):
                if len(comp) == 1:
                    n = list(comp)[0]
                    if n in endpoint_set:
                        supported_ends += 2
                    else:
                        unsupported_ends += 2
                else:
                    sg = subg.subgraph(comp)
                    for n in comp:
                        if sg.degree(n) == 1:
                            if n in endpoint_set:
                                supported_ends += 1
                            else:
                                unsupported_ends += 1
            # 不支持的端点暗示裂变事件，给予惩罚
            penalty += unsupported_ends * 120
            # 支持的端点数量应接近 2 * target_chromosomes
            expected_ends = target_chromosomes * 2
            penalty += abs(supported_ends - expected_ends) * 25

        n_cycles = 0
        for comp in nx.connected_components(subg):
            sg = subg.subgraph(comp)
            if sg.number_of_edges() >= len(comp) and len(comp) > 1:
                n_cycles += 1
        penalty += n_cycles * 300

        isolated = sum(1 for n in nodes if degree.get(n, 0) == 0)
        penalty += isolated * 200

        return score - penalty

    def _edges_to_aag(self, aag, edges, nodes, edge_weights, target_chromosomes=None):
        """将优化后的边集转换回AncestralAdjacencyGraph"""
        subg = nx.Graph()
        for n in nodes:
            subg.add_node(n)
        for u, v in edges:
            subg.add_edge(u, v)

        new_graph = nx.DiGraph()
        seen_nodes = set()

        for comp_nodes in nx.connected_components(subg):
            if len(comp_nodes) == 1:
                n = list(comp_nodes)[0]
                new_graph.add_node(n)
                seen_nodes.add(n)
                continue

            comp = subg.subgraph(comp_nodes).copy()

            # 打破环
            if comp.number_of_edges() == len(comp_nodes):
                weakest = None
                min_w = float('inf')
                for u, v in comp.edges():
                    key = (u, v) if u.hog_id < v.hog_id else (v, u)
                    w = edge_weights.get(key, 0)
                    if w < min_w:
                        min_w = w
                        weakest = (u, v)
                if weakest:
                    comp.remove_edge(*weakest)

            endpoints = [n for n in comp_nodes if comp.degree(n) == 1]
            if endpoints:
                start = endpoints[0]
            else:
                start = list(comp_nodes)[0]

            curr = start
            prev = None
            while True:
                nbrs = [n for n in comp.neighbors(curr) if n != prev]
                if not nbrs:
                    break
                nxt = nbrs[0]
                new_graph.add_edge(curr, nxt)
                prev = curr
                curr = nxt

            seen_nodes.update(comp_nodes)

        for n in nodes:
            if n not in seen_nodes:
                new_graph.add_node(n)

        for n in new_graph.nodes():
            if n in aag.graph:
                for k, v in aag.graph.nodes[n].items():
                    if k not in new_graph.nodes[n]:
                        new_graph.nodes[n][k] = v

        result = AncestralAdjacencyGraph(node_id=aag.node_id)
        result.species_set = aag.species_set
        result.hog_map = dict(aag.hog_map)
        result.graph = new_graph
        result.gene_nodes = set(new_graph.nodes())
        result.events = list(aag.events)
        for n in new_graph.nodes():
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom

        # 端粒约束后处理：丢弃碎片路径
        result = self._filter_fragments(result, aag, target_chromosomes)
        return result
