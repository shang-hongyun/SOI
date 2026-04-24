import sys
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
        'translocation',
        'reciprocal_translocation',
        'eej',
        'fission',
        'ncf',
    }

    def __init__(self, event_type, node,
                 genes_involved=None, parent_chroms=None,
                 desc='', support=None):
        if event_type not in self.EVENT_TYPES:
            logger.warning('Unknown event type: {}'.format(event_type))
        self.event_type = event_type
        self.node = node
        self.genes_involved = genes_involved or []
        self.parent_chroms = parent_chroms or []
        self.desc = desc
        self.support = support

    def __repr__(self):
        return '<{} on {}: {}>'.format(self.event_type, self.node, self.desc)

    def to_dict(self):
        return {
            'event_type': self.event_type,
            'node': self.node,
            'genes_involved': [str(g) for g in self.genes_involved],
            'parent_chroms': self.parent_chroms,
            'desc': self.desc,
            'support': self.support,
        }


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
        """主运行函数"""
        logger.info("=== 开始祖先核型重建 (AKR) ===")
        self._start_time = time.time()

        self._build_hogs()
        self._parse_ploidy_map()
        self._build_leaf_graphs()

        for node in self.tree.traverse(strategy="postorder"):
            if node.is_leaf():
                continue
            elapsed = time.time() - self._start_time
            if self.timeout > 0 and elapsed > self.timeout:
                logger.warning("Global timeout ({}s) reached at {:.1f}s, skipping remaining nodes".format(
                    self.timeout, elapsed))
                break
            self._reconstruct_node(node)

            # WGD节点：额外生成pre-WGD基因组
            if node.name in self.ploidy_map and self.ploidy_map[node.name] > 1:
                self._collapse_wgd(node)

        for i in range(self.rounds):
            logger.info('Optimization round {}'.format(i))
            self._optimize_round()

        self._detect_events_topdown()
        self._export_results()

        logger.info("=== 祖先核型重建完成 ===")
        return self.anc_graphs

    def _build_hogs(self):
        """运行 HOG 流程并建立索引"""
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
        """为每个叶物种构建 AncestralAdjacencyGraph"""
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
                    aag.hog_map[gene] = self.gene_to_all_hogs[gene.id][0]
                    mapped += 1
            self.leaf_graphs[sp] = aag
            self.anc_graphs[sp] = aag
            logger.info("Built leaf graph for {}: {} genes, {} chromosomes".format(
                sp, len(aag.gene_nodes), len(list(aag.chromosomes))))

    def _reconstruct_node(self, node):
        """
        重建一个内部节点的祖先核型。
        核心流程：
        1. 将两个子节点图逐级映射到当前节点的 HOG
        2. 合并映射后的邻接图（产生冲突边）
        3. 外类群反向映射到当前节点 HOG
        4. 基于外类群投票解决合并冲突
        5. 线性化并重建端粒结构
        6. 逐个解决小规模重排（indel / duplication / inversion）
        7. 解决 telomere-centric 大规模重排
        """
        node_id = node.name
        logger.info("Reconstructing node {}".format(node_id))
        t0 = time.time()

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

        # 打印子节点信息
        for cg in child_graphs:
            n_chrom = len(list(cg.chromosomes))
            n_genes = len(cg.gene_nodes)
            n_edges = len(list(cg.get_adjacencies(include_telomere=False)))
            logger.info("  Child {}: {} genes, {} chromosomes, {} adjacencies".format(
                cg.node_id, n_genes, n_chrom, n_edges))

        # ============================================================
        # Step 1: 将每个子节点映射到当前节点的 HOGrecord
        # ============================================================
        t1 = time.time()
        mapped_children = [self._map_to_parent_hogs(node_id, cg) for cg in child_graphs]
        for i, mc in enumerate(mapped_children):
            n_edges = len(list(mc.get_adjacencies(include_telomere=False)))
            logger.info("  Mapped child {} -> {} HOGs, {} adjacencies ({:.2f}s)".format(
                child_graphs[i].node_id, len(mc.gene_nodes), n_edges, time.time()-t1))

        # ============================================================
        # Step 2: 合并两个映射后的邻接图
        # ============================================================
        t1 = time.time()
        merged = self._merge_two_graphs(node_id, mapped_children[0], mapped_children[1])
        logger.info("  After merge: {} HOGs, {} conflict edges ({:.2f}s)".format(
            len(merged.gene_nodes),
            sum(1 for u, v, d in merged.graph.edges(data=True) if d.get('conflict')),
            time.time()-t1))

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

        # 导出重建前 GFA（原始合并结果，含冲突）
        t1 = time.time()
        merge_gfa = "{}.{}.merge.gfa".format(self.outpre, node_id)
        with open(merge_gfa, 'w') as fout:
            merged.to_gfa(fout)
        logger.info("  Exported pre-reconstruction GFA: {} ({:.2f}s)".format(merge_gfa, time.time()-t1))

        # ============================================================
        # Step 3: 外类群反向映射到当前节点 HOG（带系统发育距离权重）
        # ============================================================
        t1 = time.time()
        outgroup_graphs = self._get_outgroup_graphs(node)
        mapped_outgroups = [(self._map_outgroup_to_current_hogs(node_id, og), weight)
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
        target_chromosomes = self._estimate_target_chromosomes(node)
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
        t3 = time.time()
        merged = self._resolve_inversions(merged, mapped_children, child_graphs)
        logger.info("  Small rearrangements: indel {:.2f}s, dup {:.2f}s, inv {:.2f}s".format(
            t2-t1, t3-t2, time.time()-t3))

        # 节点超时检查：若已接近超时，跳过剩余步骤
        if self.node_timeout > 0 and (time.time() - t0) > self.node_timeout:
            logger.warning("  Node {} approaching timeout, skipping telomere rearrangements".format(node_id))
        else:
            # ============================================================
            # Step 6: Telomere-centric 大规模重排
            # ============================================================
            t1 = time.time()
            merged = self._resolve_telomere_rearrangements(node, merged, mapped_children, mapped_outgroups)
            logger.info("  Telomere rearrangements: {:.2f}s".format(time.time()-t1))

        self.anc_graphs[node_id] = merged

        # 导出重建后 GFA
        anc_gfa = "{}.{}.anc.gfa".format(self.outpre, node_id)
        with open(anc_gfa, 'w') as fout:
            merged.to_gfa(fout)
        logger.info("  Exported post-reconstruction GFA: {}".format(anc_gfa))

        # 事件分类统计与染色体数目验证
        n_chrom = len(list(merged.chromosomes))
        if merged.events:
            evt_counts = Counter(e.event_type for e in merged.events)
            evt_summary = ', '.join('{}:{}'.format(k, v) for k, v in sorted(evt_counts.items()))
            # 染色体数目验证：fission (+1), ncf (-1), eej (-1)
            n_fission = evt_counts.get('fission', 0)
            n_ncf = evt_counts.get('ncf', 0)
            n_eej = evt_counts.get('eej', 0)
            chrom_delta = n_fission - n_ncf - n_eej
            # 计算子节点染色体数（取平均）
            child_chroms = []
            for cg in child_graphs:
                child_chroms.append(len(list(cg.chromosomes)))
            avg_child_chrom = sum(child_chroms) / len(child_chroms) if child_chroms else n_chrom
            expected_delta = avg_child_chrom - n_chrom
            logger.info("Node {} final: {} chromosomes, events summary: {}".format(
                node_id, n_chrom, evt_summary))
            logger.info("  Chromosome validation: fission({}) - ncf({}) - eej({}) = delta {}".format(
                n_fission, n_ncf, n_eej, chrom_delta))
            logger.info("  Expected delta (child {} - anc {}) = {}, match={}".format(
                avg_child_chrom, n_chrom, expected_delta,
                "YES" if chrom_delta == expected_delta else "NO"))
        else:
            logger.info("Node {} final: {} chromosomes, 0 events".format(
                node_id, n_chrom))

        node_elapsed = time.time() - t0
        if self.node_timeout > 0 and node_elapsed > self.node_timeout:
            logger.warning("Node {} reconstruction took {:.1f}s, exceeding node_timeout {}s".format(
                node_id, node_elapsed, self.node_timeout))

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

    def _map_to_parent_hogs(self, node_id, graph):
        """
        将子图（叶子图或祖先图）的节点映射到指定 node_id 的 HOGrecord，
        返回一个新的 AncestralAdjacencyGraph，节点全部为 HOGrecord 对象。
        边按照原图的邻接关系复制，去除自环和端粒边。
        """
        mapping = self._map_child_node_to_hog(node_id, graph)
        mapped = AncestralAdjacencyGraph(node_id="{}_mapped".format(graph.node_id))
        mapped.species_set = set(graph.species_set)

        # 添加节点：只添加成功映射的 HOGrecord
        for rec in set(mapping.values()):
            mapped.graph.add_node(rec, hog=rec)
            mapped.gene_nodes.add(rec)
            mapped.hog_map[rec] = rec

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
        for chrom in graph.chromosomes:
            if not chrom:
                continue
            genes = [n for n in chrom if n not in graph.telomeres]
            if not genes:
                continue
            first_g = genes[0]
            last_g = genes[-1]
            left_tel = chrom[0] if chrom[0] in graph.telomeres else None
            right_tel = chrom[-1] if chrom[-1] in graph.telomeres else None
            for g in genes:
                if g not in mapping:
                    continue
                h = mapping[g]
                if g == first_g:
                    mapped.hog_endpoints[h]['left'].append((graph.node_id, left_tel))
                if g == last_g:
                    mapped.hog_endpoints[h]['right'].append((graph.node_id, right_tel))

        # 复制边：只复制两端都成功映射的邻接，且去除自环和双向重复
        seen_edges = set()
        for n1, n2 in graph.get_adjacencies(include_telomere=False):
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

    def _merge_two_graphs(self, node_id, g1, g2):
        """
        合并两个已经映射到同一层 HOG 的邻接图。
        节点取并集，边按支持度投票：
        - support=2：两个子节点都支持，高置信度保留
        - support=1：仅一个子节点支持，标记为冲突边
        同时标记 conflict 边是否涉及子节点染色体端点（端点邻接更可能是重排断点），
        以及冲突边在另一子节点中是否位于同一染色体（同一染色体更可能是祖先邻接断裂）。
        返回合并后的 AncestralAdjacencyGraph。
        """
        merged = AncestralAdjacencyGraph(node_id=node_id)
        merged.species_set = g1.species_set | g2.species_set

        # 节点并集
        all_hogs = set(g1.gene_nodes) | set(g2.gene_nodes)
        for rec in all_hogs:
            merged.graph.add_node(rec, hog=rec)
            merged.gene_nodes.add(rec)
            merged.hog_map[rec] = rec

        # 预计算每个子节点的 HOG -> 染色体索引 / 位置
        child_chrom_info = []
        for g in (g1, g2):
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
        for i, g in enumerate((g1, g2)):
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

        # 为 conflict 边计算：在另一子节点中是否位于同一染色体及距离
        same_chrom_other = {}
        dist_other = {}
        for key, count in edge_support.items():
            if count >= 2:
                continue
            child_idx = list(edge_children[key])[0]
            other_idx = 1 - child_idx
            h1, h2 = key
            other_h2c = child_chrom_info[other_idx][0]
            other_h2p = child_chrom_info[other_idx][1]
            c1 = other_h2c.get(h1)
            c2 = other_h2c.get(h2)
            if c1 is not None and c2 is not None and c1 == c2:
                same_chrom_other[key] = True
                p1 = other_h2p.get(h1, 0)
                p2 = other_h2p.get(h2, 0)
                dist_other[key] = abs(p1 - p2)
            else:
                same_chrom_other[key] = False
                dist_other[key] = 999999

        for (h1, h2), count in edge_support.items():
            is_endpoint = edge_endpoint.get((h1, h2), False)
            merged.graph.add_edge(h1, h2, support=count, conflict=(count < 2), endpoint=is_endpoint)
            merged.graph.add_edge(h2, h1, support=count, conflict=(count < 2), endpoint=is_endpoint)

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

    def _estimate_target_chromosomes(self, node):
        """
        估计目标染色体数目。
        使用直接子节点的染色体数，取最大值作为估计。
        对于根节点，同时考虑所有叶子的最大值（祖先数可能大于任何直接后代）。
        """
        child_counts = []
        for child in node.children:
            child_name = child.name
            if child_name in self.anc_graphs:
                child_counts.append(len(list(self.anc_graphs[child_name].chromosomes)))
            elif child_name in self.leaf_graphs:
                child_counts.append(len(list(self.leaf_graphs[child_name].chromosomes)))

        # 回退到后代叶子最大值
        clade_leaves = set(node.get_leaf_names())
        leaf_counts = []
        for leaf in clade_leaves:
            if leaf in self.leaf_graphs:
                leaf_counts.append(len(list(self.leaf_graphs[leaf].chromosomes)))

        if not child_counts:
            return max(leaf_counts) if leaf_counts else 1

        base = max(max(child_counts), max(leaf_counts) if leaf_counts else 0)

        # 启发式修正：当所有子节点染色体数相同且较小时，它们可能在各自分支上
        # 独立丢失了染色体，导致祖先实际染色体数更多。
        if len(child_counts) >= 2 and len(set(child_counts)) == 1:
            same_count = child_counts[0]
            if same_count <= 5:
                base += 1

        # 根节点额外修正：若所有叶子的染色体数都较少，根祖先可能更多。
        if node.is_root() and leaf_counts:
            if max(leaf_counts) <= 7:
                base += 1

        return base

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
        if target_chromosomes is not None and count_components() > target_chromosomes:
            for (u, v), weight in low_quality:
                if degree[u] >= 2 or degree[v] >= 2:
                    continue
                if find(u) == find(v):
                    continue
                selected.add((u, v))
                degree[u] += 1
                degree[v] += 1
                union(u, v)
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
        for n in new_graph.nodes():
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom
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

    def _resolve_indels(self, merged, mapped_children, mapped_outgroups):
        """
        解决 Indel（基因丢失/插入）。
        若某 HOG 仅出现在部分子节点中，且外类群也缺失该 HOG，
        则推断该 HOG 是在有它的子节点中独立插入的，应从祖先图中删除。
        注意：若某子节点本身无数据（如外类群叶节点无基因），不应计入分母。
        """
        node_id = merged.node_id
        # 统计每个 HOG 在哪些子节点中出现
        hog_to_children = defaultdict(set)
        for mc in mapped_children:
            for rec in mc.gene_nodes:
                hog_to_children[rec].add(mc.node_id)

        # 收集外类群中的 HOG（用于判断祖先是否应有），按系统发育距离加权
        # HOG权重 >= 总权重的 1/3 才视为外类群支持存在
        outgroup_hog_weights = defaultdict(float)
        for og, weight in mapped_outgroups:
            for rec in og.gene_nodes:
                outgroup_hog_weights[rec] += weight
        total_weight = sum(w for _, w in mapped_outgroups) if mapped_outgroups else 0
        weight_threshold = total_weight / 3.0 if total_weight > 0 else 0
        outgroup_hogs = {hog for hog, w in outgroup_hog_weights.items()
                         if w >= weight_threshold}

        # 只统计有数据的子节点
        children_with_data = [mc for mc in mapped_children if len(mc.gene_nodes) > 0]
        n_children_data = len(children_with_data)

        genes_to_remove = set()
        for hog_rec in list(merged.gene_nodes):
            child_present = len(hog_to_children.get(hog_rec, set()))
            outgroup_present = hog_rec in outgroup_hogs
            # 若该 HOG 在有数据的子节点中并非全部出现，且外类群也缺失，则判为 indel
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
        解决 Duplication。
        若子节点映射后一个 HOG 对应多个原始基因，记录 duplication 事件。
        祖先图本身不创建多拷贝（每个 HOG 只保留一个代表节点）。
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

    def _resolve_inversions(self, merged, mapped_children, child_graphs):
        """
        底向上 Inversion 检测已禁用：底向上策略把子节点分支上的反转误归因到当前节点，
        导致大量假阳性。Inversion 检测统一由顶向下阶段 (_detect_events_topdown) 完成，
        该阶段比较祖先与其父节点的链向，能准确识别发生在该分支上的 inversion。
        """
        return merged

    def _resolve_telomere_rearrangements(self, node, merged, mapped_children, mapped_outgroups):
        """
        Telomere-centric 大规模重排：
        EEJ, Fission, NCF, Reciprocal Translocation
        基于断点图（breakpoint graph）思想，分析 HOG 端粒邻接关系。
        所有分析均在当前节点 HOG 层级进行。
        """
        node_id = node.name

        # 染色体数目记录（用于后续验证，但不基于外类群推断fission/eej）
        anc_chrom_count = len(merged.get_chromosome_ends())
        child_chrom_counts = [len(list(mc.chromosomes)) for mc in mapped_children]
        outgroup_chrom_counts = [len(list(og.chromosomes)) for og, _ in mapped_outgroups]

        # 只在祖先染色体数显著偏离子节点范围时发出警告，不生成虚假事件
        if child_chrom_counts:
            min_child = min(child_chrom_counts)
            max_child = max(child_chrom_counts)
            if anc_chrom_count > max_child * 2:
                logger.warning("  Node {} has {} chromosomes, much higher than child range [{}-{}], possible over-fragmentation".format(
                    node_id, anc_chrom_count, min_child, max_child))

        # NCF 和 telomere 变化检测：
        # 1. 子节点染色体的两端 HOG 在 merged 中位于不同染色体 -> NCF
        # 2. 子节点端点 HOG 在 merged 中变为内部节点 -> telomere inversion 或其他端粒重排
        # 缓存 merged 染色体列表和每个HOG的位置信息
        merged_chroms_cache = list(merged.chromosomes)
        merged_hogs_by_chrom = []
        merged_hog_positions = {}  # hog -> {'chrom_idx': i, 'is_left_end': bool, 'is_right_end': bool}
        for i, chrom in enumerate(merged_chroms_cache):
            hogs = [g for g in chrom if g not in merged.telomeres]
            merged_hogs_by_chrom.append(hogs)
            for j, h in enumerate(hogs):
                pos = {'chrom_idx': i, 'is_left_end': (j == 0), 'is_right_end': (j == len(hogs) - 1)}
                merged_hog_positions[h] = pos

        for mc in mapped_children:
            for left_tel, first_g, last_g, right_tel in mc.get_chromosome_ends():
                if not first_g or not last_g:
                    continue
                h_first = first_g if first_g not in mc.telomeres else None
                h_last = last_g if last_g not in mc.telomeres else None
                if not h_first or not h_last:
                    continue
                pos_first = merged_hog_positions.get(h_first)
                pos_last = merged_hog_positions.get(h_last)

                # NCF 和 telomere inversion 检测已移至顶向下阶段（_detect_events_topdown），
                # 底向上阶段仅用于重建祖先图，不记录事件，避免将分支差异误判为节点事件。
                pass

        return merged

    def _detect_events_topdown(self):
        """
        自顶向下事件检测：将每个内部节点与其父节点比较，
        通过断点图（breakpoint graph）思想检测 EEJ、Fission、NCF、
        Translocation 和 Inversion 事件。
        """
        for node in self.tree.traverse(strategy="preorder"):
            if node.is_root() or node.is_leaf():
                continue
            node_id = node.name
            parent_id = node.up.name
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

        # 将规范化后的邻接还原为 HOGrecord 对
        lost_pairs = list(lost)
        gained_pairs = list(gained)

        # 2. 父节点端点 HOG
        par_ends = set()
        for chrom in par.chromosomes:
            genes = [n for n in chrom if n not in par.telomeres]
            if genes:
                par_ends.add(genes[0])
                par_ends.add(genes[-1])

        # 3. 分类 gained adjacencies（新连接）
        # 只有涉及父节点端点的邻接变化才报告为事件（EEJ 或 NCF）。
        for p1, p2 in gained_pairs:
            p1_end = p1 in par_ends
            p2_end = p2 in par_ends
            if p1_end and p2_end:
                etype = 'eej'
            elif p1_end or p2_end:
                etype = 'ncf'
            else:
                continue
            anc.events.append(RearrangementEvent(
                etype, node_id,
                genes_involved=[p1, p2],
                desc="{}: {} and {}".format(etype.upper(), p1.hog_id, p2.hog_id),
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

        for p1, p2 in lost_pairs:
            a1 = par_to_anc.get(p1)
            a2 = par_to_anc.get(p2)
            c1 = anc_hog_to_chrom.get(a1) if a1 else None
            c2 = anc_hog_to_chrom.get(a2) if a2 else None
            # 关键过滤：要求 a1 和 a2 在 anc 中都是端点，避免重建内部断裂误判为 fission
            if c1 is not None and c2 is not None and c1 != c2:
                if a1 in anc_chrom_ends and a2 in anc_chrom_ends:
                    anc.events.append(RearrangementEvent(
                        'fission', node_id,
                        genes_involved=[p1, p2],
                        desc="Fission: {} and {} separated".format(p1.hog_id, p2.hog_id),
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

        event_file = "{}.events.tsv".format(self.outpre)
        with open(event_file, 'w') as f:
            header = ["node", "event_type", "genes", "chroms", "desc", "support"]
            f.write('\t'.join(header) + '\n')
            for node_id, aag in self.anc_graphs.items():
                if node_id in self.leaf_graphs:
                    continue
                for ev in aag.events:
                    genes = ','.join(str(g) for g in ev.genes_involved)
                    chroms = ','.join(str(c) for c in ev.parent_chroms)
                    line = [node_id, ev.event_type, genes, chroms, ev.desc, str(ev.support)]
                    f.write('\t'.join(line) + '\n')
        logger.info("Exported rearrangement events to {}".format(event_file))

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

    # =====================
    # WGD pre-WGD 重建
    # =====================

    def _collapse_wgd(self, node):
        """
        将WGD节点的post-WGD基因组推断为pre-WGD基因组。
        核心思路：
        1. 将post-WGD HOG映射到父节点层级的HOG
        2. 基于parent HOG集合的Jaccard相似度，将post-WGD染色体配对（聚类）
        3. 对每对同源染色体，按平均位置投票确定pre-WGD的基因顺序
        4. 生成pre-WGD的AncestralAdjacencyGraph，存入self.pre_wgd_graphs
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

        # 1. 建立post-WGD HOG -> parent HOG映射
        mapping = self._map_child_node_to_hog(parent_id, post_graph)
        if not mapping:
            logger.warning("No mapping from {} to parent {} for WGD collapse".format(node_id, parent_id))
            return

        # 2. 收集每条post-WGD染色体的parent信息
        chroms = []
        for chrom in post_graph.chromosomes:
            genes = [n for n in chrom if n not in post_graph.telomeres]
            parent_seq = [mapping.get(g) for g in genes if g in mapping]
            parent_set = set(p for p in parent_seq if p is not None)
            chroms.append({
                'genes': genes,
                'parent_seq': parent_seq,
                'parent_set': parent_set,
            })

        n = len(chroms)
        if n == 0:
            return

        # 3. 配对染色体：基于parent_set Jaccard相似度
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                s1, s2 = chroms[i]['parent_set'], chroms[j]['parent_set']
                if not s1 or not s2:
                    continue
                inter = len(s1 & s2)
                union = len(s1 | s2)
                sim = inter / union if union > 0 else 0
                pairs.append((sim, i, j))
        pairs.sort(reverse=True)

        # 贪心匹配：每条染色体最多匹配一次
        matched = set()
        groups = []
        for sim, i, j in pairs:
            if i in matched or j in matched:
                continue
            if sim < 0.2:
                break
            groups.append([i, j])
            matched.add(i)
            matched.add(j)

        # 未匹配的单独成组
        for i in range(n):
            if i not in matched:
                groups.append([i])

        # 4. 对每组，推断pre-WGD染色体
        pre_chroms = []
        for group in groups:
            parent_pos = defaultdict(list)
            for chrom_idx in group:
                for pos, p in enumerate(chroms[chrom_idx]['parent_seq']):
                    if p is not None:
                        parent_pos[p].append((chrom_idx, pos))
            if not parent_pos:
                continue
            sorted_parents = sorted(
                parent_pos.keys(),
                key=lambda p: sum(pos for _, pos in parent_pos[p]) / len(parent_pos[p])
            )
            pre_chroms.append(sorted_parents)

        # 5. 构建pre-WGD图
        pre_graph = AncestralAdjacencyGraph(node_id="{}_pre".format(node_id))
        pre_graph.species_set = set(post_graph.species_set)

        # 链向投票
        strand_votes = defaultdict(lambda: {'+': 0, '-': 0})
        for chrom in post_graph.chromosomes:
            genes = [n for n in chrom if n not in post_graph.telomeres]
            for g in genes:
                if g in mapping:
                    strand = post_graph.graph.nodes[g].get('strand', '+')
                    strand_votes[mapping[g]][strand] += 1

        for chrom in pre_chroms:
            for p in chrom:
                pre_graph.graph.add_node(p, hog=p)
                pre_graph.gene_nodes.add(p)
                pre_graph.hog_map[p] = p
                consensus = '+' if strand_votes[p]['+'] >= strand_votes[p]['-'] else '-'
                pre_graph.graph.nodes[p]['strand'] = consensus
            for i in range(len(chrom) - 1):
                pre_graph.graph.add_edge(chrom[i], chrom[i + 1])
                pre_graph.graph.add_edge(chrom[i + 1], chrom[i])

        pre_graph._add_telomeres()
        self.pre_wgd_graphs[node_id] = pre_graph
        logger.info("WGD collapse for {}: {} post -> {} pre chromosomes (groups={})".format(
            node_id, n, len(pre_chroms), [len(g) for g in groups]))

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
        if target_chromosomes is not None and target_chromosomes > 0:
            target_edges = max(n_nodes - target_chromosomes, 0)
            prob += pulp.lpSum(x[key] for key in edge_weights) == target_edges, "edge_count"

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
        current_score = self._evaluate_edge_set(current_edges, edge_weights, nodes, target_chromosomes)
        best_edges = set(current_edges)
        best_score = current_score

        temperature = 100.0
        cooling_rate = 0.9995
        min_temp = 0.01
        all_edges = list(edge_weights.keys())

        iteration = 0
        max_iter = self.sa_iterations

        while temperature > min_temp and iteration < max_iter:
            iteration += 1
            op = random.choice(['remove_add', 'swap'])
            new_edges = None

            if op == 'remove_add':
                if not current_edges:
                    continue
                e_remove = random.choice(list(current_edges))
                u_r, v_r = e_remove
                candidates = []
                for e in all_edges:
                    if e in current_edges:
                        continue
                    u_a, v_a = e
                    # 计算当前度
                    deg = {n: 0 for n in nodes}
                    for e2 in current_edges:
                        deg[e2[0]] += 1
                        deg[e2[1]] += 1
                    # 移除e_remove后的度
                    deg[u_r] -= 1
                    deg[v_r] -= 1
                    # 添加e后的度
                    deg[u_a] += 1
                    deg[v_a] += 1
                    if all(deg[n] <= 2 for n in (u_a, v_a, u_r, v_r)):
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

            new_score = self._evaluate_edge_set(new_edges, edge_weights, nodes, target_chromosomes)
            delta = new_score - current_score

            if delta > 0 or random.random() < math.exp(delta / max(temperature, 0.001)):
                current_edges = new_edges
                current_score = new_score
                if current_score > best_score:
                    best_edges = set(current_edges)
                    best_score = current_score

            temperature *= cooling_rate

        logger.info("  SA: {} iters, score {:.1f} -> {:.1f} in {:.2f}s".format(
            iteration, self._evaluate_edge_set(selected_edges, edge_weights, nodes, target_chromosomes),
            best_score, time.time() - t0))

        # Phase 3: 构建结果图
        return self._edges_to_aag(aag, best_edges, nodes, edge_weights)

    def _evaluate_edge_set(self, edges, edge_weights, nodes, target_chromosomes):
        """评估边集质量：权重和 - 惩罚项"""
        score = sum(edge_weights.get(e, 0) for e in edges)

        degree = {n: 0 for n in nodes}
        for u, v in edges:
            degree[u] = degree.get(u, 0) + 1
            degree[v] = degree.get(v, 0) + 1

        penalty = 0
        for n, d in degree.items():
            if d > 2:
                penalty += (d - 2) * 500

        if target_chromosomes is not None:
            n_components = len(nodes) - len(edges)
            penalty += abs(n_components - target_chromosomes) * 200

        subg = nx.Graph()
        for n in nodes:
            subg.add_node(n)
        for u, v in edges:
            subg.add_edge(u, v)

        n_cycles = 0
        for comp in nx.connected_components(subg):
            sg = subg.subgraph(comp)
            if sg.number_of_edges() >= len(comp) and len(comp) > 1:
                n_cycles += 1
        penalty += n_cycles * 300

        isolated = sum(1 for n in nodes if degree.get(n, 0) == 0)
        penalty += isolated * 100

        return score - penalty

    def _edges_to_aag(self, aag, edges, nodes, edge_weights):
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
        for n in new_graph.nodes():
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom
        return result
