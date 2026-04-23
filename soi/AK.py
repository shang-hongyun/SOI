import sys
import itertools
import random
from collections import Counter, defaultdict, OrderedDict
from typing import Optional, Set, Dict, List, Tuple
import networkx as nx

from .mcscan import ColinearGroups, Gff, GffGraph, SyntenyGraph, GffLine
from .hog import HOG, HOGrecord
from .RunCmdsMP import logger


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
        for start in starts:
            chrom_nodes = list(self.iter_chrom(start))
            if not chrom_nodes:
                continue
            first = chrom_nodes[0]
            last = chrom_nodes[-1]
            if hasattr(first, 'chrom'):
                chrom = first.chrom
            elif hasattr(last, 'chrom'):
                chrom = last.chrom
            else:
                continue
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
        while True:
            succs = list(self.graph.successors(current))
            if not succs:
                break
            current = succs[0]
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
        gg = GffGraph()
        for n in self.graph.nodes():
            gg.add_node(n, **dict(self.graph.nodes[n]))
        for n1, n2 in self.graph.edges():
            gg.add_edge(n1, n2, **dict(self.graph[n1][n2]))
        return gg


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

        self.hog = None
        self.tree = None
        self.leaf_graphs = {}
        self.anc_graphs = {}
        self.events = defaultdict(list)

        self.hogs_by_node = defaultdict(list)
        self.gene_to_hog = {}
        self.node_gene_to_hog = defaultdict(dict)

    def run(self):
        """主运行函数"""
        logger.info("=== 开始祖先核型重建 (AKR) ===")

        self._build_hogs()
        self._build_leaf_graphs()

        for node in self.tree.traverse(strategy="postorder"):
            if node.is_leaf():
                continue
            self._reconstruct_node(node)

        for i in range(self.rounds):
            logger.info('Optimization round {}'.format(i))
            self._optimize_round()

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
            for gene in rec.genes:
                self.gene_to_hog[gene] = hog_id
                self.node_gene_to_hog[rec.node_id][gene] = hog_id

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

        # 诊断：检查 gene_to_hog 与 GFF 基因 ID 的格式匹配
        sample_hog_keys = list(self.gene_to_hog.keys())[:5]
        sample_gff_ids = []
        for line in Gff(self.gfffile):
            if line.species in self.hog.species:
                sample_gff_ids.append(line.id)
            if len(sample_gff_ids) >= 5:
                break
        logger.info("Sample gene_to_hog keys: {}".format(sample_hog_keys))
        logger.info("Sample GFF IDs: {}".format(sample_gff_ids))

        # 自动检测并修复格式不匹配
        has_pipe_hog = any('|' in k for k in sample_hog_keys)
        has_pipe_gff = any('|' in k for k in sample_gff_ids)

        if has_pipe_hog and not has_pipe_gff:
            # HOG 是 SPECIES|GENEID，GFF 是纯 GENEID -> 提取 HOG 的后半部分建立映射
            logger.info("Format mismatch detected: HOG has SPECIES|GENEID, GFF has plain GENEID. Auto-fixing...")
            self._gene_to_hog_fixed = {}
            for gid, hid in self.gene_to_hog.items():
                plain = gid.split('|', 1)[1] if '|' in gid else gid
                self._gene_to_hog_fixed[plain] = hid
            self.gene_to_hog = self._gene_to_hog_fixed
            # 同步修复 node_gene_to_hog
            for node_id, mapping in self.node_gene_to_hog.items():
                fixed = {}
                for gid, hid in mapping.items():
                    plain = gid.split('|', 1)[1] if '|' in gid else gid
                    fixed[plain] = hid
                self.node_gene_to_hog[node_id] = fixed
        elif not has_pipe_hog and has_pipe_gff:
            # HOG 是纯 GENEID，GFF 是 SPECIES|GENEID -> 为 GFF 提取后半部分匹配
            logger.info("Format mismatch detected: HOG has plain GENEID, GFF has SPECIES|GENEID. Auto-fixing...")
            self._gene_to_hog_fixed = {}
            for gid, hid in self.gene_to_hog.items():
                self._gene_to_hog_fixed[gid] = hid
            self.gene_to_hog = self._gene_to_hog_fixed
            # 同步修复 node_gene_to_hog（此处不需要转换键，因为映射保持原样）
        else:
            # 格式一致或无法判断，直接验证交集比例
            matched = sum(1 for gid in sample_gff_ids if gid in self.gene_to_hog)
            logger.info("Format check: {}/{} sample GFF IDs matched in gene_to_hog".format(
                matched, len(sample_gff_ids)))

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
            aag = AncestralAdjacencyGraph.from_gffgraph(sp, gg, add_telomere=True)
            aag.species_set = {sp}
            mapped = 0
            sample_gff_ids = []
            for i, gene in enumerate(aag.gene_nodes):
                if i < 5:
                    sample_gff_ids.append(gene.id)
                if gene.id in self.gene_to_hog:
                    aag.hog_map[gene] = self.gene_to_hog[gene.id]
                    mapped += 1
            logger.info("Species {} GFF sample IDs: {}".format(sp, sample_gff_ids))
            logger.info("Species {} hog_map: {}/{} genes mapped to {} HOGs".format(
                sp, mapped, len(aag.gene_nodes), len(set(aag.hog_map.values()))))
            self.leaf_graphs[sp] = aag
            self.anc_graphs[sp] = aag
            logger.info("Built leaf graph for {}: {} genes, {} chromosomes".format(
                sp, len(aag.gene_nodes), len(list(aag.chromosomes))))

    def _reconstruct_node(self, node):
        """重建一个内部节点的祖先核型"""
        node_id = node.name
        logger.info("Reconstructing node {}".format(node_id))

        children = node.children
        if len(children) != 2:
            logger.warning("Node {} has {} children; expecting binary tree".format(
                node_id, len(children)))

        child_graphs = []
        for child in children:
            child_id = child.name
            if child_id in self.anc_graphs:
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

        # Step A: 合并子图为祖先候选图
        anc = self._merge_child_graphs(node_id, child_graphs)
        logger.info("  After merge: {} genes, {} chromosomes, {} edges".format(
            len(anc.gene_nodes), len(list(anc.chromosomes)), len(list(anc.get_adjacencies(include_telomere=False)))))

        # 导出重建前 GFA
        merge_gfa = "{}.{}.merge.gfa".format(self.outpre, node_id)
        with open(merge_gfa, 'w') as fout:
            anc.to_gffgraph().to_gfa(fout)
        logger.info("  Exported pre-reconstruction GFA: {}".format(merge_gfa))

        # Step B: 检测小规模重排
        anc = self._resolve_small_rearrangements(node, anc, child_graphs)
        logger.info("  After small rearrangements: {} genes, {} chromosomes, {} events".format(
            len(anc.gene_nodes), len(list(anc.chromosomes)), len(anc.events)))

        # Step C: Telomere-centric 大规模重排
        anc = self._resolve_telomere_rearrangements(node, anc, child_graphs)
        logger.info("  After telomere rearrangements: {} genes, {} chromosomes, {} events".format(
            len(anc.gene_nodes), len(list(anc.chromosomes)), len(anc.events)))

        self.anc_graphs[node_id] = anc

        # 导出重建后 GFA
        anc_gfa = "{}.{}.anc.gfa".format(self.outpre, node_id)
        with open(anc_gfa, 'w') as fout:
            anc.to_gffgraph().to_gfa(fout)
        logger.info("  Exported post-reconstruction GFA: {}".format(anc_gfa))

        # 事件分类统计
        if anc.events:
            evt_counts = Counter(e.event_type for e in anc.events)
            evt_summary = ', '.join('{}:{}'.format(k, v) for k, v in sorted(evt_counts.items()))
            logger.info("Node {} final: {} chromosomes, events summary: {}".format(
                node_id, len(list(anc.chromosomes)), evt_summary))
        else:
            logger.info("Node {} final: {} chromosomes, 0 events".format(
                node_id, len(list(anc.chromosomes))))

    def _merge_child_graphs(self, node_id, child_graphs):
        """
        基于当前节点的 HOG 合并子节点图为祖先候选图。
        重建节点 N 时，使用 N 的 HOG ID（如 OG0001.N.hog0）统一映射两个子节点，
        这样共享 HOG 数应接近该节点的总 HOG 数（10000+）。
        """
        anc = AncestralAdjacencyGraph(node_id=node_id)
        anc.species_set = set()

        # 使用当前节点 node_id 的 HOG 映射（而非子图自带的 hog_map）
        node_hog_map = self.node_gene_to_hog.get(node_id, {})
        if not node_hog_map:
            logger.warning("  Node {} has no HOG mappings!".format(node_id))

        hog_to_genes = defaultdict(list)
        for cg in child_graphs:
            anc.species_set.update(cg.species_set)
            for gene in cg.gene_nodes:
                gid = getattr(gene, 'id', str(gene))
                if gid in node_hog_map:
                    hog_id = node_hog_map[gid]
                    hog_to_genes[hog_id].append((gene, cg.node_id))

        # 诊断：打印 HOG 分组统计
        logger.info("  HOG grouping (node {}): {} unique HOGs from {} child graphs".format(
            node_id, len(hog_to_genes), len(child_graphs)))
        for cg in child_graphs:
            cg_hogs = set()
            for gene in cg.gene_nodes:
                gid = getattr(gene, 'id', str(gene))
                if gid in node_hog_map:
                    cg_hogs.add(node_hog_map[gid])
            logger.info("    Child {} contributed {} HOGs".format(cg.node_id, len(cg_hogs)))

        hog_representatives = {}
        for hog_id, gene_list in hog_to_genes.items():
            reps = []
            seen_sp = set()
            for gene, sp in gene_list:
                sp_name = getattr(gene, 'species', sp)
                if sp_name not in seen_sp:
                    seen_sp.add(sp_name)
                    reps.append(gene)
            if reps:
                hog_representatives[hog_id] = reps[0]
                anc.hog_map[reps[0]] = hog_id

        adj_support = Counter()
        for cg in child_graphs:
            for chrom in cg.chromosomes:
                hog_seq = []
                for gene in chrom:
                    if gene in cg.telomeres:
                        continue
                    gid = getattr(gene, 'id', str(gene))
                    hid = node_hog_map.get(gid)
                    if hid:
                        hog_seq.append(hid)
                filtered = []
                for hid in hog_seq:
                    if not filtered or filtered[-1] != hid:
                        filtered.append(hid)
                for i in range(len(filtered) - 1):
                    h1, h2 = filtered[i], filtered[i + 1]
                    if h1 != h2:
                        key = tuple(sorted([h1, h2]))
                        adj_support[key] += 1

        for hog_id, rep in hog_representatives.items():
            anc.graph.add_node(rep)
            anc.gene_nodes.add(rep)
            if hasattr(rep, 'chrom'):
                anc.chrom_map[rep] = rep.chrom

        for (h1, h2), count in adj_support.items():
            if h1 in hog_representatives and h2 in hog_representatives:
                r1, r2 = hog_representatives[h1], hog_representatives[h2]
                anc.graph.add_edge(r1, r2, support=count)
                anc.graph.add_edge(r2, r1, support=count)

        anc = self._linearize_graph(anc)
        anc._add_telomeres()
        return anc

    def _linearize_graph(self, aag):
        """将无向连通分量贪婪串成染色体路径"""
        undirected = aag.graph.to_undirected()
        new_graph = nx.DiGraph()
        seen_edges = set()

        for cmpt in nx.connected_components(undirected):
            cmpt = list(cmpt)
            if len(cmpt) == 1:
                new_graph.add_node(cmpt[0])
                continue
            subg = undirected.subgraph(cmpt)
            ends = [n for n in cmpt if subg.degree(n) == 1]
            start = ends[0] if ends else cmpt[0]

            path = [start]
            visited = {start}
            while True:
                curr = path[-1]
                neighbors = [n for n in subg.neighbors(curr) if n not in visited]
                if not neighbors:
                    break
                nxt = neighbors[0]
                path.append(nxt)
                visited.add(nxt)

            for i in range(len(path) - 1):
                n1, n2 = path[i], path[i + 1]
                if (n1, n2) not in seen_edges:
                    new_graph.add_edge(n1, n2)
                    seen_edges.add((n1, n2))

        result = AncestralAdjacencyGraph(node_id=aag.node_id)
        result.species_set = aag.species_set
        result.hog_map = dict(aag.hog_map)
        result.graph = new_graph
        result.gene_nodes = set(new_graph.nodes())
        for n in new_graph.nodes():
            if hasattr(n, 'chrom'):
                result.chrom_map[n] = n.chrom
        return result

    def _resolve_small_rearrangements(self, node, anc, child_graphs):
        """
        识别小规模重排：
        - duplication (tandem, proximal, dispersed)
        - indel（外类群投票）
        - inversion（内部 / 端粒）
        - translocation（小规模非相互易位）
        """
        node_id = node.name
        outgroup_graphs = self._get_outgroup_graphs(node)

        # Indel：外类群一致缺失则删除
        genes_to_remove = set()
        for gene in anc.gene_nodes:
            hog_id = anc.hog_map.get(gene)
            if not hog_id:
                continue
            outgroup_total = len(outgroup_graphs)
            outgroup_presence = sum(1 for og in outgroup_graphs
                                    if hog_id in set(og.hog_map.values()))
            child_presence = sum(1 for cg in child_graphs
                                 if hog_id in set(cg.hog_map.values()))

            if outgroup_total > 0 and outgroup_presence == 0 and child_presence < len(child_graphs):
                genes_to_remove.add(gene)
                self.events[node_id].append(RearrangementEvent(
                    'indel', node_id,
                    genes_involved=[gene],
                    desc="HOG {} absent in outgroup, inferred deletion".format(hog_id),
                    support="{}/{} outgroups".format(outgroup_presence, outgroup_total)
                ))

        if genes_to_remove:
            anc.remove_nodes(genes_to_remove)
            anc = self._linearize_graph(anc)
            anc._add_telomeres()

        # Duplication
        for cg in child_graphs:
            hog_counts = Counter(cg.hog_map.values())
            for hog_id, count in hog_counts.items():
                if count <= 1:
                    continue
                copies = [g for g, h in cg.hog_map.items() if h == hog_id]
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
                self.events[node_id].append(RearrangementEvent(
                    etype, node_id,
                    genes_involved=copies,
                    desc="{}: HOG {} x{} in {}".format(etype, hog_id, count, cg.node_id),
                    support=min_dist
                ))

        # Inversion 改进：检测染色体上的链方向块（strand blocks）
        # 大片段 inversion 会表现为连续多个基因链方向与两侧区域相反
        for cg in child_graphs:
            for chrom in cg.chromosomes:
                genes = [g for g in chrom if g not in cg.telomeres and hasattr(g, 'strand')]
                if len(genes) < 3:
                    continue
                # 划分 strand block
                blocks = []
                start = 0
                curr_strand = genes[0].strand
                for i in range(1, len(genes)):
                    if genes[i].strand != curr_strand:
                        blocks.append((start, i - 1, curr_strand, genes[start:i]))
                        start = i
                        curr_strand = genes[i].strand
                blocks.append((start, len(genes) - 1, curr_strand, genes[start:]))

                if len(blocks) <= 1:
                    continue

                # 假设最大的 block 为祖先方向，其余为 inversion
                blocks_sorted = sorted(blocks, key=lambda b: len(b[3]), reverse=True)
                main_block = blocks_sorted[0]
                for start_idx, end_idx, strand, block_genes in blocks:
                    if len(block_genes) < 2:
                        continue
                    if (start_idx, end_idx, strand) == (main_block[0], main_block[1], main_block[2]):
                        continue
                    # 判断是内部 inversion 还是端粒 inversion
                    if start_idx == 0:
                        etype = 'telomere_inversion'
                    elif end_idx == len(genes) - 1:
                        etype = 'telomere_inversion'
                    else:
                        etype = 'internal_inversion'
                    hog_ids = [cg.hog_map.get(g) for g in block_genes]
                    hog_ids = [h for h in hog_ids if h]
                    self.events[node_id].append(RearrangementEvent(
                        etype, node_id,
                        genes_involved=block_genes,
                        desc="{}: {} genes (strand={}) in {} chrom, HOGs={}..{}".format(
                            etype, len(block_genes), strand, cg.node_id,
                            hog_ids[0] if hog_ids else '?',
                            hog_ids[-1] if hog_ids else '?'),
                        support=len(block_genes)
                    ))

        return anc

    def _resolve_telomere_rearrangements(self, node, anc, child_graphs):
        """
        Telomere-centric 大规模重排：
        EEJ, Fission, NCF, Reciprocal Translocation
        """
        node_id = node.name
        outgroup_graphs = self._get_outgroup_graphs(node)

        anc_ends = anc.get_chromosome_ends()
        outgroup_chrom_counts = [len(list(og.chromosomes)) for og in outgroup_graphs]
        anc_chrom_count = len(anc_ends)

        if outgroup_chrom_counts:
            median_outgroup = sorted(outgroup_chrom_counts)[len(outgroup_chrom_counts) // 2]
            if anc_chrom_count < median_outgroup:
                self.events[node_id].append(RearrangementEvent(
                    'eej', node_id,
                    desc="Chromosome number reduced from ~{} to {}".format(median_outgroup, anc_chrom_count),
                    support=outgroup_chrom_counts
                ))
            elif anc_chrom_count > median_outgroup:
                self.events[node_id].append(RearrangementEvent(
                    'fission', node_id,
                    desc="Chromosome number increased from ~{} to {}".format(median_outgroup, anc_chrom_count),
                    support=outgroup_chrom_counts
                ))

        # NCF
        for cg in child_graphs:
            for left_tel, first_g, last_g, right_tel in cg.get_chromosome_ends():
                if not first_g or not last_g:
                    continue
                h_first = cg.hog_map.get(first_g)
                h_last = cg.hog_map.get(last_g)
                if not h_first or not h_last:
                    continue
                anc_chroms = list(anc.chromosomes)
                chrom_first = None
                chrom_last = None
                for i, chrom in enumerate(anc_chroms):
                    genes = [g for g in chrom if g not in anc.telomeres]
                    hogs = [anc.hog_map.get(g) for g in genes]
                    if h_first in hogs:
                        chrom_first = i
                    if h_last in hogs:
                        chrom_last = i
                if chrom_first is not None and chrom_last is not None and chrom_first != chrom_last:
                    outgroup_separated = 0
                    for og in outgroup_graphs:
                        og_chroms = list(og.chromosomes)
                        og_first = None
                        og_last = None
                        for j, ochrom in enumerate(og_chroms):
                            og_genes = [g for g in ochrom if g not in og.telomeres]
                            og_hogs = [og.hog_map.get(g) for g in og_genes]
                            if h_first in og_hogs:
                                og_first = j
                            if h_last in og_hogs:
                                og_last = j
                        if og_first is not None and og_last is not None and og_first != og_last:
                            outgroup_separated += 1
                    self.events[node_id].append(RearrangementEvent(
                        'ncf', node_id,
                        genes_involved=[first_g, last_g],
                        desc="NCF: ends of {} chrom fused from anc chroms {} and {}".format(
                            cg.node_id, chrom_first, chrom_last),
                        support=outgroup_separated
                    ))

        # Reciprocal Translocation
        if len(child_graphs) == 2:
            cg1, cg2 = child_graphs
            ends1 = cg1.get_chromosome_ends()
            ends2 = cg2.get_chromosome_ends()

            def build_end_graph(ends, source_graph):
                g = nx.Graph()
                for left_tel, first_g, last_g, right_tel in ends:
                    h_left = source_graph.hog_map.get(first_g) if first_g else None
                    h_right = source_graph.hog_map.get(last_g) if last_g else None
                    if h_left and h_right:
                        g.add_edge(h_left, h_right)
                return g

            g1 = build_end_graph(ends1, cg1)
            g2 = build_end_graph(ends2, cg2)
            shared_ends = set(g1.nodes()) & set(g2.nodes())
            for h in shared_ends:
                partners1 = set(g1.neighbors(h))
                partners2 = set(g2.neighbors(h))
                if partners1 != partners2 and partners1 and partners2:
                    self.events[node_id].append(RearrangementEvent(
                        'reciprocal_translocation', node_id,
                        desc="End {} partners changed: {} -> {}".format(h, partners1, partners2),
                        support=1
                    ))

        return anc

    def _get_outgroup_graphs(self, node):
        """获取外类群图"""
        ingroup_leaves = set(node.get_leaf_names())
        all_leaves = set(self.tree.get_leaf_names())
        outgroup_leaves = all_leaves - ingroup_leaves
        return [self.leaf_graphs[leaf] for leaf in outgroup_leaves
                if leaf in self.leaf_graphs]

    def _optimize_round(self):
        """迭代优化：清理孤立节点"""
        for node_id, aag in self.anc_graphs.items():
            if aag.node_id in self.leaf_graphs:
                continue
            to_remove = {n for n in aag.gene_nodes if aag.graph.degree(n) == 0}
            if to_remove:
                aag.remove_nodes(to_remove)
                new_aag = self._linearize_graph(aag)
                new_aag._add_telomeres()
                self.anc_graphs[node_id] = new_aag

    def _export_results(self):
        """导出祖先核型和重排事件"""
        for node_id, aag in self.anc_graphs.items():
            if node_id in self.leaf_graphs:
                continue
            prefix = "{}.{}".format(self.outpre, node_id)
            gg = aag.to_gffgraph()
            try:
                gg.to_wgdi(prefix)
                logger.info("Exported {} to {}".format(node_id, prefix))
            except Exception as e:
                logger.warning("Failed to export {}: {}".format(node_id, e))

        event_file = "{}.events.tsv".format(self.outpre)
        with open(event_file, 'w') as f:
            header = ["node", "event_type", "genes", "chroms", "desc", "support"]
            f.write('\t'.join(header) + '\n')
            for node_id, events in self.events.items():
                for ev in events:
                    genes = ','.join(str(g) for g in ev.genes_involved)
                    chroms = ','.join(str(c) for c in ev.parent_chroms)
                    line = [node_id, ev.event_type, genes, chroms, ev.desc, str(ev.support)]
                    f.write('\t'.join(line) + '\n')
        logger.info("Exported rearrangement events to {}".format(event_file))
