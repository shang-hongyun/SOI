#!/usr/bin/env python3
"""
Phase 2 indel 测试 — 每种场景独立，复杂度递增。

不走模拟器，直接构造 ColoredGraph。
验证 find_indel_shortcuts + resolve_indels 的正确性。
"""
import sys, os, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import logging
logging.basicConfig(level=logging.WARNING)

from soi.takr_colored_graph import ColoredGraph


# ── 辅助 ─────────────────────────────────────────────────────

class MockChild:
    """最小化 child graph，只提供 add_child 需要的接口。"""
    def __init__(self, chromosomes, telomeres=None):
        self._chromosomes = chromosomes
        if telomeres is None:
            telomeres = set()
            for chrom in chromosomes:
                if chrom:
                    telomeres.add(chrom[0])
                    telomeres.add(chrom[-1])
        self.telomeres = telomeres
        self.gene_nodes = set()
        for chrom in chromosomes:
            for n in chrom:
                if n not in self.telomeres:
                    self.gene_nodes.add(n)

    @property
    def chromosomes(self):
        return self._chromosomes

    def get_adjacencies(self, include_telomere=False):
        adjs = set()
        for chrom in self._chromosomes:
            for i in range(len(chrom) - 1):
                n1, n2 = chrom[i], chrom[i + 1]
                if not include_telomere and (n1 in self.telomeres or n2 in self.telomeres):
                    continue
                adjs.add((n1, n2))
        return adjs


def T(name, side):
    """Telomere node."""
    return (name, side)


def linear_chrom(genes, chrom_name='c'):
    """构造一条染色体 [T_L, g1, g2, ..., gN, T_R]。"""
    return [T(chrom_name, 'L')] + list(genes) + [T(chrom_name, 'R')]


def build_graph(children_spec):
    """children_spec: [(child_id, [chrom1, chrom2, ...]), ...]"""
    G = ColoredGraph(hog_level='test')
    for cid, chroms in children_spec:
        G.add_child(cid, MockChild(chroms))
    return G


def events_by_type(G):
    from collections import Counter
    return dict(Counter(e.event_type for e in G.events))


# ═══════════════════════════════════════════════════════════════
#  INDEL — 基因插入/删除
# ═══════════════════════════════════════════════════════════════

class TestIndelBasic:
    """基础 indel 测试：单基因插入/删除。"""

    def test_single_gene_deletion(self):
        """单基因删除：c1 有 A-B-C-D-E，c2 有 A-B-D-E（C 被删）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 indel shortcut, got {len(shortcuts)}"
        # 检查 spanned HOG 包含 C
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned, f"Expected 'C' in spanned, got {spanned}"

    def test_single_gene_insertion(self):
        """单基因插入：c1 有 A-B-D-E，c2 有 A-B-C-D-E（C 被插入）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABDE', 'c1')]),
            ('c2', [linear_chrom('ABCDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 indel shortcut, got {len(shortcuts)}"
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned, f"Expected 'C' in spanned, got {spanned}"

    def test_resolve_removes_edge(self):
        """resolve_indels 应移除跨越边并记录事件。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABDE', 'c2')]),
        ])
        edges_before = G._graph.number_of_edges()
        G.resolve_indels()
        edges_after = G._graph.number_of_edges()
        # 应移除至少 1 条边
        assert edges_after < edges_before, \
            f"Expected edges removed: {edges_before} → {edges_after}"
        # 应记录事件
        indel_events = [e for e in G.events
                        if e.event_type in ('gene_loss', 'gene_indel', 'indel')]
        assert len(indel_events) >= 1, \
            f"Expected ≥1 indel event, got {events_by_type(G)}"


class TestIndelMultiGene:
    """多基因 indel 测试。"""

    def test_2gene_deletion(self):
        """2 基因删除：c1 有 A-B-C-D-E，c2 有 A-B-E（C,D 被删）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 shortcut, got {len(shortcuts)}"
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned and 'D' in spanned, \
            f"Expected C,D in spanned, got {spanned}"

    def test_3gene_insertion(self):
        """3 基因插入：c1 有 A-B-F，c2 有 A-B-C-D-E-F。"""
        G = build_graph([
            ('c1', [linear_chrom('ABF', 'c1')]),
            ('c2', [linear_chrom('ABCDEF', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 shortcut, got {len(shortcuts)}"
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned and 'D' in spanned and 'E' in spanned, \
            f"Expected C,D,E in spanned, got {spanned}"


class TestIndelPosition:
    """不同位置的 indel。"""

    def test_deletion_at_start(self):
        """染色体开头删除：c1 有 A-B-C-D，c2 有 B-C-D（A 被删）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCD', 'c1')]),
            ('c2', [linear_chrom('BCD', 'c2')]),
        ])
        # A 在 c1 但不在 c2 — 这是 gene_loss，不是 indel shortcut
        # indel shortcut 需要 h1,h2 都在两个孩子中
        shortcuts = G.find_indel_shortcuts()
        # 开头删除可能不被检测为 indel（A 不在 c2 中）
        # 但 resolve_indels 应该能处理
        G.resolve_indels()
        # 至少不应该崩溃
        assert True

    def test_deletion_at_end(self):
        """染色体末尾删除：c1 有 A-B-C-D，c2 有 A-B-C（D 被删）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCD', 'c1')]),
            ('c2', [linear_chrom('ABC', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        # D 在 c1 但不在 c2 — 不是 indel shortcut
        G.resolve_indels()
        assert True

    def test_internal_deletion(self):
        """内部删除：c1 有 A-B-C-D-E-F，c2 有 A-B-D-E-F（C 被删）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABDEF', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 shortcut, got {len(shortcuts)}"


class TestIndelNotRearrangement:
    """不应被误判为 indel 的重排。"""

    def test_inversion_not_indel(self):
        """倒位不应被检测为 indel。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABDCEF', 'c2')]),  # C-D-E 倒位
        ])
        shortcuts = G.find_indel_shortcuts()
        # 倒位的 spanned HOGs 在两个孩子中都存在 → 应被识别为重排
        # 不应该有 indel shortcut
        for h1, h2, child_id, spanned in shortcuts:
            # 如果有 shortcut，spanned 中的 HOG 不应全在 child 中
            child_hogs = set('ABCDEF')
            n_in_child = sum(1 for h in spanned if h in child_hogs)
            assert n_in_child <= len(spanned) * 0.5, \
                f"Rearrangement misclassified as indel: spanned={spanned}"

    def test_no_indel_when_same_genes(self):
        """两个孩子基因完全相同：不应有 indel。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABCDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 0, \
            f"Expected 0 shortcuts for identical genes, got {len(shortcuts)}"


class TestIndelMultiple:
    """多个 indel 同时存在。"""

    def test_2_deletions(self):
        """两个独立删除。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ABCEGHIJ', 'c2')]),  # D 删, F 删
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, \
            f"Expected ≥1 shortcut, got {len(shortcuts)}"
        G.resolve_indels()
        indel_events = [e for e in G.events
                        if e.event_type in ('gene_loss', 'gene_indel', 'indel')]
        assert len(indel_events) >= 1, \
            f"Expected ≥1 indel event, got {events_by_type(G)}"


# ═══════════════════════════════════════════════════════════════
#  模拟器集成测试
# ═══════════════════════════════════════════════════════════════

class TestIndelSimulation:
    """用模拟器数据验证 indel 检测。"""

    def test_sim_with_indels(self, tmp_path):
        """模拟含基因得失的数据，验证 indel 检测。"""
        from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
        from soi.AK import AKR
        from soi.evolution_simulator_ak import parse_tree
        import random

        rng = random.Random(42)
        tree = generate_tree(4, rng)
        sim = EvolutionSimulator(
            seed=42, num_chroms=3, min_genes=100, max_genes=200,
            inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
            fission_rate=0, unidir_trans_rate=0,
            gene_gain_rate=5, gene_loss_rate=5,
            tandem_dup_rate=0, dispersed_dup_rate=0,
            seg_dup_rate=0, wgd_rate=0,
        )
        orig_nw = tree.write(format=1)
        sim.run(tree, {})
        sim_dir = str(tmp_path / 'sim_indel')
        recon_dir = str(tmp_path / 'sim_indel' / 'recon')
        os.makedirs(recon_dir, exist_ok=True)
        sim.generate_outputs(sim_dir, tree, {}, orig_nw)

        akr = AKR(
            ogfile=f'{sim_dir}/ortholog_groups.txt',
            orthfiles=[f'{sim_dir}/ortholog_pairs.txt'],
            gfffile=f'{sim_dir}/all_species_gene.gff',
            sptreefile=f'{sim_dir}/species_tree.nwk',
            outpre=f'{recon_dir}/AKR',
            reconstruction_algorithm='v4_colored', min_genes=0, timeout=600,
        )
        akr._build_hogs()
        akr._build_leaf_graphs()
        tree2, _, _ = parse_tree(f'{sim_dir}/species_tree.nwk')

        # 构建合并图并检测 indel
        for node in tree2.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            for cname in children:
                mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                                  source_id=cname)
                G.add_child(cname, mapped)
            shortcuts = G.find_indel_shortcuts()
            # 有 gene_gain/gene_loss 的数据应产生 indel shortcuts
            # 但不强制要求（取决于模拟器的事件分布）
            if shortcuts:
                G.resolve_indels()
                indel_events = [e for e in G.events
                                if e.event_type in ('gene_loss', 'gene_indel', 'indel')]
                assert len(indel_events) >= 1, \
                    f"Expected ≥1 indel event after resolve, got {events_by_type(G)}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
