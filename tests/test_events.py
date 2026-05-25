#!/usr/bin/env python3
"""
事件级单元测试 — 每种事件类型独立测试，复杂度递增。

不走模拟器，直接构造 ColoredGraph。
每个测试只含一种事件类型。
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
    return [T(f'{chrom_name}', 'L')] + list(genes) + [T(f'{chrom_name}', 'R')]


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
#  INVERSION
# ═══════════════════════════════════════════════════════════════

class TestInversion:
    """倒位：A-B-C-D → A-C-B-D。

    无向图中产生 3-cycle（三角形），不是 4-cycle。
    块级检测需要 3-cycle 中 2 条唯一边来自不同孩子 + 1 条共享边。
    """

    def test_simple_6hogs(self):
        """最简倒位：6 HOGs，1 个倒位。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABDCEF', 'c2')]),  # C-D-E 段倒位
        ])
        G.resolve_all_events()
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        assert len(inv) >= 1, f"6-HOG inversion: {events_by_type(G)}"

    def test_10hogs_1inversion(self):
        """10 HOGs，1 个内部倒位。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ABCFEDGHIJ', 'c2')]),  # D-E-F 倒位
        ])
        G.resolve_all_events()
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        assert len(inv) >= 1, f"10-HOG inversion: {events_by_type(G)}"

    def test_20hogs_1inversion(self):
        """20 HOGs，1 个倒位 — 测试规模影响。"""
        genes = 'ABCDEFGHIJKLMNOPQRST'
        inverted = 'ABCDEFGKJIHLMNOPQRST'  # H-I-J-K 倒位
        G = build_graph([
            ('c1', [linear_chrom(genes, 'c1')]),
            ('c2', [linear_chrom(inverted, 'c2')]),
        ])
        G.resolve_all_events()
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        assert len(inv) >= 1, f"20-HOG inversion: {events_by_type(G)}"

    def test_2inversions(self):
        """10 HOGs，2 个独立倒位。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ACBDFEHGIJ', 'c2')]),  # A-B-C 和 E-F-G 各倒一次
        ])
        G.resolve_all_events()
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        assert len(inv) >= 2, f"2 inversions: {events_by_type(G)}"

    def test_direction_conflict_detected(self):
        """倒位必须产生方向冲突边。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABDCEF', 'c2')]),
        ])
        n_conflict = sum(1 for h1, h2 in G.shared_edges()
                         if G.edge_has_direction_conflict(h1, h2))
        assert n_conflict > 0, "Inversion should create direction conflict edges"


# ═══════════════════════════════════════════════════════════════
#  RECIPROCAL TRANSLOCATION (RT)
# ═══════════════════════════════════════════════════════════════

class TestRT:
    """RT：两条染色体交换片段。

    块级产生 4-cycle，4 种颜色（不同 child+chrom_idx）。
    """

    def test_simple_rt(self):
        """最简 RT：2 孩子各 2 条染色体，交换尾部。"""
        G = build_graph([
            ('c1', [
                linear_chrom('ABCD', 'c1_0'),
                linear_chrom('EFGH', 'c1_1'),
            ]),
            ('c2', [
                linear_chrom('ABGH', 'c2_0'),  # C-D 换成 G-H
                linear_chrom('EFCD', 'c2_1'),  # G-H 换成 C-D
            ]),
        ])
        G.resolve_all_events()
        rt = [e for e in G.events if e.event_type == 'reciprocal_translocation']
        assert len(rt) >= 1, f"Simple RT: {events_by_type(G)}"

    def test_rt_10hogs_per_chrom(self):
        """RT：每条染色体 10 HOGs。"""
        G = build_graph([
            ('c1', [
                linear_chrom('ABCDEFGHIJ', 'c1_0'),
                linear_chrom('KLMNOPQRST', 'c1_1'),
            ]),
            ('c2', [
                linear_chrom('ABCDEFPQRST', 'c2_0'),  # G-J 换成 P-Q-R-S-T
                linear_chrom('KLMNOGHIJ', 'c2_1'),    # P-T 换成 G-H-I-J
            ]),
        ])
        G.resolve_all_events()
        rt = [e for e in G.events if e.event_type == 'reciprocal_translocation']
        assert len(rt) >= 1, f"10-HOG RT: {events_by_type(G)}"


# ═══════════════════════════════════════════════════════════════
#  EEJ (End-to-End Join / Fusion)
# ═══════════════════════════════════════════════════════════════

class TestEEJ:
    """EEJ：两条染色体端对端融合。

    孩子 A: [chr1, chr2] → 孩子 B: [chr1+chr2]
    桥接检测：唯一边连接两个不同共享组件。
    """

    def test_simple_eej(self):
        """最简 EEJ：3+3 HOGs，一条边连接两条染色体。"""
        G = ColoredGraph(hog_level='test')
        # c1: 两条独立染色体（祖先态）
        G.add_child('c1', MockChild([
            linear_chrom('ABC', 'c1_0'),
            linear_chrom('DEF', 'c1_1'),
        ]))
        # c2: A-B-C-D-E-F 连成一条（融合）
        G.add_child('c2', MockChild([
            linear_chrom('ABCDEF', 'c2_0'),
        ]))
        # 手动加桥接边 C-D (c2 unique)
        # c2 有 C-D 邻接（在同一条染色体上），c1 没有
        # 实际上 c1 的 C 和 D 在不同染色体，所以 C-D 是 unique edge
        G.resolve_all_events()
        eej = [e for e in G.events if e.event_type == 'eej']
        ncf = [e for e in G.events if e.event_type == 'ncf']
        bridge = [e for e in G.events if e.event_type == 'bridge_unclassified']
        assert len(eej) + len(ncf) + len(bridge) >= 1, \
            f"EEJ: eej={len(eej)}, ncf={len(ncf)}, bridge={len(bridge)}"

    def test_eej_10hogs(self):
        """EEJ：每条染色体 10 HOGs。"""
        G = ColoredGraph(hog_level='test')
        G.add_child('c1', MockChild([
            linear_chrom('ABCDEFGHIJ', 'c1_0'),
            linear_chrom('KLMNOPQRST', 'c1_1'),
        ]))
        G.add_child('c2', MockChild([
            linear_chrom('ABCDEFGHIJKLMNOPQRST', 'c2_0'),
        ]))
        G.resolve_all_events()
        eej = [e for e in G.events if e.event_type in ('eej', 'ncf', 'bridge_unclassified')]
        assert len(eej) >= 1, f"10-HOG EEJ: {events_by_type(G)}"


# ═══════════════════════════════════════════════════════════════
#  NCF (Non-Collinear Fusion)
# ═══════════════════════════════════════════════════════════════

class TestNCF:
    """NCF：一条染色体的片段插入另一条中间。

    孩子 A: [A-B-C-D, E-F] → 孩子 B: [A-B-E-F-C-D]
    桥接检测：唯一边连接大组件和小组件。
    """

    def test_simple_ncf(self):
        """最简 NCF：6 HOGs。"""
        G = ColoredGraph(hog_level='test')
        G.add_child('c1', MockChild([
            linear_chrom('ABCD', 'c1_0'),
            linear_chrom('EF', 'c1_1'),
        ]))
        G.add_child('c2', MockChild([
            linear_chrom('ABEFCD', 'c2_0'),  # E-F 插入 B 和 C 之间
        ]))
        G.resolve_all_events()
        ncf = [e for e in G.events if e.event_type == 'ncf']
        eej = [e for e in G.events if e.event_type == 'eej']
        bridge = [e for e in G.events if e.event_type == 'bridge_unclassified']
        assert len(ncf) + len(eej) + len(bridge) >= 1, \
            f"NCF: ncf={len(ncf)}, eej={len(eej)}, bridge={len(bridge)}"


# ═══════════════════════════════════════════════════════════════
#  INDEL (gene insertion/deletion)
# ═══════════════════════════════════════════════════════════════

class TestIndel:
    """Indel：一个基因被插入或删除。

    孩子 A: A-B-C-D，孩子 B: A-C-D（B 被删除）
    """

    def test_simple_indel(self):
        """最简 indel：1 个基因插入/删除。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCD', 'c1')]),
            ('c2', [linear_chrom('ACD', 'c2')]),  # B 缺失
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, "Simple indel not detected"

    def test_indel_2genes(self):
        """2 个连续基因被删除。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ACE', 'c2')]),  # B-D 缺失
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, "2-gene indel not detected"

    def test_indel_at_start(self):
        """染色体开头的基因被删除。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCD', 'c1')]),
            ('c2', [linear_chrom('BCD', 'c2')]),  # A 缺失
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1, "Start indel not detected"


# ═══════════════════════════════════════════════════════════════
#  图完整性
# ═══════════════════════════════════════════════════════════════

class TestGraphIntegrity:
    """图的基本完整性检查。"""

    def test_no_events_linear(self):
        """两个完全相同的线性图 → 无事件。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ABCDEFGHIJ', 'c2')]),
        ])
        G.resolve_all_events()
        major = [e for e in G.events
                 if e.event_type in ('inversion', 'reciprocal_translocation',
                                     'eej', 'ncf', 'fission')]
        assert len(major) == 0, f"No events expected: {events_by_type(G)}"

    def test_path_cover_covers_all(self):
        """路径覆盖必须覆盖所有 HOG。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ABCDEFGHIJ', 'c2')]),
        ])
        G.resolve_all_events()
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        covered = set()
        for p in paths:
            covered.update(p)
        assert covered == G.all_hogs(), \
            f"Uncovered: {len(G.all_hogs() - covered)} HOGs"

    def test_chrom_count_no_events(self):
        """无事件时染色体数不变。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1'), linear_chrom('FGHIJ', 'c1b')]),
            ('c2', [linear_chrom('ABCDE', 'c2'), linear_chrom('FGHIJ', 'c2b')]),
        ])
        G.resolve_all_events()
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        assert len(paths) == 2, f"Expected 2 chroms, got {len(paths)}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
