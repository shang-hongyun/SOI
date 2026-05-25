#!/usr/bin/env python3
"""
真值驱动集成测试 — 用模拟数据验证每种事件类型。

要求：tests/benchmark_run/sim_data/ 存在（由 benchmark_rak.sh 生成）。
运行：python -m pytest tests/test_integration.py -v
"""
import sys, os, csv, pytest
from collections import defaultdict, Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import logging
logging.basicConfig(level=logging.WARNING)

SIM = os.path.join(os.path.dirname(__file__), 'benchmark_run', 'sim_data')

# ── 真值加载 ──────────────────────────────────────────────────

def load_truth_events(path):
    """{branch: [event_type, ...]}"""
    events = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            events[row['branch']].append(row['event_type'])
    return events


def load_truth_karyotypes(path):
    """{node_id: [[gene, ...], ...]}"""
    karyo, current = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current = line.lstrip('>').split('\t')[0]
                karyo[current] = []
            elif current and line:
                karyo[current].append(line.split('\t')[-1].split())
    return karyo


# ── 重建单个节点 ──────────────────────────────────────────────

def reconstruct_node(akr, tree, node_id):
    """重建指定内部节点，返回 (ColoredGraph, detected_events)"""
    from soi.takr_colored_graph import ColoredGraph

    node = None
    for n in tree.traverse():
        if n.name == node_id:
            node = n
            break
    assert node is not None, f"Node {node_id} not found in tree"

    G = ColoredGraph(hog_level=node_id)
    for child in node.children:
        cname = child.name
        if cname in akr.leaf_graphs:
            mapped = akr._map_to_parent_hogs(node_id, akr.leaf_graphs[cname], source_id=cname)
            G.add_child(cname, mapped)
        elif cname in akr.anc_graphs:
            G.add_child(cname, akr.anc_graphs[cname])

    G.resolve_all_events()
    return G


# ── Fixture: 加载 AKR + 模拟数据 ─────────────────────────────

@pytest.fixture(scope='module')
def akr_and_tree():
    """加载一次 AKR 和树，所有测试复用。"""
    from soi.AK import AKR
    from soi.evolution_simulator_ak import parse_tree

    if not os.path.exists(os.path.join(SIM, 'events.tsv')):
        pytest.skip('模拟数据不存在，先运行 bash tests/benchmark_rak.sh')

    outdir = os.path.join(os.path.dirname(__file__), 'benchmark_run', 'recon_test')
    os.makedirs(outdir, exist_ok=True)

    akr = AKR(
        ogfile=os.path.join(SIM, 'ortholog_groups.txt'),
        orthfiles=[os.path.join(SIM, 'ortholog_pairs.txt')],
        gfffile=os.path.join(SIM, 'all_species_gene.gff'),
        sptreefile=os.path.join(SIM, 'species_tree.nwk'),
        outpre=os.path.join(outdir, 'AKR'),
        reconstruction_algorithm='v4_colored',
        min_genes=0, timeout=600,
    )
    akr._build_hogs()
    akr._build_leaf_graphs()
    tree, _, _ = parse_tree(os.path.join(SIM, 'species_tree.nwk'))
    return akr, tree


@pytest.fixture(scope='module')
def truth():
    return {
        'events': load_truth_events(os.path.join(SIM, 'events.tsv')),
        'karyotypes': load_truth_karyotypes(os.path.join(SIM, 'ancestors_karyotypes.txt')),
    }


# ── 测试：染色体数 ───────────────────────────────────────────

class TestChromosomeCounts:
    """每个内部节点的重建染色体数必须等于真值。"""

    def test_N1_chrom_count(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        truth_n = len(truth['karyotypes']['N1'])
        assert len(paths) == truth_n, \
            f"N1: {len(paths)} chroms (truth {truth_n})"

    def test_N2_chrom_count(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N2')
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        truth_n = len(truth['karyotypes']['N2'])
        assert len(paths) == truth_n, \
            f"N2: {len(paths)} chroms (truth {truth_n})"


# ── 测试：RT 检测 ─────────────────────────────────────────────

class TestReciprocalTranslocation:
    """真值有 4 个 RT（N0-N1:1, N1-Sp_1:2, N2-Sp_2:1）。"""

    def test_N1_detects_rt(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        rt = [e for e in G.events if e.event_type == 'reciprocal_translocation']
        truth_rt = truth['events'].get('N1-Sp_1', []).count('reciprocal_translocation') + \
                   truth['events'].get('N1-Sp_4', []).count('reciprocal_translocation')
        assert len(rt) >= truth_rt, \
            f"N1 RT: detected {len(rt)}, truth {truth_rt}"

    def test_N2_detects_rt(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N2')
        rt = [e for e in G.events if e.event_type == 'reciprocal_translocation']
        truth_rt = truth['events'].get('N2-Sp_2', []).count('reciprocal_translocation') + \
                   truth['events'].get('N2-Sp_3', []).count('reciprocal_translocation')
        assert len(rt) >= truth_rt, \
            f"N2 RT: detected {len(rt)}, truth {truth_rt}"


# ── 测试：inversion 检测 ──────────────────────────────────────

class TestInversion:
    """真值有 3 个 inversion（N1-Sp_1:1, N1-Sp_4:1, N2-Sp_3:1）。"""

    def test_N1_direction_conflicts_exist(self, akr_and_tree):
        """N1 图中必须有方向冲突边（inversion 的前提条件）。"""
        akr, tree = akr_and_tree
        from soi.takr_colored_graph import ColoredGraph
        G = ColoredGraph(hog_level='N1')
        for cname in ['Sp_1', 'Sp_4']:
            mapped = akr._map_to_parent_hogs('N1', akr.leaf_graphs[cname], source_id=cname)
            G.add_child(cname, mapped)
        n_conflict = sum(1 for h1, h2 in G.shared_edges()
                         if G.edge_has_direction_conflict(h1, h2))
        assert n_conflict > 0, "N1 should have direction conflict edges"

    def test_N1_detects_inversion(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        truth_inv = truth['events'].get('N1-Sp_1', []).count('inversion') + \
                    truth['events'].get('N1-Sp_4', []).count('inversion')
        assert len(inv) >= truth_inv, \
            f"N1 inversion: detected {len(inv)}, truth {truth_inv}"

    def test_N2_detects_inversion(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N2')
        inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
        truth_inv = truth['events'].get('N2-Sp_2', []).count('inversion') + \
                    truth['events'].get('N2-Sp_3', []).count('inversion')
        assert len(inv) >= truth_inv, \
            f"N2 inversion: detected {len(inv)}, truth {truth_inv}"


# ── 测试：EEJ 检测 ────────────────────────────────────────────

class TestEEJ:
    """真值有 1 个 eej（N1-Sp_4:1）。"""

    def test_N1_detects_eej(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        eej = [e for e in G.events if e.event_type == 'eej']
        truth_eej = truth['events'].get('N1-Sp_1', []).count('eej') + \
                    truth['events'].get('N1-Sp_4', []).count('eej')
        assert len(eej) >= truth_eej, \
            f"N1 eej: detected {len(eej)}, truth {truth_eej}"

    def test_eej_not_over_detected(self, akr_and_tree, truth):
        """eej 不应过度检测（允许 2x 容差）。"""
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        eej = [e for e in G.events if e.event_type == 'eej']
        truth_eej = truth['events'].get('N1-Sp_1', []).count('eej') + \
                    truth['events'].get('N1-Sp_4', []).count('eej')
        assert len(eej) <= truth_eej * 2 + 1, \
            f"N1 eej over-detected: {len(eej)} vs truth {truth_eej}"


# ── 测试：NCF 检测 ────────────────────────────────────────────

class TestNCF:
    """真值有 3 个 ncf（N1-Sp_4:1, N2-Sp_2:1, N2-Sp_3:1）。"""

    def test_detects_ncf(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        total_ncf = 0
        total_truth = 0
        for node_id in ['N1', 'N2']:
            G = reconstruct_node(akr, tree, node_id)
            ncf = [e for e in G.events if e.event_type == 'ncf']
            total_ncf += len(ncf)
            for branch in [f'{node_id}-Sp_1', f'{node_id}-Sp_2',
                           f'{node_id}-Sp_3', f'{node_id}-Sp_4']:
                total_truth += truth['events'].get(branch, []).count('ncf')
        assert total_ncf >= total_truth, \
            f"NCF: detected {total_ncf}, truth {total_truth}"


# ── 测试：unidir_trans 检测 ───────────────────────────────────

class TestUnidirTrans:
    """真值有 4 个 unidir_trans（N1-Sp_4:2, N2-Sp_2:1, N2-Sp_3:1）。"""

    def test_detects_unidir_trans(self, akr_and_tree, truth):
        akr, tree = akr_and_tree
        total_ut = 0
        total_truth = 0
        for node_id in ['N1', 'N2']:
            G = reconstruct_node(akr, tree, node_id)
            ut = [e for e in G.events if e.event_type == 'unidir_trans']
            total_ut += len(ut)
            for branch in [f'{node_id}-Sp_1', f'{node_id}-Sp_2',
                           f'{node_id}-Sp_3', f'{node_id}-Sp_4']:
                total_truth += truth['events'].get(branch, []).count('unidir_trans')
        assert total_ut >= total_truth, \
            f"unidir_trans: detected {total_ut}, truth {total_truth}"


# ── 测试：图结构完整性 ────────────────────────────────────────

class TestGraphIntegrity:
    """重建图的基本完整性检查。"""

    def test_N1_no_uncovered_hogs(self, akr_and_tree):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        covered = set()
        for p in paths:
            covered.update(p)
        hog_set = G.all_hogs()
        uncovered = hog_set - covered
        assert len(uncovered) == 0, f"N1: {len(uncovered)} uncovered HOGs"

    def test_N2_no_uncovered_hogs(self, akr_and_tree):
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N2')
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        covered = set()
        for p in paths:
            covered.update(p)
        hog_set = G.all_hogs()
        uncovered = hog_set - covered
        assert len(uncovered) == 0, f"N2: {len(uncovered)} uncovered HOGs"

    def test_N1_path_sizes_reasonable(self, akr_and_tree, truth):
        """路径大小应与真值染色体大小在同数量级。"""
        akr, tree = akr_and_tree
        G = reconstruct_node(akr, tree, 'N1')
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        truth_sizes = sorted([len(c) for c in truth['karyotypes']['N1']], reverse=True)
        recon_sizes = sorted([len(p) for p in paths], reverse=True)
        # 最大路径应接近最大真值染色体
        assert recon_sizes[0] >= truth_sizes[0] * 0.8, \
            f"Max path {recon_sizes[0]} too small (truth max {truth_sizes[0]})"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
