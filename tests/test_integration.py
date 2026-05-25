#!/usr/bin/env python3
"""
真值驱动集成测试 — 用模拟数据验证每种事件类型。

要求：tests/sim_data_v2/ 存在。
运行：python -m pytest tests/test_integration.py -v
"""
import sys, os, csv, pytest
from collections import defaultdict, Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import logging
logging.basicConfig(level=logging.WARNING)

SIM = os.path.join(os.path.dirname(__file__), 'sim_data_v2')


def load_truth_events(path):
    events = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            events[row['branch']].append(row['event_type'])
    return events


def load_truth_karyotypes(path):
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


def get_leaf_adjacent_nodes(tree):
    """返回所有孩子全是叶节点的内部节点。"""
    nodes = []
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            continue
        if all(c.is_leaf() for c in node.children):
            nodes.append(node.name)
    return nodes


def reconstruct_node(akr, tree, node_id):
    from soi.takr_colored_graph import ColoredGraph
    node = None
    for n in tree.traverse():
        if n.name == node_id:
            node = n
            break
    G = ColoredGraph(hog_level=node_id)
    for child in node.children:
        cname = child.name
        if cname in akr.leaf_graphs:
            mapped = akr._map_to_parent_hogs(node_id, akr.leaf_graphs[cname], source_id=cname)
            G.add_child(cname, mapped)
    G.resolve_all_events()
    return G


@pytest.fixture(scope='module')
def akr_and_tree():
    from soi.AK import AKR
    from soi.evolution_simulator_ak import parse_tree
    if not os.path.exists(os.path.join(SIM, 'events.tsv')):
        pytest.skip('模拟数据不存在')
    outdir = os.path.join(os.path.dirname(__file__), 'benchmark_run', 'recon_test_v2')
    os.makedirs(outdir, exist_ok=True)
    akr = AKR(
        ogfile=os.path.join(SIM, 'ortholog_groups.txt'),
        orthfiles=[os.path.join(SIM, 'ortholog_pairs.txt')],
        gfffile=os.path.join(SIM, 'all_species_gene.gff'),
        sptreefile=os.path.join(SIM, 'species_tree.nwk'),
        outpre=os.path.join(outdir, 'AKR'),
        reconstruction_algorithm='v4_colored', min_genes=0, timeout=600,
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


@pytest.fixture(scope='module')
def leaf_nodes(akr_and_tree):
    _, tree = akr_and_tree
    return get_leaf_adjacent_nodes(tree)


# ── 测试：染色体数 ───────────────────────────────────────────

class TestChromosomeCounts:
    def test_chrom_counts(self, akr_and_tree, truth, leaf_nodes):
        """所有叶相邻节点的染色体数必须等于真值。"""
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            if node_id not in truth['karyotypes']:
                continue
            G = reconstruct_node(akr, tree, node_id)
            paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
            truth_n = len(truth['karyotypes'][node_id])
            assert len(paths) == truth_n, \
                f"{node_id}: {len(paths)} chroms (truth {truth_n})"


# ── 测试：RT 检测 ─────────────────────────────────────────────

class TestReciprocalTranslocation:
    def test_rt_detection(self, akr_and_tree, truth, leaf_nodes):
        """每个节点检测到的 RT 数 >= 真值。"""
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            truth_rt = 0
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                truth_rt += truth['events'].get(f'{node_id}-{child}', []).count('reciprocal_translocation')
            if truth_rt == 0:
                continue
            G = reconstruct_node(akr, tree, node_id)
            rt = [e for e in G.events if e.event_type == 'reciprocal_translocation']
            assert len(rt) >= truth_rt, \
                f"{node_id} RT: detected {len(rt)}, truth {truth_rt}"


# ── 测试：inversion 检测 ──────────────────────────────────────

class TestInversion:
    def test_direction_conflicts_exist(self, akr_and_tree, leaf_nodes, truth):
        """有 inversion 的节点必须有方向冲突边。"""
        from soi.takr_colored_graph import ColoredGraph
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            truth_inv = 0
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                truth_inv += truth['events'].get(f'{node_id}-{child}', []).count('inversion')
            if truth_inv == 0:
                continue
            G = ColoredGraph(hog_level=node_id)
            for child in tree.search(node_id).children:
                if child.is_leaf() and child.name in akr.leaf_graphs:
                    mapped = akr._map_to_parent_hogs(node_id, akr.leaf_graphs[child.name], source_id=child.name)
                    G.add_child(child.name, mapped)
            n_conflict = sum(1 for h1, h2 in G.shared_edges()
                             if G.edge_has_direction_conflict(h1, h2))
            assert n_conflict > 0, f"{node_id}: should have direction conflicts (truth {truth_inv} inversions)"

    def test_inversion_detection(self, akr_and_tree, truth, leaf_nodes):
        """每个节点检测到的 inversion 数 >= 真值。"""
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            truth_inv = 0
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                truth_inv += truth['events'].get(f'{node_id}-{child}', []).count('inversion')
            if truth_inv == 0:
                continue
            G = reconstruct_node(akr, tree, node_id)
            inv = [e for e in G.events if e.event_type in ('inversion', 'telomere_inversion')]
            assert len(inv) >= truth_inv, \
                f"{node_id} inversion: detected {len(inv)}, truth {truth_inv}"


# ── 测试：EEJ 检测 ────────────────────────────────────────────

class TestEEJ:
    def test_eej_detection(self, akr_and_tree, truth, leaf_nodes):
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            truth_eej = 0
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                truth_eej += truth['events'].get(f'{node_id}-{child}', []).count('eej')
            if truth_eej == 0:
                continue
            G = reconstruct_node(akr, tree, node_id)
            eej = [e for e in G.events if e.event_type == 'eej']
            assert len(eej) >= truth_eej, \
                f"{node_id} eej: detected {len(eej)}, truth {truth_eej}"


# ── 测试：NCF 检测 ────────────────────────────────────────────

class TestNCF:
    def test_ncf_detection(self, akr_and_tree, truth, leaf_nodes):
        akr, tree = akr_and_tree
        total_ncf, total_truth = 0, 0
        for node_id in leaf_nodes:
            G = reconstruct_node(akr, tree, node_id)
            total_ncf += len([e for e in G.events if e.event_type == 'ncf'])
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                total_truth += truth['events'].get(f'{node_id}-{child}', []).count('ncf')
        assert total_ncf >= total_truth, f"NCF: detected {total_ncf}, truth {total_truth}"


# ── 测试：unidir_trans 检测 ───────────────────────────────────

class TestUnidirTrans:
    def test_unidir_trans_detection(self, akr_and_tree, truth, leaf_nodes):
        akr, tree = akr_and_tree
        total_ut, total_truth = 0, 0
        for node_id in leaf_nodes:
            G = reconstruct_node(akr, tree, node_id)
            total_ut += len([e for e in G.events if e.event_type == 'unidir_trans'])
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                total_truth += truth['events'].get(f'{node_id}-{child}', []).count('unidir_trans')
        assert total_ut >= total_truth, f"unidir_trans: detected {total_ut}, truth {total_truth}"


# ── 测试：fission 检测 ────────────────────────────────────────

class TestFission:
    def test_fission_detection(self, akr_and_tree, truth, leaf_nodes):
        akr, tree = akr_and_tree
        total_fis, total_truth = 0, 0
        for node_id in leaf_nodes:
            G = reconstruct_node(akr, tree, node_id)
            total_fis += len([e for e in G.events if e.event_type == 'fission'])
            for child in [n.name for n in tree.search(node_id).children if n.is_leaf()]:
                total_truth += truth['events'].get(f'{node_id}-{child}', []).count('fission')
        assert total_fis >= total_truth, f"fission: detected {total_fis}, truth {total_truth}"


# ── 测试：图完整性 ────────────────────────────────────────────

class TestGraphIntegrity:
    def test_no_uncovered_hogs(self, akr_and_tree, leaf_nodes):
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            G = reconstruct_node(akr, tree, node_id)
            paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
            covered = set()
            for p in paths:
                covered.update(p)
            uncovered = G.all_hogs() - covered
            assert len(uncovered) == 0, f"{node_id}: {len(uncovered)} uncovered HOGs"

    def test_path_sizes_reasonable(self, akr_and_tree, truth, leaf_nodes):
        akr, tree = akr_and_tree
        for node_id in leaf_nodes:
            if node_id not in truth['karyotypes']:
                continue
            G = reconstruct_node(akr, tree, node_id)
            paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
            truth_sizes = sorted([len(c) for c in truth['karyotypes'][node_id]], reverse=True)
            recon_sizes = sorted([len(p) for p in paths], reverse=True)
            if recon_sizes and truth_sizes:
                assert recon_sizes[0] >= truth_sizes[0] * 0.5, \
                    f"{node_id}: max path {recon_sizes[0]} too small (truth {truth_sizes[0]})"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
