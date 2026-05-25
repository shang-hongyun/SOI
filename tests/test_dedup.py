#!/usr/bin/env python3
"""
Phase 1 dedup 测试 — 自动生成含各类 dup 的模拟数据，验证 dedup 正确性。

运行：python -m pytest tests/test_dedup.py -v
"""
import sys, os, shutil, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import logging
logging.basicConfig(level=logging.WARNING)

from soi.takr_colored_graph import ColoredGraph
from soi.AK import AKR
from soi.evolution_simulator_ak import parse_tree


# ── 模拟数据生成 ──────────────────────────────────────────────

def generate_dup_simulation(outdir, seed=42, num_species=4, num_chroms=3,
                            tandem_rate=15, dispersed_rate=8, seg_rate=5,
                            min_genes=200, max_genes=500):
    """生成只含 dup 事件（无重排）的模拟数据。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    import random
    rng = random.Random(seed)
    tree = generate_tree(num_species, rng)
    sim = EvolutionSimulator(
        seed=seed, num_chroms=num_chroms, min_genes=min_genes, max_genes=max_genes,
        inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
        fission_rate=0, unidir_trans_rate=0,
        gene_gain_rate=0, gene_loss_rate=0, frac_rate=0,
        tandem_dup_rate=tandem_rate, dispersed_dup_rate=dispersed_rate,
        seg_dup_rate=seg_rate, wgd_rate=0,
    )
    ploidy_map = {}
    orig_nw = tree.write(format=1)
    sim.run(tree, ploidy_map)
    sim.generate_outputs(outdir, tree, ploidy_map, orig_nw)
    return sim


def load_akr_and_tree(sim_dir):
    """从模拟目录加载 AKR 和树。"""
    outdir = os.path.join(sim_dir, 'recon')
    os.makedirs(outdir, exist_ok=True)
    akr = AKR(
        ogfile=os.path.join(sim_dir, 'ortholog_groups.txt'),
        orthfiles=[os.path.join(sim_dir, 'ortholog_pairs.txt')],
        gfffile=os.path.join(sim_dir, 'all_species_gene.gff'),
        sptreefile=os.path.join(sim_dir, 'species_tree.nwk'),
        outpre=os.path.join(outdir, 'AKR'),
        reconstruction_algorithm='v4_colored', min_genes=0, timeout=600,
    )
    akr._build_hogs()
    akr._build_leaf_graphs()
    tree, _, _ = parse_tree(os.path.join(sim_dir, 'species_tree.nwk'))
    return akr, tree


def count_dups_per_leaf(akr, tree):
    """统计每个叶图中的重复 HOG 数。"""
    stats = {}
    for leaf in tree.iter_leaves():
        if leaf.name not in akr.leaf_graphs:
            continue
        lg = akr.leaf_graphs[leaf.name]
        all_hogs = []
        for c in lg.chromosomes:
            for n in c:
                if n not in lg.telomeres:
                    all_hogs.append(str(n))
        from collections import Counter
        hog_counts = Counter(all_hogs)
        n_dups = sum(c - 1 for c in hog_counts.values() if c > 1)
        n_dup_genes = sum(1 for c in hog_counts.values() if c > 1)
        stats[leaf.name] = {'n_dups': n_dups, 'n_dup_genes': n_dup_genes,
                            'total_hogs': len(all_hogs)}
    return stats


# ── Fixture ────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def sim_normal(tmp_path_factory):
    """正常 dup 率模拟。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_normal'))
    generate_dup_simulation(sim_dir, seed=42)
    akr, tree = load_akr_and_tree(sim_dir)
    return akr, tree


@pytest.fixture(scope='module')
def sim_heavy(tmp_path_factory):
    """高 dup 率模拟（压力测试）。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_heavy'))
    generate_dup_simulation(sim_dir, seed=123, num_species=6, num_chroms=4,
                            tandem_rate=40, dispersed_rate=20, seg_rate=15,
                            min_genes=300, max_genes=800)
    akr, tree = load_akr_and_tree(sim_dir)
    return akr, tree


@pytest.fixture(scope='module')
def sim_stress(tmp_path_factory):
    """极高 dup 率 + 多物种（极端压力）。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_stress'))
    generate_dup_simulation(sim_dir, seed=999, num_species=8, num_chroms=5,
                            tandem_rate=80, dispersed_rate=40, seg_rate=30,
                            min_genes=500, max_genes=1500)
    akr, tree = load_akr_and_tree(sim_dir)
    return akr, tree


# ── 基础测试 ────────────────────────────────────────────────────

class TestDedupBasic:
    """正常 dup 率下的基础测试。"""

    def test_chrom_count_no_ref(self, sim_normal):
        akr, tree = sim_normal
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            before = len(list(lg.chromosomes))
            after = len(list(G._deduplicate_single_child(lg, leaf.name).chromosomes))
            assert before == after, f"{leaf.name}: {before} → {after}"

    def test_chrom_count_with_ref(self, sim_normal):
        akr, tree = sim_normal
        G = ColoredGraph(hog_level='test')
        for node in tree.traverse('postorder'):
            if node.is_leaf() or len(node.children) < 2:
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            graphs = [akr.leaf_graphs[c] for c in children]
            for i, (cname, cg) in enumerate(zip(children, graphs)):
                before = len(list(cg.chromosomes))
                refs = [graphs[j] for j in range(len(graphs)) if j != i]
                after = len(list(G._deduplicate_single_child(cg, cname, refs).chromosomes))
                assert before == after, f"{cname}: {before} → {after}"

    def test_no_duplicates_after_dedup(self, sim_normal):
        akr, tree = sim_normal
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, chrom in enumerate(deduped.chromosomes):
                hogs = [str(n) for n in chrom if n not in deduped.telomeres]
                assert len(hogs) == len(set(hogs)), \
                    f"{leaf.name} chrom{ci}: duplicates remain"

    def test_gene_count_monotone(self, sim_normal):
        akr, tree = sim_normal
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, (orig, ded) in enumerate(
                    zip(lg.chromosomes, deduped.chromosomes)):
                n_orig = sum(1 for n in orig if n not in lg.telomeres)
                n_ded = sum(1 for n in ded if n not in deduped.telomeres)
                assert n_ded <= n_orig, f"{leaf.name} chrom{ci}: {n_orig} → {n_ded}"


# ── 压力测试 ────────────────────────────────────────────────────

class TestDedupStress:
    """高 dup 率压力测试。"""

    def test_heavy_dup_chrom_count(self, sim_heavy):
        """高 dup 率：染色体数不变。"""
        akr, tree = sim_heavy
        stats = count_dups_per_leaf(akr, tree)
        total_dups = sum(s['n_dups'] for s in stats.values())
        print(f"\n  heavy: {total_dups} total dup copies across {len(stats)} leaves")

        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            before = len(list(lg.chromosomes))
            after = len(list(G._deduplicate_single_child(lg, leaf.name).chromosomes))
            assert before == after, f"{leaf.name}: {before} → {after}"

    def test_heavy_dup_no_remaining(self, sim_heavy):
        """高 dup 率：dedup 后无重复。"""
        akr, tree = sim_heavy
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, chrom in enumerate(deduped.chromosomes):
                hogs = [str(n) for n in chrom if n not in deduped.telomeres]
                assert len(hogs) == len(set(hogs)), \
                    f"{leaf.name} chrom{ci}: duplicates remain after heavy dedup"

    def test_stress_dup_chrom_count(self, sim_stress):
        """极高 dup 率 + 8 物种：染色体数不变。"""
        akr, tree = sim_stress
        stats = count_dups_per_leaf(akr, tree)
        total_dups = sum(s['n_dups'] for s in stats.values())
        print(f"\n  stress: {total_dups} total dup copies across {len(stats)} leaves")

        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            before = len(list(lg.chromosomes))
            after = len(list(G._deduplicate_single_child(lg, leaf.name).chromosomes))
            assert before == after, f"{leaf.name}: {before} → {after}"

    def test_stress_no_remaining(self, sim_stress):
        """极高 dup 率：dedup 后无重复。"""
        akr, tree = sim_stress
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, chrom in enumerate(deduped.chromosomes):
                hogs = [str(n) for n in chrom if n not in deduped.telomeres]
                assert len(hogs) == len(set(hogs)), \
                    f"{leaf.name} chrom{ci}: duplicates remain after stress dedup"

    def test_stress_with_ref(self, sim_stress):
        """极高 dup 率 + 参照图：染色体数不变。"""
        akr, tree = sim_stress
        G = ColoredGraph(hog_level='test')
        for node in tree.traverse('postorder'):
            if node.is_leaf() or len(node.children) < 2:
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            graphs = [akr.leaf_graphs[c] for c in children]
            for i, (cname, cg) in enumerate(zip(children, graphs)):
                before = len(list(cg.chromosomes))
                refs = [graphs[j] for j in range(len(graphs)) if j != i]
                after = len(list(G._deduplicate_single_child(cg, cname, refs).chromosomes))
                assert before == after, f"{cname}: {before} → {after}"


# ── 多次模拟 ────────────────────────────────────────────────────

class TestDedupMultiSeed:
    """不同种子多次模拟，验证稳定性。"""

    @pytest.mark.parametrize('seed', [1, 7, 42, 100, 256])
    def test_different_seeds(self, seed, tmp_path):
        """不同种子的模拟数据，dedup 后染色体数不变。"""
        sim_dir = str(tmp_path / f'sim_{seed}')
        generate_dup_simulation(sim_dir, seed=seed, num_species=4, num_chroms=3,
                                tandem_rate=20, dispersed_rate=10, seg_rate=8)
        akr, tree = load_akr_and_tree(sim_dir)

        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            before = len(list(lg.chromosomes))
            after = len(list(G._deduplicate_single_child(lg, leaf.name).chromosomes))
            assert before == after, f"seed={seed} {leaf.name}: {before} → {after}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
