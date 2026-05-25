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

def generate_dup_simulation(outdir, seed=42):
    """生成只含 dup 事件（无重排）的模拟数据。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    import random
    rng = random.Random(seed)
    tree = generate_tree(4, rng)
    sim = EvolutionSimulator(
        seed=seed, num_chroms=3,
        inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
        fission_rate=0, unidir_trans_rate=0,
        gene_gain_rate=0, gene_loss_rate=0, frac_rate=0,
        tandem_dup_rate=15, dispersed_dup_rate=8, seg_dup_rate=5,
        wgd_rate=0,
    )
    ploidy_map = {}
    orig_nw = tree.write(format=1)
    sim.run(tree, ploidy_map)
    sim.generate_outputs(outdir, tree, ploidy_map, orig_nw)
    return sim


@pytest.fixture(scope='module')
def sim_env(tmp_path_factory):
    """生成模拟数据，返回 (akr, tree, sim_dir)。"""
    sim_dir = str(tmp_path_factory.mktemp('sim'))
    sim = generate_dup_simulation(sim_dir, seed=42)

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
    return akr, tree, sim_dir


# ── 测试 ──────────────────────────────────────────────────────

class TestDedupSimulation:
    """用模拟数据验证 dedup：染色体数不变，dup 类型全覆盖。"""

    def test_chrom_count_preserved_no_ref(self, sim_env):
        """每个叶图 dedup（无参照）后染色体数不变。"""
        akr, tree, _ = sim_env
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            before = len(list(lg.chromosomes))
            deduped = G._deduplicate_single_child(lg, leaf.name)
            after = len(list(deduped.chromosomes))
            assert before == after, f"{leaf.name}: {before} → {after}"

    def test_chrom_count_preserved_with_ref(self, sim_env):
        """用兄弟孩子做参照，dedup 后染色体数不变。"""
        akr, tree, _ = sim_env
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
                deduped = G._deduplicate_single_child(cg, cname, ref_graphs=refs)
                after = len(list(deduped.chromosomes))
                assert before == after, f"{cname} (node {node.name}): {before} → {after}"

    def test_dup_events_exist(self, sim_env):
        """模拟数据中确实包含 dup 事件。"""
        _, _, sim_dir = sim_env
        import csv
        with open(os.path.join(sim_dir, 'events.tsv')) as f:
            events = list(csv.DictReader(f, delimiter='\t'))
        dup_types = {'tandem_dup', 'dispersed_dup', 'seg_duplication'}
        found = {e['event_type'] for e in events}
        assert dup_types & found, f"No dup events in simulation: {found}"

    def test_dedup_removes_duplicates(self, sim_env):
        """dedup 后每个 HOG 在每个染色体中只出现一次。"""
        akr, tree, _ = sim_env
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, chrom in enumerate(deduped.chromosomes):
                hogs = [n for n in chrom if n not in deduped.telomeres]
                assert len(hogs) == len(set(str(h) for h in hogs)), \
                    f"{leaf.name} chrom{ci}: duplicate HOGs remain after dedup"

    def test_dedup_preserves_gene_count_order(self, sim_env):
        """dedup 后每个染色体的基因数 ≤ 原始数（只减不增）。"""
        akr, tree, _ = sim_env
        G = ColoredGraph(hog_level='test')
        for leaf in tree.iter_leaves():
            if leaf.name not in akr.leaf_graphs:
                continue
            lg = akr.leaf_graphs[leaf.name]
            deduped = G._deduplicate_single_child(lg, leaf.name)
            for ci, (orig, dedup) in enumerate(
                    zip(lg.chromosomes, deduped.chromosomes)):
                orig_genes = [n for n in orig if n not in lg.telomeres]
                dedup_genes = [n for n in dedup if n not in deduped.telomeres]
                assert len(dedup_genes) <= len(orig_genes), \
                    f"{leaf.name} chrom{ci}: {len(orig_genes)} → {len(dedup_genes)}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
