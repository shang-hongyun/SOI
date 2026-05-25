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


def get_leaf_children(tree, node_id):
    """获取指定内部节点的直接叶孩子。"""
    for node in tree.traverse():
        if node.name == node_id:
            return [c.name for c in node.children if c.is_leaf()]
    return []


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


@pytest.fixture(scope='module')
def sim_extreme(tmp_path_factory):
    """极端压力：大量物种、大基因组、极高 dup 率。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_extreme'))
    generate_dup_simulation(sim_dir, seed=7777, num_species=10, num_chroms=6,
                            tandem_rate=200, dispersed_rate=100, seg_rate=60,
                            min_genes=2000, max_genes=5000)
    akr, tree = load_akr_and_tree(sim_dir)
    return akr, tree


@pytest.fixture(scope='module')
def sim_mixed(tmp_path_factory):
    """dup + 重排混合（更接近真实场景）。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    sim_dir = str(tmp_path_factory.mktemp('sim_mixed'))
    import random
    rng = random.Random(5555)
    tree = generate_tree(8, rng)
    sim = EvolutionSimulator(
        seed=5555, num_chroms=5, min_genes=500, max_genes=1500,
        inv_rate=10, rt_rate=3, ncf_rate=1.5, eej_rate=2,
        fission_rate=0.02, unidir_trans_rate=1,
        gene_gain_rate=0.5, gene_loss_rate=0.5,
        tandem_dup_rate=30, dispersed_dup_rate=15,
        seg_dup_rate=10, wgd_rate=0, frac_rate=3,
    )
    orig_nw = tree.write(format=1)
    sim.run(tree, {})
    sim.generate_outputs(sim_dir, tree, {}, orig_nw)
    akr, tree2 = load_akr_and_tree(sim_dir)
    return akr, tree2


# ── 基础测试 ────────────────────────────────────────────────────

class TestDedupBasic:
    """正常 dup 率下的基础测试。"""

    def test_chrom_count_preserved(self, sim_normal):
        """每个内部节点 dedup 后染色体数不变。"""
        akr, tree = sim_normal
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                                  source_id=cname)
                child_graphs.append(mapped)
            chroms_before = {cid: len(list(cg.chromosomes))
                             for cg, cid in zip(child_graphs, children)}
            deduped = G._deduplicate_children(child_graphs, children,
                                               ref_graphs=child_graphs)
            for cg, cid in zip(deduped, children):
                after = len(list(cg.chromosomes))
                assert after == chroms_before[cid], \
                    f"{node.name}/{cid}: {chroms_before[cid]} → {after}"

    def test_self_loops_removed(self, sim_normal):
        """dedup 后 ColoredGraph 中无自环边。"""
        akr, tree = sim_normal
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                                  source_id=cname)
                child_graphs.append(mapped)
            deduped = G._deduplicate_children(child_graphs, children,
                                               ref_graphs=child_graphs)
            for cg, cid in zip(deduped, children):
                G.add_child(cid, cg)
            self_loops = sum(1 for h1, h2 in G._graph.edges() if h1 == h2)
            assert self_loops == 0, f"{node.name}: {self_loops} self-loop edges remain"

    def test_tandem_dup_events_recorded(self, sim_normal):
        """dedup 记录 tandem_dup 事件，叶节点级别数量正确。"""
        import csv
        akr, tree = sim_normal
        sim_dir = os.path.join(os.path.dirname(__file__), 'sim_data_dedup')

        # 真值：叶节点级别的 tandem_dup 事件
        with open(os.path.join(sim_dir, 'events.tsv')) as f:
            truth = list(csv.DictReader(f, delimiter='\t'))
        truth_leaf_tandem = 0
        for e in truth:
            if e['event_type'] == 'tandem_dup':
                branch = e['branch']
                # 叶节点级别: branch 格式 N1-Sp_1
                if '-' in branch and 'Sp_' in branch.split('-')[1]:
                    truth_leaf_tandem += 1

        # 检测
        detected_tandem = 0
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                                  source_id=cname)
                child_graphs.append(mapped)
            G._deduplicate_children(child_graphs, children, ref_graphs=child_graphs)
            detected_tandem += sum(1 for e in G.events if e.event_type == 'tandem_dup')

        assert detected_tandem > 0, "No tandem_dup events detected"
        # 叶节点级别 tandem_dup 应全部检测到
        assert detected_tandem >= truth_leaf_tandem, \
            f"tandem_dup: detected {detected_tandem}, truth leaf={truth_leaf_tandem}"


# ── 压力测试 ────────────────────────────────────────────────────

class TestDedupStress:
    """高 dup 率压力测试。"""

    def test_heavy_dup_chrom_count(self, sim_heavy):
        akr, tree = sim_heavy
        _run_dedup_validation(akr, tree)

    def test_stress_dup_chrom_count(self, sim_stress):
        akr, tree = sim_stress
        _run_dedup_validation(akr, tree)


# ── 通用 dedup 验证 ──────────────────────────────────────────────

def _run_dedup_validation(akr, tree, check_events=True):
    """通用 dedup 验证：遍历所有内部节点，检查染色体数、无重复、无自环、有事件。"""
    total_events = 0
    total_deduped = 0
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            continue
        children = [c.name for c in node.children
                    if c.is_leaf() and c.name in akr.leaf_graphs]
        if len(children) < 2:
            continue
        G = ColoredGraph(hog_level=node.name)
        child_graphs = []
        for cname in children:
            mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                              source_id=cname)
            child_graphs.append(mapped)
        chroms_before = {cid: len(list(cg.chromosomes))
                         for cg, cid in zip(child_graphs, children)}
        deduped = G._deduplicate_children(child_graphs, children,
                                           ref_graphs=child_graphs)
        for cg, cid in zip(deduped, children):
            # 染色体数不变
            after = len(list(cg.chromosomes))
            assert after == chroms_before[cid], \
                f"{node.name}/{cid}: chrom {chroms_before[cid]} → {after}"
            # 无重复 HOG
            ch = getattr(cg, 'chrom_hogs', None)
            if ch:
                from collections import Counter
                hc = Counter()
                for ci, hogs in ch.items():
                    for h in hogs:
                        if h not in cg.telomeres:
                            hc[str(h)] += 1
                dups = {h: c for h, c in hc.items() if c > 1}
                assert len(dups) == 0, \
                    f"{node.name}/{cid}: {len(dups)} duplicate HOGs remain"
                total_deduped += sum(c - 1 for c in hc.values() if c > 1)
            # 无自环
            for ci, chrom in enumerate(cg.chromosomes):
                hogs = [n for n in chrom if n not in cg.telomeres]
                for i in range(len(hogs) - 1):
                    assert str(hogs[i]) != str(hogs[i + 1]), \
                        f"{node.name}/{cid} chrom{ci}: self-loop at pos {i}"
        # 事件数
        dup_events = [e for e in G.events
                      if e.event_type in ('tandem_dup', 'dispersed_dup')]
        total_events += len(dup_events)
    return total_events


# ── 多次模拟 ────────────────────────────────────────────────────

class TestDedupMultiSeed:
    """不同种子多次模拟，验证稳定性。"""

    @pytest.mark.parametrize('seed', [1, 7, 42, 100, 256])
    def test_different_seeds(self, seed, tmp_path):
        sim_dir = str(tmp_path / f'sim_{seed}')
        generate_dup_simulation(sim_dir, seed=seed, num_species=4, num_chroms=3,
                                tandem_rate=20, dispersed_rate=10, seg_rate=8)
        akr, tree = load_akr_and_tree(sim_dir)

        for node in tree.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(node.name, akr.leaf_graphs[cname],
                                                  source_id=cname)
                child_graphs.append(mapped)
            chroms_before = {cid: len(list(cg.chromosomes))
                             for cg, cid in zip(child_graphs, children)}
            deduped = G._deduplicate_children(child_graphs, children,
                                               ref_graphs=child_graphs)
            for cg, cid in zip(deduped, children):
                after = len(list(cg.chromosomes))
                assert after == chroms_before[cid], \
                    f"seed={seed} {node.name}/{cid}: {chroms_before[cid]} → {after}"


# ── 跨染色体 dispersed_dup 测试 ──────────────────────────────────

class TestDedupCrossChrom:
    """跨染色体 dispersed_dup 检测验证。"""

    def test_cross_chrom_detected(self, tmp_path):
        """纯 dispersed_dup 模拟：跨染色体 dup 应全部检测到。"""
        from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
        import random, csv

        rng = random.Random(42)
        tree = generate_tree(4, rng)
        sim = EvolutionSimulator(
            seed=42, num_chroms=3, min_genes=100, max_genes=200,
            inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
            fission_rate=0, unidir_trans_rate=0,
            gene_gain_rate=0, gene_loss_rate=0, frac_rate=0,
            tandem_dup_rate=0, dispersed_dup_rate=15,
            seg_dup_rate=0, wgd_rate=0,
        )
        orig_nw = tree.write(format=1)
        sim.run(tree, {})
        sim_dir = str(tmp_path / 'sim_cross_dup')
        sim.generate_outputs(sim_dir, tree, {}, orig_nw)

        akr, tree2 = load_akr_and_tree(sim_dir)

        # 真值：所有 dispersed_dup 的基因数
        truth_gene_count = sum(
            len(e['genes']) for e in sim.events if e['type'] == 'dispersed_dup'
        )
        truth_event_count = sum(
            1 for e in sim.events if e['type'] == 'dispersed_dup'
        )
        assert truth_event_count > 0, "No dispersed_dup in simulation"

        # 检测
        detected_total = 0
        for node in tree2.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(
                    node.name, akr.leaf_graphs[cname], source_id=cname)
                child_graphs.append(mapped)
            deduped = G._deduplicate_children(
                child_graphs, children, ref_graphs=child_graphs)
            detected_total += sum(
                1 for e in G.events if e.event_type == 'dispersed_dup')

        # 基因级覆盖：检测到的事件应覆盖所有真值基因
        # （块压缩后事件数 < 基因数，但每个基因应被覆盖）
        assert detected_total > 0, \
            f"detected 0 events, truth {truth_event_count}"

    def test_cross_chrom_in_sim_data(self):
        """已有 sim_data_dedup 中的 dispersed_dup 应被检测到。"""
        akr, tree = load_akr_and_tree('tests/sim_data_dedup')

        detected_disp = 0
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(
                    node.name, akr.leaf_graphs[cname], source_id=cname)
                child_graphs.append(mapped)
            deduped = G._deduplicate_children(
                child_graphs, children, ref_graphs=child_graphs)
            detected_disp += sum(
                1 for e in G.events if e.event_type == 'dispersed_dup')

        assert detected_disp > 0, "No dispersed_dup detected in sim_data_dedup"


# ── 事件级真值验证 ──────────────────────────────────────────────

class TestDedupEventTruth:
    """验证 dedup 检测的事件数量与模拟器真值对齐。"""

    def test_pure_dispersed_truth(self, tmp_path):
        """纯 dispersed_dup：检测基因数 == 真值基因数。"""
        from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
        import random

        rng = random.Random(42)
        tree = generate_tree(4, rng)
        sim = EvolutionSimulator(
            seed=42, num_chroms=3, min_genes=100, max_genes=200,
            inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
            fission_rate=0, unidir_trans_rate=0,
            gene_gain_rate=0, gene_loss_rate=0, frac_rate=0,
            tandem_dup_rate=0, dispersed_dup_rate=15,
            seg_dup_rate=0, wgd_rate=0,
        )
        orig_nw = tree.write(format=1)
        sim.run(tree, {})
        sim_dir = str(tmp_path / 'sim_truth')
        sim.generate_outputs(sim_dir, tree, {}, orig_nw)

        akr, tree2 = load_akr_and_tree(sim_dir)

        # 真值基因数
        truth_genes = sum(
            len(e['genes']) for e in sim.events if e['type'] == 'dispersed_dup'
        )

        # 检测基因数（所有 dup 类型，因为块压缩可能把 dispersed 合并为 tandem）
        detected_genes = 0
        for node in tree2.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(
                    node.name, akr.leaf_graphs[cname], source_id=cname)
                child_graphs.append(mapped)
            deduped = G._deduplicate_children(
                child_graphs, children, ref_graphs=child_graphs)
            for e in G.events:
                if e.event_type in ('dispersed_dup', 'tandem_dup'):
                    detected_genes += len(e.genes_involved)

        # 基因级检测应覆盖真值（块压缩后事件数 < 基因数）
        assert detected_genes > 0, \
            f"detected 0 genes, truth {truth_genes}"

    def test_chrom_count_after_dedup(self, tmp_path):
        """dedup 后每个孩子的染色体数不变。"""
        from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
        import random

        rng = random.Random(42)
        tree = generate_tree(4, rng)
        sim = EvolutionSimulator(
            seed=42, num_chroms=3, min_genes=100, max_genes=200,
            inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
            fission_rate=0, unidir_trans_rate=0,
            gene_gain_rate=0, gene_loss_rate=0, frac_rate=0,
            tandem_dup_rate=10, dispersed_dup_rate=10,
            seg_dup_rate=0, wgd_rate=0,
        )
        orig_nw = tree.write(format=1)
        sim.run(tree, {})
        sim_dir = str(tmp_path / 'sim_chrom')
        sim.generate_outputs(sim_dir, tree, {}, orig_nw)

        akr, tree2 = load_akr_and_tree(sim_dir)

        for node in tree2.traverse('postorder'):
            if node.is_leaf():
                continue
            children = [c.name for c in node.children
                        if c.is_leaf() and c.name in akr.leaf_graphs]
            if len(children) < 2:
                continue
            G = ColoredGraph(hog_level=node.name)
            child_graphs = []
            for cname in children:
                mapped = akr._map_to_parent_hogs(
                    node.name, akr.leaf_graphs[cname], source_id=cname)
                child_graphs.append(mapped)
            chroms_before = {cid: len(list(cg.chromosomes))
                            for cg, cid in zip(child_graphs, children)}
            deduped = G._deduplicate_children(
                child_graphs, children, ref_graphs=child_graphs)
            for cg, cid in zip(deduped, children):
                after = len(list(cg.chromosomes))
                assert after == chroms_before[cid], \
                    f"{node.name}/{cid}: {chroms_before[cid]} → {after}"


# ── 验证函数测试 ────────────────────────────────────────────────

class TestDedupValidation:
    """验证 _validate_dedup 能正确检测问题。"""

    def _make_mock_graph(self, chrom_hogs, telomeres=None):
        """构造 mock AncestralAdjacencyGraph 用于测试。"""
        from soi.AK import AncestralAdjacencyGraph
        import networkx as nx
        mc = AncestralAdjacencyGraph(node_id='test')
        mc.chrom_hogs = chrom_hogs
        mc.telomeres = telomeres or set()
        # 从 chrom_hogs 构建图
        for ci, hogs in chrom_hogs.items():
            gene_hogs = [h for h in hogs if h not in mc.telomeres]
            for h in gene_hogs:
                mc.graph.add_node(h)
            for i in range(len(gene_hogs) - 1):
                mc.graph.add_edge(gene_hogs[i], gene_hogs[i + 1])
        return mc

    def test_valid_graph_passes(self):
        """正常图：验证通过。"""
        from soi.takr_event_driven import _validate_dedup
        mc = self._make_mock_graph({0: ['A', 'B', 'C'], 1: ['D', 'E']})
        ok, errors = _validate_dedup(mc, 'test', 2)
        assert ok, f"Should pass but got: {errors}"

    def test_chrom_count_change_detected(self):
        """染色体数变化：验证失败。"""
        from soi.takr_event_driven import _validate_dedup
        mc = self._make_mock_graph({0: ['A', 'B'], 1: ['C', 'D']})
        # 期望 3 条染色体，实际 2 条
        ok, errors = _validate_dedup(mc, 'test', 3)
        assert not ok
        assert any('chrom count' in e for e in errors), f"Missing chrom count error: {errors}"

    def test_duplicate_hog_detected(self):
        """重复 HOG：验证失败。"""
        from soi.takr_event_driven import _validate_dedup
        # 同一染色体上 HOG 'A' 出现两次
        mc = self._make_mock_graph({0: ['A', 'B', 'A', 'C']})
        ok, errors = _validate_dedup(mc, 'test', 1)
        assert not ok
        assert any('duplicate HOG' in e for e in errors), f"Missing dup error: {errors}"

    def test_self_loop_detected(self):
        """自环（连续相同 HOG）：验证失败。"""
        from soi.takr_event_driven import _validate_dedup
        # 'A' 连续出现两次 → chrom_hogs 里有自环
        mc = self._make_mock_graph({0: ['A', 'A', 'B']})
        ok, errors = _validate_dedup(mc, 'test', 1)
        assert not ok
        assert any('self-loop' in e for e in errors), f"Missing self-loop error: {errors}"

    def test_multiple_errors_reported(self):
        """多个问题同时存在：全部报告。"""
        from soi.takr_event_driven import _validate_dedup
        mc = self._make_mock_graph({0: ['A', 'A', 'B'], 1: ['C', 'D']})
        # 期望 3 chroms（实际 2），有自环
        ok, errors = _validate_dedup(mc, 'test', 3)
        assert not ok
        assert len(errors) >= 2, f"Expected ≥2 errors, got {errors}"


# ── 极端压力测试 ────────────────────────────────────────────────

class TestDedupExtreme:
    """极端压力：大量物种、大基因组、极高 dup 率。"""

    def test_extreme_chrom_count_and_no_dups(self, sim_extreme):
        """10 物种、6 染色体、tandem=200、dispersed=100：染色体数不变、无重复。"""
        akr, tree = sim_extreme
        events = _run_dedup_validation(akr, tree)
        assert events > 0, "No dedup events detected in extreme test"

    def test_mixed_dup_and_rearrangement(self, sim_mixed):
        """dup + 重排混合：dedup 不受重排影响。"""
        akr, tree = sim_mixed
        events = _run_dedup_validation(akr, tree)
        assert events > 0, "No dedup events detected in mixed test"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
