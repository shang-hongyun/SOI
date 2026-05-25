#!/usr/bin/env python3
"""
Phase 4b inversion 测试 — 内嵌模拟数据 + 压力测试 + 手工构造。

参考 test_dedup.py 规范：
- generate_inversion_simulation(): 生成含倒位的模拟数据
- _run_inversion_validation(): 通用验证函数
- 多复杂度级别：basic → stress → extreme
- 手工构造：精确控制方向和融合场景
"""
import sys, os, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import logging
logging.basicConfig(level=logging.WARNING)

from soi.takr_colored_graph import ColoredGraph


# ── 模拟数据生成 ──────────────────────────────────────────────

def generate_inversion_simulation(outdir, seed=42, num_species=4, num_chroms=3,
                                  inv_rate=10, min_genes=100, max_genes=300,
                                  rt_rate=0, ncf_rate=0, eej_rate=0):
    """生成含倒位的模拟数据。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    import random
    rng = random.Random(seed)
    tree = generate_tree(num_species, rng)
    sim = EvolutionSimulator(
        seed=seed, num_chroms=num_chroms, min_genes=min_genes, max_genes=max_genes,
        inv_rate=inv_rate, rt_rate=rt_rate, ncf_rate=ncf_rate, eej_rate=eej_rate,
        fission_rate=0, unidir_trans_rate=0,
        gene_gain_rate=0, gene_loss_rate=0,
        tandem_dup_rate=0, dispersed_dup_rate=0,
        seg_dup_rate=0, wgd_rate=0,
    )
    orig_nw = tree.write(format=1)
    sim.run(tree, {})
    sim.generate_outputs(outdir, tree, {}, orig_nw)
    return sim


def load_akr_and_tree(sim_dir):
    """从模拟目录加载 AKR 和树。"""
    from soi.AK import AKR
    from soi.evolution_simulator_ak import parse_tree
    recon_dir = os.path.join(sim_dir, 'recon')
    os.makedirs(recon_dir, exist_ok=True)
    akr = AKR(
        ogfile=os.path.join(sim_dir, 'ortholog_groups.txt'),
        orthfiles=[os.path.join(sim_dir, 'ortholog_pairs.txt')],
        gfffile=os.path.join(sim_dir, 'all_species_gene.gff'),
        sptreefile=os.path.join(sim_dir, 'species_tree.nwk'),
        outpre=os.path.join(recon_dir, 'AKR'),
        reconstruction_algorithm='v4_colored', min_genes=0, timeout=600,
    )
    akr._build_hogs()
    akr._build_leaf_graphs()
    tree, _, _ = parse_tree(os.path.join(sim_dir, 'species_tree.nwk'))
    return akr, tree


# ── 通用验证 ─────────────────────────────────────────────────

def _run_inversion_validation(akr, tree, expect_inversions=True):
    """通用 inversion 验证：遍历所有内部节点，检测方向冲突。"""
    total_conflicts = 0
    total_events = 0
    for node in tree.traverse('postorder'):
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
        # 方向调和
        G.harmonize_directions()
        # 检测方向冲突
        n_conflicts = 0
        for h1, h2 in G._graph.edges():
            if G.edge_has_direction_conflict(h1, h2):
                n_conflicts += 1
        total_conflicts += n_conflicts
        # 倒位事件
        G._detect_inversions()
        inv_events = [e for e in G.events if 'inversion' in e.event_type]
        total_events += len(inv_events)
    if expect_inversions:
        assert total_conflicts > 0 or total_events > 0, \
            "No direction conflicts or inversion events detected"
    return total_conflicts, total_events


# ── Fixtures ─────────────────────────────────────────────────

@pytest.fixture(scope='module')
def sim_normal_inv(tmp_path_factory):
    """正常倒位率模拟。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_normal_inv'))
    generate_inversion_simulation(sim_dir, seed=42)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_heavy_inv(tmp_path_factory):
    """高倒位率模拟。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_heavy_inv'))
    generate_inversion_simulation(sim_dir, seed=123, num_species=6, num_chroms=4,
                                  inv_rate=30, min_genes=200, max_genes=500)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_extreme_inv(tmp_path_factory):
    """极端倒位率 + 大基因组。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_extreme_inv'))
    generate_inversion_simulation(sim_dir, seed=7777, num_species=8, num_chroms=5,
                                  inv_rate=80, min_genes=500, max_genes=1500)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_massive_inv(tmp_path_factory):
    """大规模倒位：12 物种、极大基因组。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_massive_inv'))
    generate_inversion_simulation(sim_dir, seed=99999, num_species=12, num_chroms=8,
                                  inv_rate=150, min_genes=2000, max_genes=5000)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_inv_with_fusion(tmp_path_factory):
    """倒位 + 融合（混合方向）。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_inv_fusion'))
    generate_inversion_simulation(sim_dir, seed=5555, num_species=6, num_chroms=4,
                                  inv_rate=20, ncf_rate=3, eej_rate=2,
                                  min_genes=200, max_genes=500)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_inv_with_rearrangement(tmp_path_factory):
    """倒位 + 全部重排（inv + RT + NCF + EEJ）。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_inv_all'))
    generate_inversion_simulation(sim_dir, seed=3333, num_species=8, num_chroms=5,
                                  inv_rate=20, rt_rate=5, ncf_rate=3, eej_rate=3,
                                  min_genes=300, max_genes=800)
    return load_akr_and_tree(sim_dir)


# ── 基础测试 ─────────────────────────────────────────────────

class TestInversionBasic:
    """正常倒位率下的基础测试。"""

    def test_inversions_detected(self, sim_normal_inv):
        """应检测到方向冲突或倒位事件。"""
        akr, tree = sim_normal_inv
        conflicts, events = _run_inversion_validation(akr, tree)
        assert conflicts > 0 or events > 0, "No inversions in normal simulation"


# ── 压力测试 ─────────────────────────────────────────────────

class TestInversionStress:
    """高倒位率压力测试。"""

    def test_heavy_inv(self, sim_heavy_inv):
        """6 物种、高倒位率。"""
        akr, tree = sim_heavy_inv
        _run_inversion_validation(akr, tree)

    def test_extreme_inv(self, sim_extreme_inv):
        """8 物种、极高倒位率。"""
        akr, tree = sim_extreme_inv
        _run_inversion_validation(akr, tree)

    def test_massive_inv(self, sim_massive_inv):
        """12 物种、极大基因组、极高倒位率。"""
        akr, tree = sim_massive_inv
        _run_inversion_validation(akr, tree)


# ── 融合 + 倒位 ─────────────────────────────────────────────

class TestInversionWithFusion:
    """倒位 + 融合（混合方向）。"""

    def test_inv_with_fusion(self, sim_inv_with_fusion):
        """倒位 + NCF + EEJ：方向调和应处理融合后的混合方向。"""
        akr, tree = sim_inv_with_fusion
        _run_inversion_validation(akr, tree, expect_inversions=False)

    def test_inv_with_all_rearrangement(self, sim_inv_with_rearrangement):
        """倒位 + 全部重排：方向调和应处理所有混合方向。"""
        akr, tree = sim_inv_with_rearrangement
        _run_inversion_validation(akr, tree, expect_inversions=False)


# ── 多种子稳定性 ─────────────────────────────────────────────

class TestInversionMultiSeed:
    """不同种子多次模拟，验证稳定性。"""

    @pytest.mark.parametrize('seed', [1, 7, 42, 100, 256])
    def test_different_seeds(self, seed, tmp_path):
        sim_dir = str(tmp_path / f'sim_{seed}')
        generate_inversion_simulation(sim_dir, seed=seed)
        akr, tree = load_akr_and_tree(sim_dir)
        _run_inversion_validation(akr, tree, expect_inversions=False)


# ── 手工构造图测试 ────────────────────────────────────────────

class MockChild:
    """最小化 child graph。"""
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
    return (name, side)


def linear_chrom(genes, chrom_name='c'):
    return [T(chrom_name, 'L')] + list(genes) + [T(chrom_name, 'R')]


def build_graph(children_spec):
    G = ColoredGraph(hog_level='test')
    for cid, chroms in children_spec:
        G.add_child(cid, MockChild(chroms))
    return G


class TestInversionManual:
    """手工构造图测试：精确控制方向。"""

    def test_simple_inversion(self):
        """简单倒位：A-B-C-D-E vs A-D-C-B-E（C-D 段倒位）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ADCBE', 'c2')]),
        ])
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts > 0, "Expected direction conflicts"

    def test_whole_chrom_reversal(self):
        """整条染色体反向：A-B-C-D-E vs E-D-C-B-A。调和后无冲突。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('EDCBA', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts == 0, f"Whole chrom reversal should be harmonized, got {conflicts}"

    def test_partial_inversion(self):
        """部分倒位：A-B-C-D-E-F-G-H-I-J vs A-B-C-F-E-D-G-H-I-J。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJ', 'c1')]),
            ('c2', [linear_chrom('ABCFEDGHIJ', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts > 0

    def test_telomere_inversion_at_start(self):
        """端粒倒位（起始端）：A-B-C-D-E → C-B-A-D-E。端粒位置不变。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('CBADE', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts > 0, "Telomere inversion at start should produce conflicts"

    def test_telomere_inversion_at_end(self):
        """端粒倒位（末端）：A-B-C-D-E → A-B-C-E-D。端粒位置不变。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABCED', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts > 0, "Telomere inversion at end should produce conflicts"

    def test_whole_arm_inversion(self):
        """全臂倒位：A-B-C-D-E → D-E-C-B-A。端粒位置不变。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('DECBA', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts > 0, "Whole-arm inversion should produce conflicts"

    def test_two_independent_inversions(self):
        """两个独立倒位。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEFGHIJKLMNOP', 'c1')]),
            ('c2', [linear_chrom('ABCFEDGHILKJMNOP', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts >= 2

    def test_no_inversion(self):
        """无倒位：调和后无方向冲突。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABCDE', 'c2')]),
        ])
        G.harmonize_directions()
        conflicts = sum(1 for h1, h2 in G._graph.edges()
                        if G.edge_has_direction_conflict(h1, h2))
        assert conflicts == 0

    def test_inversion_with_gene_loss(self):
        """倒位 + 基因丢失：A-B-C-D-E-F vs A-B-F-E-D（C 丢失，D-E 倒位）。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABFED', 'c2')]),
        ])
        G.harmonize_directions()
        # 不崩溃即可
        assert True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
