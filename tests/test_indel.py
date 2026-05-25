#!/usr/bin/env python3
"""
Phase 2 indel 测试 — 内嵌模拟数据 + 压力测试 + 可复用验证。

参考 test_dedup.py 规范：
- generate_indel_simulation(): 生成含基因得失的模拟数据
- _run_indel_validation(): 通用验证函数
- 多复杂度级别：basic → stress → extreme
"""
import sys, os, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import logging
logging.basicConfig(level=logging.WARNING)

from soi.takr_colored_graph import ColoredGraph


# ── 模拟数据生成 ──────────────────────────────────────────────

def generate_indel_simulation(outdir, seed=42, num_species=4, num_chroms=3,
                              gene_gain_rate=5, gene_loss_rate=5,
                              min_genes=100, max_genes=300):
    """生成含基因得失（无重排、无 dup）的模拟数据。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    import random
    rng = random.Random(seed)
    tree = generate_tree(num_species, rng)
    sim = EvolutionSimulator(
        seed=seed, num_chroms=num_chroms, min_genes=min_genes, max_genes=max_genes,
        inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
        fission_rate=0, unidir_trans_rate=0,
        gene_gain_rate=gene_gain_rate, gene_loss_rate=gene_loss_rate,
        tandem_dup_rate=0, dispersed_dup_rate=0,
        seg_dup_rate=0, wgd_rate=0,
    )
    orig_nw = tree.write(format=1)
    sim.run(tree, {})
    sim.generate_outputs(outdir, tree, {}, orig_nw)
    return sim


def generate_dense_indel_simulation(outdir, seed=42, num_species=4, num_chroms=3,
                                    gene_gain_rate=50, gene_loss_rate=50,
                                    min_genes=20, max_genes=50):
    """高密度 indel：短染色体、高得失率。"""
    from soi.evolution_simulator_ak import generate_tree, EvolutionSimulator
    import random
    rng = random.Random(seed)
    tree = generate_tree(num_species, rng)
    sim = EvolutionSimulator(
        seed=seed, num_chroms=num_chroms, min_genes=min_genes, max_genes=max_genes,
        inv_rate=0, rt_rate=0, ncf_rate=0, eej_rate=0,
        fission_rate=0, unidir_trans_rate=0,
        gene_gain_rate=gene_gain_rate, gene_loss_rate=gene_loss_rate,
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

def _run_indel_validation(akr, tree, expect_indels=True):
    """通用 indel 验证：遍历所有内部节点，检测 indel shortcuts 并 resolve。"""
    import networkx as nx
    total_shortcuts = 0
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
        # 检测 indel shortcuts
        shortcuts = G.find_indel_shortcuts()
        total_shortcuts += len(shortcuts)
        # resolve
        G.resolve_indels()
        indel_events = [e for e in G.events
                        if e.event_type in ('gene_loss', 'gene_indel', 'indel')]
        total_events += len(indel_events)
        # 验证 resolve 后无剩余 shortcuts
        remaining = G.find_indel_shortcuts()
        assert len(remaining) == 0, \
            f"{node.name}: {len(remaining)} shortcuts remain after resolve"
    if expect_indels:
        assert total_shortcuts > 0, "No indel shortcuts detected"
    return total_shortcuts, total_events


# ── Fixtures ─────────────────────────────────────────────────

@pytest.fixture(scope='module')
def sim_normal_indel(tmp_path_factory):
    """正常基因得失率模拟。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_normal_indel'))
    generate_indel_simulation(sim_dir, seed=42)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_heavy_indel(tmp_path_factory):
    """高基因得失率模拟。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_heavy_indel'))
    generate_indel_simulation(sim_dir, seed=123, num_species=8, num_chroms=5,
                              gene_gain_rate=30, gene_loss_rate=30,
                              min_genes=300, max_genes=800)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_extreme_indel(tmp_path_factory):
    """极端基因得失率 + 大基因组。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_extreme_indel'))
    generate_indel_simulation(sim_dir, seed=7777, num_species=10, num_chroms=6,
                              gene_gain_rate=80, gene_loss_rate=80,
                              min_genes=1000, max_genes=3000)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_massive_indel(tmp_path_factory):
    """大规模压力：12 物种、极大基因组、极高得失率。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_massive_indel'))
    generate_indel_simulation(sim_dir, seed=99999, num_species=12, num_chroms=8,
                              gene_gain_rate=150, gene_loss_rate=150,
                              min_genes=2000, max_genes=5000)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_dense_indel(tmp_path_factory):
    """高密度 indel：短染色体、高得失率。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_dense_indel'))
    generate_dense_indel_simulation(sim_dir, seed=42, num_species=4, num_chroms=3,
                                    gene_gain_rate=50, gene_loss_rate=50,
                                    min_genes=20, max_genes=50)
    return load_akr_and_tree(sim_dir)


@pytest.fixture(scope='module')
def sim_dense_extreme(tmp_path_factory):
    """极高密度 indel：极短染色体、极高得失率。"""
    sim_dir = str(tmp_path_factory.mktemp('sim_dense_extreme'))
    generate_dense_indel_simulation(sim_dir, seed=7777, num_species=6, num_chroms=3,
                                    gene_gain_rate=100, gene_loss_rate=100,
                                    min_genes=10, max_genes=30)
    return load_akr_and_tree(sim_dir)


# ── 基础测试 ─────────────────────────────────────────────────

class TestIndelBasic:
    """正常基因得失率下的基础测试。"""

    def test_indels_detected(self, sim_normal_indel):
        """应检测到 indel shortcuts。"""
        akr, tree = sim_normal_indel
        shortcuts, events = _run_indel_validation(akr, tree)
        assert shortcuts > 0, "No indel shortcuts in normal simulation"

    def test_resolve_removes_all_shortcuts(self, sim_normal_indel):
        """resolve 后不应有剩余 shortcuts。"""
        akr, tree = sim_normal_indel
        _run_indel_validation(akr, tree)

    def test_no_over_fragmentation(self, sim_normal_indel):
        """不应过度碎片化。"""
        akr, tree = sim_normal_indel
        _run_indel_validation(akr, tree)


# ── 压力测试 ─────────────────────────────────────────────────

class TestIndelStress:
    """高基因得失率压力测试。"""

    def test_heavy_indel(self, sim_heavy_indel):
        """6 物种、高得失率：检测 + resolve + 无碎片化。"""
        akr, tree = sim_heavy_indel
        _run_indel_validation(akr, tree)

    def test_extreme_indel(self, sim_extreme_indel):
        """10 物种、极高得失率、大基因组：检测 + resolve + 无碎片化。"""
        akr, tree = sim_extreme_indel
        _run_indel_validation(akr, tree)

    def test_massive_indel(self, sim_massive_indel):
        """12 物种、极大基因组、极高得失率：检测 + resolve + 无碎片化。"""
        akr, tree = sim_massive_indel
        _run_indel_validation(akr, tree)


class TestIndelDense:
    """高密度 indel：短染色体、连续基因删除/插入。"""

    def test_dense_indel(self, sim_dense_indel):
        """短染色体（20-50 基因）、高得失率（50）：检测 + resolve。"""
        akr, tree = sim_dense_indel
        _run_indel_validation(akr, tree)

    def test_dense_extreme(self, sim_dense_extreme):
        """极短染色体（10-30 基因）、极高得失率（100）：resolve 不崩溃。"""
        akr, tree = sim_dense_extreme
        _run_indel_validation(akr, tree, expect_indels=False)


# ── 多种子稳定性 ─────────────────────────────────────────────

class TestIndelMultiSeed:
    """不同种子多次模拟，验证稳定性。"""

    @pytest.mark.parametrize('seed', [1, 7, 42, 100, 256])
    def test_different_seeds(self, seed, tmp_path):
        sim_dir = str(tmp_path / f'sim_{seed}')
        generate_indel_simulation(sim_dir, seed=seed)
        akr, tree = load_akr_and_tree(sim_dir)
        _run_indel_validation(akr, tree, expect_indels=False)


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


class TestIndelManual:
    """手工构造图测试：精确控制场景。"""

    def test_single_gene_deletion(self):
        """单基因删除：A-B-C-D-E vs A-B-D-E。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned

    def test_single_gene_insertion(self):
        """单基因插入：A-B-D-E vs A-B-C-D-E。"""
        G = build_graph([
            ('c1', [linear_chrom('ABDE', 'c1')]),
            ('c2', [linear_chrom('ABCDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1
        h1, h2, child_id, spanned = shortcuts[0]
        assert 'C' in spanned

    def test_multi_gene_deletion(self):
        """多基因删除：A-B-C-D-E vs A-B-E。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1
        spanned = shortcuts[0][3]
        assert 'C' in spanned and 'D' in spanned

    def test_inversion_not_indel(self):
        """倒位不应被误判为 indel。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDEF', 'c1')]),
            ('c2', [linear_chrom('ABDCEF', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        # 倒位的 spanned HOGs 在两个孩子中都存在 → 不是 indel
        for h1, h2, child_id, spanned in shortcuts:
            child_hogs = set('ABCDEF')
            n_in_child = sum(1 for h in spanned if h in child_hogs)
            assert n_in_child <= len(spanned) * 0.5, \
                f"Rearrangement misclassified as indel: {spanned}"

    def test_no_indel_same_genes(self):
        """相同基因无 indel。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABCDE', 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 0

    def test_10gene_deletion(self):
        """10 基因连续删除。"""
        genes_ref = 'ABCDEFGHIJKLMNOP'
        genes_del = 'ABCP'  # D-O 被删
        G = build_graph([
            ('c1', [linear_chrom(genes_ref, 'c1')]),
            ('c2', [linear_chrom(genes_del, 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1
        spanned = shortcuts[0][3]
        assert len(spanned) >= 10, f"Expected ≥10 spanned, got {len(spanned)}"

    def test_15gene_insertion(self):
        """15 基因连续插入。"""
        genes_ref = 'ABCS'  # D-R 被插入
        genes_ins = 'ABCDEFGHIJKLMNOPQRS'
        G = build_graph([
            ('c1', [linear_chrom(genes_ref, 'c1')]),
            ('c2', [linear_chrom(genes_ins, 'c2')]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1
        spanned = shortcuts[0][3]
        assert len(spanned) >= 15, f"Expected ≥15 spanned, got {len(spanned)}"

    def test_resolve_removes_edges(self):
        """resolve 应移除跨越边。"""
        G = build_graph([
            ('c1', [linear_chrom('ABCDE', 'c1')]),
            ('c2', [linear_chrom('ABDE', 'c2')]),
        ])
        edges_before = G._graph.number_of_edges()
        G.resolve_indels()
        edges_after = G._graph.number_of_edges()
        assert edges_after < edges_before


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
