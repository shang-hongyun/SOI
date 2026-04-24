#!/usr/bin/env python3
"""Fix evolution_simulator_ak.py in one shot."""
import re

path = "evolution_simulator_ak.py"
with open(path) as f:
    content = f.read()

# 1. Fix __init__ signature and body
old_init = '''    def __init__(self, seed=42, num_chroms=8, min_genes=200, max_genes=1000,
                 inv_rate=10.0, rt_rate=3.0, ncf_rate=1.5, eej_rate=2.0,
                 fission_rate=0.02, wgd_rate=0.4,
                 gene_gain_rate=0.5, tandem_dup_rate=2.0,
                 dispersed_dup_rate=0.3, unidir_trans_rate=1.0,
                 frac_rate=5.0, scale=1.0,
                 rt_mode="random", eej_mode="random"):
        self.seed = seed
        self.rng = random.Random(seed)
        self.num_chroms = num_chroms
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.inv_rate = inv_rate * scale
        self.rt_rate = rt_rate * scale
        self.ncf_rate = ncf_rate * scale
        self.eej_rate = eej_rate * scale
        self.fission_rate = fission_rate * scale
        self.wgd_rate = wgd_rate * scale
        self.gene_gain_rate = gene_gain_rate * scale
        self.tandem_dup_rate = tandem_dup_rate * scale
        self.dispersed_dup_rate = dispersed_dup_rate * scale
        self.unidir_trans_rate = unidir_trans_rate * scale
        self.frac_rate = frac_rate
        self.rt_mode = rt_mode
        self.eej_mode = eej_mode

        self.tracker = GeneTracker()
        self.leaf_karyotypes = {}
        self.all_node_karyotypes = {}
        self.events = []
        self.wgd_map = {}  # species -> max WGD factor (for tree annotation)'''

new_init = '''    def __init__(self, seed=42, num_chroms=8, min_genes=200, max_genes=1000,
                 inv_rate=10.0, rt_rate=3.0, ncf_rate=1.5, eej_rate=2.0,
                 fission_rate=0.02, wgd_rate=0.4,
                 gene_gain_rate=0.5, gene_loss_rate=0.5, tandem_dup_rate=2.0,
                 dispersed_dup_rate=0.3, unidir_trans_rate=1.0,
                 frac_rate=5.0, scale=1.0,
                 seg_del_rate=0.5, seg_dup_rate=0.3, chromothripsis_rate=0.05,
                 rt_mode="random", eej_mode="random"):
        self.seed = seed
        self.rng = random.Random(seed)
        self.num_chroms = num_chroms
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.inv_rate = inv_rate * scale
        self.rt_rate = rt_rate * scale
        self.ncf_rate = ncf_rate * scale
        self.eej_rate = eej_rate * scale
        self.fission_rate = fission_rate * scale
        self.wgd_rate = wgd_rate * scale
        self.gene_gain_rate = gene_gain_rate * scale
        self.gene_loss_rate = gene_loss_rate * scale
        self.tandem_dup_rate = tandem_dup_rate * scale
        self.dispersed_dup_rate = dispersed_dup_rate * scale
        self.unidir_trans_rate = unidir_trans_rate * scale
        self.frac_rate = frac_rate * scale
        self.seg_del_rate = seg_del_rate * scale
        self.seg_dup_rate = seg_dup_rate * scale
        self.chromothripsis_rate = chromothripsis_rate * scale
        self.rt_mode = rt_mode
        self.eej_mode = eej_mode

        self.tracker = GeneTracker()
        self.leaf_karyotypes = {}
        self.all_node_karyotypes = {}
        self.events = []
        self.wgd_map = {}
        self.centromeres = {}'''

content = content.replace(old_init, new_init)

# 2. Fix run / _init_ancestor / _evolve / _apply_wgd
old_run = '''    def run(self, tree, ploidy_map):
        print("Initializing ancestor with {} chromosomes ({}-{} genes each)...".format(
            self.num_chroms, self.min_genes, self.max_genes))
        root_karyo = self._init_ancestor()
        n_leaves = len(tree.get_leaves())
        print("Simulating evolution along tree with {} species...".format(n_leaves))
        self._evolve(tree, root_karyo, ploidy_map)
        self._print_summary()
        return self.leaf_karyotypes

    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
        return karyo

    def _evolve(self, tree, root_karyo, ploidy_map):
        for node in tree.traverse("preorder"):
            if node.is_root():
                self.all_node_karyotypes[node.name] = root_karyo
                continue

            karyo = copy.deepcopy(self.all_node_karyotypes[node.up.name])
            name = node.name
            bl = node.dist if node.dist else 0.0

            # WGD
            if name in ploidy_map and ploidy_map[name] > 1:
                factor = ploidy_map[name]
                print("  Applying {}x genome duplication at {} (from tree annotation)".format(
                    factor, name))
                karyo = self._apply_wgd(karyo, name, factor)
                if self.frac_rate > 0:
                    self._apply_fractionation(karyo, name)
            elif self.wgd_rate > 0 and bl > 0:
                n_wgd = poisson_sample(self.rng, self.wgd_rate * bl)
                for _ in range(n_wgd):
                    factor = sample_wgd_factor(self.rng)
                    print("  Applying {}x genome duplication at {} (rate-based)".format(
                        factor, name))
                    karyo = self._apply_wgd(karyo, name, factor)
                    if self.frac_rate > 0:
                        self._apply_fractionation(karyo, name)

            if bl > 0:
                self._apply_rearrangements(karyo, name, bl)

            self.all_node_karyotypes[name] = karyo
            if node.is_leaf():
                self.leaf_karyotypes[name] = karyo
                for cid, genes in karyo.items():
                    for gid, orient in genes:
                        self.tracker.set_species(gid, name)'''

new_run = '''    def run(self, tree, ploidy_map):
        print("Initializing ancestor with {} chromosomes ({}-{} genes each)...".format(
            self.num_chroms, self.min_genes, self.max_genes))
        root_karyo, root_centros = self._init_ancestor()
        n_leaves = len(tree.get_leaves())
        print("Simulating evolution along tree with {} species...".format(n_leaves))
        self._evolve(tree, root_karyo, root_centros, ploidy_map)
        self._print_summary()
        return self.leaf_karyotypes

    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        centros = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
            centros[cid] = self.rng.randint(2, n - 2) if n > 4 else n // 2
        return karyo, centros

    def _evolve(self, tree, root_karyo, root_centros, ploidy_map):
        for node in tree.traverse("preorder"):
            if node.is_root():
                self.all_node_karyotypes[node.name] = root_karyo
                self.centromeres[node.name] = root_centros
                continue

            karyo = copy.deepcopy(self.all_node_karyotypes[node.up.name])
            centros = copy.deepcopy(self.centromeres[node.up.name])
            name = node.name
            bl = node.dist if node.dist else 0.0

            # WGD
            if name in ploidy_map and ploidy_map[name] > 1:
                factor = ploidy_map[name]
                print("  Applying {}x genome duplication at {} (from tree annotation)".format(
                    factor, name))
                karyo, centros = self._apply_wgd(karyo, centros, name, factor)
            elif self.wgd_rate > 0 and bl > 0:
                n_wgd = poisson_sample(self.rng, self.wgd_rate * bl)
                for _ in range(n_wgd):
                    factor = sample_wgd_factor(self.rng)
                    print("  Applying {}x genome duplication at {} (rate-based)".format(
                        factor, name))
                    karyo, centros = self._apply_wgd(karyo, centros, name, factor)

            if bl > 0:
                self._apply_rearrangements(karyo, centros, name, bl)

            self.all_node_karyotypes[name] = karyo
            self.centromeres[name] = centros
            if node.is_leaf():
                self.leaf_karyotypes[name] = karyo
                for cid, genes in karyo.items():
                    for gid, orient in genes:
                        self.tracker.set_species(gid, name)'''

content = content.replace(old_run, new_run)

# 3. Fix _apply_wgd
old_wgd = '''    def _apply_wgd(self, karyo, node_name, factor):
        new_k = {}
        for c, g in karyo.items():
            for copy_i in range(1, factor + 1):
                new_c = "{}_{}".format(c, copy_i)
                new_g = []
                for gid, orient in g:
                    new_gid = self.tracker.new_wgd_copy(gid, node_name, copy_i)
                    new_g.append((new_gid, orient))
                new_k[new_c] = new_g
        self.events.append({"node": node_name, "type": "WGD", "factor": factor,
                            "desc": "{}->{} chromosomes".format(len(karyo), len(new_k))})
        # Track max ploidy per species
        cur = self.wgd_map.get(node_name, 1)
        self.wgd_map[node_name] = max(cur, factor)
        return new_k'''

new_wgd = '''    def _apply_wgd(self, karyo, centros, node_name, factor):
        new_k = {}
        new_centros = {}
        for c, g in karyo.items():
            for copy_i in range(1, factor + 1):
                new_c = "{}_{}".format(c, copy_i)
                new_g = []
                for gid, orient in g:
                    new_gid = self.tracker.new_wgd_copy(gid, node_name, copy_i)
                    new_g.append((new_gid, orient))
                new_k[new_c] = new_g
                new_centros[new_c] = centros.get(c, len(new_g) // 2)
        self.events.append({"node": node_name, "type": "WGD", "factor": factor,
                            "desc": "{}->{} chromosomes".format(len(karyo), len(new_k))})
        # Track max ploidy per species
        cur = self.wgd_map.get(node_name, 1)
        self.wgd_map[node_name] = max(cur, factor)
        return new_k, new_centros'''

content = content.replace(old_wgd, new_wgd)

# 4. Fix _apply_rearrangements: add n_loss, n_sdel, n_sdup, n_cht
old_rearr = '''    def _apply_rearrangements(self, karyo, centros, node_name, branch_length):
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)
        n_frac = poisson_sample(self.rng, self.frac_rate * branch_length)

        event_list = (
            ["inv"] * n_inv + ["rt"] * n_rt + ["ncf"] * n_ncf +
            ["eej"] * n_eej + ["fis"] * n_fis + ["gain"] * n_gain +
            ["loss"] * n_loss + ["tdup"] * n_tdup + ["ddup"] * n_ddup +
            ["utrans"] * n_utrans + ["frac"] * n_frac +
            ["segdel"] * n_sdel + ["segdup"] * n_sdup + ["cht"] * n_cht
        )'''

new_rearr = '''    def _apply_rearrangements(self, karyo, centros, node_name, branch_length):
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_loss = poisson_sample(self.rng, self.gene_loss_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)
        n_frac = poisson_sample(self.rng, self.frac_rate * branch_length)
        n_sdel = poisson_sample(self.rng, self.seg_del_rate * branch_length)
        n_sdup = poisson_sample(self.rng, self.seg_dup_rate * branch_length)
        n_cht = poisson_sample(self.rng, self.chromothripsis_rate * branch_length)

        event_list = (
            ["inv"] * n_inv + ["rt"] * n_rt + ["ncf"] * n_ncf +
            ["eej"] * n_eej + ["fis"] * n_fis + ["gain"] * n_gain +
            ["loss"] * n_loss + ["tdup"] * n_tdup + ["ddup"] * n_ddup +
            ["utrans"] * n_utrans + ["frac"] * n_frac +
            ["segdel"] * n_sdel + ["segdup"] * n_sdup + ["cht"] * n_cht
        )'''

content = content.replace(old_rearr, new_rearr)

# 5. Fix xmain: add new params
old_xmain = '''        unidir_trans_rate=kargs.get('unidir_trans_rate', 1.0),
        frac_rate=kargs.get('frac_rate', 5.0),
        scale=kargs.get('scale', 1.0),
    )'''

new_xmain = '''        unidir_trans_rate=kargs.get('unidir_trans_rate', 1.0),
        frac_rate=kargs.get('frac_rate', 5.0),
        seg_del_rate=kargs.get('seg_del_rate', 0.5),
        seg_dup_rate=kargs.get('seg_dup_rate', 0.3),
        chromothripsis_rate=kargs.get('chromothripsis_rate', 0.05),
        scale=kargs.get('scale', 1.0),
    )'''

content = content.replace(old_xmain, new_xmain)

# 6. Remove duplicate _apply_chromothripsis if present
dup_pattern = r'(    def _apply_chromothripsis\(self, karyo, node_name, centros=None\):.*?return True\n)\n    def _apply_chromothripsis\(self, karyo, node_name, centros=None\):.*?return True\n'
match = re.search(dup_pattern, content, re.DOTALL)
if match:
    content = content[:match.start()] + match.group(1) + content[match.end():]

with open(path, 'w') as f:
    f.write(content)

print("Fixes applied.")
