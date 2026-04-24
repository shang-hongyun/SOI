import re

with open('evolution_simulator_ak.py') as f:
    content = f.read()

# Fix 1: Add gene_loss_rate and other params to __init__
content = content.replace(
    'gene_gain_rate=0.5, tandem_dup_rate=2.0,',
    'gene_gain_rate=0.5, gene_loss_rate=0.5, tandem_dup_rate=2.0,'
)
content = content.replace(
    'self.gene_gain_rate = gene_gain_rate * scale\n        self.tandem_dup_rate',
    'self.gene_gain_rate = gene_gain_rate * scale\n        self.gene_loss_rate = gene_loss_rate * scale\n        self.tandem_dup_rate'
)

# Fix 2: Add centromeres
if 'self.centromeres = {}' not in content:
    content = content.replace(
        'self.wgd_map = {}  # species -> max WGD factor (for tree annotation)',
        'self.wgd_map = {}\n        self.centromeres = {}'
    )

# Fix 3: Fix run/_init_ancestor/_evolve/_apply_wgd
content = content.replace(
    'root_karyo = self._init_ancestor()',
    'root_karyo, root_centros = self._init_ancestor()'
)
content = content.replace(
    'self._evolve(tree, root_karyo, ploidy_map)',
    'self._evolve(tree, root_karyo, root_centros, ploidy_map)'
)

# _init_ancestor
old_init_anc = '''    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
        return karyo'''
new_init_anc = '''    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        centros = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
            centros[cid] = self.rng.randint(2, n - 2) if n > 4 else n // 2
        return karyo, centros'''
content = content.replace(old_init_anc, new_init_anc)

# _evolve - need to handle carefully
old_evolve = '''    def _evolve(self, tree, root_karyo, ploidy_map):
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

new_evolve = '''    def _evolve(self, tree, root_karyo, root_centros, ploidy_map):
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
content = content.replace(old_evolve, new_evolve)

# _apply_wgd
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

with open('evolution_simulator_ak.py', 'w') as f:
    f.write(content)
print("Done")
