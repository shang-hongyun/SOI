#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Rearrangement Evolution Simulator (AK.py output format)
==================================================================

Simulates realistic chromosome rearrangements along a phylogenetic tree,
producing output compatible with AK.py (Ancestral Karyotype Reconstruction).

Rearrangement events modeled:
  - Inversion
  - RT (Reciprocal Translocation, may be unbalanced)
  - NCF (Non-Centric Fusion)
  - EEJ (End-End Joining)
  - Fission (very low rate, matching real biology)
  - WGD / Polyploidy (rate-based; 2x most common, 3x less, etc.)
  - Fractionation (gene loss after WGD)
  - Tandem duplication (mostly single gene, some contiguous multi-gene)
  - Dispersed duplication (rare)
  - Gene gain (de novo, species-specific)
  - Unidirectional translocation (small segment moves one-way)

Output files (AK.py input + validation):
  1. ortholog_groups.txt   -- SOG file
  2. ortholog_pairs.txt    -- pairwise ortholog relationships
  3. all_species_gene.gff  -- leaf species gene order
  4. species_tree.nwk      -- Newick tree with branch lengths
  5. ancestors_gene.gff    -- ALL node gene orders (for validation)
  6. true_hogs.tsv         -- true HOG hierarchy (for validation)
  7. events.tsv            -- event log

Usage as subcommand:
    soi sim -t tree.nwk -o output_dir [options]

Tree format:
    Newick with branch lengths.  WGD is now rate-based (no [p=N] needed).
    Example: ((A:0.15, B:0.12):0.08, (C:0.10, (D:0.18, E:0.14):0.06):0.25);
"""

from __future__ import print_function
import sys
import os
import copy
import math
import random
import string
import argparse
import re
from collections import defaultdict, Counter

from ete3 import Tree
import networkx as nx


# ============================================================================
# Helper functions
# ============================================================================

Gene = tuple  # (gene_id: str, orientation: str '+'/'-')
Karyotype = dict  # {chrom_id: [Gene, ...]}


def flip_orient(o):
    return "+" if o == "-" else "-"


def reverse_segment(genes):
    return [(g, flip_orient(o)) for g, o in reversed(genes)]


def get_chr_labels(n):
    base = list(string.ascii_uppercase)
    if n <= 26:
        return base[:n]
    labels = list(base)
    for a in base:
        for b in base:
            labels.append(a + b)
            if len(labels) >= n:
                return labels
    return labels


def poisson_sample(rng, lam):
    if lam <= 0:
        return 0
    if lam < 30:
        L = math.exp(-lam)
        k = 0
        p = 1.0
        while True:
            k += 1
            p *= rng.random()
            if p <= L:
                return k - 1
    else:
        val = rng.gauss(lam, math.sqrt(lam))
        return max(0, int(round(val)))


def sample_wgd_factor(rng):
    """Sample WGD factor: P(2) > P(3) > P(4) > ...
    P(n) = 1/2^(n-1)  for n >= 2, normalized."""
    r = rng.random()
    # P(2) = 1/2, P(3) = 1/4, P(4) = 1/8, ...
    n = 2
    cum = 0.5
    while r > cum:
        n += 1
        cum += 0.5 ** (n - 1)
    return n


def parse_tree(tree_file, wgd_rate=0.0):
    """Parse tree file. If wgd_rate > 0, WGD events are sampled along
    branches (overriding any [p=N] annotations).  If wgd_rate == 0,
    [p=N] annotations in the tree are used."""
    with open(tree_file) as f:
        orig_nw = f.read().strip()

    ploidy_map = {}
    if wgd_rate == 0.0:
        # Use [p=N] annotations from tree
        for m in re.finditer(r'([\w. -]+)\[p=(\d+)\]', orig_nw):
            name = m.group(1).strip()
            ploidy = int(m.group(2))
            ploidy_map[name] = ploidy

    clean_nw = re.sub(r'\[p=\d+\]', '', orig_nw)
    tree = None
    for fmt in [1, 0, 2, 3, 8]:
        try:
            tree = Tree(clean_nw, format=fmt)
            break
        except Exception:
            continue
    if tree is None:
        raise ValueError("Cannot parse tree file: " + tree_file)

    # Collect existing names to avoid collisions
    existing_names = set(node.name for node in tree.traverse() if node.name)
    counter = 0
    for node in tree.traverse("preorder"):
        if not node.name or node.name.strip() == "":
            while True:
                candidate = "N{}".format(counter)
                counter += 1
                if candidate not in existing_names:
                    node.name = candidate
                    existing_names.add(candidate)
                    break
    return tree, ploidy_map, orig_nw


def number_tree_nodes(tree):
    """Number internal tree nodes as N0, N1, ... (matching AK.py's number_nodes)."""
    i = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        if node.name and re.match(r'^N\d+$', node.name):
            continue
        old = node.name
        node.name = "N{}".format(i)
        i += 1
    internal_nodes = [n for n in tree.traverse("postorder") if not n.is_leaf()]
    name_map = {}
    for idx, node in enumerate(internal_nodes):
        old_name = node.name
        new_name = "N{}".format(idx)
        name_map[old_name] = new_name
        node.name = new_name
    return name_map


# ============================================================================
# GeneTracker -- tracks gene ancestry, WGD copies, duplications
# Dot-notation naming: homologous genes share a prefix.
#   Ancestral gene: g1, g2, ...
#   WGD copy:       g1.1, g1.2, ...
#   Tandem dup:     g1.3, g1.4, ...
#   Dispersed dup:  g1.5, ...
#   Gene gain:      gn1, gn2, ... (own root, no homologs)
# ============================================================================

class GeneTracker:
    def __init__(self):
        self._anc_counter = 0      # for ancestral genes
        self._gain_counter = 0     # for gene gains
        self._homolog_seq = {}     # anc_gene -> next suffix index
        self.ancestor_genes = []       # ordered list of root ancestral gene IDs
        self.gene_ancestor = {}        # gene_id -> ultimate ancestral gene ID
        self.gene_species = {}         # gene_id -> species (for leaf genes)
        self.gene_wgd_copies = {}      # gene_id -> {wgd_node_name: copy_index}
        self.gene_is_gain = {}         # gene_id -> bool
        self.gene_dup_of = {}          # gene_id -> original gene_id
        self.gene_dup_type = {}        # gene_id -> 'tandem'|'dispersed'|None
        self.gain_genes = []           # list of gene gain IDs

    def _next_homolog_id(self, anc_gene):
        """Get the next dot-suffix ID for a homolog of anc_gene."""
        idx = self._homolog_seq.get(anc_gene, 1)
        self._homolog_seq[anc_gene] = idx + 1
        return "{}.{}".format(anc_gene, idx)

    def new_ancestral_gene(self):
        self._anc_counter += 1
        gid = "g{}".format(self._anc_counter)
        self.gene_ancestor[gid] = gid
        self.ancestor_genes.append(gid)
        self.gene_wgd_copies[gid] = {}
        self.gene_is_gain[gid] = False
        return gid

    def new_wgd_copy(self, original_gid, wgd_node, copy_index, species=None):
        """Create a new gene as WGD copy of original_gid."""
        anc = self.gene_ancestor[original_gid]
        gid = self._next_homolog_id(anc)
        self.gene_ancestor[gid] = anc
        if species is not None:
            self.gene_species[gid] = species
        wgd = dict(self.gene_wgd_copies.get(original_gid, {}))
        wgd[wgd_node] = copy_index
        self.gene_wgd_copies[gid] = wgd
        self.gene_is_gain[gid] = False
        return gid

    def new_duplicate_gene(self, original_gid, dup_type, species=None):
        """Create a duplicate gene (tandem or dispersed)."""
        anc = self.gene_ancestor[original_gid]
        gid = self._next_homolog_id(anc)
        self.gene_ancestor[gid] = anc
        if species is not None:
            self.gene_species[gid] = species
        self.gene_wgd_copies[gid] = dict(self.gene_wgd_copies.get(original_gid, {}))
        self.gene_is_gain[gid] = False
        self.gene_dup_of[gid] = original_gid
        self.gene_dup_type[gid] = dup_type
        return gid

    def new_gain_gene(self, species):
        """Create a completely new gene (gene gain)."""
        self._gain_counter += 1
        gid = "gn{}".format(self._gain_counter)
        self.gene_ancestor[gid] = gid  # own ancestor
        self.gene_species[gid] = species
        self.gene_wgd_copies[gid] = {}
        self.gene_is_gain[gid] = True
        self.gain_genes.append(gid)
        return gid

    def get_ancestor(self, gid):
        return self.gene_ancestor.get(gid, gid)

    def set_species(self, gid, species):
        self.gene_species[gid] = species

    def total_genes(self):
        total = self._anc_counter + self._gain_counter
        total += sum(v - 1 for v in self._homolog_seq.values())
        return total


# ============================================================================
# EvolutionSimulator
# ============================================================================

class EvolutionSimulator:
    def __init__(self, seed=42, num_chroms=8, min_genes=200, max_genes=1000,
                 inv_rate=10.0, rt_rate=3.0, ncf_rate=1.5, eej_rate=2.0,
                 fission_rate=0.05, wgd_rate=0.05,
                 gene_gain_rate=0.5, tandem_dup_rate=2.0,
                 dispersed_dup_rate=0.3, unidir_trans_rate=1.0,
                 frac_rate=5.0, scale=1.0,
                 rt_mode="random", eej_mode="random"):
        # Apply scale factor to all rates
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
        self.frac_rate = frac_rate  # frac_rate is probability per gene, NOT scaled
        self.rt_mode = rt_mode
        self.eej_mode = eej_mode

        self.tracker = GeneTracker()
        self.leaf_karyotypes = {}
        self.all_node_karyotypes = {}   # ALL nodes (internal + leaves)
        self.events = []

    def run(self, tree, ploidy_map):
        print("Initializing ancestor with {} chromosomes ({}-{} genes each)...".format(
            self.num_chroms, self.min_genes, self.max_genes))
        root_karyo = self._init_ancestor()
        n_leaves = len(tree.get_leaves())
        print("Simulating evolution along tree with {} species...".format(n_leaves))
        self._evolve(tree, root_karyo, ploidy_map)
        self._print_summary()
        return self.leaf_karyotypes

    # ------------------------------------------------------------------
    # Ancestor initialization
    # ------------------------------------------------------------------
    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
        return karyo

    # ------------------------------------------------------------------
    # Core evolution loop
    # ------------------------------------------------------------------
    def _evolve(self, tree, root_karyo, ploidy_map):
        for node in tree.traverse("preorder"):
            if node.is_root():
                self.all_node_karyotypes[node.name] = root_karyo
                continue

            karyo = copy.deepcopy(self.all_node_karyotypes[node.up.name])
            name = node.name
            bl = node.dist if node.dist else 0.0

            # WGD / Polyploidy (rate-based or from ploidy_map)
            if name in ploidy_map and ploidy_map[name] > 1:
                factor = ploidy_map[name]
                print("  Applying {}x genome duplication at {} (from tree annotation)".format(
                    factor, name))
                karyo = self._apply_wgd(karyo, name, factor)
                if self.frac_rate > 0:
                    self._apply_fractionation(karyo, name)
            elif self.wgd_rate > 0 and bl > 0:
                # Rate-based WGD: Poisson number of WGD events on this branch
                n_wgd = poisson_sample(self.rng, self.wgd_rate * bl)
                for _ in range(n_wgd):
                    factor = sample_wgd_factor(self.rng)
                    print("  Applying {}x genome duplication at {} (rate-based)".format(
                        factor, name))
                    karyo = self._apply_wgd(karyo, name, factor)
                    if self.frac_rate > 0:
                        self._apply_fractionation(karyo, name)

            # Chromosome rearrangements scaled by branch length
            if bl > 0:
                self._apply_rearrangements(karyo, name, bl)

            self.all_node_karyotypes[name] = karyo
            if node.is_leaf():
                self.leaf_karyotypes[name] = karyo
                # Register species for all genes in this leaf
                for cid, genes in karyo.items():
                    for gid, orient in genes:
                        self.tracker.set_species(gid, name)

    # ------------------------------------------------------------------
    # WGD / Polyploidy
    # ------------------------------------------------------------------
    def _apply_wgd(self, karyo, node_name, factor):
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
        return new_k

    # ------------------------------------------------------------------
    # Fractionation
    # ------------------------------------------------------------------
    def _apply_fractionation(self, karyo, node_name):
        lost = 0
        for cid in list(karyo.keys()):
            retained = []
            for gid, orient in karyo[cid]:
                if self.rng.random() >= self.frac_rate:
                    retained.append((gid, orient))
                else:
                    self.tracker.gene_ancestor.pop(gid, None)
                    self.tracker.gene_species.pop(gid, None)
                    self.tracker.gene_wgd_copies.pop(gid, None)
                    lost += 1
            if retained:
                karyo[cid] = retained
            else:
                del karyo[cid]
        if lost > 0:
            self.events.append({"node": node_name, "type": "fractionation",
                                "lost": lost, "rate": self.frac_rate})

    # ------------------------------------------------------------------
    # Rearrangement dispatch
    # ------------------------------------------------------------------
    def _apply_rearrangements(self, karyo, node_name, branch_length):
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)

        event_list = (
            ["inv"] * n_inv + ["rt"] * n_rt + ["ncf"] * n_ncf +
            ["eej"] * n_eej + ["fis"] * n_fis + ["gain"] * n_gain +
            ["tdup"] * n_tdup + ["ddup"] * n_ddup + ["utrans"] * n_utrans
        )
        if not event_list:
            return
        self.rng.shuffle(event_list)

        applied = defaultdict(int)
        for e in event_list:
            if e == "inv" and self._apply_inversion(karyo, node_name):
                applied["inv"] += 1
            elif e == "rt" and self._apply_rt(karyo, node_name):
                applied["rt"] += 1
            elif e == "ncf" and self._apply_ncf(karyo, node_name):
                applied["ncf"] += 1
            elif e == "eej" and self._apply_eej(karyo, node_name):
                applied["eej"] += 1
            elif e == "fis" and self._apply_fission(karyo, node_name):
                applied["fis"] += 1
            elif e == "gain" and self._apply_gene_gain(karyo, node_name):
                applied["gain"] += 1
            elif e == "tdup" and self._apply_tandem_dup(karyo, node_name):
                applied["tdup"] += 1
            elif e == "ddup" and self._apply_dispersed_dup(karyo, node_name):
                applied["ddup"] += 1
            elif e == "utrans" and self._apply_unidir_trans(karyo, node_name):
                applied["utrans"] += 1

        self.events.append({"node": node_name, "type": "rearrangements",
                            "branch_length": branch_length,
                            "sampled": {"inv": n_inv, "rt": n_rt, "ncf": n_ncf,
                                        "eej": n_eej, "fis": n_fis, "gain": n_gain,
                                        "tdup": n_tdup, "ddup": n_ddup, "utrans": n_utrans},
                            "applied": dict(applied)})

    # ------------------------------------------------------------------
    # Inversion
    # ------------------------------------------------------------------
    def _apply_inversion(self, karyo, node_name):
        cids = [c for c in karyo if len(karyo[c]) >= 2]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        p1, p2 = sorted(self.rng.sample(range(len(genes)), 2))
        if p1 == p2:
            return False
        karyo[cid] = genes[:p1] + reverse_segment(genes[p1:p2]) + genes[p2:]
        self.events.append({"node": node_name, "type": "inversion",
                            "chrom": cid, "pos": "{}-{}".format(p1, p2)})
        return True

    # ------------------------------------------------------------------
    # RT (Reciprocal Translocation, may be unbalanced)
    # ------------------------------------------------------------------
    def _apply_rt(self, karyo, node_name):
        cids = [c for c in karyo if len(karyo[c]) >= 4]
        if len(cids) < 2:
            return False
        c1, c2 = self.rng.sample(cids, 2)
        g1, g2 = karyo[c1], karyo[c2]
        b1 = self.rng.randint(2, len(g1) - 2)
        b2 = self.rng.randint(2, len(g2) - 2)
        aL, aR = g1[:b1], g1[b1:]
        bL, bR = g2[:b2], g2[b2:]

        mode = self.rt_mode
        if mode == "random":
            mode = self.rng.choice(["tailswap", "alternate"])

        if mode == "alternate":
            karyo[c1] = aL + reverse_segment(bL)
            karyo[c2] = reverse_segment(aR) + bR
        else:
            karyo[c1] = aL + bR
            karyo[c2] = bL + aR

        self.events.append({"node": node_name, "type": "RT",
                            "chroms": [c1, c2],
                            "pos": "{}:{}|{}:{}".format(c1, b1, c2, b2),
                            "mode": mode})
        return True

    # ------------------------------------------------------------------
    # NCF (Non-Centric Fusion)
    # ------------------------------------------------------------------
    def _apply_ncf(self, karyo, node_name):
        if len(karyo) < 2:
            return False
        keys = list(karyo.keys())
        recipient_cid, donor_cid = self.rng.sample(keys, 2)
        recipient = karyo[recipient_cid]
        donor = list(karyo[donor_cid])
        if len(recipient) < 2:
            return False
        insert_pos = self.rng.randint(1, len(recipient) - 1)
        if self.rng.choice([True, False]):
            donor = reverse_segment(donor)
        karyo[recipient_cid] = recipient[:insert_pos] + donor + recipient[insert_pos:]
        del karyo[donor_cid]
        self.events.append({"node": node_name, "type": "NCF",
                            "chroms": [recipient_cid, donor_cid],
                            "pos": "insert {} into {}:{}".format(donor_cid, recipient_cid, insert_pos)})
        return True

    # ------------------------------------------------------------------
    # EEJ (End-End Joining)
    # ------------------------------------------------------------------
    def _apply_eej(self, karyo, node_name):
        if len(karyo) < 2:
            return False
        c1, c2 = self.rng.sample(list(karyo.keys()), 2)
        g1, g2 = list(karyo[c1]), list(karyo[c2])
        e1 = self.rng.choice(["H", "T"])
        e2 = self.rng.choice(["H", "T"])
        if e1 == "H":
            g1 = reverse_segment(g1)
        if e2 == "T":
            g2 = reverse_segment(g2)
        karyo[c1] = g1 + g2
        del karyo[c2]
        self.events.append({"node": node_name, "type": "EEJ",
                            "chroms": [c1, c2], "mode": "{}{}".format(e1, e2)})
        return True

    # ------------------------------------------------------------------
    # Fission (very low rate)
    # ------------------------------------------------------------------
    def _apply_fission(self, karyo, node_name):
        cids = [c for c in karyo if len(karyo[c]) >= 4]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        pos = self.rng.randint(2, len(genes) - 2)
        new_id = cid + "f"
        while new_id in karyo:
            new_id += "f"
        karyo[new_id] = genes[pos:]
        karyo[cid] = genes[:pos]
        self.events.append({"node": node_name, "type": "fission",
                            "chroms": [cid, new_id], "pos": "split {} at {}".format(cid, pos)})
        return True

    # ------------------------------------------------------------------
    # Gene gain (de novo)
    # ------------------------------------------------------------------
    def _apply_gene_gain(self, karyo, node_name):
        if not karyo:
            return False
        cid = self.rng.choice(list(karyo.keys()))
        genes = karyo[cid]
        new_gid = self.tracker.new_gain_gene(species=None)
        orient = self.rng.choice(["+", "-"])
        pos = self.rng.randint(0, len(genes))
        genes.insert(pos, (new_gid, orient))
        self.events.append({"node": node_name, "type": "gene_gain",
                            "chrom": cid, "gene": new_gid, "pos": pos})
        return True

    # ------------------------------------------------------------------
    # Tandem duplication
    # ------------------------------------------------------------------
    def _apply_tandem_dup(self, karyo, node_name):
        cids = [c for c in karyo if len(karyo[c]) >= 1]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        if len(genes) <= 2 or self.rng.random() < 0.75:
            seg_len = 1
        elif self.rng.random() < 0.8:
            seg_len = 2
        else:
            seg_len = min(3, len(genes))
        start = self.rng.randint(0, len(genes) - seg_len)
        segment = genes[start:start + seg_len]
        dup_segment = []
        for gid, orient in segment:
            new_gid = self.tracker.new_duplicate_gene(gid, "tandem")
            dup_segment.append((new_gid, orient))
        insert_pos = start + seg_len
        genes[insert_pos:insert_pos] = dup_segment
        self.events.append({"node": node_name, "type": "tandem_dup",
                            "chrom": cid, "pos": start,
                            "len": seg_len, "genes": [g for g, _ in segment]})
        return True

    # ------------------------------------------------------------------
    # Dispersed duplication
    # ------------------------------------------------------------------
    def _apply_dispersed_dup(self, karyo, node_name):
        if not karyo:
            return False
        src_cid = self.rng.choice(list(karyo.keys()))
        src_genes = karyo[src_cid]
        if not src_genes:
            return False
        src_idx = self.rng.randint(0, len(src_genes) - 1)
        src_gid, src_orient = src_genes[src_idx]
        new_gid = self.tracker.new_duplicate_gene(src_gid, "dispersed")
        tgt_cid = self.rng.choice(list(karyo.keys()))
        tgt_genes = karyo[tgt_cid]
        tgt_pos = self.rng.randint(0, len(tgt_genes))
        new_orient = src_orient if self.rng.random() < 0.5 else flip_orient(src_orient)
        tgt_genes.insert(tgt_pos, (new_gid, new_orient))
        self.events.append({"node": node_name, "type": "dispersed_dup",
                            "source": src_gid, "source_chrom": src_cid,
                            "target_chrom": tgt_cid, "target_pos": tgt_pos,
                            "new_gene": new_gid})
        return True

    # ------------------------------------------------------------------
    # Unidirectional translocation (small segment, one-way)
    # ------------------------------------------------------------------
    def _apply_unidir_trans(self, karyo, node_name):
        cids = [c for c in karyo if len(karyo[c]) >= 3]
        if len(cids) < 2:
            return False
        src_cid, tgt_cid = self.rng.sample(cids, 2)
        src_genes = karyo[src_cid]
        seg_len = self.rng.randint(1, min(3, len(src_genes) - 2))
        start = self.rng.randint(1, len(src_genes) - seg_len - 1)
        segment = src_genes[start:start + seg_len]
        if self.rng.random() < 0.5:
            segment = reverse_segment(segment)
        del src_genes[start:start + seg_len]
        tgt_genes = karyo[tgt_cid]
        tgt_pos = self.rng.randint(0, len(tgt_genes))
        tgt_genes[tgt_pos:tgt_pos] = segment
        self.events.append({"node": node_name, "type": "unidir_trans",
                            "source_chrom": src_cid, "target_chrom": tgt_cid,
                            "seg_len": seg_len, "pos": start})
        return True

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def _print_summary(self):
        print("\n=== Simulation Summary ===")
        print("Ancestral genes: {}".format(len(self.tracker.ancestor_genes)))
        print("Total genes created: {}".format(self.tracker.total_genes()))
        print("Gene gains: {}".format(len(self.tracker.gain_genes)))
        print("Species simulated: {}".format(len(self.leaf_karyotypes)))
        for sp in sorted(self.leaf_karyotypes):
            karyo = self.leaf_karyotypes[sp]
            n_chr = len(karyo)
            n_genes = sum(len(g) for g in karyo.values())
            print("  {}: {} chromosomes, {} genes".format(sp, n_chr, n_genes))
        event_counts = defaultdict(int)
        for e in self.events:
            event_counts[e["type"]] += 1
        print("Events:")
        for etype, count in sorted(event_counts.items()):
            print("  {}: {}".format(etype, count))

    # ==================================================================
    # Output generation
    # ==================================================================
    def generate_outputs(self, outdir, tree, ploidy_map, orig_nw):
        os.makedirs(outdir, exist_ok=True)
        species = sorted(self.leaf_karyotypes.keys())
        print("\nGenerating output files in {}/...".format(outdir))

        # Register species for all genes in leaf karyotypes
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    self.tracker.set_species(gid, sp)

        # Write tree first -- this also renumbers internal nodes
        self._write_tree(outdir, orig_nw, tree)

        self._write_leaf_gff(outdir, species)
        self._write_ancestor_gff(outdir, tree)
        self._write_og(outdir, species)
        self._write_orthologs(outdir, species)
        self._write_events(outdir)
        self._write_true_hogs(outdir, tree, species)
        self._write_ancestor_karyotypes(outdir, tree)

        print("Done! Output files:")
        for f in sorted(os.listdir(outdir)):
            print("  {}/{}".format(outdir, f))

    # ------------------------------------------------------------------
    # Leaf GFF
    # ------------------------------------------------------------------
    def _write_leaf_gff(self, outdir, species):
        path = os.path.join(outdir, "all_species_gene.gff")
        with open(path, 'w') as f:
            for sp in species:
                karyo = self.leaf_karyotypes[sp]
                sorted_cids = self._sorted_chr_ids(karyo)
                chrom_idx = 0
                for cid in sorted_cids:
                    chrom_idx += 1
                    for i, (gid, orient) in enumerate(karyo[cid]):
                        start = (i + 1) * 100 - 99
                        end = (i + 1) * 100
                        f.write("Chr{}\t{}|{}\t{}\t{}\t{}\n".format(
                            chrom_idx, sp, gid, start, end, orient))

    # ------------------------------------------------------------------
    # Ancestor GFF (ALL nodes for validation)
    # ------------------------------------------------------------------
    def _write_ancestor_gff(self, outdir, tree):
        path = os.path.join(outdir, "ancestors_gene.gff")
        with open(path, 'w') as f:
            for node in tree.traverse("preorder"):
                if node.name not in self.all_node_karyotypes:
                    continue
                karyo = self.all_node_karyotypes[node.name]
                sorted_cids = self._sorted_chr_ids(karyo)
                chrom_idx = 0
                for cid in sorted_cids:
                    chrom_idx += 1
                    for i, (gid, orient) in enumerate(karyo[cid]):
                        start = (i + 1) * 100 - 99
                        end = (i + 1) * 100
                        f.write("Chr{}\t{}|{}\t{}\t{}\t{}\n".format(
                            chrom_idx, node.name, gid, start, end, orient))

    # ------------------------------------------------------------------
    # Ancestor karyotypes (human-readable)
    # ------------------------------------------------------------------
    def _write_ancestor_karyotypes(self, outdir, tree):
        path = os.path.join(outdir, "ancestors_karyotypes.txt")
        with open(path, 'w') as f:
            for node in tree.traverse("preorder"):
                if node.name not in self.all_node_karyotypes:
                    continue
                karyo = self.all_node_karyotypes[node.name]
                f.write(">{}\t{} chromosomes\n".format(node.name, len(karyo)))
                sorted_cids = self._sorted_chr_ids(karyo)
                for cid in sorted_cids:
                    genes = karyo[cid]
                    gene_str = " ".join(
                        ("-" + g if o == "-" else g) for g, o in genes)
                    f.write("{}\t{}\n".format(cid, gene_str))

    # ------------------------------------------------------------------
    # OG / SOG output
    # ------------------------------------------------------------------
    def _write_og(self, outdir, species):
        """SOG file: group all genes by their ultimate ancestral gene.
        Gene gains get their own SOGs."""
        anc_to_genes = defaultdict(list)
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    anc = self.tracker.get_ancestor(gid)
                    full_name = "{}|{}".format(sp, gid)
                    anc_to_genes[anc].append(full_name)

        path = os.path.join(outdir, "ortholog_groups.txt")
        with open(path, 'w') as f:
            for i, (anc, gene_sp_list) in enumerate(sorted(anc_to_genes.items()), 1):
                gene_names = sorted(gene_sp_list)
                f.write("SOG{}: {}\n".format(i, " ".join(gene_names)))

    # ------------------------------------------------------------------
    # Ortholog pairs (only true orthologs, not cross-WGD-copy paralogs)
    # ------------------------------------------------------------------
    def _write_orthologs(self, outdir, species):
        """Output true ortholog pairs only."""
        sp_gene_groups = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_gene_groups[sp][(anc, wgd_key)].append(full_name)

        path = os.path.join(outdir, "ortholog_pairs.txt")
        with open(path, 'w') as f:
            sp_list = sorted(species)
            for i in range(len(sp_list)):
                for j in range(i + 1, len(sp_list)):
                    sp1, sp2 = sp_list[i], sp_list[j]
                    groups1 = sp_gene_groups[sp1]
                    groups2 = sp_gene_groups[sp2]
                    common_keys = set(groups1.keys()) & set(groups2.keys())
                    for key in common_keys:
                        for g1 in sorted(groups1[key]):
                            for g2 in sorted(groups2[key]):
                                f.write("{}\t{}\n".format(g1, g2))

            self._write_cross_wgd_orthologs(f, species, sp_gene_groups)

    def _write_cross_wgd_orthologs(self, fout, species, sp_gene_groups):
        """Write ortholog pairs between species at different WGD levels."""
        sp_anc_genes = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_genes[sp][anc].append((full_name, wgd_key))

        sp_list = sorted(species)
        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                common_ancs = set(sp_anc_genes[sp1].keys()) & set(sp_anc_genes[sp2].keys())
                for anc in common_ancs:
                    genes1 = sp_anc_genes[sp1][anc]
                    genes2 = sp_anc_genes[sp2][anc]
                    for g1, wk1 in genes1:
                        for g2, wk2 in genes2:
                            if wk1 == wk2:
                                continue
                            shared_wgd_nodes = set(k for k, _ in wk1) & set(k for k, _ in wk2)
                            compatible = all(
                                dict(wk1).get(n) == dict(wk2).get(n)
                                for n in shared_wgd_nodes
                            )
                            if compatible:
                                fout.write("{}\t{}\n".format(
                                    min(g1, g2), max(g1, g2)))

    # ------------------------------------------------------------------
    # Tree output
    # ------------------------------------------------------------------
    def _write_tree(self, outdir, orig_nw, tree):
        """Write tree file compatible with AK.py's number_nodes."""
        path = os.path.join(outdir, "species_tree.nwk")
        nw = orig_nw.strip()
        if not nw.endswith(';'):
            nw += ';'
        # Remove internal node names: )Name: → ):
        nw = re.sub(
            r'\)([A-Za-z_]\w*):',
            r'):',
            nw)
        with open(path, 'w') as f:
            f.write(nw + '\n')

        self._renumber_internal_nodes(tree)

    # ------------------------------------------------------------------
    # Renumber internal nodes (matching AK.py's number_nodes)
    # ------------------------------------------------------------------
    def _renumber_internal_nodes(self, tree):
        """Renumber internal nodes to match AK.py's number_nodes output."""
        old_to_new = {}
        i = 0
        for node in tree.traverse():  # default = preorder
            if node.is_leaf():
                continue
            old_name = node.name
            new_name = "N{}".format(i)
            old_to_new[old_name] = new_name
            node.name = new_name
            i += 1
        new_karyotypes = {}
        for old_name, karyo in self.all_node_karyotypes.items():
            new_name = old_to_new.get(old_name, old_name)
            new_karyotypes[new_name] = karyo
        self.all_node_karyotypes = new_karyotypes
        for evt in self.events:
            if evt.get("node") in old_to_new:
                evt["node"] = old_to_new[evt["node"]]

    # ------------------------------------------------------------------
    # Events log
    # ------------------------------------------------------------------
    def _write_events(self, outdir, path="events.tsv"):
        evt_path = os.path.join(outdir, path)
        with open(evt_path, 'w') as f:
            f.write("node\tevent_type\tdetails\n")
            for e in self.events:
                details = {k: v for k, v in e.items() if k not in ("node", "type")}
                f.write("{}\t{}\t{}\n".format(e['node'], e['type'], details))

    # ------------------------------------------------------------------
    # True HOG output (for validation)
    # ------------------------------------------------------------------
    def _write_true_hogs(self, outdir, tree, species):
        """Build true HOG hierarchy using ortholog graph + tree structure."""
        tree_copy = tree.copy("deepcopy")
        name_map = number_tree_nodes(tree_copy)

        G = nx.Graph()
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    full_name = "{}|{}".format(sp, gid)
                    G.add_node(full_name)

        sp_anc_genes = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_genes[sp][(anc, wgd_key)].append(full_name)

        sp_list = sorted(species)
        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                groups1 = sp_anc_genes[sp1]
                groups2 = sp_anc_genes[sp2]
                common_keys = set(groups1.keys()) & set(groups2.keys())
                for key in common_keys:
                    for g1 in groups1[key]:
                        for g2 in groups2[key]:
                            G.add_edge(g1, g2)

        sp_anc_only = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_only[sp][anc].append((full_name, wgd_key))

        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                common_ancs = set(sp_anc_only[sp1].keys()) & set(sp_anc_only[sp2].keys())
                for anc in common_ancs:
                    genes1 = sp_anc_only[sp1][anc]
                    genes2 = sp_anc_only[sp2][anc]
                    for g1, wk1 in genes1:
                        for g2, wk2 in genes2:
                            if wk1 == wk2:
                                continue
                            shared = set(k for k, _ in wk1) & set(k for k, _ in wk2)
                            if all(dict(wk1).get(n) == dict(wk2).get(n) for n in shared):
                                G.add_edge(g1, g2)

        og_genes = defaultdict(list)
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    anc = self.tracker.get_ancestor(gid)
                    full_name = "{}|{}".format(sp, gid)
                    og_genes[anc].append(full_name)

        for gid in self.tracker.gain_genes:
            sp = self.tracker.gene_species.get(gid, "unknown")
            full_name = "{}|{}".format(sp, gid)
            og_genes[gid].append(full_name)

        all_hogs = []
        node_species = {}
        for node in tree_copy.traverse("postorder"):
            if node.is_leaf():
                node_species[node.name] = {node.name}
            else:
                node_species[node.name] = set()
                for child in node.children:
                    node_species[node.name] |= node_species[child.name]

        og_counter = 0
        for anc in sorted(og_genes.keys()):
            og_counter += 1
            og_id = "SOG{}".format(og_counter)
            og_gene_set = set(og_genes[anc])

            node_hogs = {}
            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                sp_set = node_species[nid]
                node_genes = set()
                for g in og_gene_set:
                    sp = g.split("|")[0]
                    if sp in sp_set:
                        node_genes.add(g)
                if not node_genes:
                    continue
                subG = G.subgraph(node_genes)
                ccs = list(nx.connected_components(subG))
                hog_list = []
                for idx, cc in enumerate(sorted(ccs, key=lambda x: sorted(x)[0] if x else "")):
                    hog_id = "{}.{}.hog{}".format(og_id, nid, idx)
                    hog_list.append((hog_id, cc))
                node_hogs[nid] = hog_list

            hog_parent = {}
            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                if nid not in node_hogs:
                    continue
                if node.up is None:
                    for hog_id, cc in node_hogs[nid]:
                        hog_parent[hog_id] = "Root"
                else:
                    parent_nid = node.up.name
                    if parent_nid not in node_hogs:
                        for hog_id, cc in node_hogs[nid]:
                            hog_parent[hog_id] = "Root"
                    else:
                        for hog_id, cc in node_hogs[nid]:
                            found_parent = False
                            for parent_hog_id, parent_cc in node_hogs[parent_nid]:
                                if cc.issubset(parent_cc):
                                    hog_parent[hog_id] = parent_hog_id
                                    found_parent = True
                                    break
                            if not found_parent:
                                hog_parent[hog_id] = "Root"

            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                if nid not in node_hogs:
                    continue
                for hog_id, cc in node_hogs[nid]:
                    all_hogs.append({
                        "hog_id": hog_id,
                        "og_id": og_id,
                        "node_id": nid,
                        "parent": hog_parent.get(hog_id, "Root"),
                        "genes": sorted(cc)
                    })

        path = os.path.join(outdir, "true_hogs.tsv")
        with open(path, 'w') as f:
            f.write("{}\n".format("\t".join(["HOG", "OG", "Tree_Node", "Parent", "Genes"])))
            for hog in all_hogs:
                gene_str = " ".join(hog["genes"])
                f.write("{}\n".format("\t".join([
                    hog["hog_id"], hog["og_id"], hog["node_id"],
                    hog["parent"], gene_str
                ])))

    # ------------------------------------------------------------------
    # Helper: sort chromosome IDs
    # ------------------------------------------------------------------
    def _sorted_chr_ids(self, karyo):
        def sort_key(c):
            m = re.search(r'(\d+)', c)
            return (0, int(m.group(1))) if m else (1, c)
        return sorted(karyo.keys(), key=sort_key)


# ============================================================================
# CLI argument definition (for options.py integration)
# ============================================================================

def sim_args(parser):
    """Define arguments for the sim subcommand."""
    parser.add_argument('-t', '-tree', required=True, type=str,
                        dest='tree', metavar='FILE',
                        help='Newick tree with branch lengths [required]')
    parser.add_argument('-o', '-outdir', required=True, type=str,
                        dest='outdir', metavar='DIR',
                        help='Output directory [required]')

    g_anc = parser.add_argument_group('Ancestor genome')
    g_anc.add_argument('--num-chroms', type=int, default=8,
                       help='Number of ancestral chromosomes [default=%(default)s]')
    g_anc.add_argument('--min-genes', type=int, default=200,
                       help='Min genes per ancestral chromosome [default=%(default)s]')
    g_anc.add_argument('--max-genes', type=int, default=1000,
                       help='Max genes per ancestral chromosome [default=%(default)s]')

    g_rate = parser.add_argument_group('Rates (per unit branch length)')
    g_rate.add_argument('--inv-rate', type=float, default=10.0,
                        help='Inversion rate [default=%(default)s]')
    g_rate.add_argument('--rt-rate', type=float, default=3.0,
                        help='Reciprocal translocation rate [default=%(default)s]')
    g_rate.add_argument('--ncf-rate', type=float, default=1.5,
                        help='Non-centric fusion rate [default=%(default)s]')
    g_rate.add_argument('--eej-rate', type=float, default=2.0,
                        help='End-end joining rate [default=%(default)s]')
    g_rate.add_argument('--fission-rate', type=float, default=0.05,
                        help='Fission rate (very low) [default=%(default)s]')
    g_rate.add_argument('--wgd-rate', type=float, default=0.05,
                        help='WGD rate (2x most common, 3x less, etc.) [default=%(default)s]')
    g_rate.add_argument('--gene-gain-rate', type=float, default=0.5,
                        help='Gene gain rate [default=%(default)s]')
    g_rate.add_argument('--tandem-dup-rate', type=float, default=2.0,
                        help='Tandem duplication rate [default=%(default)s]')
    g_rate.add_argument('--dispersed-dup-rate', type=float, default=0.3,
                        help='Dispersed duplication rate [default=%(default)s]')
    g_rate.add_argument('--unidir-trans-rate', type=float, default=1.0,
                        help='Unidirectional translocation rate [default=%(default)s]')
    g_rate.add_argument('--frac-rate', type=float, default=5.0,
                        help='Fractionation rate after WGD (probability per gene) [default=%(default)s]')
    g_rate.add_argument('--scale', type=float, default=1.0,
                        help='Global rate scaling factor [default=%(default)s]')

    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed [default=%(default)s]')


def xmain(**kargs):
    """Main entry point for the sim subcommand."""
    tree, ploidy_map, orig_nw = parse_tree(kargs['tree'], kargs.get('wgd_rate', 0.05))
    if ploidy_map:
        print("Ploidy annotations found: {}".format(ploidy_map))
    if kargs.get('wgd_rate', 0.05) > 0:
        print("WGD rate: {} (rate-based, factor sampled: 2x most common)".format(
            kargs['wgd_rate']))

    sim = EvolutionSimulator(
        seed=kargs.get('seed', 42),
        num_chroms=kargs.get('num_chroms', 8),
        min_genes=kargs.get('min_genes', 200),
        max_genes=kargs.get('max_genes', 1000),
        inv_rate=kargs.get('inv_rate', 10.0),
        rt_rate=kargs.get('rt_rate', 3.0),
        ncf_rate=kargs.get('ncf_rate', 1.5),
        eej_rate=kargs.get('eej_rate', 2.0),
        fission_rate=kargs.get('fission_rate', 0.05),
        wgd_rate=kargs.get('wgd_rate', 0.05),
        gene_gain_rate=kargs.get('gene_gain_rate', 0.5),
        tandem_dup_rate=kargs.get('tandem_dup_rate', 2.0),
        dispersed_dup_rate=kargs.get('dispersed_dup_rate', 0.3),
        unidir_trans_rate=kargs.get('unidir_trans_rate', 1.0),
        frac_rate=kargs.get('frac_rate', 5.0),
        scale=kargs.get('scale', 1.0),
    )

    sim.run(tree, ploidy_map)
    sim.generate_outputs(kargs['outdir'], tree, ploidy_map, orig_nw)


# ============================================================================
# Standalone CLI (for direct execution)
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Chromosome Rearrangement Evolution Simulator (AK.py format)",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sim_args(parser)
    args = parser.parse_args()
    xmain(**vars(args))


if __name__ == "__main__":
    main()
