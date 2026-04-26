#!/usr/bin/env python3
"""Run v2 reconstruction directly, bypassing AKR._build_hogs() pipeline.

Usage: python3.11 run_v2.py --sim-dir tests/sim_data/ --out tests/v2_test/
"""

import argparse
import logging
import os
import sys
from collections import defaultdict

# Setup
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim-dir', required=True, help='Simulator output directory')
    parser.add_argument('--out', required=True, help='Output directory')
    parser.add_argument('--min-hogs', type=int, default=3)
    args = parser.parse_args()

    sim_dir = os.path.abspath(args.sim_dir)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    ogfile = os.path.join(sim_dir, 'ortholog_groups.txt')
    gfffile = os.path.join(sim_dir, 'all_species_gene.gff')
    treefile = os.path.join(sim_dir, 'species_tree.nwk')

    # Build AKR instance with leaf graphs pre-loaded
    from soi.AK import AKR
    from soi.takr_event_driven import reconstruct_event_driven_v2

    akr = AKR(
        ogfile=ogfile,
        gfffile=gfffile,
        sptreefile=treefile,
        outpre=os.path.join(out_dir, 'AKR'),
        use_v4=False, use_v3=False,
        min_genes=0, timeout=600,
    )

    # Build HOG tree directly (skip mcscan pipeline)
    logger.info("Building HOG data structure...")
    _build_hogs_direct(akr)
    logger.info("  HOGs per node: %s", {k: len(v) for k, v in akr.hogs_by_node.items()})

    # Build leaf graphs
    logger.info("Building leaf graphs...")
    akr._build_leaf_graphs()
    logger.info("  Leaf graphs: %s", list(akr.leaf_graphs.keys()))

    # Run v2 reconstruction
    logger.info("Running v2 reconstruction...")
    anc_graphs = reconstruct_event_driven_v2(akr, min_hogs=args.min_hogs)

    # Export events
    export_events(akr, os.path.join(out_dir, 'AKR.events.tsv'))
    logger.info("Done! Ancestors: %s", list(anc_graphs.keys()))


def _build_hogs_direct(akr):
    """Build HOG tree directly from ortholog groups file."""
    from ete3 import Tree
    from collections import defaultdict

    # Parse tree
    tree = Tree(akr.sptreefile, format=1)
    from soi.tree import number_nodes
    tree = number_nodes(tree)
    akr.tree = tree

    # Initialize HOG storage
    akr.hogs_by_node = defaultdict(list)
    akr.gene_to_hog = {}
    akr.node_gene_to_hog = defaultdict(dict)
    akr.gene_to_all_hogs = defaultdict(list)
    akr.hog_species = defaultdict(set)
    akr.hog_parent = {}
    akr.hog_children = defaultdict(list)
    akr._hog_node_cache = {}

    # Read ortholog groups
    from .OrthoFinder import OrthoMCLGroup
    all_hogs = {}
    for og in OrthoMCLGroup(akr.ogfile):
        og_id = og.ogid
        og_genes = og.genes
        og_species = set(og.species)

        # For simulated data, all genes in a SOG map to the same ancestor HOG
        # Create a simple HOG structure
        hog_id = og_id
        node_id = _find_hog_node(tree, og_species)

        rec = SimpleHOG(hog_id, og_genes, node_id, og_species)
        all_hogs[hog_id] = rec
        akr.hogs_by_node[node_id].append(rec)

        for gene in og_genes:
            akr.gene_to_hog[gene] = hog_id
            akr.node_gene_to_hog[node_id][gene] = hog_id
            akr.gene_to_all_hogs[gene].append(hog_id)
            if '|' in gene:
                akr.hog_species[hog_id].add(gene.split('|')[0])
            else:
                akr.hog_species[hog_id].add('unknown')

    akr.hog = SimpleHOGWrapper(all_hogs)
    logger.info("  Indexed %d genes into %d HOGs",
                len(akr.gene_to_hog), len(all_hogs))


def _find_hog_node(tree, species_set):
    """Find the lowest node that contains all given species."""
    species = list(species_set)
    if len(species) <= 1:
        return species[0] if species else 'root'
    # Find MRCA
    nodes = [tree.search_nodes(name=sp)[0] for sp in species]
    mrca = tree.get_common_ancestor(nodes)
    return mrca.name if mrca else 'root'


class SimpleHOG:
    """Minimal HOG record for testing."""
    def __init__(self, hog_id, genes, node_id, species):
        self.hog_id = hog_id
        self.genes = genes
        self.node_id = node_id
        self.species = species


class SimpleHOGWrapper:
    """Minimal HOG object that has an all_hogs dict."""
    def __init__(self, all_hogs):
        self.all_hogs = all_hogs
        self.tree = None


def export_events(akr, outpath):
    """Export events to TSV."""
    import csv
    events = []
    for node_id, aag in akr.anc_graphs.items():
        for e in getattr(aag, 'events', []):
            events.append(e)

    with open(outpath, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['branch', 'event_type', 'genes', 'chroms', 'desc', 'support'])
        for e in events:
            branch = getattr(e, 'branch', akr.anc_graphs[0].node_id if akr.anc_graphs else '')
            genes = ','.join(str(g) for g in e.genes_involved) if e.genes_involved else ''
            writer.writerow([
                branch,
                e.event_type,
                genes,
                '',
                e.desc,
                e.support,
            ])


if __name__ == '__main__':
    main()
