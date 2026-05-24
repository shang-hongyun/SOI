#!/usr/bin/env python3
"""用模拟真值逐步验证重建管线。直接比较图结构，不做 benchmark 评估。"""
import sys, os, csv, logging
from collections import defaultdict, Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
logging.basicConfig(level=logging.WARNING)

from soi.takr_colored_graph import ColoredGraph
from soi.AK import AKR
from soi.evolution_simulator_ak import parse_tree

SIM = os.path.join(os.path.dirname(__file__), 'benchmark_run', 'sim_data')


def load_truth_karyotypes(path):
    karyo = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current = line.lstrip('>').split('\t')[0]
                karyo[current] = []
            elif current and line:
                chrom = line.split('\t')[-1].split()
                karyo[current].append(chrom)
    return karyo


def load_truth_events(path):
    events = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            events[row['branch']].append({
                'type': row['event_type'],
                'genes': row.get('genes_involved', ''),
            })
    return events


def gene_adjacencies(chroms):
    adjs = set()
    for chrom in chroms:
        for i in range(len(chrom) - 1):
            a, b = chrom[i].lstrip('-+'), chrom[i+1].lstrip('-+')
            if a and b:
                adjs.add((min(a, b), max(a, b)))
    return adjs


def test_N1_step_by_step():
    print("=" * 60)
    print("  STEP-BY-STEP VERIFICATION: N1")
    print("=" * 60)

    truth_karyo = load_truth_karyotypes(os.path.join(SIM, 'ancestors_karyotypes.txt'))
    truth_events = load_truth_events(os.path.join(SIM, 'events.tsv'))

    n1_truth = truth_karyo['N1']
    n1_truth_adjs = gene_adjacencies(n1_truth)
    n1_chrom_sizes = [len(c) for c in n1_truth]
    print(f"\n[TRUTH] N1: {len(n1_truth)} chroms, sizes={n1_chrom_sizes}")
    print(f"[TRUTH] N1 adjacencies: {len(n1_truth_adjs)}")

    n1_events = truth_events.get('N1-Sp_1', []) + truth_events.get('N1-Sp_4', [])
    event_types = Counter(e['type'] for e in n1_events)
    print(f"[TRUTH] N1 events (Sp_1+Sp_4): {dict(event_types)}")

    # 加载 AKR
    outdir = 'tests/benchmark_run/recon_diag4'
    os.makedirs(outdir, exist_ok=True)
    akr = AKR(
        ogfile=os.path.join(SIM, 'ortholog_groups.txt'),
        orthfiles=[os.path.join(SIM, 'ortholog_pairs.txt')],
        gfffile=os.path.join(SIM, 'all_species_gene.gff'),
        sptreefile=os.path.join(SIM, 'species_tree.nwk'),
        outpre=os.path.join(outdir, 'AKR'),
        reconstruction_algorithm='v4_colored',
        min_genes=0, timeout=600,
    )
    akr._build_hogs()
    tree, _, _ = parse_tree(os.path.join(SIM, 'species_tree.nwk'))
    akr._build_leaf_graphs()

    # 关键步骤：用 _map_to_parent_hogs 映射到 HOG 级别
    G = ColoredGraph(hog_level='N1')
    hog_level = 'N1'
    for cname in ['Sp_1', 'Sp_4']:
        if cname in akr.leaf_graphs:
            mapped = akr._map_to_parent_hogs(hog_level, akr.leaf_graphs[cname], source_id=cname)
            n = G.add_child(cname, mapped)
            print(f"[ADD] {cname}: {n} chroms (HOG-mapped)")

    print(f"\n[BUILT] N1 ColoredGraph: {G.node_count()} nodes, {G.edge_count()} edges")
    print(f"[BUILT] Shared: {len(G.shared_edges())}, Unique: {len(G.unique_edges())}")

    hog_set = G.all_hogs()
    print(f"[BUILT] HOG count: {len(hog_set)}")

    # Phase 2: indels
    print(f"\n--- Phase 2: Indel Resolution ---")
    n_before = G.edge_count()
    G.resolve_indels()
    n_after = G.edge_count()
    n_events = len(G.events)
    print(f"  Edges: {n_before} → {n_after} (removed {n_before - n_after})")
    print(f"  Events: {n_events}")
    spanning = G.find_indel_shortcuts()
    print(f"  Postcondition: {len(spanning)} spanning edges remain {'WARN' if spanning else 'OK'}")

    # Phase 3: blocks
    print(f"\n--- Phase 3: Block Compression ---")
    G._build_synteny_blocks()
    G._compress_to_block_level()
    multi = [b for b, hogs in G._blocks.items() if len(hogs) >= 2]
    single = [b for b, hogs in G._blocks.items() if len(hogs) == 1]
    print(f"  Blocks: {len(G._blocks)} (multi={len(multi)}, singleton={len(single)})")

    bg = G._block_graph
    n_unique = sum(1 for _, _, d in bg.edges(data=True) if len(d.get('colors', set())) == 1)
    n_shared = sum(1 for _, _, d in bg.edges(data=True) if len(d.get('colors', set())) > 1)
    print(f"  Block graph: {bg.number_of_nodes()} nodes, {bg.number_of_edges()} edges ({n_unique} unique, {n_shared} shared)")

    import networkx as nx
    shared_bg = nx.Graph()
    for b1, b2, data in bg.edges(data=True):
        if len(data.get('colors', set())) > 1:
            shared_bg.add_edge(b1, b2)
    components = list(nx.connected_components(shared_bg))
    print(f"  Shared components: {len(components)} (truth chroms: {len(n1_truth)})")

    # Phase 4c-4e: structural
    print(f"\n--- Phase 4c-4e: Structural Events ---")
    G._save_original_shared_components()
    n_before_e = len(G.events)
    G._resolve_block_structural_events()
    struct_events = G.events[n_before_e:]
    struct_types = Counter(e.event_type for e in struct_events)
    print(f"  Events: {len(struct_events)} {dict(struct_types)}")

    # Phase 4f: bridge
    print(f"\n--- Phase 4f: Bridge Events ---")
    n_before_e = len(G.events)
    G._resolve_block_bridge_events()
    bridge_events = G.events[n_before_e:]
    bridge_types = Counter(e.event_type for e in bridge_events)
    print(f"  Events: {len(bridge_events)} {dict(bridge_types)}")

    # Phase 5: path cover
    print(f"\n--- Phase 5: Path Cover ---")
    paths = G.path_cover()
    path_sizes = sorted([len(p) for p in paths], reverse=True)
    print(f"  Paths: {len(paths)} (truth: {len(n1_truth)})")
    print(f"  Path sizes (top10): {path_sizes[:10]}")
    covered = set()
    for p in paths:
        covered.update(p)
    uncovered = len(hog_set) - len(covered)
    print(f"  Uncovered HOGs: {uncovered}")

    # 最终对比
    print(f"\n{'=' * 60}")
    print(f"  RESULT: N1 truth={len(n1_truth)} chroms, reconstructed={len(paths)} chroms")
    if len(paths) == len(n1_truth):
        print(f"  OK CHROMOSOME COUNT MATCH")
    else:
        print(f"  FAIL CHROMOSOME COUNT MISMATCH (diff={len(paths) - len(n1_truth)})")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    test_N1_step_by_step()
