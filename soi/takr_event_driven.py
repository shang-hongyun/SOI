#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""takr_event_driven.py - Event-driven ancestral karyotype reconstruction.

Core algorithm for TAKR v4. Called by AK.run().

Algorithm:
  1. Postorder traverse tree
  2. For each ancestor node:
     a. Map children to ancestor HOG level
     b. Build consensus graph: shared adjacencies = ancestral state
     c. Detect small-scale events (inversion, unidir_trans) per child vs consensus
     d. Detect large-scale events (EEJ, NCF, fission, RT/URT) per child vs consensus
  3. Consensus = ancestor. Events are the differences per child branch.

Key insight: The consensus of children IS the ancestor.
Events are differences between each child and the consensus.
"""

import itertools
import logging
import os
import time
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

logger = logging.getLogger(__name__)


# =========================================================================
#  Helpers
# =========================================================================

def _graph_stats(graph):
    """Return (nodes, edges, connected_components). Auto-detect directed/undirected."""
    n = graph.number_of_nodes()
    e = graph.number_of_edges()
    if graph.is_directed():
        cc = nx.number_weakly_connected_components(graph)
    else:
        cc = nx.number_connected_components(graph)
    return (n, e, cc)


def _validate_dedup(mc, cid, expected_chrom_count):
    """Validate graph after dedup. Returns (ok, errors)."""
    errors = []
    ch = getattr(mc, 'chrom_hogs', None)
    if not ch:
        errors.append("no chrom_hogs found")
        return (False, errors)
    n_chrom = len(ch)
    if n_chrom != expected_chrom_count:
        errors.append(f"chrom count changed {expected_chrom_count} → {n_chrom}")
    # No duplicate HOGs (within each chromosome)
    for ci, hogs in ch.items():
        gene_hogs = [str(h) for h in hogs if h not in mc.telomeres]
        seen = set()
        for h in gene_hogs:
            if h in seen:
                errors.append(f"chrom{ci}: duplicate HOG {h}")
            seen.add(h)
    # No self-loops (consecutive identical HOGs)
    for ci, hogs in ch.items():
        gene_hogs = [str(h) for h in hogs if h not in mc.telomeres]
        for i in range(len(gene_hogs) - 1):
            if gene_hogs[i] == gene_hogs[i + 1]:
                errors.append(f"chrom{ci}: self-loop at pos {i} ({gene_hogs[i]})")
    return (len(errors) == 0, errors)


def _log_child_events(events, phase_label, child_id=None):
    """Unified event log format, grouped by child."""
    child_events = defaultdict(list)
    for e in events:
        parts = e.branch.split('-')
        cid = parts[-1] if parts else 'unknown'
        if child_id and cid != child_id:
            continue
        child_events[cid].append(e)

    for cid in sorted(child_events.keys()):
        evts = child_events[cid]
        type_counts = defaultdict(int)
        type_lens = defaultdict(list)
        for e in evts:
            type_counts[e.event_type] += 1
            type_lens[e.event_type].append(len(e.genes_involved))
        parts = []
        for etype in sorted(type_counts.keys()):
            count = type_counts[etype]
            lens = type_lens[etype]
            if len(lens) == 1:
                parts.append(f"{etype}={count} (len {lens[0]})")
            else:
                parts.append(f"{etype}={count} (len {min(lens)}-{max(lens)})")
        logger.info("  [%s] %s events: %s", phase_label, cid, ", ".join(parts))


# =========================================================================
#  ReconstructorV2 class
# =========================================================================

class ReconstructorV2:
    """Orchestrates event-driven ancestral karyotype reconstruction.

    Encapsulates: tree traversal, WGD collapse, outgroup computation,
    dedup/merge, event resolution, and visualization.
    """

    def __init__(self, akr, min_hogs=3, gfa_debug=False):
        self.akr = akr
        self.min_hogs = min_hogs
        self.gfa_debug = gfa_debug
        self.outgroup_leaves_map: Dict[str, List[str]] = {}

    # ── Entry point ──────────────────────────────────────────────────

    def run(self) -> Dict:
        """Main entry point. Returns {node_id: AncestralAdjacencyGraph}."""
        logger.info("=== Event-driven reconstruction v2 (ColoredGraph) ===")
        t0 = time.time()

        self.outgroup_leaves_map = self._build_outgroup_map()

        for node in self.akr.tree.traverse(strategy="postorder"):
            self._process_node(node)

        logger.info("=== Event-driven reconstruction v2 done (%.1fs) ===",
                     time.time() - t0)
        return self.akr.anc_graphs

    # ── Outgroup precomputation ──────────────────────────────────────

    def _build_outgroup_map(self) -> Dict[str, List[str]]:
        """For each internal node, find outgroup (sibling subtree leaves)."""
        outgroup_map = {}
        for node in self.akr.tree.traverse(strategy="postorder"):
            if node.is_leaf() or not node.up:
                continue
            node_id = node.name
            sibling_leaves = []
            for sibling in node.up.children:
                if sibling.name == node_id:
                    continue
                for leaf in sibling.iter_leaves():
                    sibling_leaves.append(leaf.name)
            if sibling_leaves:
                outgroup_map[node_id] = sibling_leaves
        return outgroup_map

    # ── Tree traversal: per-node dispatch ────────────────────────────

    def _process_node(self, node):
        """Process one tree node — dispatch by leaf/internal/polyploid."""
        node_id = node.name
        ploidy = self.akr.ploidy_map.get(node_id, 1)

        if node.is_leaf():
            if ploidy > 1 and node_id in self.akr.leaf_graphs \
                    and node_id not in self.akr.pre_wgd_graphs:
                self._collapse_leaf_wgd(node, node_id, ploidy)
            return

        # Internal node — full reconstruction
        self._reconstruct_internal_node(node, node_id, ploidy)

    def _collapse_leaf_wgd(self, node, node_id, ploidy):
        """Collapse polyploid leaf to pre-WGD: HOG map, build ColoredGraph, collapse."""
        from .takr_colored_graph import ColoredGraph

        leaf_nodes = self.akr.tree.search_nodes(name=node_id)
        if not leaf_nodes or not leaf_nodes[0].up:
            return
        parent_id = leaf_nodes[0].up.name
        mapped = self.akr._map_to_parent_hogs(
            parent_id, self.akr.leaf_graphs[node_id], source_id=node_id)
        G_leaf = ColoredGraph(hog_level=node_id)
        G_leaf.add_child(node_id, mapped)
        pre_anc, _ = self._wgd_collapse(G_leaf, node_id, ploidy, parent_id)
        self.akr.pre_wgd_graphs[node_id] = pre_anc

    # ── Internal node reconstruction ─────────────────────────────────

    def _reconstruct_internal_node(self, node, node_id, ploidy):
        """Full reconstruction for one internal node."""
        from .takr_colored_graph import ColoredGraph

        logger.info("Reconstructing node %s [v2 ColoredGraph]", node_id)
        t_node = time.time()

        # 1. Collect child graphs
        child_graphs, child_source_ids = self._collect_child_graphs(node)
        if len(child_graphs) < 2:
            return

        self._log_children(child_graphs, child_source_ids)

        # 2. Map to ancestor HOG level
        mapped_children = [
            self.akr._map_to_parent_hogs(node_id, cg, source_id=cid)
            for cg, cid in zip(child_graphs, child_source_ids)
        ]

        # 2b. 输出子图 GFA（dedup 之前，block 压缩后）
        if self.gfa_debug:
            self._output_child_gfas(mapped_children, child_source_ids, node_id)

        # 3. Phase 1: dedup + merge
        G = ColoredGraph(hog_level=node_id)
        deduped_children = self._phase1_dedup_and_merge(
            G, mapped_children, child_source_ids, node_id)

        # 4. Viz: raw graph
        self._viz_raw_graph(G, node_id)
        # 5. Viz: child paths
        self._viz_child_paths(deduped_children, child_source_ids, node_id)

        # 6. Outgroup info
        outgroup_adjacency, outgroup_hogs = self._collect_outgroup_info(
            node_id, child_source_ids)

        # 7. Event resolution
        G.resolve_all_events(
            outgroups=outgroup_hogs,
            outgroup_adjacency=outgroup_adjacency,
            min_hogs=self.min_hogs,
            gfa_debug=self.gfa_debug,
            gfa_prefix=self._gfa_prefix(node_id),
        )

        # 8. Viz: resolved graph
        self._viz_resolved_graph(G, node_id)

        # 9. Store ancestor
        ancestor = G.to_ancestral_graph()
        self.akr.anc_graphs[node_id] = ancestor
        n_chrom = len(list(ancestor.chromosomes))
        logger.info("  Done: %d chroms, %d events (%.1fs)",
                     n_chrom, len(ancestor.events), time.time() - t_node)

        # 10. WGD collapse if polyploid
        if ploidy > 1:
            parent_name = node.up.name if node.up else node_id
            pre_anc, events = self._wgd_collapse(G, node_id, ploidy, parent_name)
            self.akr.pre_wgd_graphs[node_id] = pre_anc
            virtual_branch = f"{node_id}_preWGD-{node_id}"
            for e in pre_anc.events:
                e.branch = virtual_branch

    # ── Child graph collection ───────────────────────────────────────

    def _collect_child_graphs(self, node) -> Tuple[List, List]:
        """Collect child graphs for an internal node."""
        child_graphs, child_source_ids = [], []
        for child in node.children:
            cid = child.name
            if cid in self.akr.pre_wgd_graphs:
                child_graphs.append(self.akr.pre_wgd_graphs[cid])
                child_source_ids.append(f"{cid}_pre")
            elif child.is_leaf():
                if cid in self.akr.leaf_graphs:
                    child_graphs.append(self.akr.leaf_graphs[cid])
                    child_source_ids.append(cid)
            elif cid in self.akr.anc_graphs:
                child_graphs.append(self.akr.anc_graphs[cid])
                child_source_ids.append(cid)
        return child_graphs, child_source_ids

    def _log_children(self, child_graphs, child_source_ids):
        """Log child chromosome counts."""
        info = []
        for cg, cid in zip(child_graphs, child_source_ids):
            info.append(f"{cid}={len(list(cg.chromosomes))}")
        logger.info("  Children: %s", ", ".join(info))

    # ── Phase 1: dedup + merge ───────────────────────────────────────

    def _phase1_dedup_and_merge(self, G, mapped_children, child_source_ids, node_id):
        """Per-child dedup, validate, merge into ColoredGraph."""
        pre_dedup_chroms = {}
        for mc, cid in zip(mapped_children, child_source_ids):
            n_chrom = len(list(mc.chromosomes))
            pre_dedup_chroms[cid] = n_chrom
            n_nodes, n_edges, n_cc = _graph_stats(mc.graph)
            logger.info("  [Phase 1] %s: %d chroms, %d nodes, %d edges, %d cc",
                        cid, n_chrom, n_nodes, n_edges, n_cc)

        deduped = G._deduplicate_children(mapped_children, child_source_ids,
                                           ref_graphs=mapped_children)

        for mc, cid in zip(deduped, child_source_ids):
            mc.rebuild_edges_from_chrom_hogs()

            n_chrom = len(list(mc.chromosomes))
            n_nodes, n_edges, n_cc = _graph_stats(mc.graph)
            logger.info("  [Phase 1] %s deduped: %d chroms, %d nodes, %d edges, %d cc",
                        cid, n_chrom, n_nodes, n_edges, n_cc)

            ok, errors = _validate_dedup(mc, cid, pre_dedup_chroms[cid])
            if ok:
                logger.info("  [Phase 1] %s: dedup validated (chrom count ok, no dups, no self-loops)", cid)
            else:
                for err in errors:
                    logger.error("  [Phase 1] %s: %s", cid, err)

            G.add_child(cid, mc)

        # Dedup events
        dedup_events = [e for e in G.events
                        if e.event_type in ('tandem_dup', 'dispersed_dup')]
        if dedup_events:
            _log_child_events(dedup_events, "Phase 1")

        n_nodes, n_edges, n_cc = _graph_stats(G)
        logger.info("  [Phase 1] merged: %d nodes, %d edges, %d cc",
                     n_nodes, n_edges, n_cc)

        return deduped

    # ── Outgroup info ────────────────────────────────────────────────

    def _collect_outgroup_info(self, node_id, child_source_ids):
        """Collect outgroup adjacency and HOG info for event polarity."""
        outgroup_adjacency = None
        outgroup_hogs = {}

        if node_id not in self.outgroup_leaves_map:
            return outgroup_adjacency, outgroup_hogs

        og_leaves = self.outgroup_leaves_map[node_id]
        # Find parent HOG level
        node = self.akr.tree.search_nodes(name=node_id)
        if not node or not node[0].up:
            return outgroup_adjacency, outgroup_hogs
        parent_hog_level = node[0].up.name

        # Collect outgroup adjacencies
        og_adj_counts = defaultdict(int)
        n_og_leaves = 0
        logger.debug("  [outgroup] %s: sibling leaves %s -> mapping to %s",
                     node_id, og_leaves, parent_hog_level)

        for leaf_name in og_leaves:
            if leaf_name not in self.akr.leaf_graphs:
                continue
            try:
                mapped = self.akr._map_to_parent_hogs(
                    parent_hog_level,
                    self.akr.leaf_graphs[leaf_name],
                    source_id=leaf_name)
                n_og_leaves += 1
                leaf_adj = set()
                for h1, h2 in mapped.get_adjacencies(include_telomere=False):
                    h1_id = h1.hog_id if hasattr(h1, 'hog_id') else str(h1)
                    h2_id = h2.hog_id if hasattr(h2, 'hog_id') else str(h2)
                    key = (h1_id, h2_id) if h1_id < h2_id else (h2_id, h1_id)
                    leaf_adj.add(key)
                for key in leaf_adj:
                    og_adj_counts[key] += 1
                logger.debug("  [outgroup]   %s: %d adjacencies at %s level",
                             leaf_name, len(leaf_adj), parent_hog_level)
            except Exception as e:
                logger.debug("  [outgroup]   %s: skip (%s)", leaf_name, e)

        if og_adj_counts and n_og_leaves >= 2:
            outgroup_adjacency = {k for k, cnt in og_adj_counts.items()
                                  if cnt >= n_og_leaves}
            logger.info("  [outgroup] %s: %d/%d adjacencies conserved in all %d leaves at %s level",
                        node_id, len(outgroup_adjacency), len(og_adj_counts),
                        n_og_leaves, parent_hog_level)

        return outgroup_adjacency, outgroup_hogs

    # ── WGD collapse ─────────────────────────────────────────────────

    def _wgd_collapse(self, post_graph, node_id, ploidy, parent_hog_level):
        """WGD collapse: pair post-WGD chromosomes by HOG similarity → pre-WGD."""
        from .takr_colored_graph import ColoredGraph, _merge_chromosome_paths

        logger.info("Reconstructing node %s_preWGD [v2 ColoredGraph]", node_id)
        t0 = time.time()

        child_chroms = post_graph._child_chromosomes.get(node_id, [])
        n_orig = len(child_chroms)
        logger.info("  Children: %s=%d", node_id, n_orig)

        if ploidy <= 1 or n_orig < ploidy:
            pre_anc = self.akr._map_to_parent_hogs(
                parent_hog_level, post_graph, source_id=f"{node_id}_pre")
            pre_anc.node_id = f"{node_id}_pre"
            logger.info("  Done: no WGD (ploidy=%d), 0 chroms (%.1fs)",
                         ploidy, time.time() - t0)
            return pre_anc, []

        # Pair chromosomes by HOG Jaccard similarity
        chrom_hogs = [set(p) for p in child_chroms]
        n = len(chrom_hogs)
        pairs = []
        for i, j in itertools.combinations(range(n), 2):
            inter = len(chrom_hogs[i] & chrom_hogs[j])
            union = len(chrom_hogs[i] | chrom_hogs[j])
            jaccard = inter / union if union > 0 else 0.0
            if inter > 0:
                pairs.append((jaccard, inter, i, j))

        pairs.sort(reverse=True)
        paired = set()
        target = n // ploidy
        pre_chroms = []
        for jaccard, inter, i, j in pairs:
            if i in paired or j in paired:
                continue
            if len(pre_chroms) >= target:
                break
            pre_chroms.append((child_chroms[i], child_chroms[j]))
            paired.add(i)
            paired.add(j)

        for i in range(n):
            if i not in paired:
                pre_chroms.append((child_chroms[i],))

        # Build pre-WGD graph
        pre_G = ColoredGraph(hog_level=f"{node_id}_preWGD")
        for chrom_idx, cp in enumerate(pre_chroms):
            merged = _merge_chromosome_paths(cp[0], cp[1]) if len(cp) == 2 else list(cp[0])
            for k in range(len(merged) - 1):
                pre_G.add_edge(merged[k], merged[k + 1],
                               f"pre_chr{chrom_idx}", chrom_idx)

        pre_G.resolve_all_events(
            outgroups=None, min_hogs=self.min_hogs,
            gfa_debug=self.gfa_debug,
            gfa_prefix=self._gfa_prefix(node_id) + "_preWGD")

        all_events = post_graph.events + pre_G.events
        pre_anc = pre_G.to_ancestral_graph()
        pre_anc.node_id = f"{node_id}_pre"
        n_pre = len(list(pre_anc.chromosomes))

        if parent_hog_level != node_id:
            saved_events = list(all_events)
            pre_anc = self.akr._map_to_parent_hogs(
                parent_hog_level, pre_anc, source_id=f"{node_id}_pre")
            pre_anc.events = saved_events

        pre_anc.node_id = f"{node_id}_pre"
        logger.info("  Done: %d -> %d chroms, %d events (%.1fs)",
                     n_orig, n_pre, len(all_events), time.time() - t0)
        return pre_anc, all_events

    # ── Visualization ────────────────────────────────────────────────

    def _viz_prefix(self, node_id):
        """Visualization output prefix for a node."""
        viz_dir = os.path.dirname(self.akr.outpre) if hasattr(self.akr, 'outpre') else '.'
        viz_base = os.path.basename(self.akr.outpre) if hasattr(self.akr, 'outpre') else 'AKR'
        return os.path.join(viz_dir, f'{viz_base}.{node_id}')

    def _gfa_prefix(self, node_id):
        """GFA debug file prefix for a node."""
        return self._viz_prefix(node_id)

    def _output_child_gfas(self, children, child_source_ids, node_id):
        """输出每个子图的 block 级 GFA（dedup 之前）。"""
        from .takr_colored_graph import ColoredGraph
        import networkx as nx
        for mc, cid in zip(children, child_source_ids):
            try:
                tmp = ColoredGraph(hog_level=node_id)
                tmp.add_child(cid, mc)

                # Debug: HOG 图连通分量
                hog_components = list(nx.weakly_connected_components(tmp))
                logger.info("  [gfa:child %s] HOG graph: %d nodes, %d edges, %d cc",
                            cid, tmp.node_count(), tmp.edge_count(), len(hog_components))

                # Debug: 度数分布 + per-chromosome 分支/端粒
                degs = dict(tmp.degree())
                max_deg = max(degs.values()) if degs else 0
                branch_count = sum(1 for d in degs.values() if d > 2)
                tel_count = sum(1 for n in tmp.nodes
                                if tmp.nodes[n].get('telomere'))
                logger.info("  [gfa:child %s] degree: max=%d, branch(d>2)=%d, telomere_HOGs=%d",
                            cid, max_deg, branch_count, tel_count)

                # Per-chromosome
                chr_info = defaultdict(lambda: {'nodes':0, 'branch':0, 'tel':0})
                for n in tmp.nodes:
                    srcs = tmp.nodes[n].get('sources', set())
                    if not srcs: continue
                    ch = next(iter(srcs))[1]
                    chr_info[ch]['nodes'] += 1
                    if degs.get(n, 0) > 2:
                        chr_info[ch]['branch'] += 1
                    if tmp.nodes[n].get('telomere'):
                        chr_info[ch]['tel'] += 1
                for ch in sorted(chr_info.keys()):
                    ci = chr_info[ch]
                    logger.info("  [gfa:child %s]   chr %s: nodes=%d, branch=%d, tel=%d",
                                cid, ch, ci['nodes'], ci['branch'], ci['tel'])

                # Debug: 每染色体 HOG 边数（压缩前）
                chr_edges = defaultdict(int)
                for h1, h2 in tmp.edges():
                    srcs = tmp.nodes[h1].get('sources', set())
                    if srcs:
                        ch = next(iter(srcs))[1]
                        chr_edges[ch] += 1
                for ch in sorted(chr_info.keys()):
                    logger.info("  [gfa:child %s]   chr %s: nodes=%d, hog_edges=%d, branch=%d, tel=%d",
                                cid, ch, chr_info[ch]['nodes'], chr_edges.get(ch, 0),
                                chr_info[ch]['branch'], chr_info[ch]['tel'])

                # HOG 级 GFA（压缩前对照）
                hog_path = f'{self._viz_prefix(node_id)}.child_{cid}.hog.gfa'
                with open(hog_path, 'w') as fout:
                    fout.write(f"H\ttype:child_hog\tparent:{node_id}\tchild:{cid}\t"
                               f"nodes:{tmp.node_count()}\tedges:{tmp.edge_count()}\n")
                    tmp.to_gfa(fout, use_blocks=False)
                logger.info("  [gfa] child HOG %s -> %s", cid, hog_path)

                # Debug: 每染色体块数
                tmp._ensure_blocks()
                # 重建块间边：沿有向图线性行走，保证非端粒 block 有 in/out
                tmp._compress_to_block_level()
                chr_blocks = defaultdict(list)
                for bid, hogs in tmp._blocks.items():
                    # 取第一个 HOG 的 sources 中的 chrom
                    chrom = '?'
                    for h in hogs[:1]:
                        srcs = tmp.nodes[h].get('sources', set())
                        if srcs:
                            chrom = next(iter(srcs))[1]
                            break
                    chr_blocks[chrom].append((bid, len(hogs)))

                # 孤立 block（block graph 中 degree=0）
                bg = getattr(tmp, '_block_graph', None)
                if bg:
                    isolated = [n for n in bg.nodes if bg.degree(n) == 0]
                    if isolated:
                        logger.warning("  [gfa:child %s] %d isolated blocks in block_graph: %s",
                                       cid, len(isolated), isolated[:5])

                for chrom in sorted(chr_blocks.keys()):
                    blks = chr_blocks[chrom]
                    sizes = [s for _, s in blks]
                    logger.info("  [gfa:child %s]   chr %s: %d blocks, sizes %d-%d, cc_in_hog=%d",
                                cid, chrom, len(blks),
                                min(sizes), max(sizes),
                                sum(1 for comp in hog_components
                                    if any(h in comp for h in tmp.nodes
                                           if tmp.nodes[h].get('sources', set())
                                           and next(iter(tmp.nodes[h].get('sources', set())))[1] == chrom)))

                # 诊断：per-chromosome block 连通性
                if bg:
                    for chrom in sorted(chr_blocks.keys()):
                        blks = chr_blocks[chrom]
                        for bid, _ in blks:
                            if not bg.has_node(bid):
                                continue
                            nbrs = list(bg.neighbors(bid))
                            # 找邻居的染色体
                            nbr_chroms = []
                            for nb in nbrs:
                                nb_hogs = tmp._blocks.get(nb, [])
                                nb_chrom = '?'
                                for h in nb_hogs[:1]:
                                    srcs = tmp.nodes[h].get('sources', set())
                                    if srcs:
                                        nb_chrom = next(iter(srcs))[1]
                                        break
                                nbr_chroms.append(nb_chrom)
                            # 只报告有跨染色体连接或度≠2 的内部 block
                            is_tel = bool(tmp._blocks.get(bid) and 
                                         tmp.nodes[tmp._blocks[bid][0]].get('telomere'))
                            cross = [nc for nc in nbr_chroms if nc != chrom]
                            if cross or (not is_tel and len(nbrs) != 2) or (is_tel and len(nbrs) > 2):
                                logger.warning("  [gfa:child %s]   BLOCK %s (chr=%s, size=%d, tel=%s): "
                                               "deg=%d, cross_chr=%s, nbrs=%s",
                                               cid, bid, chrom, len(tmp._blocks.get(bid, [])),
                                               is_tel, len(nbrs), cross,
                                               [(nb, nc) for nb, nc in zip(nbrs, nbr_chroms)])

                path = f'{self._viz_prefix(node_id)}.child_{cid}.gfa'
                with open(path, 'w') as fout:
                    fout.write(f"H\ttype:child\tparent:{node_id}\tchild:{cid}\t"
                               f"nodes:{tmp.node_count()}\tedges:{tmp.edge_count()}\n")
                    tmp.to_gfa(fout)
                logger.info("  [gfa] child %s -> %s", cid, path)

                # 诊断：验证每个非端粒 block 都有 successors 和 predecessors
                broken = tmp.check_block_degrees()
                if broken:
                    by_chr = defaultdict(list)
                    for bid, in_d, out_d, ch, sz in broken:
                        by_chr[ch].append((bid, in_d, out_d, sz))
                    logger.warning("  [gfa:child %s] %d non-tel BLOCKS missing in/out edges:",
                                   cid, len(broken))
                    for ch in sorted(by_chr.keys()):
                        items = by_chr[ch]
                        examples = [(bid, f"in={id_},out={od},sz={sz}")
                                   for bid, id_, od, sz in items[:5]]
                        logger.warning("  [gfa:child %s]   chr %s: %d broken blocks: %s",
                                       cid, ch, len(items), examples)
                else:
                    logger.info("  [gfa:child %s] all non-tel BLOCKS pass in/out check ✓",
                                cid)
            except Exception as e:
                logger.exception("  [gfa] child %s failed: %s", cid, e)

    def _viz_child_paths(self, children, child_source_ids, node_id):
        """Draw per-child block graphs before merging."""
        from .takr_colored_graph import ColoredGraph
        try:
            for mc, cid in zip(children, child_source_ids):
                outpath = f'{self._viz_prefix(node_id)}.child_{cid}.png'
                try:
                    tmp = ColoredGraph(hog_level=node_id)
                    tmp.add_child(cid, mc)
                    tmp._build_synteny_blocks()
                    tmp._compress_to_block_level()
                    tmp.draw_block_graph(
                        outpath,
                        title=f'Child {cid} ({node_id}): block graph')
                except Exception as e:
                    logger.debug("  [viz] child %s block graph failed: %s", cid, e)
        except Exception as e:
            logger.debug("  [viz] child paths skipped: %s", e)

    def _viz_raw_graph(self, G, node_id):
        """Draw raw merged block graph before event resolution."""
        try:
            G._build_synteny_blocks()
            n_blocks = len(getattr(G, '_blocks', []))
            if n_blocks <= 200:
                G.draw_block_graph(
                    f'{self._viz_prefix(node_id)}.raw_block_graph.png',
                    title=f'Raw Block Graph: {node_id} (before resolution)')
            else:
                logger.info("  [viz] skipping raw graph: %d blocks (too many)", n_blocks)
        except Exception as e:
            logger.debug("  [viz] raw graph skipped: %s", e)

    def _viz_resolved_graph(self, G, node_id):
        """Draw resolved block graph + adjacency heatmap."""
        try:
            n_blocks = len(getattr(G, '_blocks', []))
            if n_blocks <= 200:
                G.draw_block_graph(
                    f'{self._viz_prefix(node_id)}.block_graph.png',
                    title=f'Block Graph: {node_id}')
                G.draw_adjacency_heatmap(
                    f'{self._viz_prefix(node_id)}.adj_heatmap.png',
                    title=f'Adjacency Matrix: {node_id}')
            else:
                logger.info("  [viz] skipping resolved graph: %d blocks (too many)", n_blocks)
        except Exception as e:
            logger.debug("  [viz] skipped: %s", e)


# =========================================================================
#  Main entry point (delegates to ReconstructorV2)
# =========================================================================

def reconstruct_event_driven_v2(akr, min_hogs=3, gfa_debug=False):
    """Event-driven reconstruction v2 — ColoredGraph pipeline.

    Thin wrapper around ReconstructorV2 for backward compatibility.
    """
    recon = ReconstructorV2(akr, min_hogs=min_hogs, gfa_debug=gfa_debug)
    return recon.run()


