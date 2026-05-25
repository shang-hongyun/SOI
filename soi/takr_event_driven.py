#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""takr_event_driven.py - Event-driven ancestral karyotype reconstruction.

Core algorithm for TAKR v4. Called by AK.run() when use_v4=True.

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

import copy
import logging
import os
import time
from collections import defaultdict
from typing import Dict, List, Optional, Set

import networkx as nx

logger = logging.getLogger(__name__)


# =========================================================================
#  Helpers
# =========================================================================

def _graph_stats(graph):
    """返回 (nodes, edges, connected_components)。自动区分有向/无向图。"""
    n = graph.number_of_nodes()
    e = graph.number_of_edges()
    if graph.is_directed():
        cc = nx.number_weakly_connected_components(graph)
    else:
        cc = nx.number_connected_components(graph)
    return (n, e, cc)


def _validate_dedup(mc, cid, expected_chrom_count):
    """验证 dedup 后的图状态。返回 (ok, errors)。"""
    errors = []
    ch = getattr(mc, 'chrom_hogs', None)
    if not ch:
        errors.append("no chrom_hogs found")
        return (False, errors)
    # 1. 染色体数不变
    n_chrom = len(ch)
    if n_chrom != expected_chrom_count:
        errors.append(f"chrom count changed {expected_chrom_count} → {n_chrom}")
    # 2. 无重复 HOG（每条染色体内）
    for ci, hogs in ch.items():
        gene_hogs = [str(h) for h in hogs if h not in mc.telomeres]
        seen = set()
        for h in gene_hogs:
            if h in seen:
                errors.append(f"chrom{ci}: duplicate HOG {h}")
            seen.add(h)
    # 3. 无自环（连续相同 HOG）
    for ci, hogs in ch.items():
        gene_hogs = [str(h) for h in hogs if h not in mc.telomeres]
        for i in range(len(gene_hogs) - 1):
            if gene_hogs[i] == gene_hogs[i + 1]:
                errors.append(f"chrom{ci}: self-loop at pos {i} ({gene_hogs[i]})")
    return (len(errors) == 0, errors)


# =========================================================================
#  Visualization helpers
# =========================================================================

def _draw_child_paths(child_graph, child_id, node_id, outpath):
    """Draw a child's block-level graph using the same logic as merged block graph."""
    from soi.takr_colored_graph import ColoredGraph
    import logging as _logging
    _logger = _logging.getLogger(__name__)

    try:
        # Build a temporary ColoredGraph with just this one child
        tmp = ColoredGraph(hog_level=node_id)
        tmp.add_child(child_id, child_graph)
        tmp._build_synteny_blocks()
        tmp._compress_to_block_level()
        tmp.draw_block_graph(outpath,
                             title=f'Child {child_id} ({node_id}): block graph')
    except Exception as e:
        _logger.debug("  [viz] child %s block graph failed: %s", child_id, e)


# =========================================================================
#  Main entry point
# =========================================================================

def reconstruct_event_driven(akr, min_hogs=3):
    """Event-driven reconstruction — main entry (v1: consensus ancestor + branch detection).

    Args:
        akr: AKR instance (tree, leaf_graphs, hog already built)
        min_hogs: Minimum HOGs for a valid event (noise filter)

    Returns:
        {node_name: AncestralAdjacencyGraph}
    """
    logger.info("=== Event-driven reconstruction (v4) ===")
    t0 = time.time()

    og_graphs_cache = {}

    for node in akr.tree.traverse(strategy="postorder"):
        if node.is_leaf():
            continue

        node_id = node.name
        logger.info("Reconstructing node %s [v4 event-driven]", node_id)
        t_node = time.time()

        # Collect child graphs (check for pre-WGD graphs first)
        child_graphs, child_source_ids = [], []
        for child in node.children:
            cid = child.name
            # Use pre-WGD graph if available (post-WGD already collapsed)
            if cid in akr.pre_wgd_graphs:
                child_graphs.append(akr.pre_wgd_graphs[cid])
                child_source_ids.append(cid)
            elif child.is_leaf():
                if cid in akr.leaf_graphs:
                    child_graphs.append(akr.leaf_graphs[cid])
                    child_source_ids.append(cid)
            elif cid in akr.anc_graphs:
                child_graphs.append(akr.anc_graphs[cid])
                child_source_ids.append(cid)

        if len(child_graphs) < 2:
            continue

        # Map children to ancestor HOG level
        hog_level = node_id
        mapped_children = []
        for cg, cid in zip(child_graphs, child_source_ids):
            mc = akr._map_to_parent_hogs(hog_level, cg, source_id=cid)
            mapped_children.append(mc)

        # Build consensus ancestor
        ancestor = _build_consensus_ancestor(node_id, mapped_children, child_source_ids)

        # Build outgroup adjacency set for polarity voting
        og_adj_set = set()
        og_weighted = set()
        if node_id in og_graphs_cache:
            total_og_weight = sum(w for _, w in og_graphs_cache[node_id])
            og_threshold = total_og_weight / 3.0 if total_og_weight > 0 else 0
            for og_graph, weight in og_graphs_cache[node_id]:
                for h1, h2 in og_graph.get_adjacencies(include_telomere=False):
                    key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
                    og_adj_set.add(key)
                    if weight >= og_threshold:
                        og_weighted.add(key)

        # Detect events on each child's branch
        all_events = []
        for ci, (mc, src) in enumerate(zip(mapped_children, child_source_ids)):
            branch = "%s-%s" % (node_id, src)
            events = _detect_branch_events(
                ancestor, mc, branch, src, min_hogs,
                mapped_children, child_source_ids, og_adj_set)
            all_events.extend(events)

        ancestor.events = all_events
        akr.anc_graphs[node_id] = ancestor

        # If this is a WGD node, collapse to pre-WGD and detect events
        ploidy = akr.ploidy_map.get(node_id, 1)
        if ploidy > 1:
            _handle_wgd_node(akr, node, ploidy)

        n_chrom = len(list(ancestor.chromosomes))
        logger.info("  Done: %d chroms, %d events (%.1fs)",
                     n_chrom, len(all_events), time.time() - t_node)

    logger.info("=== Event-driven reconstruction done (%.1fs) ====",
                time.time() - t0)
    return akr.anc_graphs

def reconstruct_event_driven_v2(akr, min_hogs=3):
    """Event-driven reconstruction v2 -- ColoredGraph: colored graph + cycle detection + path cover."""
    from .takr_colored_graph import ColoredGraph

    def _wgd_collapse(post_graph, node_id, ploidy, parent_hog_level, min_hogs=3):
        """WGD collapse: pair post-WGD chromosomes by HOG similarity → pre-WGD.

        post_graph: ColoredGraph (node's own genome, already at the correct HOG level).
        parent_hog_level: 父节点的 HOG 层级（结果映射到此层级）。
        """
        logger.info("Reconstructing node %s_preWGD [v2 ColoredGraph]", node_id)
        t_collapse = time.time()

        # 桥接检测：跨染色体边 = post-WGD 融合/嵌入
        n_before = len(post_graph.events)
        post_graph.resolve_bridge_events()
        n_bridges = len(post_graph.events) - n_before
        logger.info("  [colored] bridges: %d events", n_bridges)

        paths = post_graph.path_cover()
        n_post = len(paths)
        logger.info("  Children: %s=%d", node_id, n_post)

        # If ploidy==1 or fewer chroms than ploidy, no collapse needed
        if ploidy <= 1 or n_post < ploidy:
            anc = post_graph.to_ancestral_graph()
            anc.node_id = "{}_pre".format(node_id)
            logger.info("  Done: %d -> %d chroms, 0 events (%.1fs)",
                         n_post, n_post, time.time() - t_collapse)
            # Map to parent HOG level
            if parent_hog_level != node_id:
                saved_events = list(anc.events)
                anc = akr._map_to_parent_hogs(parent_hog_level, anc,
                                               source_id="{}_pre".format(node_id))
                anc.events = saved_events
            anc.node_id = "{}_pre".format(node_id)
            return anc, []

        # Build chromosome-level similarity matrix
        chrom_hogs = [set(p) for p in paths]
        n = len(chrom_hogs)
        import itertools
        pairs = []
        for i, j in itertools.combinations(range(n), 2):
            inter = len(chrom_hogs[i] & chrom_hogs[j])
            union = len(chrom_hogs[i] | chrom_hogs[j])
            jaccard = inter / union if union > 0 else 0.0
            if inter > 0:
                pairs.append((jaccard, inter, i, j))

        # Greedy matching: pair highest similarity chromosomes
        pairs.sort(reverse=True)
        paired = set()
        target = n // ploidy
        pre_chroms = []
        for jaccard, inter, i, j in pairs:
            if i in paired or j in paired:
                continue
            if len(pre_chroms) >= target:
                break
            pre_chroms.append((paths[i], paths[j]))
            paired.add(i)
            paired.add(j)

        # Unpaired chromosomes stay as singles
        for i in range(n):
            if i not in paired:
                pre_chroms.append((paths[i],))

        # Build pre-WGD graph from merged chromosomes
        from .takr_colored_graph import _merge_chromosome_paths
        pre_G = ColoredGraph(hog_level="{}_preWGD".format(node_id))
        for chrom_idx, cp in enumerate(pre_chroms):
            if len(cp) == 2:
                merged = _merge_chromosome_paths(cp[0], cp[1])
            else:
                merged = list(cp[0])
            for i in range(len(merged) - 1):
                pre_G.add_edge(merged[i], merged[i + 1],
                               "pre_chr{}".format(chrom_idx), chrom_idx)

        pre_G.resolve_all_events(outgroups=None, min_hogs=min_hogs)
        all_events = post_graph.events + pre_G.events
        pre_anc = pre_G.to_ancestral_graph()
        pre_anc.node_id = "{}_pre".format(node_id)
        n_pre = len(list(pre_anc.chromosomes))
        # Map to parent HOG level
        if parent_hog_level != node_id:
            saved_events = list(all_events)
            pre_anc = akr._map_to_parent_hogs(parent_hog_level, pre_anc,
                                               source_id="{}_pre".format(node_id))
            pre_anc.events = saved_events
        pre_anc.node_id = "{}_pre".format(node_id)
        logger.info("  Done: %d -> %d chroms, %d events (%.1fs)",
                     n_post, n_pre, len(all_events), time.time() - t_collapse)
        return pre_anc, all_events

    logger.info("=== Event-driven reconstruction v2 (ColoredGraph) ===")
    t0 = time.time()
    og_graphs_cache = {}

    # Precompute: for each internal node, find its outgroup (sibling subtree leaves)
    # outgroup_leaves_map: node_id -> [leaf_name, ...]
    outgroup_leaves_map = {}
    for node in akr.tree.traverse(strategy="postorder"):
        if node.is_leaf():
            continue
        node_id = node.name
        if not node.up:
            # Root node has no outgroup
            continue
        # Find the sibling: another child of parent that is not this node
        sibling_leaves = []
        for sibling in node.up.children:
            if sibling.name == node_id:
                continue
            # Collect all leaves under sibling
            for leaf in sibling.iter_leaves():
                sibling_leaves.append(leaf.name)
        if sibling_leaves:
            outgroup_leaves_map[node_id] = sibling_leaves

    for node in akr.tree.traverse(strategy="postorder"):
        node_id = node.name
        ploidy = akr.ploidy_map.get(node_id, 1)
        if node.is_leaf() and ploidy > 1 and node_id in akr.leaf_graphs:
            if node_id in akr.pre_wgd_graphs:
                continue
            # Polyploid leaf: collapse to pre-WGD before parent uses it
            leaf_node_s = akr.tree.search_nodes(name=node_id)
            if not leaf_node_s or not leaf_node_s[0].up:
                continue
            parent_id = leaf_node_s[0].up.name
            mapped = akr._map_to_parent_hogs(parent_id, akr.leaf_graphs[node_id], source_id=node_id)
            G_leaf = ColoredGraph(hog_level=node_id)
            G_leaf.add_child(node_id, mapped)
            pre_anc, _ = _wgd_collapse(G_leaf, node_id, ploidy, parent_id)
            akr.pre_wgd_graphs[node_id] = pre_anc
            continue
        elif node.is_leaf():
            continue
        # Internal node
        logger.info("Reconstructing node %s [v2 ColoredGraph]", node_id)
        t_node = time.time()
        child_graphs, child_source_ids = [], []
        for child in node.children:
            cid = child.name
            if cid in akr.pre_wgd_graphs:
                child_graphs.append(akr.pre_wgd_graphs[cid])
                child_source_ids.append(cid)
            elif child.is_leaf():
                if cid in akr.leaf_graphs:
                    child_graphs.append(akr.leaf_graphs[cid])
                    child_source_ids.append(cid)
            elif cid in akr.anc_graphs:
                child_graphs.append(akr.anc_graphs[cid])
                child_source_ids.append(cid)
        if len(child_graphs) < 2:
            continue
        # Log children info
        child_info = []
        for cg, cid in zip(child_graphs, child_source_ids):
            n_ch = len(list(cg.chromosomes))
            child_info.append("{}={}".format(cid, n_ch))
        logger.info("  Children: %s", ", ".join(child_info))
        hog_level = node_id
        mapped_children = []
        for cg, cid in zip(child_graphs, child_source_ids):
            mc = akr._map_to_parent_hogs(hog_level, cg, source_id=cid)
            mapped_children.append(mc)

        # Phase 1: 每孩子内部 deduplication (tandem/dispersed/proximal/seg_dup)
        G = ColoredGraph(hog_level=node_id)
        pre_dedup_chroms = {}
        for mc, cid in zip(mapped_children, child_source_ids):
            n_chrom = len(list(mc.chromosomes))
            pre_dedup_chroms[cid] = n_chrom
            n_nodes, n_edges, n_cc = _graph_stats(mc.graph)
            logger.info("  [Phase 1] %s: %d chroms, %d nodes, %d edges, %d cc",
                        cid, n_chrom, n_nodes, n_edges, n_cc)
        deduped_children = G._deduplicate_children(mapped_children, child_source_ids,
                                                       ref_graphs=mapped_children)
        for mc, cid in zip(deduped_children, child_source_ids):
            # 从 dedup 后的 chrom_hogs 重建图边（简化图）
            mc.rebuild_edges_from_chrom_hogs()

            n_chrom = len(list(mc.chromosomes))
            n_nodes, n_edges, n_cc = _graph_stats(mc.graph)
            logger.info("  [Phase 1] %s deduped: %d chroms, %d nodes, %d edges, %d cc",
                        cid, n_chrom, n_nodes, n_edges, n_cc)

            # 该孩子的事件汇总
            child_events = [e for e in G.events
                            if e.event_type in ('tandem_dup', 'dispersed_dup')
                            and cid in e.branch]
            if child_events:
                tandem = [e for e in child_events if e.event_type == 'tandem_dup']
                disp = [e for e in child_events if e.event_type == 'dispersed_dup']
                parts = []
                if tandem:
                    lens = [len(e.genes_involved) for e in tandem]
                    parts.append(f"tandem_dup={len(tandem)} (len {min(lens)}-{max(lens)})")
                if disp:
                    lens = [len(e.genes_involved) for e in disp]
                    parts.append(f"dispersed_dup={len(disp)} (len {min(lens)}-{max(lens)})")
                logger.info("  [Phase 1] %s events: %s", cid, ", ".join(parts))

            # ── dedup 后验证（重建后） ──
            ok, errors = _validate_dedup(mc, cid, pre_dedup_chroms[cid])
            if ok:
                logger.info("  [Phase 1] %s: dedup validated (chrom count ok, no dups, no self-loops)", cid)
            else:
                for err in errors:
                    logger.error("  [Phase 1] %s: %s", cid, err)

            G.add_child(cid, mc)

        # 合图后统计
        n_merged_nodes, n_merged_edges, n_merged_cc = _graph_stats(G._graph)
        logger.info("  [Phase 1] merged: %d nodes, %d edges, %d cc",
                    n_merged_nodes, n_merged_edges, n_merged_cc)
        # Visualization: per-child chromosome paths (before merging)
        try:
            viz_dir = os.path.dirname(akr.outpre) if hasattr(akr, 'outpre') else '.'
            viz_base = os.path.basename(akr.outpre) if hasattr(akr, 'outpre') else 'AKR'
            for mc, cid in zip(deduped_children, child_source_ids):
                outpath = os.path.join(viz_dir, f'{viz_base}.{node_id}.child_{cid}.png')
                _draw_child_paths(mc, cid, node_id, outpath)
        except Exception as e:
            logger.debug("  [viz] child paths skipped: %s", e)

        # Visualization: raw merged graph (before event resolution)
        # 跳过块数太多的情况（graphviz 渲染太慢）
        try:
            viz_dir = os.path.dirname(akr.outpre) if hasattr(akr, 'outpre') else '.'
            viz_base = os.path.basename(akr.outpre) if hasattr(akr, 'outpre') else 'AKR'
            G._build_synteny_blocks()
            n_blocks = len(G._blocks) if hasattr(G, '_blocks') else 0
            if n_blocks <= 200:
                G.draw_block_graph(
                    os.path.join(viz_dir, f'{viz_base}.{node_id}.raw_block_graph.png'),
                    title=f'Raw Block Graph: {node_id} (before resolution)')
            else:
                logger.info("  [viz] skipping raw graph: %d blocks (too many)", n_blocks)
        except Exception as e:
            logger.debug("  [viz] raw graph skipped: %s", e)

        # === 收集外类群邻接信息 ===
        # 外类群在 parent HOG level, 用于 bridge 事件极性判定
        outgroup_adjacency = None
        if node_id in outgroup_leaves_map and node.up:
            parent_hog_level = node.up.name
            og_adj_counts = defaultdict(int)  # key -> number of outgroup leaves with it
            og_leaves = outgroup_leaves_map[node_id]
            n_og_leaves = 0
            logger.debug("  [outgroup] %s: sibling leaves %s -> mapping to %s",
                         node_id, og_leaves, parent_hog_level)
            for leaf_name in og_leaves:
                if leaf_name in akr.leaf_graphs:
                    try:
                        mapped = akr._map_to_parent_hogs(
                            parent_hog_level,
                            akr.leaf_graphs[leaf_name],
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
            # 只有所有 outgroup 物种都有的邻接才算祖先态
            if og_adj_counts and n_og_leaves >= 2:
                outgroup_adjacency = {k for k, cnt in og_adj_counts.items()
                                      if cnt >= n_og_leaves}
                logger.info("  [outgroup] %s: %d/%d adjacencies conserved in all %d leaves at %s level",
                            node_id, len(outgroup_adjacency), len(og_adj_counts),
                            n_og_leaves, parent_hog_level)

        outgroup_hogs = {}
        if node_id in og_graphs_cache:
            for og_graph, _weight in og_graphs_cache[node_id]:
                for cid in child_source_ids:
                    if cid not in outgroup_hogs:
                        outgroup_hogs[cid] = set()
                    for h1, h2 in og_graph.get_adjacencies(include_telomere=False):
                        outgroup_hogs[cid].add(h1)
                        outgroup_hogs[cid].add(h2)
        G.resolve_all_events(outgroups=outgroup_hogs,
                             outgroup_adjacency=outgroup_adjacency,
                             min_hogs=min_hogs)

        # Visualization: block graph + adjacency heatmap
        try:
            viz_dir = os.path.dirname(akr.outpre) if hasattr(akr, 'outpre') else '.'
            viz_base = os.path.basename(akr.outpre) if hasattr(akr, 'outpre') else 'AKR'
            n_blocks = len(G._blocks) if hasattr(G, '_blocks') else 0
            if n_blocks <= 200:
                G.draw_block_graph(
                    os.path.join(viz_dir, f'{viz_base}.{node_id}.block_graph.png'),
                    title=f'Block Graph: {node_id}')
                G.draw_adjacency_heatmap(
                    os.path.join(viz_dir, f'{viz_base}.{node_id}.adj_heatmap.png'),
                    title=f'Adjacency Matrix: {node_id}')
            else:
                logger.info("  [viz] skipping resolved graph: %d blocks (too many)", n_blocks)
        except Exception as e:
            logger.debug("  [viz] skipped: %s", e)

        ancestor = G.to_ancestral_graph()
        akr.anc_graphs[node_id] = ancestor
        n_chrom = len(list(ancestor.chromosomes))
        logger.info("  Done: %d chroms, %d events (%.1fs)",
                     n_chrom, len(ancestor.events), time.time() - t_node)
        if ploidy > 1:
            parent_node = node.up
            parent_name = parent_node.name if parent_node else node_id
            pre_anc, events = _wgd_collapse(G, node_id, ploidy, parent_name)
            akr.pre_wgd_graphs[node_id] = pre_anc
            virtual_branch = "{}_preWGD-{}".format(node_id, node_id)
            for e in pre_anc.events:
                e.branch = virtual_branch
    logger.info("=== Event-driven reconstruction v2 done (%.1fs) ===",
                time.time() - t0)
    return akr.anc_graphs


# =========================================================================
#  Consensus ancestor graph
# =========================================================================


# =========================================================================
#  Consensus ancestor graph
# =========================================================================

def _build_consensus_ancestor(node_id, mapped_children, child_source_ids):
    """Build consensus ancestor graph from shared child adjacencies.

    The consensus = ancestral state. Adjacencies shared by >=2 children
    are ancestral; unique adjacencies are derived events on that child's branch.
    """
    from .AK import AncestralAdjacencyGraph

    n_children = len(mapped_children)
    ancestor = AncestralAdjacencyGraph(node_id=node_id)

    # Collect all HOGs
    all_hogs = set()
    for mc in mapped_children:
        all_hogs.update(mc.gene_nodes)
    for h in all_hogs:
        ancestor.graph.add_node(h)
        ancestor.gene_nodes.add(h)

    # Build adjacency support counts
    child_adj_sets = []
    for mc in mapped_children:
        adjs = set()
        for h1, h2 in mc.get_adjacencies(include_telomere=False):
            key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
            adjs.add(key)
        child_adj_sets.append(adjs)

    # Consensus: adjacencies shared by >= ceil(n_children/2) children
    # This is the ancestral state estimate (majority rule)
    adj_counts = defaultdict(int)
    for adjs in child_adj_sets:
        for key in adjs:
            adj_counts[key] += 1

    threshold = max(n_children // 2, 1)
    for (h1, h2), count in adj_counts.items():
        if count >= threshold:
            ancestor.graph.add_edge(h1, h2, support=count)

    # Add telomeres
    ancestor._add_telomeres()
    return ancestor


# =========================================================================
#  Branch-level event detection
# =========================================================================

def _detect_branch_events(ancestor, mc, branch, src, min_hogs,
                          all_mapped_children=None, all_child_ids=None,
                          og_adj_set=None):
    """Detect all events on one child branch: ancestor vs child (mc).

    Events detected:
    - inversion (internal or telomere)
    - unidirectional translocation
    - EEJ (end-end join)
    - fission
    - NCF (nested chromosome fusion)
    - RT/URT (reciprocal / unbalanced reciprocal translocation)
    """
    from .takr_events import TAKREvent

    events = []
    child_adj = mc.get_adjacencies(include_telomere=False)
    ancest_adj = ancestor.get_adjacencies(include_telomere=False)

    # Normalize
    child_adj_set = set()
    for h1, h2 in child_adj:
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        child_adj_set.add(key)

    ancest_adj_set = set()
    for h1, h2 in ancest_adj:
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        ancest_adj_set.add(key)

    # Gene-level events: compare HOG copy numbers across children
    # After WGD or duplication, HOGs may have multiple copies.
    # If another child has MORE copies of a HOG than this child → loss
    # If this child has MORE copies than others → gain/duplication
    consensus_hogs = set(ancestor.gene_nodes)
    child_hogs = set(mc.gene_nodes)

    if hasattr(mc, 'hog_map') and mc.hog_map and all_mapped_children:
        # Build copy counts per HOG per child
        child_copy_counts = []
        for mc2 in all_mapped_children:
            counts = defaultdict(int)
            if hasattr(mc2, 'hog_map') and mc2.hog_map:
                for gene, hog in mc2.hog_map.items():
                    counts[hog] += 1
            child_copy_counts.append(counts)

        # For each HOG, find the max copy count across all children = ancestral state
        max_copies = defaultdict(int)
        for counts in child_copy_counts:
            for hog, n in counts.items():
                if n > max_copies[hog]:
                    max_copies[hog] = n

        # This child's copy count
        this_idx = all_child_ids.index(src) if src in all_child_ids else -1
        if this_idx >= 0:
            this_counts = child_copy_counts[this_idx]
            for hog, anc_n in max_copies.items():
                child_n = this_counts.get(hog, 0)
                if child_n < anc_n:
                    # This child has fewer copies → fractionation/gene loss
                    events.append(TAKREvent(
                        event_type='fractionation',
                        branch=branch,
                        genes_involved=[hog],
                        desc="fractionation: %s copies %d->%d in %s" % (
                            hog.hog_id, anc_n, child_n, src),
                        support=anc_n - child_n,
                    ))

    # Gain: HOG present in child but absent in consensus
    gained_hogs = child_hogs - consensus_hogs
    for h in gained_hogs:
        events.append(TAKREvent(
            event_type='gene_gain',
            branch=branch,
            genes_involved=[h],
            desc="gene_gain: %s gained in %s" % (h.hog_id, src),
            support=1,
        ))

    # Inversion detection via adjacency comparison
    inversions = _detect_inversions(
        ancestor, ancest_adj_set, mc, child_adj_set, branch, src, min_hogs)
    events.extend(inversions)

    # Unidirectional translocation
    ut_events = _detect_unidir_trans(
        ancestor, ancest_adj_set, mc, child_adj_set, branch, src)
    events.extend(ut_events)

    # EEJ
    eej = _detect_eej(ancestor, mc, branch, src, og_adj_set)
    events.extend(eej)

    # Fission: internal adjacency in ancestor → chromosome endpoints in child
    fission = _detect_fission(ancestor, ancest_adj_set, mc, branch, src, og_adj_set)
    events.extend(fission)

    # NCF
    ncf = _detect_ncf(ancestor, mc, branch, src)
    events.extend(ncf)

    # RT/URT: use sibling edges with outgroup polarity
    sibling_adj_set = set()
    if all_mapped_children:
        for ci2, mc2 in enumerate(all_mapped_children):
            if all_child_ids and all_child_ids[ci2] == src:
                continue
            for h1, h2 in mc2.get_adjacencies(include_telomere=False):
                key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
                sibling_adj_set.add(key)
    rt = _detect_rt(ancestor, sibling_adj_set, mc, child_adj_set, branch, src, og_adj_set)
    events.extend(rt)

    return events


def _detect_inversions(ancestor, ancest_adj_set, mc, child_adj_set,
                       branch, src, min_hogs):
    """Detect inversions by comparing child HOG order against consensus.

    Key insight: Inversion means a block of HOGs appears reversed
    in the child relative to the ancestor. The ancestral order is
    determined from consensus adjacency relationships:

    For each adjacency A-B in the consensus, A comes before B in
    the ancestral order (assuming no other constraints). We build
    the ancestral order by topological sorting of consensus adjacencies.

    Then compare child's HOG order against this ancestral order
    to find reversed segments.
    """
    from .takr_events import TAKREvent

    events = []

    # Build the consensus adjacency graph as a simple direction map
    # Adjacency edges give us (h1, h2) pairs where h1 is linked to h2
    # For a simple path graph, we can extract the order by following edges
    consensus_adjs = list(ancestor.get_adjacencies(include_telomere=False))

    # Build adjacency-based ordering: start from nodes with in-degree 0
    in_degree = {}
    out_degree = {}
    adj_set = set()
    for h1, h2 in consensus_adjs:
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        adj_set.add(key)
        in_degree[h2] = in_degree.get(h2, 0) + 1
        in_degree.setdefault(h1, 0)
        out_degree[h1] = out_degree.get(h1, 0) + 1
        out_degree.setdefault(h2, 0)

    # Find start nodes (in-degree 0)
    starts = [h for h in in_degree if in_degree[h] == 0]

    # Build ancestral order by walking consensus graph
    anc_order = []
    visited = set()
    for start in starts:
        curr = start
        while curr and curr not in visited:
            if curr not in ancestor.telomeres:
                anc_order.append(curr)
            visited.add(curr)
            # Find the next HOG (outgoing edge)
            found_next = False
            for h1, h2 in consensus_adjs:
                if h1 == curr and h2 not in visited:
                    curr = h2
                    found_next = True
                    break
                if h2 == curr and h1 not in visited:
                    curr = h1
                    found_next = True
                    break
            if not found_next:
                break

    anc_idx = {h: i for i, h in enumerate(anc_order)}

    # For each child chromosome, map HOGs to ancestor order
    for chrom in mc.chromosomes:
        hogs = [h for h in chrom if h not in mc.telomeres]
        if len(hogs) < min_hogs:
            continue

        # Map positions in ancestor
        positions = [anc_idx.get(h, -1) for h in hogs]
        valid = [(i, positions[i]) for i in range(len(positions)) if positions[i] >= 0]
        if len(valid) < min_hogs:
            continue

        # Find the longest reversed segment
        # A reversed segment means positions in child are decreasing
        # when they should be increasing (in ancestor order)

        # Simple approach: find all decreasing consecutive pairs
        # and merge them into blocks
        rev_blocks = []
        i = 0
        while i < len(valid) - 1:
            j = i
            # Look for decreasing subsequence
            while j < len(valid) - 1 and valid[j][1] > valid[j + 1][1]:
                j += 1
            if j - i >= min_hogs - 1:  # at least min_hogs HOGs
                start_idx = valid[i][0]
                end_idx = valid[j][0]
                rev_hogs = hogs[start_idx:end_idx + 1]
                # Verify: check the whole segment is reversed
                seg_pos = [p for _, p in valid[i:j + 1]]
                if seg_pos == sorted(seg_pos, reverse=True) and seg_pos != sorted(seg_pos):
                    rev_blocks.append(rev_hogs)
            i = max(i + 1, j)
            if i >= len(valid):
                break

        # Also do sliding window for long contiguous segments
        if not rev_blocks:
            hog_pos_in_anc = [(h, anc_idx.get(h, -1)) for h in hogs if h in anc_idx]
            positions_list = [p for h, p in hog_pos_in_anc]
            if len(positions_list) >= min_hogs:
                sorted_pos = sorted(positions_list)
                for seg_len in range(len(positions_list), min_hogs - 1, -1):
                    found = False
                    for start in range(len(positions_list) - seg_len + 1):
                        seg = positions_list[start:start + seg_len]
                        if seg == sorted(seg, reverse=True) and seg != sorted(seg):
                            rev_hogs = [h for h, p in hog_pos_in_anc[start:start + seg_len]]
                            rev_blocks.append(rev_hogs)
                            found = True
                            break
                    if found:
                        break

        # Deduplicate and add events
        seen_blocks = set()
        for block_hogs in rev_blocks:
            block_set = frozenset(block_hogs)
            if block_set in seen_blocks or len(block_hogs) < min_hogs:
                continue
            seen_blocks.add(block_set)

            # Classify telomere vs internal
            is_tel = False
            if hogs:
                if block_hogs[0] == hogs[0] or block_hogs[-1] == hogs[-1]:
                    is_tel = True

            inv_type = 'telomere_inversion' if is_tel else 'internal_inversion'
            events.append(TAKREvent(
                event_type=inv_type,
                branch=branch,
                genes_involved=block_hogs,
                desc="%s: %d HOGs [%s..%s] in %s" % (
                    inv_type, len(block_hogs),
                    block_hogs[0].hog_id, block_hogs[-1].hog_id, src),
                support=len(block_hogs),
            ))

    return events


def _detect_unidir_trans(ancestor, ancest_adj_set, mc, child_adj_set,
                         branch, src):
    """Detect unidirectional translocations.

    An adjacency unique to the child where the two HOGs are on
    different chromosomes in the ancestor → segment moved.
    """
    from .takr_events import TAKREvent

    events = []
    mc_only = child_adj_set - ancest_adj_set

    def _hog_to_chrom(aag):
        chrom_of = {}
        for ci, chrom in enumerate(aag.chromosomes):
            for h in chrom:
                if h not in aag.telomeres:
                    chrom_of[h] = ci
        return chrom_of

    anc_chrom = _hog_to_chrom(ancestor)

    for (h1, h2) in mc_only:
        if h1 not in anc_chrom or h2 not in anc_chrom:
            continue
        if anc_chrom[h1] == anc_chrom[h2]:
            continue  # same chromosome, not a translocation

        events.append(TAKREvent(
            event_type='unidir_trans',
            branch=branch,
            genes_involved=[h1, h2],
            desc="unidir_trans: %s-%s cross-chrom in %s" % (
                h1.hog_id, h2.hog_id, src),
            support=1,
        ))

    return events


def _detect_eej(ancestor, mc, branch, src, og_adj_set=None):
    """Detect EEJ with outgroup polarity voting.

    EEJ = two separate chromosomes in the ancestor fused in the child.

    Detection: For each child adjacency (h1, h2), check if h1 and h2 are
    on DIFFERENT chromosomes in the consensus (separate connected components).

    Polarity (outgroup): If outgroup also has h1,h2 separate → this child's
    fusion is derived (true EEJ on this branch). If outgroup also has them
    fused → ancestral state is fused, this branch maintained it.
    """
    from .takr_events import TAKREvent

    events = []

    # Build consensus connected components
    consensus_graph = nx.Graph()
    for h1, h2 in ancestor.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        consensus_graph.add_edge(h1, h2)

    hog_to_comp = {}
    for ci, comp in enumerate(nx.connected_components(consensus_graph)):
        for h in comp:
            hog_to_comp[h] = ci

    for h1, h2 in mc.get_adjacencies(include_telomere=False):
        c1 = hog_to_comp.get(h1)
        c2 = hog_to_comp.get(h2)
        if c1 is not None and c2 is not None and c1 != c2:
            # Different consensus components → potential EEJ
            # Outgroup polarity check: confirm outgroup also has them separate
            key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
            if og_adj_set and key in og_adj_set:
                # Outgroup also has this cross-component adjacency
                # → ancestral state IS fused, skip (not EEJ on this branch)
                continue

            events.append(TAKREvent(
                event_type='eej',
                branch=branch,
                genes_involved=[h1, h2],
                desc="EEJ: %s-%s (comp %d+%d) in %s" % (
                    h1.hog_id, h2.hog_id, c1, c2, src),
                support=1,
            ))

    return events


def _detect_fission(ancestor, ancest_adj_set, mc, branch, src, og_adj_set=None):
    """Detect fission with outgroup polarity voting.

    Fission = an internal adjacency in the ancestor becomes two chromosome
    endpoints in the child.

    Polarity: if outgroup also has the adjacency broken → ancestral state is
    split, not fission. If outgroup maintains the adjacency → fission confirmed.
    """
    from .takr_events import TAKREvent

    events = []
    mc_adj = mc.get_adjacencies(include_telomere=False)
    mc_adj_set = set()
    for h1, h2 in mc_adj:
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        mc_adj_set.add(key)

    ancest_only = ancest_adj_set - mc_adj_set

    # Find child endpoints
    mc_ends = set()
    for chrom in mc.chromosomes:
        genes = [n for n in chrom if n not in mc.telomeres]
        if len(genes) >= 2:
            mc_ends.add(genes[0])
            mc_ends.add(genes[-1])

    for (h1, h2) in ancest_only:
        if h1 in mc_ends and h2 in mc_ends:
            # Both endpoints now split → potential fission
            # Outgroup check: if outgroup maintains this adjacency → fission confirmed
            key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
            if og_adj_set and key not in og_adj_set:
                # Outgroup also lost this adjacency → ancestral state is split, skip
                continue

            events.append(TAKREvent(
                event_type='fission',
                branch=branch,
                genes_involved=[h1, h2],
                desc="fission: %s-%s split in %s" % (h1.hog_id, h2.hog_id, src),
                support=1,
            ))

    return events


def _detect_rt(ancestor, ancest_adj_set, mc, child_adj_set, branch, src, og_adj_set=None):
    """Detect RT/URT via cross-connection pattern with outgroup polarity.

    RT: A-B + C-D → A-D + C-B

    Sibling comparison finds the cross pattern. Outgroup determines which
    side is derived: the child whose pattern DIFFERS from the outgroup
    is the branch where RT occurred.
    """
    from .takr_events import TAKREvent

    events = []
    ancest_only = ancest_adj_set - child_adj_set
    mc_only = child_adj_set - ancest_adj_set

    if len(ancest_only) < 2 or len(mc_only) < 2:
        return events

    for (a, b) in ancest_only:
        for (c, d) in ancest_only:
            if a == c and b == d:
                continue
            cross1 = (a, d) if a.hog_id < d.hog_id else (d, a)
            cross2 = (c, b) if c.hog_id < b.hog_id else (b, c)
            if cross1 in mc_only and cross2 in mc_only:
                # Cross-connection pattern found
                # Outgroup polarity check:
                # If outgroup has the SAME pattern as THIS child → ancestral,
                # the SIBLING has the derived pattern → RT is NOT on this branch
                if og_adj_set:
                    mc_cross = {cross1, cross2}
                    if mc_cross & og_adj_set:
                        # This child's pattern matches outgroup = ancestral
                        continue

                # Check for URT: one breakpoint at telomere
                mc_ends = set()
                for chrom in mc.chromosomes:
                    genes = [n for n in chrom if n not in mc.telomeres]
                    if genes:
                        mc_ends.add(genes[0])
                        mc_ends.add(genes[-1])
                has_tel = sum(1 for h in (a, b, c, d) if h in mc_ends)
                rt_type = 'unbalanced_reciprocal_translocation' if has_tel else 'reciprocal_translocation'

                events.append(TAKREvent(
                    event_type=rt_type,
                    branch=branch,
                    genes_involved=[a, b, c, d],
                    desc="reciprocal_translocation: %s,%s,%s,%s cross in %s" % (
                        a.hog_id, b.hog_id, c.hog_id, d.hog_id, src),
                    support=2,
                ))

    return events


def _detect_ncf(ancestor, mc, branch, src):
    """Detect NCF: one chromosome fully embedded inside another."""
    from .takr_events import TAKREvent

    events = []

    def _chrom_gene_sets(aag):
        return [set(h for h in chrom if h not in aag.telomeres)
                for chrom in aag.chromosomes]

    anc_sets = _chrom_gene_sets(ancestor)
    mc_sets = _chrom_gene_sets(mc)

    for mi, m_set in enumerate(anc_sets):
        for ci, c_set in enumerate(mc_sets):
            if not m_set or not c_set:
                continue
            if m_set.issubset(c_set) and len(m_set) < len(c_set):
                for mj, mj_set in enumerate(anc_sets):
                    if mi != mj and (mj_set & c_set):
                        rest = c_set - m_set
                        if rest.issubset(mj_set) and len(rest) >= 3:
                            events.append(TAKREvent(
                                event_type='ncf',
                                branch=branch,
                                genes_involved=list(m_set | mj_set),
                                desc="NCF: chr%d into chr%d in %s" % (mi, mj, src),
                                support=1,
                            ))

    return events


# =========================================================================
#  Phase 3: WGD pre→post detection
# =========================================================================

def _collapse_polyploid_event_driven(akr, leaf_id):
    """Collapse a polyploid leaf to pre-WGD using event-driven approach.

    Instead of CP-SAT (v3), split the polyploid graph into two subgenomes
    by finding sister chromosome pairs via shared HOGs, then treat them
    as two children for the standard event-driven pipeline.

    Returns pre-WGD graph or None.
    """
    from .AK import AncestralAdjacencyGraph

    post_graph = akr.leaf_graphs.get(leaf_id)
    if post_graph is None:
        return None

    ploidy = akr.ploidy_map.get(leaf_id, 2)
    if ploidy < 2:
        return None

    # Map to parent HOG level
    parent_id = None
    for node in akr.tree.traverse():
        if node.name == leaf_id and not node.is_root():
            parent_id = node.up.name
            break
    if parent_id is None:
        return None

    hog_level = parent_id
    mapped = akr._map_to_parent_hogs(hog_level, post_graph, source_id=leaf_id)

    # Build HOG copy → chromosome mapping
    # hog_copies[hog] = list of (chrom_idx, position, original_gene)
    hog_copies = {}
    if hasattr(mapped, 'hog_map') and mapped.hog_map:
        # Build chromosome index for each HOG record
        hog_record_to_chrom = {}
        for ci, chrom in enumerate(mapped.chromosomes):
            for h in chrom:
                if h not in mapped.telomeres:
                    hog_record_to_chrom[h] = ci

        for gene, hog in mapped.hog_map.items():
            if hog not in hog_copies:
                hog_copies[hog] = []
            chrom = hog_record_to_chrom.get(hog, -1)
            hog_copies[hog].append((chrom, gene))

    # For HOGs with 2+ copies, build chromosome pairing graph
    # chrom_pairs[chrA] = {chrB: count} where chrA and chrB share WGD pairs
    chrom_pairs = {}
    for hog, copies in hog_copies.items():
        if len(copies) >= 2:
            chroms = list(set(c for c, _ in copies))
            for i in range(len(chroms)):
                for j in range(i + 1, len(chroms)):
                    c1, c2 = chroms[i], chroms[j]
                    chrom_pairs.setdefault(c1, {})
                    chrom_pairs.setdefault(c2, {})
                    chrom_pairs[c1][c2] = chrom_pairs[c1].get(c2, 0) + 1
                    chrom_pairs[c2][c1] = chrom_pairs[c2].get(c1, 0) + 1

    # Greedy pairing: assign chromosomes to subgenomes
    # For each pair (chrA, chrB) with shared HOGs, put chrA in SG1, chrB in SG2
    sg1_chroms = set()
    sg2_chroms = set()
    assigned = set()

    all_chroms = sorted(chrom_pairs.keys(), key=lambda c: -sum(chrom_pairs[c].values()))

    for c in all_chroms:
        if c in assigned:
            continue
        # Find best partner
        partners = chrom_pairs.get(c, {})
        unassigned_partners = [(p, w) for p, w in sorted(partners.items(), key=lambda x: -x[1])
                               if p not in assigned]
        if unassigned_partners:
            best_p, _ = unassigned_partners[0]
            sg1_chroms.add(c)
            sg2_chroms.add(best_p)
            assigned.add(c)
            assigned.add(best_p)
        else:
            # No partner found, put in SG1
            sg1_chroms.add(c)
            assigned.add(c)

    # Any remaining unassigned chroms go to SG1
    for ci in range(len(list(mapped.chromosomes))):
        if ci not in assigned:
            sg1_chroms.add(ci)

    # Split the mapped graph into two subgenomes
    def _extract_subgenome(mapped, chrom_set):
        """Build a subgenome graph containing only specified chromosomes."""
        sg = AncestralAdjacencyGraph(node_id="%s_sg" % leaf_id)
        for ci, chrom in enumerate(mapped.chromosomes):
            if ci not in chrom_set:
                continue
            for h in chrom:
                if h not in mapped.telomeres:
                    sg.graph.add_node(h)
                    sg.gene_nodes.add(h)
            # Add edges within this chromosome
            for i in range(len(chrom) - 1):
                h1, h2 = chrom[i], chrom[i + 1]
                if h1 not in mapped.telomeres and h2 not in mapped.telomeres:
                    sg.graph.add_edge(h1, h2)
        sg._add_telomeres()
        return sg

    sg1 = _extract_subgenome(mapped, sg1_chroms)
    sg2 = _extract_subgenome(mapped, sg2_chroms)

    if len(sg1.gene_nodes) < 3 or len(sg2.gene_nodes) < 3:
        logger.warning("  Subgenomes too small for %s, falling back to v3", leaf_id)
        akr._collapse_polyploid_leaf(leaf_id)
        return akr.pre_wgd_graphs.get(leaf_id)

    # Run event-driven: build consensus from two subgenomes
    pre_wgd = _build_consensus_ancestor("pre-WGD " + leaf_id, [sg1, sg2],
                                         [leaf_id + "_sg1", leaf_id + "_sg2"])

    # Detect events between subgenomes (fractionation)
    sg1_adj = set()
    for h1, h2 in sg1.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        sg1_adj.add(key)
    sg2_adj = set()
    for h1, h2 in sg2.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        sg2_adj.add(key)

    events = []
    branch1 = "%s_preWGD-%s_sg1" % (leaf_id, leaf_id)
    branch2 = "%s_preWGD-%s_sg2" % (leaf_id, leaf_id)

    # Detect WGD event
    from .takr_events import TAKREvent
    events.append(TAKREvent(
        event_type='WGD',
        branch="%s_preWGD-%s" % (leaf_id, leaf_id),
        desc="WGD %dx at %s" % (ploidy, leaf_id),
        support=ploidy,
    ))

    # Fractionation: HOGs in one subgenome but not the other
    sg1_hogs = set(sg1.gene_nodes)
    sg2_hogs = set(sg2.gene_nodes)
    lost_in_sg1 = sg2_hogs - sg1_hogs
    for h in lost_in_sg1:
        events.append(TAKREvent(
            event_type='fractionation',
            branch=branch1,
            genes_involved=[h],
            desc="fractionation: %s lost in sg1" % h.hog_id,
            support=1,
        ))
    lost_in_sg2 = sg1_hogs - sg2_hogs
    for h in lost_in_sg2:
        events.append(TAKREvent(
            event_type='fractionation',
            branch=branch2,
            genes_involved=[h],
            desc="fractionation: %s lost in sg2" % h.hog_id,
            support=1,
        ))

    pre_wgd.events = events
    n_post = len(list(mapped.chromosomes))
    n_pre = len(list(pre_wgd.chromosomes))
    logger.info("  Polyploid leaf %s: %d post→%d pre chroms (event-driven)",
                leaf_id, n_post, n_pre)

    return pre_wgd, events

def _handle_wgd_node(akr, node, ploidy):
    """Handle WGD node: collapse post-WGD to pre-WGD, detect events.

    1. Get the post-WGD ancestor graph (reconstructed from children)
    2. Collapse it to pre-WGD (reduce chromosome count by ploidy)
    3. Detect events between pre-WGD and post-WGD
    4. Store events on virtual branch 'node_preWGD→node'
    """
    from .takr_events import TAKREvent, branch_id

    node_id = node.name
    post_graph = akr.anc_graphs.get(node_id)
    if post_graph is None:
        return

    # Collapse post-WGD to pre-WGD using v3's method (reuse)
    akr._collapse_wgd_v3(node)
    pre_graph = akr.pre_wgd_graphs.get(node_id)
    if pre_graph is None:
        logger.warning("  WGD collapse failed for %s", node_id)
        return

    # Map both to the same HOG level for comparison
    hog_level = node_id
    pre_mapped = akr._map_to_parent_hogs(hog_level, pre_graph, source_id=node_id + "_pre")
    post_mapped = akr._map_to_parent_hogs(hog_level, post_graph, source_id=node_id)

    # Detect events between pre-WGD and post-WGD on virtual branch
    virtual_branch = "%s_preWGD-%s" % (node_id, node_id)
    events = []

    # 1. WGD event
    events.append(TAKREvent(
        event_type='WGD',
        branch=virtual_branch,
        genes_involved=[],
        desc="WGD %dx at %s" % (ploidy, node_id),
        support=ploidy,
    ))

    # 2. Structural event detection (same algorithm as branch-level)
    ancest_adj = set()
    for h1, h2 in pre_mapped.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        ancest_adj.add(key)

    post_adj = set()
    for h1, h2 in post_mapped.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        post_adj.add(key)

    wgd_events = _detect_branch_events(
        pre_mapped, post_mapped, virtual_branch, node_id + "_post", 3)

    # For WGD, also detect fractionation (gene loss between pre and post)
    pre_hogs = set(pre_mapped.gene_nodes)
    post_hogs = set(post_mapped.gene_nodes)
    lost_hogs = pre_hogs - post_hogs
    if lost_hogs:
        wgd_events.append(TAKREvent(
            event_type='fractionation',
            branch=virtual_branch,
            genes_involved=list(lost_hogs),
            desc="fractionation: %d HOGs lost post-WGD in %s" % (len(lost_hogs), node_id),
            support=len(lost_hogs),
        ))

    events.extend(wgd_events)

    # Store events
    if pre_graph.events is None:
        pre_graph.events = []
    pre_graph.events.extend(events)

    logger.info("  WGD %s: %d post→%d pre chroms, %d events on %s",
                node_id, len(list(post_graph.chromosomes)),
                len(list(pre_graph.chromosomes)), len(events), virtual_branch)


# =========================================================================
#  Obsolete / placeholder functions kept for reference
# =========================================================================

def _merge_graphs(aag1, aag2):
    """Merge two graphs (placeholder for future use)."""
    pass


__all__ = ['reconstruct_event_driven']
