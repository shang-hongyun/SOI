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


def _log_child_events(events, phase_label, child_id=None):
    """统一的事件 log 格式。按孩子分组输出。

    Args:
        events: TAKREvent 列表
        phase_label: 阶段标签（如 "Phase 1", "Phase 4a"）
        child_id: 如果指定，只输出该孩子的事件；否则按孩子分组
    """
    from collections import defaultdict
    # 按孩子分组
    child_events = defaultdict(list)
    for e in events:
        # 从 branch 提取 child_id: "HOG-node-Sp_1" → "Sp_1"
        parts = e.branch.split('-')
        cid = parts[-1] if parts else 'unknown'
        if child_id and cid != child_id:
            continue
        child_events[cid].append(e)

    for cid in sorted(child_events.keys()):
        evts = child_events[cid]
        # 按类型分组
        type_counts = defaultdict(int)
        type_lens = defaultdict(list)
        for e in evts:
            type_counts[e.event_type] += 1
            type_lens[e.event_type].append(len(e.genes_involved))
        # 格式化
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

def reconstruct_event_driven_v2(akr, min_hogs=3):
    """Event-driven reconstruction v2 -- ColoredGraph: colored graph + cycle detection + path cover."""
    from .takr_colored_graph import ColoredGraph

    def _wgd_collapse(post_graph, node_id, ploidy, parent_hog_level, min_hogs=3):
        """WGD collapse: pair post-WGD chromosomes by HOG similarity → pre-WGD."""
        logger.info("Reconstructing node %s_preWGD [v2 ColoredGraph]", node_id)
        t_collapse = time.time()

        # 用原始染色体路径（不 path_cover）
        child_chroms = post_graph._child_chromosomes.get(node_id, [])
        n_orig = len(child_chroms)
        logger.info("  Children: %s=%d", node_id, n_orig)

        if ploidy <= 1 or n_orig < ploidy:
            # 没有 WGD，直接映射到 parent HOG level
            pre_anc = akr._map_to_parent_hogs(parent_hog_level, post_graph,
                                               source_id="{}_pre".format(node_id))
            pre_anc.node_id = "{}_pre".format(node_id)
            logger.info("  Done: no WGD (ploidy=%d), 0 chroms (%.1fs)",
                         ploidy, time.time() - t_collapse)
            return pre_anc, []

        # 按原始染色体配对
        import itertools
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

        # 构建 pre-WGD 图
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

        if parent_hog_level != node_id:
            saved_events = list(all_events)
            pre_anc = akr._map_to_parent_hogs(parent_hog_level, pre_anc,
                                               source_id="{}_pre".format(node_id))
            pre_anc.events = saved_events
        pre_anc.node_id = "{}_pre".format(node_id)
        logger.info("  Done: %d -> %d chroms, %d events (%.1fs)",
                     n_orig, n_pre, len(all_events), time.time() - t_collapse)
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
                child_source_ids.append("{}_pre".format(cid))
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

            # ── dedup 后验证（重建后） ──
            ok, errors = _validate_dedup(mc, cid, pre_dedup_chroms[cid])
            if ok:
                logger.info("  [Phase 1] %s: dedup validated (chrom count ok, no dups, no self-loops)", cid)
            else:
                for err in errors:
                    logger.error("  [Phase 1] %s: %s", cid, err)

            G.add_child(cid, mc)

        # dedup 事件汇总（统一格式）
        dedup_events = [e for e in G.events
                        if e.event_type in ('tandem_dup', 'dispersed_dup')]
        if dedup_events:
            _log_child_events(dedup_events, "Phase 1")

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

