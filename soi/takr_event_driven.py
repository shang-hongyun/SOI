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
import time
from collections import defaultdict
from typing import Dict, List, Optional

import networkx as nx

logger = logging.getLogger(__name__)


# =========================================================================
#  Main entry point
# =========================================================================

def reconstruct_event_driven(akr, min_hogs=3):
    """Event-driven reconstruction — main entry.

    Args:
        akr: AKR instance (tree, leaf_graphs, hog already built)
        min_hogs: Minimum HOGs for a valid event (noise filter)

    Returns:
        {node_name: AncestralAdjacencyGraph}
    """
    logger.info("=== Event-driven reconstruction (v4) ===")
    t0 = time.time()

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

        # Detect events on each child's branch
        all_events = []
        for ci, (mc, src) in enumerate(zip(mapped_children, child_source_ids)):
            branch = "%s-%s" % (node_id, src)
            events = _detect_branch_events(
                ancestor, mc, branch, src, min_hogs,
                mapped_children, child_source_ids)
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

    logger.info("=== Event-driven reconstruction done (%.1fs) ===",
                time.time() - t0)
    return akr.anc_graphs


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

    # Consensus: adjacencies in >= ceil(n_children/2) children
    threshold = max(n_children // 2, 1)
    adj_counts = defaultdict(int)
    for adjs in child_adj_sets:
        for key in adjs:
            adj_counts[key] += 1

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
                          all_mapped_children=None, all_child_ids=None):
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

    # EEJ: telomeres adjacent in ancestor but NOT in child
    eej = _detect_eej(ancestor, mc, branch, src)
    events.extend(eej)

    # Fission: internal adjacency in ancestor → chromosome endpoints in child
    fission = _detect_fission(ancestor, ancest_adj_set, mc, branch, src)
    events.extend(fission)

    # NCF
    ncf = _detect_ncf(ancestor, mc, branch, src)
    events.extend(ncf)

    # RT/URT
    rt = _detect_rt(ancestor, ancest_adj_set, mc, child_adj_set, branch, src)
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


def _detect_eej(ancestor, mc, branch, src):
    """Detect EEJ via connected-component comparison.

    Key insight: EEJ fuses two chromosomes. If two HOGs are on DIFFERENT
    chromosomes in the consensus (separate connected components), but
    ADJACENT in the child, then the two chromosomes were fused in the
    child's lineage = EEJ on this branch.

    This avoids dependency on telomere markers entirely — uses graph
    connectivity instead.
    """
    from .takr_events import TAKREvent

    events = []

    # Build consensus connected components from adjacency graph
    # Two HOGs are in the same component if connected by shared adjacencies
    consensus_graph = nx.Graph()
    for h1, h2 in ancestor.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        consensus_graph.add_edge(h1, h2)

    # Map each HOG to its component ID
    hog_to_comp = {}
    for ci, comp in enumerate(nx.connected_components(consensus_graph)):
        for h in comp:
            hog_to_comp[h] = ci

    # For each child adjacency, check if the two HOGs are in
    # DIFFERENT consensus components → EEJ
    for h1, h2 in mc.get_adjacencies(include_telomere=False):
        key = (h1, h2) if h1.hog_id < h2.hog_id else (h2, h1)
        c1 = hog_to_comp.get(h1)
        c2 = hog_to_comp.get(h2)
        if c1 is not None and c2 is not None and c1 != c2:
            # These two HOGs are from different ancestral chromosomes
            # but adjacent in the child → EEJ
            events.append(TAKREvent(
                event_type='eej',
                branch=branch,
                genes_involved=[h1, h2],
                desc="EEJ: %s-%s (comp %d+%d) in %s" % (
                    h1.hog_id, h2.hog_id, c1, c2, src),
                support=1,
            ))

    return events


def _detect_fission(ancestor, ancest_adj_set, mc, branch, src):
    """Detect fission: ancestral internal adjacency → child endpoints."""
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
            events.append(TAKREvent(
                event_type='fission',
                branch=branch,
                genes_involved=[h1, h2],
                desc="fission: %s-%s split in %s" % (h1.hog_id, h2.hog_id, src),
                support=1,
            ))

    return events


def _detect_rt(ancestor, ancest_adj_set, mc, child_adj_set, branch, src):
    """Detect RT/URT via cross-connection pattern.

    RT: A-B + C-D → A-D + C-B
    Detection: (A,B) in ancestor, (C,D) in ancestor
               (A,D) in child,   (C,B) in child
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
            # Check for cross in child
            cross1 = (a, d) if a.hog_id < d.hog_id else (d, a)
            cross2 = (c, b) if c.hog_id < b.hog_id else (b, c)
            if cross1 in mc_only and cross2 in mc_only:
                # Determine URT: one breakpoint at telomere
                mc_ends = set()
                for chrom in mc.chromosomes:
                    genes = [n for n in chrom if n not in mc.telomeres]
                    if genes:
                        mc_ends.add(genes[0])
                        mc_ends.add(genes[-1])

                has_tel = sum(1 for h in (a, b, c, d) if h in mc_ends)
                rt_type = 'unbalanced_reciprocal_translocation' if has_tel else 'reciprocal_translocation'
                hog_ids = ','.join(str(h.hog_id) for h in (a, b, c, d))
                events.append(TAKREvent(
                    event_type=rt_type,
                    branch=branch,
                    genes_involved=[a, b, c, d],
                    desc="%s: %s cross in %s" % (rt_type, hog_ids, src),
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
