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

        # Collect child graphs
        child_graphs, child_source_ids = [], []
        for child in node.children:
            cid = child.name
            if child.is_leaf():
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
                ancestor, mc, branch, src, min_hogs)
            all_events.extend(events)

        ancestor.events = all_events
        akr.anc_graphs[node_id] = ancestor

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

def _detect_branch_events(ancestor, mc, branch, src, min_hogs):
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
    """Detect inversions via paired breakpoints.

    An inversion of block X..Y means:
    - Edges F-X and Y-G (involving block boundaries) are broken
    - Edges F-Y and X-G appear instead
    - The block X..Y appears reversed in the child

    Detection:
    1. Find all ancest_only edges
    2. Group them by child chrom
    3. For each chrom, check if boundary edges form consecutive breakpoints
    4. The HOGs between consecutive breakpoints = inversion block
    """
    from .takr_events import TAKREvent

    events = []
    ancest_only = ancest_adj_set - child_adj_set
    if not ancest_only:
        return events

    # Build chromosome HOG lists (child)
    mc_chrom_hogs = []
    for chrom in mc.chromosomes:
        hog_list = [h for h in chrom if h not in mc.telomeres]
        if hog_list:
            mc_chrom_hogs.append(hog_list)

    # Build child position map: hog -> (chrom_idx, position)
    mc_pos = {}
    for ci, hogs in enumerate(mc_chrom_hogs):
        for j, h in enumerate(hogs):
            mc_pos[h] = (ci, j)

    # For each ancest_only edge, collect breakpoints per child chromosome
    for (h1, h2) in ancest_only:
        if h1 not in mc_pos or h2 not in mc_pos:
            continue
        c1, p1 = mc_pos[h1]
        c2, p2 = mc_pos[h2]
        if c1 != c2:
            continue  # different chroms = not an inversion (more like translocation)

        # h1-h2 were adjacent in ancestor, now separated in child
        # The HOGs between them in ancestor = inversion candidate
        # Find positions in ancestor
        anc_hogs = []
        for chrom in ancestor.chromosomes:
            anc_hogs = [h for h in chrom if h not in ancestor.telomeres]
            break

        anc_idx = {h: i for i, h in enumerate(anc_hogs)}
        if h1 not in anc_idx or h2 not in anc_idx:
            continue

        ai1, ai2 = anc_idx[h1], anc_idx[h2]
        start, end = min(ai1, ai2), max(ai1, ai2)
        block = anc_hogs[start:end + 1]

        # The inversion block is the INTERNAL part (between h1 and h2)
        if len(block) < 2:
            continue

        # Check if block is reversed in child by comparing order
        mc_order = []
        seen = set()
        mc_hog_set = set(mc_chrom_hogs[c1])
        for h in block:
            if h in mc_hog_set:
                pos = mc_pos[h][1]
                mc_order.append((h, pos))
                seen.add(h)

        # If less than min_hogs survived, skip
        if len(mc_order) < min_hogs and len(mc_order) < len(block):
            continue

        # Check if MC order is reversed relative to ancestor
        # For a real inversion, the order is fully or partially reversed
        ancestor_order = [pos for h in block if h in mc_hog_set]
        mc_order_only = [p for h, p in mc_order]
        if not mc_order_only:
            continue

        # An inversion means mc_order is different from ancestor order
        # Simplest check: first HOG in ancestor is NOT first in child
        first_anc = block[0]
        is_reversed = False
        if first_anc in mc_hog_set:
            mc_first_pos = mc_pos.get(first_anc, (0, -1))[1]
            # Check if any ancestor-before-first HOG appears AFTER first_anc in child
            for h in block[1:min(len(block), 5)]:
                if h in mc_hog_set and mc_pos[h][1] < mc_first_pos:
                    is_reversed = True
                    break

        if not is_reversed:
            continue

        # Classify: telomere_inversion if block endpoint at telomere
        is_telomere = False
        if mc_chrom_hogs[c1]:
            first_gene = mc_chrom_hogs[c1][0]
            last_gene = mc_chrom_hogs[c1][-1]
            if block[0] == first_gene or block[-1] == last_gene:
                is_telomere = True

        # Also check if block endpoints reach telomeres
        if p1 == 0 or p2 == len(mc_chrom_hogs[c1]) - 1:
            is_telomere = True

        inv_type = 'telomere_inversion' if is_telomere else 'internal_inversion'

        # Deduplicate: check if this block was already detected
        block_set = set(block)
        is_dup = False
        for ev in events:
            if ev.event_type.startswith('inversion') and set(ev.genes_involved) == block_set:
                is_dup = True
                break
        if is_dup:
            continue

        events.append(TAKREvent(
            event_type=inv_type,
            branch=branch,
            genes_involved=block,
            desc="%s: %d HOGs [%s..%s] in %s" % (
                inv_type, len(block), block[0].hog_id, block[-1].hog_id, src),
            support=len(block),
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
    """Detect EEJ: telomeres adjacent in ancestor but not in child."""
    from .takr_events import TAKREvent

    events = []
    anc_tel_adjs = set(ancestor.get_telomere_adjacencies())
    mc_tel_adjs = set(mc.get_telomere_adjacencies())

    lost = anc_tel_adjs - mc_tel_adjs
    for (t1, t2) in lost:
        events.append(TAKREvent(
            event_type='eej',
            branch=branch,
            genes_involved=[t1, t2],
            desc="EEJ: %s-%s in %s" % (t1, t2, src),
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
#  Obsolete / placeholder functions kept for reference
# =========================================================================

def _merge_graphs(aag1, aag2):
    """Merge two graphs (placeholder for future use)."""
    pass


__all__ = ['reconstruct_event_driven']
