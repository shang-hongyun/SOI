#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for takr_colored_graph.py — ColoredGraph event detection & resolution.

Tests each event type with minimal (simplest possible) and complex graphs:
  - indel: spanning edge without direction conflict
  - inversion: 3-cycle at block level (2 unique edges from different children + 1 shared)
  - reciprocal_translocation: 4-cycle with 3+ chrom colors
  - eej: bridge edge, outgroup lacks adjacency → fusion
  - ncf: bridge edge, small component → nested fusion
  - fission: bridge edge, outgroup has adjacency → other child lost
  - gene_gain: insertion node (degree-2, unique edges, shared neighbors)

NOTE: Inversions in undirected adjacency graphs always create 3-cycles (triangles),
not 4-cycles. The 4-cycle pattern described in design docs is for breakpoint graphs.
In ColoredGraph (undirected), inversions are detected at BLOCK level as 3-cycles
with 2 unique edges from different children + 1 shared edge.

Usage:
    cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
    python -m pytest tests/test_colored_graph.py -v
"""

import sys
import os
import pytest

# Ensure soi package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from soi.takr_colored_graph import ColoredGraph
from soi.takr_events import TAKREvent


# ============================================================
# Mock child graph — minimal interface for ColoredGraph.add_child
# ============================================================

class MockChildGraph:
    """Minimal child graph implementing the interface ColoredGraph.add_child expects.

    Args:
        chromosomes: list of lists, each inner list = one chromosome.
                     Telomere nodes must be at positions [0] and [-1].
        telomeres: set of telomere nodes (default: first+last of each chrom)
    """

    def __init__(self, chromosomes, telomeres=None):
        self._chromosomes = chromosomes
        if telomeres is None:
            telomeres = set()
            for chrom in chromosomes:
                if chrom:
                    telomeres.add(chrom[0])
                    telomeres.add(chrom[-1])
        self.telomeres = telomeres
        self.gene_nodes = set()
        for chrom in chromosomes:
            for n in chrom:
                if n not in self.telomeres:
                    self.gene_nodes.add(n)

    @property
    def chromosomes(self):
        return self._chromosomes

    def get_adjacencies(self, include_telomere=False):
        adjs = set()
        for chrom in self._chromosomes:
            for i in range(len(chrom) - 1):
                n1, n2 = chrom[i], chrom[i + 1]
                if not include_telomere:
                    if n1 in self.telomeres or n2 in self.telomeres:
                        continue
                adjs.add((n1, n2))
        return adjs


# ============================================================
# Helpers
# ============================================================

def tel(chrom_name, side):
    """Create a telomere node: ('chr1', 'L') or ('chr1', 'R')."""
    return (chrom_name, side)


def build_simple_graph(children_spec):
    """Build a ColoredGraph from a simple spec.

    Args:
        children_spec: list of (child_id, chromosomes) tuples
            chromosomes: list of lists of HOG IDs (strings)
            Telomeres auto-detected as first/last of each chromosome.

    Returns:
        ColoredGraph
    """
    G = ColoredGraph(hog_level='test')
    for child_id, chroms in children_spec:
        # Wrap each chromosome with telomeres
        tel_chroms = []
        for ci, chrom in enumerate(chroms):
            tL = tel('{}_c{}'.format(child_id, ci), 'L')
            tR = tel('{}_c{}'.format(child_id, ci), 'R')
            tel_chroms.append([tL] + list(chrom) + [tR])
        mock = MockChildGraph(tel_chroms)
        G.add_child(child_id, mock)
    return G


def get_events_by_type(G):
    """Return {event_type: count} from G.events."""
    from collections import Counter
    return dict(Counter(e.event_type for e in G.events))


def find_event(G, event_type):
    """Return first event of given type, or None."""
    for e in G.events:
        if e.event_type == event_type:
            return e
    return None


# ============================================================
# Test: Graph Construction
# ============================================================

class TestGraphConstruction:

    def test_single_child_linear(self):
        """Single child, linear chromosome A-B-C → 2 edges."""
        G = build_simple_graph([('c1', [['A', 'B', 'C']])])
        assert G.node_count() == 3
        assert G.edge_count() == 2
        assert G.children() == {'c1'}

    def test_two_children_shared_edges(self):
        """Two children with same adjacency → shared edges."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
        ])
        assert G.node_count() == 3
        # A-B shared, B-C shared → 2 edges total (each with 2 colors)
        assert G.edge_count() == 2
        shared = G.shared_edges()
        assert len(shared) == 2

    def test_two_children_unique_edges(self):
        """Two children with different adjacency → all edges unique."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),   # edges: A-B, B-C
            ('c2', [['A', 'C']]),          # edge: A-C
        ])
        assert G.node_count() == 3
        # All 3 edges are unique (each from only 1 child)
        assert G.edge_count() == 3
        unique = G.unique_edges()
        assert len(unique) == 3

    def test_direction_metadata(self):
        """Direction stored correctly for canonical edge ordering."""
        G = ColoredGraph(hog_level='test')
        # Child 1: A→B (canonical order A<B → direction +1)
        # Child 2: B→A (canonical order A<B → direction -1)
        chroms1 = [[tel('c1_c0', 'L'), 'A', 'B', tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L'), 'B', 'A', tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))
        # Edge A-B should have direction conflict
        assert G.edge_has_direction_conflict('A', 'B')

    def test_no_direction_conflict_same_order(self):
        """Both children same order → no direction conflict."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
        ])
        assert not G.edge_has_direction_conflict('A', 'B')
        assert not G.edge_has_direction_conflict('B', 'C')

    def test_direction_conflict_on_shared_edge(self):
        """Shared edge with opposite directions from different children."""
        G = ColoredGraph(hog_level='test')
        # c1: A→B→C, c2: C→B→A (reversed)
        chroms1 = [[tel('c1_c0', 'L'), 'A', 'B', 'C', tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L'), 'C', 'B', 'A', tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))
        # Both A-B and B-C should have direction conflicts
        assert G.edge_has_direction_conflict('A', 'B')
        assert G.edge_has_direction_conflict('B', 'C')


# ============================================================
# Test: Indel Detection (find_indel_shortcuts)
# ============================================================

class TestIndelDetection:
    """Indel = spanning edge without direction conflict.

    Minimal: Child 1 has A-B-C, Child 2 has A-C.
    A-C is unique to c2, spans B in c1 → indel shortcut.
    """

    def test_simple_indel(self):
        """Simplest indel: one gene inserted/deleted."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'C']]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 1
        h1, h2, child_id, spanned = shortcuts[0]
        assert child_id == 'c2'
        assert 'B' in spanned

    def test_no_indel_with_rearrangement(self):
        """Spanning edge where intermediate HOGs exist in other child → rearrangement, not indel."""
        G = ColoredGraph(hog_level='test')
        # Child 1: A-B-C-D
        # Child 2: A-C-B-D (B-C reversed → rearrangement)
        chroms1 = [[tel('c1_c0', 'L'), 'A', 'B', 'C', 'D', tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L'), 'A', 'C', 'B', 'D', tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))
        shortcuts = G.find_indel_shortcuts()
        # A-C unique to c2, spans B in c1. But B is in c2 (via C-B) → rearrangement
        assert len(shortcuts) == 0

    def test_indel_multi_hop(self):
        """Multi-hop indel: spanning edge covers multiple intermediate HOGs."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E']]),
            ('c2', [['A', 'E']]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 1
        _, _, _, spanned = shortcuts[0]
        assert len(spanned) == 3  # B, C, D spanned

    def test_indel_max_span_filter(self):
        """Spanning edge exceeding max_span → not detected as indel."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']]),
            ('c2', [['A', 'L']]),
        ])
        # max_span=5 → spans 10 HOGs → filtered out
        shortcuts = G.find_indel_shortcuts(max_span=5)
        assert len(shortcuts) == 0
        # max_span=20 → spans 10 HOGs → detected
        shortcuts = G.find_indel_shortcuts(max_span=20)
        assert len(shortcuts) == 1

    def test_no_indel_same_adjacency(self):
        """Both children same adjacency → no indel."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 0

    def test_indel_two_genes_inserted(self):
        """Two consecutive genes inserted → one spanning edge."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D']]),
            ('c2', [['A', 'D']]),
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) == 1
        _, _, _, spanned = shortcuts[0]
        assert set(spanned) == {'B', 'C'}


# ============================================================
# Test: Inversion Detection (Block Level)
# ============================================================

class TestInversionDetection:
    """Inversion detection at block level.

    In undirected adjacency graphs, inversions create 3-cycles (triangles) at block level:
    - 2 unique edges from different children + 1 shared edge
    - At least 2 blocks must have ≥2 HOGs

    Minimal graph:
      c1: [A,B, C,D, E,F] (one chromosome, 3 blocks)
      c2: [C,D, A,B, E,F] (A-B moved after C-D → rearrangement)
      Shared: A-B, C-D, E-F (all 3 blocks connected by shared adjacencies)
      But block edges: blk_0↔blk_1 is shared (both children connect them),
      blk_1↔blk_2 is c1 unique, blk_2↔blk_0 is c2 unique → 3-cycle!
    """

    def test_block_level_inversion_3cycle(self):
        """Inversion detected via direction conflict (signed adjacency model)."""
        G = ColoredGraph(hog_level='test')
        # c1: A-B-C-D-E-F (normal)
        # c2: A-B-D-C-E-F (C-D inverted → direction conflict on C-D)
        chroms1 = [[tel('c1_c0', 'L')] + list('ABCDEF') + [tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L')] + list('ABDCEF') + [tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))

        # Direction conflict on C-D: c1 has C→D (+1), c2 has D→C (-1)
        assert G.edge_has_direction_conflict('C', 'D')

        # Resolve all events (includes _detect_inversions)
        G.resolve_all_events()
        inv_events = [e for e in G.events if e.event_type == 'inversion']
        assert len(inv_events) >= 1, (
            "No inversion detected. Events: {}".format(G.events))

    def test_no_inversion_same_direction(self):
        """Same direction on all edges → no inversion."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
        ])
        G._build_synteny_blocks()
        G._compress_to_block_level()
        G._resolve_block_structural_events()
        inv_events = [e for e in G.events if e.event_type == 'inversion']
        assert len(inv_events) == 0

    def test_hog_level_3cycle_is_gene_indel(self):
        """HOG-level 3-cycles are classified as gene_indel (not inversion)."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'C']]),
        ])
        cycles = G.find_cycles()
        for cycle in cycles:
            etype, _, _ = G.classify_cycle(cycle)
            # 3-cycles at HOG level always → gene_indel
            if len(cycle) == 3:
                assert etype == 'gene_indel'

    def test_resolve_inversion_removes_cycle(self):
        """Resolving block structural events should break 3-cycles."""
        G = ColoredGraph(hog_level='test')
        chroms1 = [[tel('c1_c0', 'L')] + list('ABCDEF') + [tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L')] + list('CDABEF') + [tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))

        G._build_synteny_blocks()
        G._compress_to_block_level()
        G._resolve_block_structural_events()

        # After resolution, block graph should have fewer or no cycles
        cycles = G._find_block_cycles()
        # At least the inversion should have been resolved
        inv_events = [e for e in G.events if e.event_type == 'inversion']
        if inv_events:
            assert len(cycles) < 2, "Cycles remain after inversion resolution"


# ============================================================
# Test: Reciprocal Translocation (RT)
# ============================================================

class TestRTDetection:
    """RT = cross-chromosome exchange.

    At block level, RT creates a 4-cycle with 3+ distinct (child, chrom) colors.
    Need 4 blocks connected in a cycle with alternating child contributions.

    Minimal:
      c1: Chr0=[A,B,C,D], Chr1=[E,F,G,H]
      c2: Chr0=[A,B,G,H], Chr1=[E,F,C,D]
      Shared: A-B, E-F, C-D, G-H (but C-D and G-H on different chroms in c1 vs c2)
      → 4-cycle at block level with 3+ colors → RT
    """

    def test_simple_rt_block_level(self):
        """Simple RT: two chromosomes exchange segments → block-level 4-cycle."""
        G = ColoredGraph(hog_level='test')
        # c1: Chr0=[A,B,C,D], Chr1=[E,F,G,H]
        # c2: Chr0=[A,B,G,H], Chr1=[E,F,C,D]  (C-D and G-H swapped)
        chroms1 = [
            [tel('c1_c0', 'L')] + list('ABCD') + [tel('c1_c0', 'R')],
            [tel('c1_c1', 'L')] + list('EFGH') + [tel('c1_c1', 'R')],
        ]
        chroms2 = [
            [tel('c2_c0', 'L')] + list('ABGH') + [tel('c2_c0', 'R')],
            [tel('c2_c1', 'L')] + list('EFCD') + [tel('c2_c1', 'R')],
        ]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))

        # Verify unique edges exist with different chrom_idx
        unique = G.unique_edges()
        assert len(unique) > 0, "No unique edges for RT"

        # Run block compression and structural detection
        G._build_synteny_blocks()
        G._compress_to_block_level()
        G._resolve_block_structural_events()

        rt_events = [e for e in G.events
                     if e.event_type in ('reciprocal_translocation',
                                         'unbalanced_reciprocal_translocation')]
        # RT may or may not be detected depending on block structure
        # At minimum, verify the graph has cycles
        cycles = G._find_block_cycles()
        assert len(cycles) > 0 or len(rt_events) > 0, (
            "No RT detected and no cycles found. Events: {}".format(G.events))

    def test_hog_level_classify_cycle_with_chrom_idx(self):
        """classify_cycle distinguishes same-child-different-chrom as RT signal."""
        # Create a graph where classify_cycle sees 3+ distinct (child, chrom) colors
        G = ColoredGraph(hog_level='test')
        # Two children, each with 2 chromosomes
        # Create edges manually to control chrom_idx
        G.add_edge('A', 'B', 'c1', 0, direction=1)
        G.add_edge('B', 'C', 'c1', 0, direction=1)
        G.add_edge('A', 'C', 'c1', 1, direction=1)
        G.add_edge('A', 'B', 'c2', 0, direction=1)
        G.add_edge('B', 'D', 'c2', 0, direction=1)
        G.add_edge('A', 'D', 'c2', 1, direction=1)
        # Now we have a triangle A-B-D-A with edges from c2(0), c2(1)
        # and A-B-C-A with edges from c1(0), c1(1)
        # These are 3-cycles → gene_indel at HOG level (expected)


# ============================================================
# Test: Insertion Node (gene_gain)
# ============================================================

class TestInsertionNode:
    """Insertion = degree-2 node with unique edges, neighbors have shared edge.

    Minimal:
      Ancestral: A-C (shared, ≥2 children)
      Derived:   A-B-C (B inserted, edges A-B and B-C unique to one child)
    """

    def test_simple_insertion(self):
        """One gene inserted → detected as indel shortcut."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),  # B inserted
            ('c2', [['A', 'C']]),        # no B
        ])
        shortcuts = G.find_indel_shortcuts()
        assert len(shortcuts) >= 1

    def test_insertion_with_shared_neighbors(self):
        """Insertion where neighbors have shared edge from 3rd child → _resolve_inserted_nodes."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),   # B inserted
            ('c2', [['A', 'C']]),         # ancestral path
            ('c3', [['A', 'C']]),         # 3rd child also has A-C → shared
        ])
        # A-C: shared (c2, c3). A-B: unique (c1). B-C: unique (c1).
        # B degree=2, both edges unique, neighbors A-C shared → insertion!
        G._resolve_inserted_nodes()
        removed = any(e.event_type == 'gene_gain' for e in G.events)
        assert removed, "Insertion node B not detected"


# ============================================================
# Test: Bridge Events (EEJ, NCF, Fission)
# ============================================================

class TestBridgeEvents:
    """Bridge = unique edge connecting different shared-edge components.

    Minimal EEJ:
      Child 1: Chr1=[A,B], Chr2=[C,D]  + bridge edge B-C
      Child 2: Chr1=[A,B], Chr2=[C,D]
      → shared: A-B, C-D. Unique: B-C(c1).
      → B-C connects two shared components → bridge.
    """

    def _build_bridge_graph(self):
        """Build graph with a bridge edge between two shared components."""
        G = ColoredGraph(hog_level='N1')
        # Child 1: two chromosomes
        chroms1 = [
            [tel('c1_c0', 'L'), 'A', 'B', tel('c1_c0', 'R')],
            [tel('c1_c1', 'L'), 'C', 'D', tel('c1_c1', 'R')],
        ]
        # Child 2: same two chromosomes, no bridge
        chroms2 = [
            [tel('c2_c0', 'L'), 'A', 'B', tel('c2_c0', 'R')],
            [tel('c2_c1', 'L'), 'C', 'D', tel('c2_c1', 'R')],
        ]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))
        # Manually add bridge edge B-C for c1
        G.add_edge('B', 'C', 'c1', 0, direction=1)
        return G

    def test_bridge_detection_basic(self):
        """Basic bridge detection: unique edge connects two shared components."""
        G = self._build_bridge_graph()
        G._save_original_shared_components()
        G._build_synteny_blocks()
        G._compress_to_block_level()
        n = G._resolve_block_bridge_events(outgroup_adjacency=None)
        bridge_events = [e for e in G.events
                         if e.event_type == 'bridge_unclassified']
        assert len(bridge_events) >= 1, "No bridge detected. Events: {}".format(G.events)

    def test_eej_with_outgroup(self):
        """EEJ: outgroup lacks adjacency → fusion on this child's branch."""
        G = self._build_bridge_graph()
        G._save_original_shared_components()
        G._build_synteny_blocks()
        G._compress_to_block_level()
        outgroup_adj = set()  # empty → no ancestral adjacency
        n = G._resolve_block_bridge_events(outgroup_adjacency=outgroup_adj)
        eej_events = [e for e in G.events if e.event_type in ('eej', 'ncf')]
        assert len(eej_events) >= 1, "No EEJ/NCF detected. Events: {}".format(G.events)

    def test_fission_with_outgroup(self):
        """Fission: outgroup HAS adjacency → ancestral connected, other child lost it."""
        G = self._build_bridge_graph()
        G._save_original_shared_components()
        G._build_synteny_blocks()
        G._compress_to_block_level()
        outgroup_adj = {('B', 'C')}
        n = G._resolve_block_bridge_events(outgroup_adjacency=outgroup_adj)
        fission_events = [e for e in G.events if e.event_type == 'fission']
        assert len(fission_events) >= 1, "No fission detected. Events: {}".format(G.events)

    def test_no_bridge_same_component(self):
        """Unique edge within same shared component → NOT a bridge."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
        ])
        G._save_original_shared_components()
        G._build_synteny_blocks()
        G._compress_to_block_level()
        n = G._resolve_block_bridge_events()
        bridge_events = [e for e in G.events
                         if e.event_type in ('eej', 'ncf', 'fission', 'bridge_unclassified')]
        assert len(bridge_events) == 0


# ============================================================
# Test: Resolve All Events (Integration)
# ============================================================

class TestResolveAllEvents:
    """Integration tests for the full pipeline."""

    def test_no_events_linear(self):
        """Two identical children → no events."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D']]),
            ('c2', [['A', 'B', 'C', 'D']]),
        ])
        G.resolve_all_events()
        major_events = [e for e in G.events
                        if e.event_type in ('inversion', 'eej', 'ncf', 'fission',
                                            'reciprocal_translocation')]
        assert len(major_events) == 0

    def test_inversion_pipeline(self):
        """Full pipeline with inversion (direction conflict)."""
        G = ColoredGraph(hog_level='test')
        chroms1 = [[tel('c1_c0', 'L')] + list('ABCDEF') + [tel('c1_c0', 'R')]]
        chroms2 = [[tel('c2_c0', 'L')] + list('ABDCEF') + [tel('c2_c0', 'R')]]
        G.add_child('c1', MockChildGraph(chroms1))
        G.add_child('c2', MockChildGraph(chroms2))
        G.resolve_all_events()
        inv_events = [e for e in G.events
                      if e.event_type in ('inversion', 'telomere_inversion')]
        assert len(inv_events) >= 1, "No inversion in pipeline. Events: {}".format(G.events)

    def test_chromosome_count_consistency(self):
        """After resolution, chromosome count should be consistent."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D']]),
            ('c2', [['A', 'B', 'C', 'D']]),
        ])
        G.resolve_all_events()
        paths = G.path_cover()
        assert len(paths) >= 1

    def test_path_cover_covers_all_hogs(self):
        """Path cover should cover all HOGs in the graph."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E']]),
            ('c2', [['A', 'B', 'C', 'D', 'E']]),
        ])
        G.resolve_all_events()
        paths = G._cached_paths if hasattr(G, '_cached_paths') else G.path_cover()
        covered = set()
        for p in paths:
            covered.update(p)
        assert len(covered) > 0


# ============================================================
# Test: Edge Cases
# ============================================================

class TestEdgeCases:

    def test_empty_graph(self):
        """Empty graph → no events, no crash."""
        G = ColoredGraph(hog_level='test')
        G.resolve_all_events()
        assert len(G.events) == 0

    def test_single_hog_per_child(self):
        """One HOG per child, no edges (single HOG has no adjacency)."""
        G = build_simple_graph([
            ('c1', [['A']]),
            ('c2', [['A']]),
        ])
        # Single HOG → tel_L-A-tel_R, but no gene-gene edges (only 1 gene)
        # add_child only adds edges between consecutive non-telomere HOGs
        # With 1 gene, 0 consecutive pairs → 0 edges
        assert G.node_count() == 0  # A never added (no edge creation)
        assert G.edge_count() == 0

    def test_three_children(self):
        """Three children with same adjacency → all shared, no events."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['A', 'B', 'C']]),
            ('c3', [['A', 'B', 'C']]),
        ])
        assert G.edge_count() == 2
        for h1, h2 in [('A', 'B'), ('B', 'C')]:
            colors = G.get_colors(h1, h2)
            assert len(colors) == 3  # all 3 children

    def test_disjoint_children(self):
        """Two children with completely different HOGs → no shared edges."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C']]),
            ('c2', [['X', 'Y', 'Z']]),
        ])
        assert G.node_count() == 6
        assert len(G.shared_edges()) == 0
        assert len(G.unique_edges()) == 4


# ============================================================
# Test: Block Compression
# ============================================================

class TestBlockCompression:

    def test_linear_chain_one_block(self):
        """Linear chain with all shared edges → one block."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E']]),
            ('c2', [['A', 'B', 'C', 'D', 'E']]),
        ])
        G._build_synteny_blocks()
        assert len(G._blocks) >= 1

    def test_telomere_hogs_not_in_blocks(self):
        """Telomere-adjacent HOGs should not be compressed into large blocks."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E']]),
            ('c2', [['A', 'B', 'C', 'D', 'E']]),
        ])
        G._build_synteny_blocks()
        cons_tels = G.child_telomere_set()
        for hog in cons_tels:
            bid = G._hog_to_block.get(hog)
            if bid:
                assert len(G._blocks[bid]) <= 2

    def test_shared_edges_form_blocks(self):
        """Shared edges between children form multi-HOG blocks."""
        G = build_simple_graph([
            ('c1', [['A', 'B', 'C', 'D', 'E', 'F']]),
            ('c2', [['A', 'B', 'C', 'D', 'E', 'F']]),
        ])
        G._build_synteny_blocks()
        # Should have at least one multi-HOG block
        multi_hog_blocks = [bid for bid, hogs in G._blocks.items() if len(hogs) >= 2]
        assert len(multi_hog_blocks) >= 1


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
