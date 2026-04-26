"""
Chromosome Path Cover: Telomere-centric ancestral karyotype reconstruction
using OR-Tools CP-SAT solver.

Problem formulation:
- Input: weighted undirected graph G=(V,E), telomere HOGs T, target chromosomes k
- Output: k disjoint paths covering V'⊆V (fragments discarded)
- Constraints: degree≤2, path endpoints prefer T, fragment length≥threshold
- Objective: max edge weight + telomere endpoint reward - penalties

Key design decisions (v1.1):
- Only create CP-SAT variables for actual edges (avoid O(n²) explosion)
- OnlyEnforceIf only accepts BoolVar, not linear expressions
- Objective includes reward for telomere HOGs as endpoints (not just penalty)
- Fragment discard penalty integrated into objective function
"""

import logging
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
import networkx as nx

try:
    from ortools.sat.python import cp_model
    ORTOOLS_AVAILABLE = True
except ImportError:
    ORTOOLS_AVAILABLE = False

from .RunCmdsMP import logger


class ChromosomePathCover:
    """
    Solve constrained path cover for ancestral karyotype reconstruction.
    Uses OR-Tools CP-SAT for exact solving with fallback to greedy.
    """

    def __init__(self,
                 edge_weights: Dict[Tuple[Any, Any], float],
                 support_counts: Dict[Tuple[Any, Any], int],
                 telomere_hogs: Set[Any],
                 target_chromosomes: int,
                 n_children: int = 2,
                 reward_telomere_endpoint: float = 80.0,
                 penalty_non_telomere_endpoint: float = 100.0,
                 penalty_fragment_discard: float = 50.0,
                 min_chromosome_genes: int = 3,
                 cpsat_timeout: int = 60,
                 cpsat_workers: int = 8,
                 telomere_weights: Optional[Dict[Any, float]] = None):
        """
        Initialize path cover solver.

        Parameters:
            edge_weights: Dict[(hog_u, hog_v), weight]
            support_counts: Dict[(hog_u, hog_v), n_children_supporting]
            telomere_hogs: Set of HOGs that are telomere candidates
            target_chromosomes: Target number of chromosomes (paths)
            n_children: Number of child graphs (for support normalization)
            reward_telomere_endpoint: Reward for telomere HOG as endpoint
            penalty_non_telomere_endpoint: Penalty for non-telomere HOG as endpoint
            penalty_fragment_discard: Penalty per discarded fragment
            min_chromosome_genes: Minimum genes to qualify as chromosome (not fragment)
            cpsat_timeout: CP-SAT timeout in seconds
            cpsat_workers: Number of parallel workers
            telomere_weights: Dict of HOG -> weight multiplier for telomere reward.
                HOGs at multiple chromosome endpoints get higher weight.
                If None, all telomere_hogs get weight 1.0.
        """
        self.edge_weights = edge_weights
        self.support_counts = support_counts
        self.telomere_hogs = telomere_hogs
        self.target_chromosomes = target_chromosomes
        self.n_children = n_children
        self.reward_telomere = reward_telomere_endpoint
        self.penalty_non_telomere = penalty_non_telomere_endpoint
        self.penalty_fragment = penalty_fragment_discard
        self.min_chromosome_genes = min_chromosome_genes
        self.cpsat_timeout = cpsat_timeout
        self.cpsat_workers = cpsat_workers

        # Build node list and index mapping from actual edges + telomere HOGs
        self.nodes = set()
        for u, v in edge_weights:
            self.nodes.add(u)
            self.nodes.add(v)
        # Include telomere HOGs even if they have no edges (isolated endpoints)
        self.nodes |= telomere_hogs
        self.node_list = list(self.nodes)
        self.n_nodes = len(self.node_list)
        self.node_idx = {n: i for i, n in enumerate(self.node_list)}

        # Precompute telomere flags and weights
        self.is_telomere = [n in telomere_hogs for n in self.node_list]
        self.telomere_weight = []
        for n in self.node_list:
            if telomere_weights and n in telomere_weights:
                self.telomere_weight.append(telomere_weights[n])
            elif n in telomere_hogs:
                self.telomere_weight.append(1.0)
            else:
                self.telomere_weight.append(0.0)

        # Build normalized edge index: only for actual edges
        # key = (min_idx, max_idx) to avoid duplicates
        self.edge_list = []  # List of (idx_i, idx_j, weight, support)
        self.edge_var_map = {}  # (min_idx, max_idx) -> position in edge_list
        for (u, v), w in edge_weights.items():
            i, j = self.node_idx[u], self.node_idx[v]
            key = (min(i, j), max(i, j))
            if key not in self.edge_var_map:
                support = support_counts.get((u, v), 1)
                self.edge_var_map[key] = len(self.edge_list)
                self.edge_list.append((key[0], key[1], w, support))

        logger.info("  CPC init: {} nodes, {} edges, target={} chromosomes, "
                    "{} telomere HOGs".format(
                        self.n_nodes, len(self.edge_list),
                        target_chromosomes, len(telomere_hogs)))

    def solve(self) -> Tuple[Optional[List[List[Any]]], Dict[str, Any]]:
        """
        Solve path cover using CP-SAT.

        Returns:
            paths: List of chromosome paths (List[HOGs]), or None if infeasible
            stats: Solver statistics dict
        """
        if not ORTOOLS_AVAILABLE:
            logger.warning("  OR-Tools not available, cannot solve")
            return None, {'status': 'NO_ORTOOLS'}

        if self.n_nodes == 0:
            return [], {'status': 'EMPTY', 'n_paths': 0}

        if self.n_nodes <= self.target_chromosomes:
            # Each node is its own chromosome
            paths = [[n] for n in self.node_list]
            return paths, {'status': 'TRIVIAL', 'n_paths': len(paths)}

        try:
            return self._solve_cpsat()
        except Exception as e:
            logger.warning("  CP-SAT failed: {}".format(e))
            return None, {'status': 'ERROR', 'error': str(e)}

    def _solve_cpsat(self) -> Tuple[List[List[Any]], Dict[str, Any]]:
        """CP-SAT core solver."""
        model = cp_model.CpModel()
        n = self.n_nodes

        # ===== Variables =====
        # x[e] = 1 if edge e is selected (only for actual edges)
        x = []
        for e_idx, (i, j, w, s) in enumerate(self.edge_list):
            x.append(model.NewBoolVar('x_{}_{}'.format(i, j)))

        # y[i] = degree of node i (0, 1, or 2)
        y = [model.NewIntVar(0, 2, 'y_{}'.format(i)) for i in range(n)]

        # z[i] = 1 if node i is an endpoint (degree <= 1)
        z = [model.NewBoolVar('z_{}'.format(i)) for i in range(n)]

        # sn[i] = 1 if node i is selected (degree >= 1)
        sn = [model.NewBoolVar('sn_{}'.format(i)) for i in range(n)]

        # ===== Constraints =====

        # 1. Degree consistency: y[i] = sum of incident selected edges
        for i in range(n):
            incident = []
            for e_idx, (ei, ej, w, s) in enumerate(self.edge_list):
                if ei == i or ej == i:
                    incident.append(x[e_idx])
            if incident:
                model.Add(y[i] == sum(incident))
            else:
                model.Add(y[i] == 0)

        # 2. Endpoint definition: z[i]=1 ↔ y[i]<=1
        #    z=1 → y<=1 ;  z=0 → y>=2
        for i in range(n):
            model.Add(y[i] <= 1).OnlyEnforceIf(z[i])
            model.Add(y[i] >= 2).OnlyEnforceIf(z[i].Not())

        # 3. Selected node: sn[i]=1 ↔ y[i]>=1
        #    sn=1 → y>=1 ;  sn=0 → y==0
        for i in range(n):
            model.Add(y[i] >= 1).OnlyEnforceIf(sn[i])
            model.Add(y[i] == 0).OnlyEnforceIf(sn[i].Not())

        # 4. Path count: sum(sn) - sum(x) = target_chromosomes
        #    For a collection of disjoint paths: #paths = #nodes_in_paths - #edges
        total_edges = sum(x)
        total_selected = sum(sn)
        model.Add(total_selected - total_edges == self.target_chromosomes)

        # ===== Objective =====
        # Scale factor for converting float weights to integers
        SCALE = 10

        # 1. Edge profit
        edge_profit = 0
        for e_idx, (i, j, w, s) in enumerate(self.edge_list):
            edge_profit += int(w * SCALE) * x[e_idx]

        # 2. Telomere endpoint reward + non-telomere endpoint penalty
        #    +reward * weight for z[i]=1 & is_telomere[i] (weight by endpoint count)
        #    -penalty for z[i]=1 & !is_telomere[i]
        endpoint_obj = 0
        for i in range(n):
            if self.is_telomere[i]:
                w = self.telomere_weight[i]
                endpoint_obj += int(w * self.reward_telomere * SCALE) * z[i]
            else:
                endpoint_obj -= int(self.penalty_non_telomere * SCALE) * z[i]

        # 3. Fragment discard penalty: penalize unselected nodes
        fragment_penalty = 0
        for i in range(n):
            fragment_penalty -= int(self.penalty_fragment * SCALE) * (1 - sn[i])

        # Maximize: profit + rewards - penalties
        model.Maximize(edge_profit + endpoint_obj + fragment_penalty)

        # ===== Solve =====
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = self.cpsat_timeout
        solver.parameters.num_search_workers = self.cpsat_workers
        solver.parameters.log_search_progress = False

        logger.info("  CP-SAT solving: {} nodes, {} edge vars, target={} paths, "
                    "timeout={}s".format(
                        n, len(x), self.target_chromosomes, self.cpsat_timeout))

        status = solver.Solve(model)

        # ===== Extract results =====
        if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
            status_name = solver.StatusName(status)
            logger.warning("  CP-SAT infeasible: {}".format(status_name))
            return None, {
                'status': status_name,
                'objective': None,
                'solve_time': solver.WallTime(),
            }

        # Extract selected edges
        selected_edges = []
        for e_idx, (i, j, w, s) in enumerate(self.edge_list):
            if solver.Value(x[e_idx]) > 0.5:
                selected_edges.append((self.node_list[i], self.node_list[j]))

        # Convert to paths
        paths = self._edges_to_paths(selected_edges)

        # Filter fragments
        final_paths, discarded = self._filter_fragments(paths)

        stats = {
            'status': solver.StatusName(status),
            'objective': solver.ObjectiveValue() / SCALE,
            'n_paths_total': len(paths),
            'n_paths_kept': len(final_paths),
            'n_fragments_discarded': len(discarded),
            'n_fragments_genes': sum(len(p) for p in discarded),
            'solve_time': solver.WallTime(),
            'nodes': n,
            'edges_selected': len(selected_edges),
        }

        logger.info("  CP-SAT solved: {}, {} paths, {} fragments discarded, "
                    "time={:.2f}s".format(
                        stats['status'], stats['n_paths_kept'],
                        stats['n_fragments_discarded'], stats['solve_time']))

        return final_paths, stats

    def _edges_to_paths(self, edges: List[Tuple[Any, Any]]) -> List[List[Any]]:
        """Convert edge set to list of paths."""
        if not edges:
            return []

        g = nx.Graph()
        g.add_edges_from(edges)

        paths = []
        for comp_nodes in nx.connected_components(g):
            if len(comp_nodes) == 1:
                paths.append([list(comp_nodes)[0]])
                continue

            sg = g.subgraph(comp_nodes).copy()

            # Break cycles if any (path count constraint suppresses pure cycles,
            # but a cycle + separate path can still satisfy target)
            if sg.number_of_edges() == len(comp_nodes):
                # Find lowest weight edge to break
                min_w = float('inf')
                break_edge = None
                for u, v in sg.edges():
                    key = (u, v)
                    if key not in self.edge_weights:
                        key = (v, u)
                    w = self.edge_weights.get(key, 0)
                    if w < min_w:
                        min_w = w
                        break_edge = (u, v)

                if break_edge:
                    sg.remove_edge(*break_edge)

            # Find endpoints (degree 1)
            endpoints = [n for n in comp_nodes if sg.degree(n) == 1]

            if len(endpoints) == 2:
                start, end = endpoints
                path = self._trace_path(sg, start, end)
                paths.append(path)
            elif len(endpoints) == 0:
                # Cycle that couldn't be broken - extract longest path
                path = self._extract_longest_path(sg)
                paths.append(path)
            else:
                # Branching (shouldn't happen with degree≤2) - extract longest path
                logger.warning("  Unexpected branching: {} endpoints".format(
                    len(endpoints)))
                path = self._extract_longest_path(sg)
                paths.append(path)

        return paths

    def _trace_path(self, g: nx.Graph, start: Any, end: Any) -> List[Any]:
        """Trace path from start to end in graph."""
        path = [start]
        current = start
        prev = None

        while current != end:
            neighbors = [n for n in g.neighbors(current) if n != prev]
            if not neighbors:
                break
            nxt = neighbors[0]
            path.append(nxt)
            prev = current
            current = nxt

        return path

    def _extract_longest_path(self, g: nx.Graph) -> List[Any]:
        """Extract longest simple path from graph (for degenerate cases)."""
        nodes = list(g.nodes())
        if not nodes:
            return []

        # BFS from arbitrary node to find one endpoint
        start = nodes[0]
        lengths = nx.single_source_shortest_path_length(g, start)
        far = max(lengths, key=lengths.get)

        # BFS from far to find other endpoint
        lengths2 = nx.single_source_shortest_path_length(g, far)
        far2 = max(lengths2, key=lengths2.get)

        # Get path
        try:
            path = nx.shortest_path(g, far, far2)
            return path
        except nx.NetworkXNoPath:
            return [far]

    def _filter_fragments(self, paths: List[List[Any]]) -> Tuple[List[List[Any]], List[List[Any]]]:
        """
        Filter out fragment paths (< min_chromosome_genes).

        Returns:
            kept_paths: Paths that qualify as chromosomes
            discarded: Fragment paths
        """
        kept = []
        discarded = []

        for path in paths:
            if len(path) >= self.min_chromosome_genes:
                kept.append(path)
            else:
                discarded.append(path)

        # If too few paths, supplement from discarded (longest first)
        if len(kept) < self.target_chromosomes and discarded:
            discarded.sort(key=len, reverse=True)
            while len(kept) < self.target_chromosomes and discarded:
                kept.append(discarded.pop(0))

        return kept, discarded


# ===== Greedy fallback (for when CP-SAT fails) =====

class GreedyPathCover:
    """
    Greedy path cover as fallback when CP-SAT is unavailable or infeasible.
    Kruskal-style edge selection with telomere preference.
    """

    def __init__(self,
                 edge_weights: Dict[Tuple[Any, Any], float],
                 telomere_hogs: Set[Any],
                 target_chromosomes: int,
                 min_chromosome_genes: int = 3,
                 telomere_weights: Optional[Dict[Any, float]] = None):
        self.edge_weights = edge_weights
        self.telomere_hogs = telomere_hogs
        self.target_chromosomes = target_chromosomes
        self.min_chromosome_genes = min_chromosome_genes
        self.telomere_weights = telomere_weights or {}

        self.nodes = set()
        for u, v in edge_weights:
            self.nodes.add(u)
            self.nodes.add(v)

    def solve(self) -> Tuple[List[List[Any]], Dict[str, Any]]:
        """Greedy edge selection with telomere preference."""
        if not self.nodes:
            return [], {'status': 'EMPTY'}

        # Sort edges by weight descending; telomere-involved edges get a bonus
        TEL_BONUS = 10.0
        def edge_sort_key(item):
            (u, v), w = item
            w_u = self.telomere_weights.get(u, 1.0 if u in self.telomere_hogs else 0.0)
            w_v = self.telomere_weights.get(v, 1.0 if v in self.telomere_hogs else 0.0)
            bonus = TEL_BONUS * max(w_u, w_v)
            return -(w + bonus)

        sorted_edges = sorted(self.edge_weights.items(), key=edge_sort_key)

        # Union-Find
        parent = {n: n for n in self.nodes}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        degree = {n: 0 for n in self.nodes}
        selected = []

        def count_components():
            return len({find(n) for n in self.nodes})

        def can_add(u, v):
            if degree[u] >= 2 or degree[v] >= 2:
                return False
            if find(u) == find(v):
                return False  # Would form cycle
            return True

        initial_components = count_components()

        # Phase 1: Add edges until component count reaches target
        for (u, v), w in sorted_edges:
            if can_add(u, v):
                selected.append((u, v))
                degree[u] += 1
                degree[v] += 1
                union(u, v)

                if count_components() <= self.target_chromosomes:
                    break

        # Convert to paths using a minimal ChromosomePathCover for helper methods
        cpc = ChromosomePathCover({}, {}, self.telomere_hogs,
                                  self.target_chromosomes)
        cpc.edge_weights = self.edge_weights  # Needed for _edges_to_paths
        paths = cpc._edges_to_paths(selected)
        final_paths, discarded = cpc._filter_fragments(paths)

        return final_paths, {
            'status': 'GREEDY',
            'n_paths_kept': len(final_paths),
            'n_fragments_discarded': len(discarded),
        }
