"""
Telomere Constraint Manager for TAKR v3.

Determines which HOGs should preferentially be placed at path endpoints
based on their observed endpoint status in child genomes.

Key design:
- Tracks total endpoint observations per HOG (not just source_ids)
- HOGs at multiple chromosome endpoints get higher confidence/weight
- Overlapping endpoints (same HOG at multiple homologous chromosomes)
  are recorded as high-confidence ancestral telomere sites
"""

from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict

from .RunCmdsMP import logger


class TelomereConstraint:
    """
    Manages telomere constraints for path cover reconstruction.

    Uses hog_endpoints (from _map_to_parent_hogs / _merge_child_graphs) to
    determine which HOGs were observed as chromosome endpoints in child genomes.
    These become soft constraints (prefer_endpoint) in the path cover solver.
    """

    def __init__(self, hog_endpoints: Dict[Any, Dict],
                 n_children: int = 2,
                 ploidy_context: Optional[int] = None):
        """
        Parameters:
            hog_endpoints: HOG -> {'left': [(source_id, tel_node), ...],
                                    'right': [(source_id, tel_node), ...]}
            n_children: Number of child graphs
            ploidy_context: WGD ploidy (for pre-WGD reconstruction).
                None for non-WGD nodes.
        """
        self.hog_endpoints = hog_endpoints
        self.n_children = n_children
        self.ploidy_context = ploidy_context

        # Precompute: total endpoint observations per HOG
        # (counts every chromosome endpoint, even if same HOG at multiple)
        self._endpoint_counts = {}  # hog -> total count (left + right)
        for hog, ends in hog_endpoints.items():
            n_left = len(ends.get('left', []))
            n_right = len(ends.get('right', []))
            self._endpoint_counts[hog] = n_left + n_right

    def is_telomere_hog(self, hog: Any) -> bool:
        """Whether this HOG is a chromosome endpoint in ANY observation."""
        return self._endpoint_counts.get(hog, 0) > 0

    def prefer_endpoint(self, hog: Any) -> bool:
        """
        Whether this HOG should preferentially be a path endpoint (soft constraint).

        Rules:
        - Non-WGD context: HOG has any endpoint observation -> prefer
        - WGD context (ploidy_context set): HOG has endpoint observation -> prefer
          (weight handles confidence; overlapping = high confidence)
        """
        return self.is_telomere_hog(hog)

    def endpoint_count(self, hog: Any) -> int:
        """Total number of chromosome endpoint observations for this HOG."""
        return self._endpoint_counts.get(hog, 0)

    def endpoint_weight(self, hog: Any) -> float:
        """
        Weight for telomere endpoint reward in CP-SAT objective.

        Proportional to number of endpoint observations:
        - 1 observation: weight 1.0 (single chromosome endpoint)
        - 2 observations: weight 2.0 (e.g., WGD both homologs share endpoint)
        - etc.

        Overlapping endpoints (same HOG at multiple homologous chromosomes)
        get higher weight as high-confidence ancestral telomere sites.
        """
        count = self._endpoint_counts.get(hog, 0)
        if count == 0:
            return 0.0
        return float(count)

    def get_telomere_hogs(self) -> Set[Any]:
        """All HOGs that should preferentially be endpoints."""
        return {h for h in self.hog_endpoints if self.prefer_endpoint(h)}

    def get_telomere_weights(self) -> Dict[Any, float]:
        """Return dict of HOG -> weight for telomere endpoint reward."""
        return {h: self.endpoint_weight(h) for h in self.get_telomere_hogs()}

    def summary(self) -> str:
        """Return summary string for logging."""
        n_total = len(self.hog_endpoints)
        n_telomere = sum(1 for h in self.hog_endpoints if self.is_telomere_hog(h))
        n_prefer = sum(1 for h in self.hog_endpoints if self.prefer_endpoint(h))
        total_obs = sum(self._endpoint_counts.values())
        return ("TelomereConstraint: {} HOGs, {} any-endpoint, "
                "{} prefer-endpoint, {} total endpoint obs, n_children={}".format(
                    n_total, n_telomere, n_prefer, total_obs, self.n_children))
