"""
TAKR Validator: Chromosome-event consistency checker for TAKR v3.

Validates that reconstructed karyotype satisfies:
1. Chromosome count vs event count balance (per-child)
2. Telomere constraint satisfaction
3. Path structure integrity
"""

from typing import Dict, List, Set, Tuple, Optional, Any
from collections import Counter, defaultdict

from .RunCmdsMP import logger


class TAKRValidator:
    """Validates reconstruction results for correctness."""

    def validate_chrom_event_balance(self,
                                     anc_chrom_count: int,
                                     events: List,
                                     child_chrom_counts: List[int],
                                     child_source_ids: Optional[List[str]] = None,
                                     ploidy_context: Optional[int] = None,
                                     is_pre_wgd: bool = False,
                                     n_indel_paths_removed: int = 0
                                     ) -> Tuple[bool, int, int]:
        """
        Validate chromosome count vs event count consistency.

        Per-child formula (events on child branch):
            child = anc - fission + ncf + eej
            i.e.  anc = child + fission - ncf - eej

        Rearrangement delta (from anc to child):
            delta = fission - ncf - eej
            (fission: +1, ncf: -1, eej: -1, translocation: 0)

        Parameters:
            anc_chrom_count: Number of chromosomes in reconstructed ancestor
            events: List of RearrangementEvent objects
            child_chrom_counts: List of chromosome counts per child
            child_source_ids: List of child source IDs (for per-child breakdown)
            ploidy_context: WGD ploidy factor (for pre-WGD nodes)
            is_pre_wgd: Whether this is a pre-WGD reconstruction
            n_indel_paths_removed: Number of paths entirely removed by indel
                (each removal reduces chromosome count by 1)

        Returns:
            (is_consistent, expected_delta, actual_delta)
        """
        # Global event counts
        evt_counts = Counter(e.event_type for e in events)
        n_fission = evt_counts.get('fission', 0)
        n_ncf = evt_counts.get('ncf', 0)
        n_eej = evt_counts.get('eej', 0)
        n_translocation = evt_counts.get('translocation', 0)

        # Rearrangement delta: fission +1, ncf -1, eej -1, translocation 0
        rearrangement_delta = n_fission - n_ncf - n_eej

        if is_pre_wgd and ploidy_context and ploidy_context > 1:
            # Pre-WGD: post = anc * ploidy + rearrangement_delta
            if child_chrom_counts:
                avg_child = sum(child_chrom_counts) / len(child_chrom_counts)
                actual_delta = rearrangement_delta
                expected_delta = int(round(avg_child)) - anc_chrom_count * ploidy_context
            else:
                expected_delta = 0
                actual_delta = rearrangement_delta
        else:
            # Normal speciation
            # Account for indel path removal: each removed path = -1 chromosome
            adjusted_anc = anc_chrom_count - n_indel_paths_removed
            actual_delta = rearrangement_delta
            if child_chrom_counts:
                avg_child = sum(child_chrom_counts) / len(child_chrom_counts)
                expected_delta = int(round(avg_child)) - adjusted_anc
            else:
                expected_delta = 0

        # Allow tolerance of ±2 (reconstruction is approximate)
        is_consistent = abs(actual_delta - expected_delta) <= 2

        if not is_consistent:
            logger.warning(
                "  Chrom-event imbalance: fission({}) - ncf({}) - eej({}) + translocation({}) = {}, "
                "expected_delta={}, actual_delta={}, anc_chrom={}, indel_paths_removed={}{}".format(
                    n_fission, n_ncf, n_eej, n_translocation, rearrangement_delta,
                    expected_delta, actual_delta, anc_chrom_count, n_indel_paths_removed,
                    " (pre-WGD x{})".format(ploidy_context)
                    if is_pre_wgd and ploidy_context else ""))

        # Per-child breakdown (if child_source_ids available)
        if child_source_ids and len(child_chrom_counts) > 1:
            per_child_events = defaultdict(Counter)
            for e in events:
                src = getattr(e, 'child_source', None)
                if src and e.event_type in ('fission', 'ncf', 'eej', 'translocation'):
                    per_child_events[src][e.event_type] += 1

            for i, (sid, child_chrom) in enumerate(
                    zip(child_source_ids, child_chrom_counts)):
                ce = per_child_events.get(sid, Counter())
                cf = ce.get('fission', 0)
                cn = ce.get('ncf', 0)
                ce_ej = ce.get('eej', 0)
                ct = ce.get('translocation', 0)
                child_delta = cf - cn - ce_ej
                child_expected = child_chrom - anc_chrom_count
                if abs(child_delta - child_expected) > 2:
                    logger.info(
                        "  Per-child [{}]: fission({}) - ncf({}) - eej({}) "
                        "= {}, expected_delta={} (child_chrom={}, anc_chrom={})".format(
                            sid, cf, cn, ce_ej, child_delta,
                            child_expected, child_chrom, anc_chrom_count))

        return is_consistent, expected_delta, actual_delta

    def validate_telomere_constraint(self,
                                     paths: List[List[Any]],
                                     tc: 'TelomereConstraint'
                                     ) -> List[str]:
        """
        Validate telomere constraint satisfaction.

        Check: for each HOG that prefer_endpoint(hog)=True, verify it
        is at a path endpoint (first or last position in a path).

        Returns:
            List of violation descriptions (empty = all satisfied)
        """
        violations = []
        # Build set of HOGs at path endpoints
        endpoint_hogs = set()
        for path in paths:
            if path:
                endpoint_hogs.add(path[0])
                if len(path) > 1:
                    endpoint_hogs.add(path[-1])

        # Check each prefer-endpoint HOG
        for hog in tc.get_telomere_hogs():
            if hog not in endpoint_hogs:
                # HOG should prefer endpoint but is internal
                # This is a soft constraint violation, not hard
                violations.append(
                    "HOG {} prefer-endpoint but is internal (support={})".format(
                        getattr(hog, 'hog_id', str(hog)),
                        tc.endpoint_support_count(hog)))

        if violations:
            logger.info("Telomere constraint: {} soft violations out of {} "
                        "prefer-endpoint HOGs".format(
                            len(violations), len(tc.get_telomere_hogs())))

        return violations

    def validate_path_structure(self, paths: List[List[Any]]) -> List[str]:
        """
        Validate structural properties of paths.

        Check:
        - No duplicate HOGs across paths
        - No empty paths
        - All paths have degree <= 2 internally (implied by path structure)

        Returns:
            List of violation descriptions (empty = all satisfied)
        """
        violations = []
        seen = set()
        for i, path in enumerate(paths):
            if not path:
                violations.append("Path {} is empty".format(i))
                continue
            for hog in path:
                if hog in seen:
                    violations.append(
                        "HOG {} appears in multiple paths".format(
                            getattr(hog, 'hog_id', str(hog))))
                seen.add(hog)
        return violations
