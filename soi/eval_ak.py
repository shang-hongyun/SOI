#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eval_ak.py - TAKR event evaluation module.

Evaluate detected rearrangement events against truth events.
Works at branch-level, matching events by type and gene overlap.

Usage:
    from soi.eval_ak import load_events, evaluate_branches, print_summary

    truth = load_events('path/to/truth/events.tsv')
    detected = load_events('path/to/detected/AKR.events.tsv')
    results, global_m = evaluate_branches(truth, detected)
    print_summary(results, global_m)
"""

import csv
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def load_events(tsv_path: str, parent_of: Optional[Dict[str, str]] = None,
                source: str = 'auto') -> Dict[str, List[dict]]:
    """Load events from TSV file, return dict of {branch: [event_dicts]}.

    Supports three formats (auto-detected by header):
    1. Unified:      branch, event_type, [genes], [chroms], desc, support, [child_source]
    2. Old AKR:      node, event_type, [genes], [chroms], desc, support
    3. Old Sim:      node, event_type, details (Python dict string)

    Args:
        tsv_path: Path to events TSV file
        parent_of: {child: parent} mapping for inferring branch from node
        source: 'truth', 'detected', or 'auto'
            'detected' applies AKR_ALIASES to canonicalize type names

    Returns:
        {branch_id: [event_dict, ...]}
    """
    if not os.path.exists(tsv_path):
        print(f"WARNING: events file not found: {tsv_path}")
        return {}

    from .takr_events import canonicalize_event_type, AKR_ALIASES

    events_by_branch = defaultdict(list)

    with open(tsv_path) as f:
        first_line = f.readline().strip()
        dialect = csv.Sniffer().sniff(first_line)
        f.seek(0)
        reader = csv.DictReader(f, dialect=dialect)

        if not reader.fieldnames:
            return {}

        headers = [h.strip().lower() for h in reader.fieldnames]
        if 'branch' in headers:
            fmt = 'unified'
        elif 'node' in headers and 'details' in headers:
            fmt = 'old_sim'
        elif 'node' in headers:
            fmt = 'old_akr'
        else:
            fmt = 'unified'

        for row in reader:
            row = {k.strip().lower(): v.strip() for k, v in row.items()}

            if fmt == 'old_sim':
                # Old simulator format: node, event_type, details (Python dict string)
                node = row.get('node', '')
                event_type = row.get('event_type', '')
                details_str = row.get('details', '{}')
                try:
                    details = eval(details_str) if details_str else {}
                except Exception:
                    details = {}
                branch = _infer_branch(node, parent_of)
                genes = details.get('genes', details.get('gene', ''))
                if isinstance(genes, list):
                    genes = ','.join(str(g) for g in genes)
                chroms = details.get('chroms', details.get('chrom', ''))
                if isinstance(chroms, list):
                    chroms = ','.join(str(x) for x in chroms)
                desc = details.get('desc', details.get('pos', ''))
                support = details.get('support', 1)

            elif fmt == 'old_akr':
                # Old AKR format: node, event_type, [genes], [chroms], desc, support
                node = row.get('node', '')
                event_type = row.get('event_type', '')
                branch = _infer_branch(node, parent_of)
                genes = row.get('genes', '')
                chroms = row.get('chroms', '')
                desc = row.get('desc', '')
                try:
                    support = int(row.get('support', 1))
                except (ValueError, TypeError):
                    support = 1

            else:
                # Unified format: branch, event_type, [genes], [chroms], desc, support
                branch = row.get('branch', '')
                event_type = row.get('event_type', '')
                genes = row.get('genes', '')
                chroms = row.get('chroms', '')
                desc = row.get('desc', '')
                try:
                    support = int(row.get('support', 1))
                except (ValueError, TypeError):
                    support = 1

            # Skip summary events if any slipped through
            if event_type in ('rearrangements',):
                continue

            # Canonicalize event type
            if source == 'detected' and event_type in AKR_ALIASES:
                # For AKR-detected events: apply AKR_ALIASES before canonicalization
                event_type = AKR_ALIASES[event_type]
            event_type = canonicalize_event_type(event_type)

            if not branch:
                continue

            events_by_branch[branch].append({
                'branch': branch,
                'event_type': event_type,
                'genes': genes,
                'chroms': chroms,
                'desc': desc,
                'support': support,
            })

    return dict(events_by_branch)


def _infer_branch(node: str, parent_of: Optional[Dict[str, str]]) -> str:
    """Infer branch identifier from node name using parent_of mapping."""
    if not node:
        return ''
    if parent_of and node in parent_of:
        parent = parent_of[node]
        if parent:
            return f"{parent}-{node}"
    return node


def match_events_branch(
    truth_events: List[dict],
    detected_events: List[dict],
    min_jaccard: float = 0.3,
    match_mode: str = 'genes',
) -> dict:
    """Match detected events to truth events within a single branch.

    Two modes:
    - 'genes' (default): Greedy matching by event_type + highest gene Jaccard
    - 'type_only': Count-based matching, ignore gene IDs (for cross-namespace)

    Args:
        truth_events: List of truth event dicts for one branch
        detected_events: List of detected event dicts for one branch
        min_jaccard: Minimum Jaccard similarity for gene set match
        match_mode: 'genes' or 'type_only'

    Returns:
        dict with keys: tp, fp, fn, matched, false_positives, false_negatives
    """
    truth_by_type = defaultdict(list)
    for e in truth_events:
        truth_by_type[e['event_type']].append(e)
    det_by_type = defaultdict(list)
    for e in detected_events:
        det_by_type[e['event_type']].append(e)

    all_types = set(list(truth_by_type.keys()) + list(det_by_type.keys()))
    matched_pairs = []
    false_positives = []
    false_negatives = []

    for etype in sorted(all_types):
        t_list = list(truth_by_type.get(etype, []))
        d_list = list(det_by_type.get(etype, []))

        if match_mode == 'type_only':
            # Count-based: TP = min(truth, detected)
            tp_count = min(len(t_list), len(d_list))
            for i in range(tp_count):
                matched_pairs.append((t_list[i], d_list[i]))
            for i in range(tp_count, len(t_list)):
                false_negatives.append(t_list[i])
            for i in range(tp_count, len(d_list)):
                false_positives.append(d_list[i])
        else:
            # Gene Jaccard matching
            used_t, used_d = set(), set()

            for di, d in enumerate(d_list):
                best_j, best_ti = 0, -1
                for ti, t in enumerate(t_list):
                    if ti in used_t:
                        continue
                    j = _genes_jaccard(t['genes'], d['genes'])
                    if j > best_j:
                        best_j, best_ti = j, ti
                if best_j >= min_jaccard:
                    used_t.add(best_ti)
                    used_d.add(di)
                    matched_pairs.append((t_list[best_ti], d))

            for ti, t in enumerate(t_list):
                if ti not in used_t:
                    false_negatives.append(t)
            for di, d in enumerate(d_list):
                if di not in used_d:
                    false_positives.append(d)

    return {
        'tp': len(matched_pairs),
        'fp': len(false_positives),
        'fn': len(false_negatives),
        'matched': matched_pairs,
        'false_positives': false_positives,
        'false_negatives': false_negatives,
    }


def _genes_jaccard(genes_a: str, genes_b: str) -> float:
    """Jaccard similarity of two gene/HOG sets."""
    set_a = set(genes_a.split(',')) if genes_a else set()
    set_b = set(genes_b.split(',')) if genes_b else set()
    # Filter empty strings
    set_a = {g for g in set_a if g}
    set_b = {g for g in set_b if g}
    if not set_a and not set_b:
        return 1.0
    intersection = set_a & set_b
    union = set_a | set_b
    return len(intersection) / len(union) if union else 0.0


def calculate_metrics(tp: int, fp: int, fn: int) -> dict:
    """Calculate Precision, Recall, F1 from TP/FP/FN counts."""
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return {
        'TP': tp, 'FP': fp, 'FN': fn,
        'Precision': round(prec, 4),
        'Recall': round(rec, 4),
        'F1': round(f1, 4),
    }


def evaluate_branches(
    truth_by_branch: Dict[str, List[dict]],
    detected_by_branch: Dict[str, List[dict]],
    match_mode: str = 'genes',
) -> Tuple[Dict, Dict]:
    """Evaluate all branches, compute per-branch and global metrics.

    Args:
        truth_by_branch: {branch: [event_dicts]} from load_events
        detected_by_branch: {branch: [event_dicts]} from load_events
        match_mode: 'genes' (Jaccard) or 'type_only' (count-based)

    Returns:
        (results_by_branch, global_metrics)
    """
    all_branches = set(list(truth_by_branch.keys()) +
                       list(detected_by_branch.keys()))

    results = {}
    total_tp, total_fp, total_fn = 0, 0, 0

    for branch in sorted(all_branches):
        truth_e = truth_by_branch.get(branch, [])
        det_e = detected_by_branch.get(branch, [])
        match = match_events_branch(truth_e, det_e, match_mode=match_mode)

        # Per-type breakdown
        by_type = {}
        all_types = set()
        for e in truth_e:
            all_types.add(e['event_type'])
        for e in det_e:
            all_types.add(e['event_type'])

        for etype in sorted(all_types):
            tp_t = sum(1 for t, d in match['matched']
                       if t['event_type'] == etype)
            fp_t = sum(1 for e in match['false_positives']
                       if e['event_type'] == etype)
            fn_t = sum(1 for e in match['false_negatives']
                       if e['event_type'] == etype)
            if tp_t or fp_t or fn_t:
                by_type[etype] = calculate_metrics(tp_t, fp_t, fn_t)

        results[branch] = {
            'total_truth': len(truth_e),
            'total_detected': len(det_e),
            'match': match,
            'metrics': calculate_metrics(match['tp'], match['fp'], match['fn']),
            'by_type': by_type,
        }

        total_tp += match['tp']
        total_fp += match['fp']
        total_fn += match['fn']

    global_micro = calculate_metrics(total_tp, total_fp, total_fn)

    return results, {'micro_F1': global_micro}


def print_summary(results: Dict, global_metrics: Dict):
    """Print evaluation summary to stdout."""
    m = global_metrics['micro_F1']
    print("=== Evaluation Report ===")
    print(f"Micro F1: {m['F1']:.3f}  (TP={m['TP']}  FP={m['FP']}  FN={m['FN']})")
    print(f"  Precision: {m['Precision']:.3f}  Recall: {m['Recall']:.3f}")
    print()

    for branch, r in sorted(results.items()):
        m2 = r['metrics']
        print(f"  {branch}")
        print(f"    Total: truth={r['total_truth']} detected={r['total_detected']}  "
              f"F1={m2['F1']:.3f} (P={m2['Precision']:.3f} R={m2['Recall']:.3f})")
        for etype, bt in sorted(r.get('by_type', {}).items()):
            total_t = bt['TP'] + bt['FN']
            total_d = bt['TP'] + bt['FP']
            print(f"    {etype:<30} truth={total_t:<3} detected={total_d:<3}  "
                  f"F1={bt['F1']:.3f}")

    print()
    print(f"Global Micro F1: {m['F1']:.3f}")
    return m


def generate_report(results: Dict, global_metrics: Dict, outpath: str):
    """Generate evaluation report TSV file."""
    with open(outpath, 'w') as f:
        # Global metrics
        m = global_metrics['micro_F1']
        f.write("# Global Metrics\n")
        f.write(f"# Micro F1\t{m['F1']}\n")
        f.write(f"# Micro Precision\t{m['Precision']}\n")
        f.write(f"# Micro Recall\t{m['Recall']}\n")
        f.write(f"# Total TP\t{m['TP']}\n")
        f.write(f"# Total FP\t{m['FP']}\n")
        f.write(f"# Total FN\t{m['FN']}\n")
        f.write("\n")

        # Per-branch x event_type table
        f.write("branch\tevent_type\ttruth\tdetected\tTP\tFP\tFN\tPrecision\tRecall\tF1\n")
        for branch in sorted(results.keys()):
            r = results[branch]
            for etype in sorted(r.get('by_type', {}).keys()):
                bt = r['by_type'][etype]
                total_t = bt['TP'] + bt['FN']
                total_d = bt['TP'] + bt['FP']
                f.write(f"{branch}\t{etype}\t{total_t}\t{total_d}\t"
                        f"{bt['TP']}\t{bt['FP']}\t{bt['FN']}\t"
                        f"{bt['Precision']:.4f}\t{bt['Recall']:.4f}\t{bt['F1']:.4f}\n")

            # Branch summary
            m2 = r['metrics']
            total_t = r['total_truth']
            total_d = r['total_detected']
            f.write(f"{branch}\t__total__\t{total_t}\t{total_d}\t"
                    f"{m2['TP']}\t{m2['FP']}\t{m2['FN']}\t"
                    f"{m2['Precision']:.4f}\t{m2['Recall']:.4f}\t{m2['F1']:.4f}\n")

    print(f"Report written to {outpath}")


# ===========================================================================
#  Chromosome count comparison
# ===========================================================================

def compare_chrom_counts(truth_karyotype_file, lens_dir=None, akr_instance=None, outpre=None):
    """Compare truth vs reconstructed chromosome counts.

    Args:
        truth_karyotype_file: Path to ancestors_karyotypes.txt (simulator output)
        lens_dir: Directory containing AKR.*.lens files, or
        akr_instance: AKR object (for anc_graphs), or
        outpre: AKR output prefix (e.g. 'tests/AKR')
    """
    # Load truth counts
    truth_counts = {}
    with open(truth_karyotype_file) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                node = parts[0][1:]
                count = int(parts[1].split()[0])
                truth_counts[node] = count

    # Load reconstructed counts
    recon_counts = {}
    if lens_dir:
        import glob, os
        for f in glob.glob(os.path.join(lens_dir, '*.lens')):
            fname = os.path.basename(f)
            node = fname.replace('.lens', '')
            # Remove prefix (e.g., 'AKR.N0' -> 'N0')
            if '.' in node:
                node = node.split('.', 1)[1]
            with open(f) as fh:
                recon_counts[node] = sum(1 for line in fh if line.strip())
    elif akr_instance:
        for node, aag in akr_instance.anc_graphs.items():
            recon_counts[node] = len(list(aag.chromosomes))
    elif outpre:
        import glob, os
        for f in glob.glob(f'{outpre}.*.lens'):
            fname = os.path.basename(f)
            node = fname.replace('.lens', '')
            if '.' in node:
                node = node.split('.', 1)[1]
            with open(f) as fh:
                recon_counts[node] = sum(1 for line in fh if line.strip())

    # Also load GFF for leaf counts if available
    gff_counts = {}
    if lens_dir:
        gff_path = os.path.join(os.path.dirname(lens_dir.rstrip('/')), 'all_species_gene.gff')
        if os.path.exists(gff_path):
            from collections import defaultdict as dd
            chroms = dd(set)
            with open(gff_path) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        chroms[parts[0]].add(parts[4])
            for sp, c in chroms.items():
                if sp not in truth_counts:
                    truth_counts[sp] = len(c)

    # Print comparison
    all_nodes = sorted(set(list(truth_counts.keys()) + list(recon_counts.keys())))
    print()
    print("Chromosome Counts")
    print(f"{'Node':<15} {'Truth':<8} {'Recon':<8} {'Match':<6} {'Type':<10}")
    print("-" * 55)
    for node in all_nodes:
        t = truth_counts.get(node, '?')
        r = recon_counts.get(node, 'N/A')
        match = "✅" if t == r and r != 'N/A' else "❌" if r != 'N/A' else "?"
        ntype = 'leaf' if node.startswith(('Sp', 'Species')) else 'internal'
        print(f"{node:<15} {str(t):<8} {str(r):<8} {match:<6} {ntype:<10}")

    return truth_counts, recon_counts


# ===========================================================================
#  Detailed event comparison
# ===========================================================================

def _event_scale(event_type):
    """Categorize event as large-scale, small-scale, or gene-level."""
    from .takr_events import CHROM_COUNT_EVENTS, CHROM_REARRANGEMENTS
    from .takr_events import canonicalize_event_type

    et = canonicalize_event_type(event_type)

    # Large-scale: changes chromosome count
    if et in CHROM_COUNT_EVENTS:
        return 'large_scale'
    # Large-scale: changes structure but not count
    if et in CHROM_REARRANGEMENTS:
        return 'small_scale'
    # Gene-level
    return 'gene_level'


def print_event_comparison(truth_by_branch, detected_by_branch):
    """Print detailed truth vs detected comparison per branch,
    categorized by event scale (large_scale, small_scale, gene_level),
    with per-type event listings."""
    from collections import defaultdict
    from .takr_events import canonicalize_event_type

    all_branches = sorted(set(list(truth_by_branch.keys()) +
                              list(detected_by_branch.keys())))

    # Categorization titles
    scale_titles = {
        'large_scale': 'Large-Scale (chrom count change)',
        'small_scale': 'Small-Scale (chrom structure change)',
        'gene_level':  'Gene-Level',
    }

    for branch in all_branches:
        print()
        print(f"{'=' * 72}")
        print(f"  Branch: {branch}")
        print(f"{'=' * 72}")

        t_events = truth_by_branch.get(branch, [])
        d_events = detected_by_branch.get(branch, [])

        # Categorize
        by_scale = {'large_scale': {'T': [], 'D': []},
                    'small_scale': {'T': [], 'D': []},
                    'gene_level':  {'T': [], 'D': []}}
        for e in t_events:
            scale = _event_scale(e['event_type'])
            by_scale[scale]['T'].append(e)
        for e in d_events:
            scale = _event_scale(e['event_type'])
            by_scale[scale]['D'].append(e)

        for scale in ('large_scale', 'small_scale', 'gene_level'):
            t_list = by_scale[scale]['T']
            d_list = by_scale[scale]['D']
            if not t_list and not d_list:
                continue

            title = scale_titles[scale]
            print(f"\n  [{title}]")
            print(f"  {'─' * 66}")
            print(f"  Truth ({len(t_list)})    Detected ({len(d_list)})")

            # Group by type
            t_by_type = defaultdict(list)
            for e in t_list:
                t_by_type[e['event_type']].append(e)
            d_by_type = defaultdict(list)
            for e in d_list:
                d_by_type[e['event_type']].append(e)

            all_types = sorted(set(list(t_by_type.keys()) + list(d_by_type.keys())))
            for etype in all_types:
                t_count = len(t_by_type[etype])
                d_count = len(d_by_type[etype])
                match = "✅" if t_count == d_count else (
                        "≈" if 0 < min(t_count, d_count) / max(t_count, d_count) >= 0.5
                        else "❌" if t_count == 0 or d_count == 0 else "⚠")
                print(f"    {etype:<30} T={t_count:<4} D={d_count:<4}  {match}")

                # Show first few event descriptions for each type
                if t_count > 0 and d_count > 0:
                    print(f"      Truth samples:")
                    for e in t_by_type[etype][:3]:
                        desc = e.get('desc', '')
                        print(f"        {desc[:80]}")
                    print(f"      Detected samples:")
                    for e in d_by_type[etype][:3]:
                        desc = e.get('desc', '')
                        print(f"        {desc[:80]}")

        # Summary for this branch
        print()
        total_t, total_d = len(t_events), len(d_events)
        match = "EXACT" if total_t == total_d and total_t > 0 else (
                "≈" if total_t > 0 and total_d > 0 else "FP" if total_d > 0 else "FN" if total_t > 0 else "-")
        print(f"  Branch total: truth={total_t}  detected={total_d}  [{match}]")


__all__ = [
    'load_events', 'match_events_branch', 'evaluate_branches',
    'calculate_metrics', 'print_summary', 'generate_report',
    'compare_chrom_counts', 'print_event_comparison',
]
