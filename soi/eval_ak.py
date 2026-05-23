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
            row = {k.strip().lower(): (v.strip() if v else '') for k, v in row.items()}
            ancestors = ''  # default, populated in unified format

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
                # Unified format: branch, event_type, [genes], [ancestors], [chroms], desc, support
                branch = row.get('branch', '')
                event_type = row.get('event_type', '')
                genes = row.get('genes', '')
                ancestors = row.get('ancestors', '')
                chroms = row.get('chroms', '')
                desc = row.get('desc', '')
                try:
                    support = int(row.get('support', 1))
                except (ValueError, TypeError):
                    support = 1
                # If branch is a node name (no '-'), infer from parent_of
                if parent_of and '-' not in branch and branch in parent_of:
                    branch = _infer_branch(branch, parent_of)

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
                'ancestors': ancestors,
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

    Three modes:
    - 'genes' (default): Greedy matching by event_type + highest gene Jaccard
    - 'type_only': Count-based matching, ignore gene IDs (for cross-namespace)
    - 'ancestors': Greedy matching by event_type + ancestor ID Jaccard
      (works when truth has 'ancestors' column from simulator and detected
       has 'genes' column containing HOG IDs that map to ancestors)

    Args:
        truth_events: List of truth event dicts for one branch
        detected_events: List of detected event dicts for one branch
        min_jaccard: Minimum Jaccard similarity for gene/ancestor set match
        match_mode: 'genes', 'type_only', or 'ancestors'

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
        elif match_mode == 'ancestors':
            # Ancestor ID Jaccard matching
            # Truth events have 'ancestors' field; detected events use 'genes' (HOG IDs)
            # We match by comparing ancestor ID sets
            used_t, used_d = set(), set()
            for di, d in enumerate(d_list):
                best_j, best_ti = 0, -1
                for ti, t in enumerate(t_list):
                    if ti in used_t:
                        continue
                    j = _ancestors_jaccard(t, d)
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


def _ancestors_jaccard(truth_event: dict, detected_event: dict) -> float:
    """Jaccard similarity using ancestor IDs.

    For truth events: use the 'ancestors' field (comma-separated ancestor IDs).
    For detected events: extract ancestor-level identifiers from 'genes' field.
    Detected genes are HOG IDs like 'SOG79.N1.hog0' — the ancestor ID is
    the OG part ('SOG79') which corresponds to the simulator's ancestor group.
    """
    # Truth: use ancestors field directly
    truth_ancs = truth_event.get('ancestors', '')
    if not truth_ancs:
        # Fallback: use genes field
        truth_ancs = truth_event.get('genes', '')
    set_t = set(truth_ancs.split(',')) if truth_ancs else set()
    set_t = {g.strip() for g in set_t if g.strip()}

    # Detected: extract ancestor-level IDs from HOG IDs
    det_genes = detected_event.get('genes', '')
    set_d = set()
    if det_genes:
        for g in det_genes.split(','):
            g = g.strip()
            if not g:
                continue
            # HOG IDs are like 'SOG79.N1.hog0' → ancestor group = 'SOG79'
            # or could be plain ancestor IDs
            if '.' in g:
                og_part = g.split('.')[0]
                set_d.add(og_part)
            else:
                set_d.add(g)

    if not set_t and not set_d:
        return 1.0
    intersection = set_t & set_d
    union = set_t | set_d
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
    """Print evaluation summary as tables, split by event scale."""
    m = global_metrics['micro_F1']
    print()
    print(f"Micro F1: {m['F1']:.3f}  (TP={m['TP']}  FP={m['FP']}  FN={m['FN']})")
    print(f"  Precision: {m['Precision']:.3f}  Recall: {m['Recall']:.3f}")
    print()

    # Group branches by scale
    for scale_name, scale_label in [('large_scale', '[Large-Scale]'),
                                     ('small_scale', '[Small-Scale]'),
                                     ('gene_level', '[Gene-Level]')]:
        rows = []
        for branch, r in sorted(results.items()):
            for etype, bt in sorted(r.get('by_type', {}).items()):
                if _event_scale(etype) == scale_name:
                    total_t = bt['TP'] + bt['FN']
                    total_d = bt['TP'] + bt['FP']
                    rows.append((branch, etype, total_t, total_d,
                                 bt['TP'], bt['FP'], bt['FN'],
                                 bt['Precision'], bt['Recall'], bt['F1']))

        if not rows:
            continue

        print(f"  {scale_label}")
        print(f"  {'Branch':<15} {'Event Type':<30} {'Truth':>5} {'Det':>5} {'TP':>4} {'FP':>4} {'FN':>4} {'Prec':>6} {'Rec':>6} {'F1':>6}")
        print(f"  {'-' * 85}")
        for row in rows:
            branch, etype, t, d, tp, fp, fn, prec, rec, f1 = row
            print(f"  {branch:<15} {etype:<30} {t:>5} {d:>5} {tp:>4} {fp:>4} {fn:>4} {prec:>6.3f} {rec:>6.3f} {f1:>6.3f}")
        print()


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
                # Only count whole genomes, not individual chroms
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

    # Print comparison — only show nodes in tree or lens files
    tree_nodes = set()
    if akr_instance and hasattr(akr_instance, 'tree') and akr_instance.tree:
        for n in akr_instance.tree.traverse():
            tree_nodes.add(n.name)
    elif outpre or lens_dir:
        tree_nodes = set(recon_counts.keys())

    all_nodes = sorted(set(list(truth_counts.keys()) + list(recon_counts.keys())))
    print()
    print("Chromosome Counts")
    print(f"{'Node':<15} {'Truth':<8} {'Recon':<8} {'Match':<6} {'Type':<10}")
    print("-" * 55)
    for node in all_nodes:
        # Skip individual chromosome entries (e.g., "Sp_1|1")
        if '|' in node and node not in tree_nodes:
            continue
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
    """Categorize event as large-scale, small-scale, or gene-level.
    
    Large-scale: EEJ, NCF, fission, RT, URT, WGD (Phase 2 events)
    Small-scale: inversion, unidir_trans (Phase 1 structural events)
    Gene-level:  fractionation, duplication, gain/loss (no structural change)
    """
    from .takr_events import canonicalize_event_type

    et = canonicalize_event_type(event_type)

    LARGE_SCALE = {'eej', 'ncf', 'fission', 'reciprocal_translocation',
                   'unbalanced_reciprocal_translocation', 'WGD'}
    SMALL_SCALE = {'inversion', 'internal_inversion', 'telomere_inversion',
                   'unidir_trans'}

    if et in LARGE_SCALE:
        return 'large_scale'
    if et in SMALL_SCALE:
        return 'small_scale'
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


# ===========================================================================
#  CLI entry point
# ===========================================================================

def main():
    """Command-line entry: run full evaluation and print summary.

    Usage:
        python3 -m soi.eval_ak \\
            --truth events.tsv \\
            --detected AKR.events.tsv \\
            --tree species_tree.nwk \\
            --karyotype ancestors_karyotypes.txt \\
            --lens-dir . \\
            [--report eval_report.tsv]

    Or minimal (no chrom counts):
        python3 -m soi.eval_ak --truth events.tsv --detected AKR.events.tsv --tree species_tree.nwk
    """
    import argparse
    parser = argparse.ArgumentParser(
        description='TAKR v4 Evaluation: compare truth vs detected events')
    parser.add_argument('--truth', required=True, help='Truth events.tsv (simulator output)')
    parser.add_argument('--detected', required=True, help='Detected AKR.events.tsv')
    parser.add_argument('--tree', required=True, help='Species tree .nwk file')
    parser.add_argument('--karyotype', help='ancestors_karyotypes.txt (for chrom counts)')
    parser.add_argument('--lens-dir', help='Directory with AKR.*.lens files')
    parser.add_argument('--report', help='Output path for eval report TSV')
    parser.add_argument('--match-mode', choices=['type_only', 'genes', 'ancestors'],
                        default='type_only',
                        help='Matching mode: type_only (count), genes (Jaccard), ancestors (HOG-level)')
    parser.add_argument('--detailed', action='store_true',
                        help='Print per-event detailed comparison')
    args = parser.parse_args()

    from .takr_events import parse_tree

    tree, parent_of, _ = parse_tree(args.tree)
    truth = load_events(args.truth, parent_of=parent_of)
    detected = load_events(args.detected, parent_of=parent_of, source='detected')

    print("=" * 72)
    print("  TAKR v4 Evaluation")
    print("=" * 72)

    # Chromosome counts
    if args.karyotype:
        compare_chrom_counts(args.karyotype, lens_dir=args.lens_dir)

    # Event evaluation
    common = set(truth.keys()) & set(detected.keys())
    print()
    print(f"Common branches: {sorted(common)}")
    print()

    results, gm = evaluate_branches(
        {b: truth[b] for b in common},
        {b: detected[b] for b in common},
        match_mode=args.match_mode)
    print_summary(results, gm)

    # Detailed comparison (optional)
    if args.detailed:
        print("\n" + "=" * 72)
        print("  Detailed Event Comparison")
        print("=" * 72)
        print_event_comparison(truth, detected)

    # Report file
    if args.report:
        generate_report(results, gm, args.report)


# ===========================================================================
#  Adjacency-level evaluation (gene_ancestor_map based)
# ===========================================================================

def load_gene_ancestor_map(tsv_path: str) -> Dict[str, str]:
    """Load gene_ancestor_map.tsv, return {gene_id: ancestor_id}.

    For WGD copies (e.g. g2.1), maps to the original ancestor (g2).
    For each (gene, species) pair, only the latest entry is kept.
    """
    mapping = {}
    if not os.path.exists(tsv_path):
        return mapping
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('gene', '').strip()
            ancestor = row.get('ancestor', '').strip()
            if gene and ancestor:
                mapping[gene] = ancestor
    return mapping


def load_karyotype_adjacencies(tsv_path: str) -> Dict[str, set]:
    """Load ancestors_karyotypes.txt, return {node: set_of_adjacent_gene_pairs}.

    Each adjacency is a frozenset({gene_a, gene_b}) for two consecutive genes
    on the same chromosome.  Using froset makes order irrelevant.
    """
    adjacencies = {}  # node -> set of frozenset pairs
    if not os.path.exists(tsv_path):
        return adjacencies

    current_node = None
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_node = line[1:].split('\t')[0]
                adjacencies[current_node] = set()
                continue
            if current_node is None:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            genes_str = parts[1]
            genes = genes_str.split()
            for i in range(len(genes) - 1):
                g1 = genes[i].lstrip('-')
                g2 = genes[i + 1].lstrip('-')
                adjacencies[current_node].add(frozenset({g1, g2}))
    return adjacencies


def gene_adjacencies_to_ancestor(
    gene_adj: set,
    gene_to_anc: Dict[str, str],
) -> set:
    """Convert gene-level adjacencies to ancestor-level adjacencies.

    Maps each gene to its ancestor ID, then creates ancestor adjacency pairs.
    Adjacencies where both genes map to the same ancestor (tandem dup) are dropped.
    """
    anc_adj = set()
    for pair in gene_adj:
        g1, g2 = tuple(pair)
        a1 = gene_to_anc.get(g1, g1)
        a2 = gene_to_anc.get(g2, g2)
        if a1 == a2:
            continue  # same ancestor — tandem dup or fractionation artifact
        anc_adj.add(frozenset({a1, a2}))
    return anc_adj


def load_recon_adjacencies_from_gfa(gfa_dir: str, prefix: str = 'AKR') -> Dict[str, set]:
    """Load reconstructed adjacencies from AKR output GFA files.

    Looks for files like AKR.N0.anc.gfa, AKR.N1.anc.gfa etc.
    Returns {node: set of frozenset HOG-ID pairs}.
    """
    import glob
    adjacencies = {}
    pattern = os.path.join(gfa_dir, f'{prefix}.*.anc.gfa')
    for gfa_path in sorted(glob.glob(pattern)):
        fname = os.path.basename(gfa_path)
        # AKR.N0.anc.gfa -> N0
        node = fname.replace(f'{prefix}.', '').replace('.anc.gfa', '')
        adj = set()
        with open(gfa_path) as f:
            for line in f:
                if line.startswith('L'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        h1, h2 = parts[1], parts[3]
                        adj.add(frozenset({h1, h2}))
        adjacencies[node] = adj
    return adjacencies


def evaluate_adjacency(
    truth_karyo_file: str,
    gene_ancestor_map_file: str,
    recon_gfa_dir: str = None,
    recon_adjacencies: Dict[str, set] = None,
    recon_prefix: str = 'AKR',
) -> dict:
    """Evaluate reconstruction by comparing ancestor-level adjacencies.

    Pipeline:
    1. Load truth karyotype (gene-level adjacencies per node)
    2. Map gene adjacencies to ancestor adjacencies using gene_ancestor_map
    3. Load reconstructed adjacencies (HOG-level, from GFA files or dict)
    4. Compare ancestor adjacency sets per node → precision/recall/F1

    Args:
        truth_karyo_file: ancestors_karyotypes.txt from simulator
        gene_ancestor_map_file: gene_ancestor_map.tsv from simulator
        recon_gfa_dir: Directory with AKR.*.anc.gfa files
        recon_adjacencies: Pre-loaded {node: set of frozenset pairs} (HOG-level)
        recon_prefix: File prefix for GFA files (default 'AKR')

    Returns:
        {node: {TP, FP, FN, Precision, Recall, F1}, 'global': {...}}
    """
    # 1. Load truth
    truth_gene_adj = load_karyotype_adjacencies(truth_karyo_file)
    gene_to_anc = load_gene_ancestor_map(gene_ancestor_map_file)

    # 2. Convert truth to ancestor level
    truth_anc_adj = {}
    for node, adj in truth_gene_adj.items():
        truth_anc_adj[node] = gene_adjacencies_to_ancestor(adj, gene_to_anc)

    # 3. Load reconstructed adjacencies
    if recon_adjacencies is None:
        recon_adjacencies = load_recon_adjacencies_from_gfa(
            recon_gfa_dir, prefix=recon_prefix)

    # 4. Compare per node
    results = {}
    total_tp, total_fp, total_fn = 0, 0, 0

    all_nodes = sorted(set(list(truth_anc_adj.keys()) +
                           list(recon_adjacencies.keys())))

    for node in all_nodes:
        t_adj = truth_anc_adj.get(node, set())
        r_adj = recon_adjacencies.get(node, set())

        tp = len(t_adj & r_adj)
        fp = len(r_adj - t_adj)
        fn = len(t_adj - r_adj)

        m = calculate_metrics(tp, fp, fn)
        results[node] = {
            'truth_adj': len(t_adj),
            'recon_adj': len(r_adj),
            **m,
        }
        total_tp += tp
        total_fp += fp
        total_fn += fn

    results['global'] = calculate_metrics(total_tp, total_fp, total_fn)
    return results


def print_adjacency_summary(results: dict):
    """Print adjacency-level evaluation summary."""
    g = results.get('global', {})
    print()
    print("Adjacency-Level Evaluation (ancestor HOG adjacencies)")
    print(f"  Micro F1: {g.get('F1', 0):.3f}  "
          f"(TP={g.get('TP', 0)}  FP={g.get('FP', 0)}  FN={g.get('FN', 0)})")
    print(f"  Precision: {g.get('Precision', 0):.3f}  "
          f"Recall: {g.get('Recall', 0):.3f}")
    print()
    print(f"  {'Node':<15} {'Truth':>6} {'Recon':>6} {'TP':>5} {'FP':>5} "
          f"{'FN':>5} {'Prec':>7} {'Rec':>7} {'F1':>7}")
    print(f"  {'-' * 68}")
    for node in sorted(results.keys()):
        if node == 'global':
            continue
        r = results[node]
        print(f"  {node:<15} {r['truth_adj']:>6} {r['recon_adj']:>6} "
              f"{r['TP']:>5} {r['FP']:>5} {r['FN']:>5} "
              f"{r['Precision']:>7.3f} {r['Recall']:>7.3f} {r['F1']:>7.3f}")
    print()


if __name__ == '__main__':
    main()
