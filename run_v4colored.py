#!/usr/bin/env python3
"""Run v4_colored reconstruction and report results."""
import sys, os, logging
from collections import Counter

BASE = '/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex'
sys.path.insert(0, BASE)
logging.basicConfig(level=logging.WARNING)

sim_dir = sys.argv[1]
out_dir = sys.argv[2]

from soi.AK import AKR
akr = AKR(
    ogfile=f'{sim_dir}/ortholog_groups.txt',
    orthfiles=f'{sim_dir}/ortholog_pairs.txt',
    gfffile=f'{sim_dir}/all_species_gene.gff',
    sptreefile=f'{sim_dir}/species_tree.nwk',
    outpre=f'{out_dir}/AKR',
    reconstruction_algorithm='v4_colored',
    min_genes=0, timeout=600,
)
akr.run()

# Truth karyotypes
truth_chroms = {}
try:
    with open(f'{sim_dir}/ancestors_karyotypes.txt') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                name = parts[0][1:]
                val = parts[1].split()[0] if len(parts) > 1 else '0'
                try:
                    truth_chroms[name] = int(val)
                except ValueError:
                    truth_chroms[name] = val
except FileNotFoundError:
    pass

# Results
for name in sorted(akr.anc_graphs.keys()):
    if not name.startswith('N'):
        continue
    aag = akr.anc_graphs[name]
    chs = sum(1 for _ in aag.chromosomes)
    evs = len(getattr(aag, 'events', []))
    typ = Counter(e.event_type for e in aag.events)
    truth = truth_chroms.get(name, '?')
    print(f'RESULT|{name}|truth={truth}|recon={chs}|events={evs}|types={dict(typ)}')
