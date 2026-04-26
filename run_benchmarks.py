#!/usr/bin/env python3
"""Run benchmark tests with different complexity levels (simplified)."""
import subprocess, sys, os, logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

BASE = '/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex'
os.chdir(BASE)

CONFIGS = [
    ("baseline", "No WGD, moderate rates",
     "--num-species 4 --num-chroms 6 --min-genes 60 --max-genes 180 "
     "--wgd-rate 0.0 --inv-rate 3.0 --rt-rate 1.0 --ncf-rate 0.5 "
     "--eej-rate 0.5 --fission-rate 0.01 --seed 42"),
    
    ("wgd", "With WGD",
     "--num-species 4 --num-chroms 6 --min-genes 60 --max-genes 180 "
     "--wgd-rate 0.5 --inv-rate 3.0 --rt-rate 1.0 --ncf-rate 0.5 "
     "--eej-rate 0.5 --fission-rate 0.01 --seed 42"),
    
    ("complex", "High rates, 6 species",
     "--num-species 6 --num-chroms 8 --min-genes 80 --max-genes 250 "
     "--wgd-rate 0.3 --inv-rate 8.0 --rt-rate 4.0 --ncf-rate 2.0 "
     "--eej-rate 3.0 --fission-rate 0.05 --seed 42"),
]

for tag, desc, sim_args in CONFIGS:
    logger.info("=" * 60)
    logger.info("CONFIG: %s - %s", tag, desc)
    logger.info("=" * 60)
    
    sim_dir = f"tests/bench_{tag}_sim"
    out_dir = f"tests/bench_{tag}_out"
    os.makedirs(sim_dir, exist_ok=True)
    
    # Phase 1: Simulate
    logger.info("[Phase 1] Simulating...")
    cmd = f"python3.11 -m soi.evolution_simulator_ak {sim_args} -o {sim_dir}/"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    # Print summary
    for line in result.stdout.split('\n'):
        if 'chromosomes' in line or 'Events' in line or 'WGD' in line or 'Summary' in line or line.strip().startswith(('  ', 'Ancestral')):
            logger.info("  %s", line.strip())
    
    # Phase 2: Reconstruct
    logger.info("[Phase 2] Reconstructing (v4_colored)...")
    os.makedirs(out_dir, exist_ok=True)
    
    py_code = f"""
import sys, logging
logging.basicConfig(level=logging.WARNING)
sys.path.insert(0, '{BASE}')
from soi.AK import AKR
akr = AKR(
    ogfile='{sim_dir}/ortholog_groups.txt',
    orthfiles='{sim_dir}/ortholog_pairs.txt',
    gfffile='{sim_dir}/all_species_gene.gff',
    sptreefile='{sim_dir}/species_tree.nwk',
    outpre='{out_dir}/AKR',
    reconstruction_algorithm='v4_colored',
    min_genes=0, timeout=600,
)
akr.run()
for name in sorted(akr.anc_graphs.keys()):
    if not name.startswith('N'): continue
    aag = akr.anc_graphs[name]
    chs = 0
    for _ in aag.chromosomes:
        chs += 1
    evs = len(getattr(aag, 'events', []))
    from collections import Counter
    typ = Counter(e.event_type for e in aag.events)
    print(f'NODE|{{name}}|chroms={{chs}}|events={{evs}}|types={{dict(typ)}}')
"""
    r = subprocess.run(f'python3.11 -c "{py_code}"', shell=True, capture_output=True, text=True, timeout=300, cwd=BASE)
    for line in r.stdout.split('\n'):
        if line.startswith('NODE|'):
            parts = line.strip().split('|')
            logger.info("  %s", " | ".join(parts[1:]))
        elif line.strip():
            logger.debug(line.strip())
    
    # Phase 3: Check results
    logger.info("[Phase 3] Summary")
    # Print karyotype comparison if available
    if os.path.exists(f'{sim_dir}/ancestors_karyotypes.txt'):
        truth_chroms = {}
        with open(f'{sim_dir}/ancestors_karyotypes.txt') as f:
            for line in f:
                if line.startswith('>'):
                    parts = line.strip().split('\t')
                    name = parts[0][1:]
                    chs = int(parts[1].split()[0]) if len(parts) > 1 else 0
                    truth_chroms[name] = chs
        logger.info("  Truth karyotypes: %s", truth_chroms)
    
    logger.info("")
