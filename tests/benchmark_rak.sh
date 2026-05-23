#!/bin/bash
# benchmark_rak.sh — TAKR v4 基准测试：模拟 → 重建 → 评估
#
# 用法:
#   bash tests/benchmark_rak.sh [--wgd] [--seed N] [--species N] [--chroms N]
#
# 示例:
#   bash tests/benchmark_rak.sh                    # 无WGD, 默认参数
#   bash tests/benchmark_rak.sh --wgd              # 带WGD
#   bash tests/benchmark_rak.sh --seed 99 --species 6 --chroms 8

set -euo pipefail

# ── 参数解析 ──────────────────────────────────────────────────────────
SEED=42
SPECIES=4
CHROMS=5
WGD=false
OUTDIR="tests/benchmark_run"

while [[ $# -gt 0 ]]; do
    case $1 in
        --wgd)      WGD=true; shift ;;
        --seed)     SEED=$2; shift 2 ;;
        --species)  SPECIES=$2; shift 2 ;;
        --chroms)   CHROMS=$2; shift 2 ;;
        --outdir)   OUTDIR=$2; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--wgd] [--seed N] [--species N] [--chroms N] [--outdir DIR]"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

export PYTHONPATH="$PROJECT_DIR:${PYTHONPATH:-}"

SIM_DIR="$OUTDIR/sim_data"
RECON_DIR="$OUTDIR/recon"
EVAL_DIR="$OUTDIR/eval"

echo "============================================================"
echo "  TAKR v4 Benchmark"
echo "============================================================"
echo "  Species: $SPECIES  Chroms: $CHROMS  Seed: $SEED  WGD: $WGD"
echo "  Output:  $OUTDIR"
echo "============================================================"
echo ""

# ── 清理 ──────────────────────────────────────────────────────────────
rm -rf "$OUTDIR"
mkdir -p "$SIM_DIR" "$RECON_DIR" "$EVAL_DIR"

# ── Step 1: 模拟 ─────────────────────────────────────────────────────
echo "[1/3] Simulating evolution..."

SIM_ARGS=(
    --num-species "$SPECIES"
    --num-chroms "$CHROMS"
    --inv-rate 5.0
    --rt-rate 3.0
    --ncf-rate 1.5
    --eej-rate 2.0
    --fission-rate 0.02
    --seed "$SEED"
    -o "$SIM_DIR"
)

if [ "$WGD" = true ]; then
    # 用树注释强制 WGD
    TREE_FILE="$SIM_DIR/species_tree_input.nwk"
    # 生成带 [p=2] 注释的树 (内部节点 N1 加 WGD)
    python3.11 -c "
from soi.evolution_simulator_ak import parse_tree
import sys
# 先生成默认树
tree, _, _ = parse_tree('$SIM_DIR/species_tree.nwk') if False else (None, None, None)
# 直接写一个带 WGD 的树
with open('$TREE_FILE', 'w') as f:
    n = $SPECIES
    # 简单二叉树，第一个内部节点加 WGD
    leaves = [f'Sp_{i}' for i in range(1, n+1)]
    # (Sp_1,Sp_2)N1[p=2],(Sp_3,Sp_4)N2...
    mid = n // 2
    left = ','.join(leaves[:mid])
    right = ','.join(leaves[mid:])
    f.write(f'(({left})N1[p=2]:0.1,({right})N2:0.1);\n')
" 2>/dev/null || true
    if [ -f "$TREE_FILE" ]; then
        SIM_ARGS+=(-t "$TREE_FILE")
    fi
    SIM_ARGS+=(--wgd-rate 0.0)  # 禁止随机WGD，只用树注释
else
    SIM_ARGS+=(--wgd-rate 0.0)
fi

python3.11 -m soi.evolution_simulator_ak "${SIM_ARGS[@]}" 2>&1 | tail -5
echo ""

# 检查模拟输出
if [ ! -f "$SIM_DIR/events.tsv" ]; then
    echo "ERROR: Simulation failed — no events.tsv generated"
    exit 1
fi

echo "  Truth events: $(tail -n +2 "$SIM_DIR/events.tsv" | wc -l)"
echo "  Truth karyotypes: $(grep '^>' "$SIM_DIR/ancestors_karyotypes.txt" | wc -l) nodes"
echo ""

# ── Step 2: 重建 ─────────────────────────────────────────────────────
echo "[2/3] Reconstructing ancestral karyotypes..."

python3.11 -c "
from soi.AK import AKR
akr = AKR(
    ogfile='$SIM_DIR/ortholog_groups.txt',
    orthfiles=['$SIM_DIR/ortholog_pairs.txt'],
    gfffile='$SIM_DIR/all_species_gene.gff',
    sptreefile='$SIM_DIR/species_tree.nwk',
    outpre='$RECON_DIR/AKR',
    reconstruction_algorithm='v4_colored',
    min_genes=0,
    timeout=600,
)
akr.run()
print(f'Done: {len(akr.events)} events detected')
" 2>&1 | tail -10
echo ""

if [ ! -f "$RECON_DIR/AKR.events.tsv" ]; then
    echo "ERROR: Reconstruction failed — no AKR.events.tsv generated"
    exit 1
fi

echo "  Detected events: $(tail -n +2 "$RECON_DIR/AKR.events.tsv" | wc -l)"
echo ""

# ── Step 3: 评估 ─────────────────────────────────────────────────────
echo "[3/3] Evaluating..."

# 3a. 事件级评估 (type_only 模式)
echo "--- Event-Level (type_only) ---"
soi rakeval \
    --truth "$SIM_DIR/events.tsv" \
    --detected "$RECON_DIR/AKR.events.tsv" \
    --tree "$SIM_DIR/species_tree.nwk" \
    --match-mode type_only \
    --report "$EVAL_DIR/eval_type_only.tsv" \
    2>&1

echo ""

# 3b. 事件级评估 (ancestors 模式)
echo "--- Event-Level (ancestors) ---"
soi rakeval \
    --truth "$SIM_DIR/events.tsv" \
    --detected "$RECON_DIR/AKR.events.tsv" \
    --tree "$SIM_DIR/species_tree.nwk" \
    --match-mode ancestors \
    --report "$EVAL_DIR/eval_ancestors.tsv" \
    2>&1

echo ""

# 3c. 染色体数比较
echo "--- Chromosome Counts ---"
soi rakeval \
    --truth "$SIM_DIR/events.tsv" \
    --detected "$RECON_DIR/AKR.events.tsv" \
    --tree "$SIM_DIR/species_tree.nwk" \
    --karyotype "$SIM_DIR/ancestors_karyotypes.txt" \
    --lens-dir "$RECON_DIR" \
    2>&1 | grep -A 30 "Chromosome Counts"

echo ""

# 3d. 共线性评估 (inter/intra-chromosomal)
echo "--- Synteny Evaluation ---"
soi rakeval \
    --truth "$SIM_DIR/events.tsv" \
    --detected "$RECON_DIR/AKR.events.tsv" \
    --tree "$SIM_DIR/species_tree.nwk" \
    --karyotype "$SIM_DIR/ancestors_karyotypes.txt" \
    --gene-map "$SIM_DIR/gene_ancestor_map.tsv" \
    --gfa-dir "$RECON_DIR" \
    --og-file "$SIM_DIR/ortholog_groups.txt" \
    --synteny-eval \
    2>&1 | grep -A 30 "Synteny"

echo ""
echo "============================================================"
echo "  Benchmark complete. Results in: $OUTDIR/"
echo "============================================================"
