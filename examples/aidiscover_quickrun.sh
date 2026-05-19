#!/usr/bin/env bash
# AI Discovery — end-to-end example targeting EGFR T790M/L858R
#
# Pre-reqs:
#   1. SchrodingerLite installed (`bash install.sh`)
#   2. ANTHROPIC_API_KEY set in your shell:
#        export ANTHROPIC_API_KEY='sk-ant-…'
set -e

ROOT="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
PROJ="$ROOT/projects/egfr_demo"
mkdir -p "$PROJ"
cd "$PROJ"

echo "▶ 1/3  Preparing EGFR T790M/L858R receptor (PDB 4ZAU)…"
prepwizard 4ZAU --extract-ligand auto --keep-hetero -o "$PROJ/receptor"

# auto-detect: find whichever ref_ligand_*.pdb prepwizard just wrote
REF_LIGAND=$(ls "$PROJ/receptor"/ref_ligand_*.pdb 2>/dev/null | head -1)
if [ -z "$REF_LIGAND" ]; then
  echo "ERROR: prepwizard did not produce a reference ligand." >&2
  exit 1
fi
echo "    reference ligand: $REF_LIGAND"

echo "▶ 2/3  Running AI discovery — 3 rounds × 15 molecules…"
aidiscover \
  --target "$ROOT/examples/egfr_aidiscover.txt" \
  --receptor "$PROJ/receptor/receptor.pdbqt" \
  --ref-ligand "$REF_LIGAND" \
  --rounds 3 \
  --per-round 15 \
  --top-k 5 \
  --exhaustiveness 8 \
  -o "$PROJ/aidiscovery"

echo "▶ 3/3  Opening top pose in Maestro viewer…"
# Use find + read to handle paths with spaces correctly
TOP_SDF=""
while IFS= read -r f; do
  TOP_SDF="$f"
  break
done < <(find "$PROJ/aidiscovery/rounds" -name "*_poses.sdf" 2>/dev/null | sort | tail -1)

if [ -n "$TOP_SDF" ]; then
  maestro "$PROJ/receptor/protein_prepared.pdb" "$TOP_SDF" || \
    echo "  (viewer step skipped — results are still saved)"
fi

echo ""
echo "✓ DONE."
echo "  Ranked CSV:    $PROJ/aidiscovery/final_ranked.csv"
echo "  Med-chem note: $PROJ/aidiscovery/summary.txt"
echo "  Per-round:     $PROJ/aidiscovery/rounds/"
