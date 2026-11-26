#!/bin/bash
# check_foam_overlaps.sh: the two overlap checks used for the final foam GDML (both 0 overlaps).
#   1 Geant4  foam_ovl <gdml> <key> <npts>: G4PVPlacement::CheckOverlaps on every world daughter whose name contains
#             <key> (default PassiveTrackerFoam = the 80 foam physvols), <npts> random points on each volume's
#             surface (default 20000), tolerance 0.01 mm. Warnings only (no abort, unlike GAPSSim).
#             Thin overlaps can be missed from one side: also run with the neighbours as key, e.g. PassiveModule,
#             PassiveTrackerTopWindow, PassiveTrackerBotWindow, PassiveTrackerModuleFrame.
#   2 ROOT    check_shape_overlap.C: TGeo CheckOverlaps(0.01) in Geant4 units (mm) on ALL volumes, prints total count.
# Rebuild foam_ovl (after a Geant4 change):
#   /home/jhyang/gcc/gcc-install/bin/g++ -O2 -o foam_ovl foam_ovl.cc $(geant4-config --cflags) $(geant4-config --libs) -lxerces-c
# Usage: ./check_foam_overlaps.sh <gdml> <outdir> [key] [npts]
set -e
[ $# -ge 2 ] || { echo "usage: $0 <gdml> <outdir> [key] [npts]"; exit 1; }
G=$(realpath "$1"); OUT=$(realpath -m "$2"); KEY=${3:-PassiveTrackerFoam}; NPTS=${4:-20000}
D=$(dirname "$(realpath "$0")")
mkdir -p "$OUT"
export LD_LIBRARY_PATH=/home/jhyang/gcc/gcc-install/lib64:$LD_LIBRARY_PATH
source /home/jhyang/geant4-install_0817_gcc13/bin/geant4.sh
source /home/jhyang/root_install_gcc13/bin/thisroot.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -v 4500000
"$D/foam_ovl" "$G" "$KEY" "$NPTS" > "$OUT/g4_ovl_$KEY.log" 2>&1
echo "G4  : $(tail -1 "$OUT/g4_ovl_$KEY.log")   (Overlap is detected: $(grep -c 'Overlap is detected' "$OUT/g4_ovl_$KEY.log" || true))"
cd "$D"; root -l -b -q "check_shape_overlap.C(\"$G\")" > "$OUT/root_ovl.log" 2>&1
echo "ROOT: $(grep '^overlaps' "$OUT/root_ovl.log")"
