#!/bin/bash
# build_foam_gdml.sh: rebuild the FINAL foam GDML from the original inputs.
#   /data/data/jhyang/SimpleDet/v26.08-JHY/simulation/gaps-assembly_foam.gdml  (~9 min, < 0.5 GB RAM, 1 core)
#
# Chain (each step reads the previous file):
#   1 make_foam_gdml_v2.py  original gaps-assembly.gdml + MEC-301 mesh + 401 CAD slabs -> v2f (simplified foam, 80 physvols)
#   2 fix_strips.py         rotate the 8 Si(Li) strip shapes to the hotfix/reco layout
#   3 fix_rib_layers.py     rib layers -> outer band only
#   4 add_rib_posts.py      + MEC-301 posts (section of v1 mesh minus windows), middle slab = CAD 11.7-83.6 mm,
#                           own 401 rib solid.  Reads input/gaps-assembly_foam_v1_mec301.gdml (v1 = original MEC-301
#                           mesh placed in all 9 gaps; the posts are cut from it) and uses gdml_mesh.py (GDML reader).
# Everything is found relative to this folder; only the two SimpleDet files below (SIM, MESH301) are outside.
# Usage: ./build_foam_gdml.sh <outdir>      (final file = <outdir>/gaps-assembly_foam.gdml; never writes into simulation/)
set -e
[ $# -eq 1 ] || { echo "usage: $0 <outdir>"; exit 1; }
OUT=$(realpath -m "$1"); mkdir -p "$OUT"
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
ulimit -v 4500000
S=$(dirname "$(realpath "$0")")
SIM=/data/data/jhyang/SimpleDet/v26.08-JHY/simulation
MESH301=/data/data/jhyang/SimpleDet/SimpleDet/simulation/resources/parts/gdml/tracker/trk-styrofoam-sec1-lay2.gdml
cd $S
python3 -B -E -s make_foam_gdml_v2.py $SIM/gaps-assembly.gdml $MESH301 $S/input/slabs_mockup.npz $OUT/step1_v2f.gdml $OUT/step1_pieces.txt > $OUT/step1.log 2>&1
python3 -B -I fix_strips.py $OUT/step1_v2f.gdml $OUT/step2_stripfix.gdml > $OUT/step2.log 2>&1
python3 -B -I fix_rib_layers.py $OUT/step2_stripfix.gdml $OUT/step3_ribband.gdml > $OUT/step3.log 2>&1
python3 -B -E -s add_rib_posts.py $OUT/step3_ribband.gdml $OUT/gaps-assembly_foam.gdml > $OUT/step4.log 2>&1
rm -f $OUT/step1_v2f.gdml $OUT/step2_stripfix.gdml $OUT/step3_ribband.gdml
md5sum $OUT/gaps-assembly_foam.gdml $SIM/gaps-assembly_foam.gdml
