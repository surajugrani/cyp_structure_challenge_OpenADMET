set -uo pipefail
module load schrodinger/2025-3
cd analyses/12_metal_constraint/work/8SPD
LIG=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/analyses/05_full_pdb_redocking/structures/8SPD/8SPD_ligprep.maegz

echo "########## A: USE_CONS fe:3  (Acceptor = built-in feature index 3) ##########"
cat > tA.in <<EOF
GRIDFILE   8SPD_mcons-gridgen.zip
LIGANDFILE   $LIG
PRECISION   SP
POSES_PER_LIG   10
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
USE_CONS   fe:3
NREQUIRED_CONS   1
EOF
$SCHRODINGER/glide -NOJOBID -WAIT tA.in 2>&1 | grep -iE "fatal|error|feature" | head -3
[ -f tA_pv.mae ] && echo "  A: POSES OK" || echo "  A: no poses"

echo "########## B: inline FEATURE section generated from glidefeature.def ##########"
$SCHRODINGER/run python3 - <<'PY'
import os
f=os.environ["SCHRODINGER"]+"/mmshare-v7.1/data/glidefeature.def"
lines=[l.rstrip("\n") for l in open(f)]
# locate the Acceptor block: "3  Acceptor 3  34  0"
i=next(k for k,l in enumerate(lines) if l.split()[1:2]==["Acceptor"])
npat=int(lines[i].split()[3])
pats=[]; j=i+1
for _ in range(npat):
    parts=lines[j].split()
    smarts=parts[0]; natoms=int(parts[1])
    atoms=[lines[j+1+k].strip() for k in range(natoms)]
    pats.append((smarts,",".join(atoms)))
    j+=1+natoms
with open("featblock.txt","w") as o:
    o.write("[FEATURE:1]\n")
    for n,(s,a) in enumerate(pats,1):
        o.write(f"   PATTERN{n}   {s} {a} include\n")
print(f"  generated {len(pats)} Acceptor patterns")
PY
cat > tB.in <<EOF
GRIDFILE   8SPD_mcons-gridgen.zip
LIGANDFILE   $LIG
PRECISION   SP
POSES_PER_LIG   10
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer

[CONSTRAINT_GROUP:1]
   USE_CONS   fe:1
   NREQUIRED_CONS   1
EOF
cat featblock.txt >> tB.in
$SCHRODINGER/glide -NOJOBID -WAIT tB.in 2>&1 | grep -iE "fatal|error|feature" | head -3
[ -f tB_pv.mae ] && echo "  B: POSES OK" || echo "  B: no poses"
