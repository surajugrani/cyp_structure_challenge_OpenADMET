"""Analysis 33 - ICM redock of all 20 ligands into submission 29's receptors.

Submission 29 (0.4638) is the best result of the project: for each ligand, whichever of
the three ESMFold2 runs scored best under Glide `inplace`. Those are the best receptors
we have. ICM has never been docked into an ESMFold2 receptor - analysis 17 used Boltz-2
and AF3 receptors, and analysis 25 ranked those poses with Glide `mininplace`, the mode
submission 29 later identified as the confound behind every failed Glide selection here.

So both halves of the earlier negative result have a defect this run removes: better
receptors, and whatever picks from these poses will not be a `mininplace` gscore.

This step only GENERATES the poses. 10 per ligand, nothing selected and nothing
packaged. How a submission gets built from them is decided after they exist.

Prospective by construction, exactly as analysis 17: the ligand is built in 3D from
SMILES with `obabel --gen3d`, never from the co-folded coordinates, so no pose
information leaks into the docking. The co-folded ligand is used only to define the
binding site for the ICM project.

Settings match analyses 06 and 17 so the numbers stay comparable: confs=10, effort=3.,
seed=42.

Recorded per pose, since the blinded set has no ground truth:
  * RMSD to the co-folded pose - how far ICM moved from what submission 29 ships
  * Fe-ligand distance - analysis 16's type I / type II diagnostic
  * ICM Score and RTCNN

What to expect, written down first. Analysis 17 measured ICM redocking into predicted
receptors at 4 of 15 ligands coordinating the iron against co-folding's 11-12, and
analysis 25 found the oracle pose - the one of ten closest to the co-folded reference -
sat a median 1.69 A away, with neither ICM Score nor Glide gscore able to find it
(2 of 20 within 2 A each). The useful question for these receptors is whether that
oracle gap narrows. If the best available pose is still 1.7 A from the co-folded one,
no selector can rescue this route; if it closes to well under 1 A, selection becomes
worth attempting.

ONE COMPLEX AT A TIME, serial, in a single job - ICM licence limit.
Resumable: anything already in the CSV is skipped on resubmission.
"""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "analyses", "17_icm_redock_predicted", "scripts",
                   "01_icm_redock.py")

spec = importlib.util.spec_from_file_location("a17", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

m.WORK = os.path.join(HERE, "structures")
m.RESULTS = os.path.join(HERE, "results")
m.RECEPTORS = {
    "sub29": os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                          "submission"),
}

if __name__ == "__main__":
    sys.exit(m.main())
