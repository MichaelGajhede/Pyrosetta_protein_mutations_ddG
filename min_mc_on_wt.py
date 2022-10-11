from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import *
import argparse
import sys
import os

input_pdb = "1bv1.pdb"

def minimize_Energy(pose):
    #Minimization#
    min_mover.apply(pose)

    #Trial_mover define#
    kT=1
    mc=MonteCarlo(pose,scorefxnRelax,kT)
    mc.boltzmann(pose)
    mc.recover_low(pose)

    trial_mover = TrialMover(min_mover,mc)
    #Monte Carlo#
    for i in range (100):
        trial_mover.apply(pose)
    
    return


ddG_parser = argparse.ArgumentParser(description='Minimize wt structure before mutations')

ddG_parser.add_argument('--pdb',
                    metavar='PDB filename',
                    type=str,
                    help='pdb file with wt protein structure')

args = ddG_parser.parse_args()

input_pdb = args.pdb

init("-ignore_unrecognized_res 1 -ex1 -ex2 -flip_HNQ -relax:cartesian -nstruct 20 -crystal_refine -optimization:default_max_cycles 4000")

testPose= Pose()
testPose = pose_from_pdb(input_pdb)

scorefxnRelax = pyrosetta.create_score_function("ref2015_cart")

# Energy minimization
min_mover = MinMover() #define a Mover in type of MinMover
mm=MoveMap()
mm.set_bb(True)
mm.set_chi(True)
min_mover.movemap(mm)
min_mover.score_function(scorefxnRelax)
min_mover.min_type("dfpmin")
min_mover.tolerance(0.01)
print(min_mover)
s0=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s0)

minimize_Energy(testPose)
s1=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s1)

# Relax the structure for later modificiation#
from pyrosetta.rosetta.protocols.relax import FastRelax
scorefxnRelax = pyrosetta.create_score_function("ref2015_cart")
relax = pyrosetta.rosetta.protocols.relax.FastRelax()
relax.constrain_relax_to_start_coords(True)
relax.coord_constrain_sidechains(True)
relax.ramp_down_constraints(False)
relax.cartesian(True)
relax.min_type("lbfgs_armijo_nonmonotone")#for non-Cartesian scorefunctions use'dfpmin_armijo_nonmonotone'
relax.set_scorefxn(scorefxnRelax)

relax.apply(testPose)
s2=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s2)

minimize_Energy(testPose)
s3=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s3)
testPose.dump_pdb(input_pdb[:-4] + '_min_mc_Fastrelaxed.pdb' )
