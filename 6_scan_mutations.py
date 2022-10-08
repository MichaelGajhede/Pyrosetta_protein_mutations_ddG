from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import *


from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols import minimization_packing
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose

import pandas as pd



def relax_it(wt_pose,resi_to_mutate, dump_cutoff, pdb_filename, s0, NEIGHBOUR_RADIUS):

    
    ddG = []
    ddG_resi = []
    
    AA=['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
    
    for i in AA:
        
        mp = packed_pose.to_pose(wt_pose)
        seq = mp.sequence()
        wt_resn = seq[resi_to_mutate - 1]
        mutate_residue(mp,resi_to_mutate,i)
        scorefxnRelax = pyrosetta.create_score_function("ref2015_cart")
        relax_local = pyrosetta.rosetta.protocols.relax.FastRelax()
        relax_local.constrain_relax_to_start_coords(True)
        relax_local.coord_constrain_sidechains(True)
        relax_local.ramp_down_constraints(False)
        relax_local.cartesian(True)
        relax_local.min_type("lbfgs_armijo_nonmonotone")#for non-Cartesian scorefunctions use'dfpmin_armijo_nonmonotone'
        relax_local.set_scorefxn(scorefxnRelax)
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(resi_to_mutate)
        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(resi_sele)
        cc_sele.threshold(NEIGHBOUR_RADIUS)
        cc_vec = cc_sele.apply(mp)
    
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=cc_vec)
        movemap.set_chi(allow_chi=cc_vec)

        relax_local.set_movemap(movemap)
    
    
        tf = rosetta.core.pack.task.TaskFactory()
        tf.push_back(rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())

    
        prevent_repacking_rlt = rosetta.core.pack.task.operation.PreventRepackingRLT()
        #True indicates here that we are flipping the selection. Taken from tutorial!!
        prevent_subset_repacking = rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, cc_sele, True)
        tf.push_back(prevent_subset_repacking)
    
        relax_local.set_task_factory(tf)
        
        relax_local.apply(mp) #relax after minimization
        s1=scorefxnRelax(mp)
        if s1-s0 < dump_cutoff:
            mp.dump_pdb(pdb_filename + '_' + str(resi_to_mutate) + '_' + i + '.pdb' )
        #minimize_Energy(mp) 
        s1=scorefxnRelax(mp)
        ddG_resi = [resi_to_mutate, wt_resn, i, s1]
        ddG.append(ddG_resi)
        
        
    return (ddG)

import argparse
import sys
import os
import multiprocessing
   
ddG_parser = argparse.ArgumentParser(description='Calculate ddG for all possible single mutations in protein with known structure')

ddG_parser.add_argument('--pdb',
                    metavar='PDB filename',
                    type=str,
                    help='pdb file with wt protein structure')

ddG_parser.add_argument('--cores',
                    metavar='Number of cores',
                    type=int,
                    help='Number of cores to be used for multiprocessing',
                    nargs='?',
                    default=multiprocessing.cpu_count())
ddG_parser.add_argument('--energy_dump_cutoff',
                    metavar='ddG kcal/mole',
                    type=float,
                    help='ddG cutoff for dump og mutant pdb in kcal/mole',
                    nargs='?',
                    default=-10.)
ddG_parser.add_argument('--radius',
                    metavar='neighbour radius in Å',
                    type=float,
                    help='Mutation neighbour radius in Å for repacking and energy minimization',
                    nargs='?',
                    default=4.)


args = ddG_parser.parse_args()

input_pdb = args.pdb

if not os.path.isfile(input_pdb):
    print('The file specified does not exist')
    sys.exit()

n_cores = args.cores

dump_cutoff = args.energy_dump_cutoff

print("Number of cores to be used: " + str(n_cores)) 

NEIGHBOUR_RADIUS = args.radius

init("-ignore_unrecognized_res 1 -ex1 -ex2 -flip_HNQ -relax:cartesian -nstruct 20 -crystal_refine -optimization:default_max_cycles 200")

testPose= Pose()
testPose = pose_from_pdb(input_pdb)
#scorefxnDDG=get_fa_scorefxn()

'''
#Optional energy minimization
min_mover = MinMover() #define a Mover in type of MinMover
mm=MoveMap()
mm.set_bb_true_range(28,36)
min_mover.movemap(mm)
min_mover.score_function(scorefxnDDG)
#min_mover.min_type("dfpmin")
min_mover.tolerance(0.01)
print(min_mover)

def minimize_Energy(pose):
    #Minimization#
    min_mover.apply(pose)

    #Trial_mover define#
    kT=1
    mc=MonteCarlo(pose,scorefxnDDG,kT)
    mc.boltzmann(pose)
    mc.recover_low(pose)

    trial_mover = TrialMover(min_mover,mc)
    #Monte Carlo#
    for i in range (100):
        trial_mover.apply(pose)
    
    return

'''


#Firstly relax the structure for later modificiation#
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
s0=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s0)
testPose.dump_pdb(input_pdb[:-4] + '_Fastrelaxed.pdb' )
#ddG = relax_it(testPose, 52)

from multiprocessing import Pool
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import itertools

if __name__ == '__main__':
    with Pool(n_cores) as p:
        work = [ (packed_pose.to_packed(testPose), i, dump_cutoff, input_pdb[:-4], s0, NEIGHBOUR_RADIUS) for i in range(1, len(testPose.residues) + 1) ]
        ddG = p.starmap(relax_it, work)

    ddG_flat = [item for sublist in ddG for item in sublist] #flatten list i.e. remove multiprocessing 20 item lists
    print(ddG_flat)
        
    for i in range(len(ddG_flat)):
        print(ddG_flat[i])
        ddG_flat[i].append(s0)
        ddG_flat[i].append(ddG_flat[i][3]-s0)
        print(ddG_flat[i])
        
    df = pd.DataFrame(ddG_flat, columns = ['resi', 'wt_resn', 'resn', 'dG_mut', 'dG_wt', 'ddG' ] )
    df.to_csv(input_pdb[:-4] + '_mut.csv')
