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

init("-ignore_unrecognized_res 1 -ex1 -ex2 -flip_HNQ -relax:cartesian -nstruct 20 -crystal_refine -optimization:default_max_cycles 200")


def relax_it(wt_pose,resi_to_mutate):

    
    ddG = []
    ddG_resi = []
    
    AA=['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
    
    for i in AA:
        
        
        mp = packed_pose.to_pose(wt_pose)
        seq = mp.sequence()
        wt_resn = seq[resi_to_mutate]
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
        cc_sele.threshold(4)
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
        #mp.dump_pdb('1BV1_repack' + '_' + str(resi_to_mutate) + '_' + i + '.pdb' )
        #minimize_Energy(mp) 
        s1=scorefxnRelax(mp)
        ddG_resi = [resi_to_mutate, wt_resn, i, s1]
        ddG.append(ddG_resi)
        
        
    return (ddG)
        


testPose= Pose()
testPose = pose_from_pdb("1BV1_relaxed.pdb")
#scorefxnDDG=get_fa_scorefxn()

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

#relax.apply(testPose)
s0=scorefxnRelax(testPose) #Record the energy score after cartesian_relax
print(s0)

#ddG = relax_it(testPose, 52)

from multiprocessing import Pool
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import itertools

if __name__ == '__main__':
    with Pool() as p:
        work = [ (packed_pose.to_packed(testPose), i) for i in range(1, len(testPose.residues) + 1) ]
        ddG = p.starmap(relax_it, work)

    ddG_flat = [item for sublist in ddG for item in sublist] #flatten list i.e. remove multiprocessing 20 item lists
    print(ddG_flat)
        
    for i in range(len(ddG_flat)):
        print(ddG_flat[i])
        ddG_flat[i].append(s0)
        ddG_flat[i].append(ddG_flat[i][3]-s0)
        print(ddG_flat[i])
        
    df = pd.DataFrame(ddG_flat, columns = ['resi', 'wt_resn', 'resn', 'dG_mut', 'dG_wt', 'ddG' ] )
    df.to_csv('1BV1_mut.csv')
