# Pyrosetta_protein_mutations_ddG

Title says it all. Requires pyrosetta serial installation. Based on snippets from other peoples code, code that I could not get to work on my systems.

Workflow:

1. Minimize wt structure with script with script min_mc_on_wt.py

python min_mc_on_wt.py -h

usage: min_mc_on_wt.py [-h] [--pdb PDB filename]

Minimize wt structure before mutations

options:
  -h, --help          show this help message and exit
  --pdb PDB filename  pdb file with wt protein structure


2. Calculate ddG for all single mutations

python Pyrosetta_protein_mutations_ddG % python 6_scan_mutations.py -h

usage: 6_scan_mutations.py [-h] [--pdb PDB filename] [--cores [Number of cores]] [--energy_dump_cutoff [ddG kcal/mole]] [--radius [neighbour radius in Å]]

Calculate ddG for all possible single mutations in protein with known structure

options:
  -h, --help            show this help message and exit
  --pdb PDB filename    pdb file with wt protein structure
  --cores [Number of cores]
                        Number of cores to be used for multiprocessing
  --energy_dump_cutoff [ddG kcal/mole]
                        ddG cutoff for dump og mutant pdb in kcal/mole
  --radius [neighbour radius in Å]
                        Mutation neighbour radius in Å for repacking and energy minimization

Use of slurm script ddG_slurm.sh recommended for larger structures 
 
 3. Generate heatmap with mutation_heatmap.ipynb:
 
 Jupyter lab notebook to produce heatmap of mutations
 
 4. Visualize mutations on 3D structure with mutations_3D.pml:
 
 Visualize lowest energy mutations in pymol
