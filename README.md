# Pyrosetta_protein_mutations_ddG

Title says it all. Requires pyrosetta serial installation. Based on snippets from other peoples code, code that I could not get to work on my systems.


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
 
 
 mutation_heatmap.ipynb:
 
 Jupyter lab notebook to produce heatmap of mutations
 
 mutations_3D:
 
 Visualize lowest energy mutations in pymol
