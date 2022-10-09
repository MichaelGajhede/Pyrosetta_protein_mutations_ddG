#!/bin/bash
#SBATCH --partition=standard                         # Name of the queue, standard or priority
#SBATCH --nodes=50                                    # Number of nodes
#SBATCH --nodelist=node name                          # If you want the job to run at a specific node
#SBATCH --ntasks-per-node=1                           # Number of tasks per node
#SBATCH --cpus-per-task=1                             # Number of cores per task
#SBATCH --mem=200G                                     # The amount of memory that should be reserved
#SBATCH --mail-type=ALL                               # When should Slurm send you an email (other options: BEGIN, END, FAIL)
#SBATCH --mail-user=mig@sund.ku.dk                    # The email address that Slurm should use
#SBATCH --output=ddG.l.                               # File to redirect the standard output to
#SBATCH --job-name jobname                            # If not set, Slurm will use the name of this script

source /opt/software/miniconda/4.12.0/etc/profile.d/conda.sh 
conda activate pyrosetta

python 6_scan_mutations.py --pdb  AF-Q9UGL1-F1-model_v3.clean.pdb --cores 50