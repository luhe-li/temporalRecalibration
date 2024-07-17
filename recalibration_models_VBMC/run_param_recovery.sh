#!/bin/bash
#
#SBATCH --job-name=RecalModel4
#SBATCH -a 63-100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

param_recovery

EOF
