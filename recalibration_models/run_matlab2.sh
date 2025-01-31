#!/bin/bash
#
#SBATCH --job-name=RecalModel2
#SBATCH -a 1-9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=21
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

fit_recal_model_VBMC(2,1)

EOF
