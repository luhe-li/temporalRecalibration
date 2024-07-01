#!/bin/bash
#
#SBATCH --job-name=FitAtheoreticalModels
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=21
#SBATCH --mem=64GB
#SBATCH --time=15:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

fit_atheo_model_VBMC(1,1)

EOF
