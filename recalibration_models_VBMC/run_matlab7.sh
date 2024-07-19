#!/bin/bash
#
#SBATCH --job-name=RecalModel7
#SBATCH -a 1-3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=41
#SBATCH --mem=64GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

fit_recal_model_VBMC(7,1)

EOF
