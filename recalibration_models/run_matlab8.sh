#!/bin/bash
#
#SBATCH --job-name=Recal8
#SBATCH -a 1-9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=21
#SBATCH --mem=32GB
#SBATCH --time=30:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

fit_recal_model_VBMC(8,1)

EOF
