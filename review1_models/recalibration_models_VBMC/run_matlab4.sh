#!/bin/bash
#
#SBATCH --job-name=RecalModel4
#SBATCH -a 1-4
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

fit_recal_model(4,1)

EOF
