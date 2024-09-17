#!/bin/bash
#
#SBATCH --job-name=Model1Recovery
#SBATCH -a 1-100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

model_recovery_s3(1)

EOF
