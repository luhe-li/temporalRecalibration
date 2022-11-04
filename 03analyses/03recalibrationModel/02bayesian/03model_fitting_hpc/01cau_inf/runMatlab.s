#!/bin/bash
#
#SBATCH --job-name=MatlabJobe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2020b
#change this (cd in prince)
#cd $/scratch/ll3981/project1/CI_model_fit

if [ "$SLURM_JOBTMP" == ""]; then
	export SLURM_JOBTMP=/state/partition1/$USER/$$
	mkdir -p $SLURM_JOBTMP
fi

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)


echo
echo "Hostname: $(hostname)"
echo "Job start: $(date)"

cat a1_modelFitting.m | srun matlab -nodisplay

