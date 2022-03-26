#!/bin/sh

#SBATCH --job-name=drgs
#SBATCH --output=drgs.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --partition=thin-shared
#SBATCH --time=23:59:00		# 
#SBATCH --mem=1MB
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=marcelinofuentes@gmail.com 
#SBATCH --array=33-37

srun $HOME/code/gnr/bin/gnr ${SLURM_ARRAY_TASK_ID}
