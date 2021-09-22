#!/bin/bash
#SBATCH -A ctbp-levine
#SBATCH --partition=ctbp-levine
#SBATCH --job-name=pf91JOBSEQ
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --export=ALL

module load MATLAB/2015a

matlab -nodisplay -r main


