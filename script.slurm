#!/bin/bash -l
#SBATCH --cluster="genius"
#SBATCH --nodes="1"
#SBATCH --job-name="Reference"
#SBATCH --mail-user="alexander.hoogsteyn@kuleuven.be"
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH --time="0:59:00"
#SBATCH --ntasks="1"
#SBATCH --cpus-per-task=16
#SBATCH --account="lp_elect_gen_modeling"
#SBATCH --partition="batch"
#SBATCH --error="%x_%A.e"

cd $VSC_DATA/GreenIndustrialProduction/demos

srun $VSC_DATA/julia-1.8.5/bin/julia --threads 16 rollingHorizonStochastic.jl

