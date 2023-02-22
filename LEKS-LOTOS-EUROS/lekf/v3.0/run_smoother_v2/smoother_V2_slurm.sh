#!/bin/bash
#SBATCH --job-name=SO2
#SBATCH --partition=thin
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=128Gb
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/home/slr/lekf/v3.0/run_smoother_v2
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load 2020
module load 2021
module load gompi/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load FlexiBLAS/3.0.4-GCC-10.3.0
module load UDUNITS/2.2.28-GCCcore-10.3.0
module use /projects/0/tm5meteo/admin/modulefiles
module load tm5/default
module load Anaconda3/2021.05
module load udunits/2.2.26_shut-up
module load makedepf90/2.8.8

source activate VU
#python3 -u main_smoother_V4.py
python3 -u SO2_Convertion_V1.py

