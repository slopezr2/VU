#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test
#SBATCH -o result_%x_%j.out      # File to which STDOUT will be written
#SBATCH -e result_%x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.lopezrestrepo@vu.nl

 
module load 2019 
module load eb 
module load intel/2018b 
module load netCDF-Fortran/4.4.4-intel-2018b 
module load Anaconda3/2018.12 
module use /home/slr/admin/modulefiles 
module load tm5 
module load makedepf90


cd /scratch/shared/slr/projects/PIPP/Test_V4/run/

./lotos-euros.x lotos-euros.rc
