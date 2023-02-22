#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=6-23:59:00
#SBATCH --job-name=op_1
#SBATCH -o result_%x_%j.out      # File to which STDOUT will be written
#SBATCH -e result_%x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ayarceb@eafit.edu.co

module load python/3.6.5_miniconda-4.5.1
module load gcc/5.4.0
module load netcdf-fortran/4.4.3_gcc-5.4.0
module load udunits/2.2.26_gcc-5.4.0

#cd /shared-dirs/grimmat-bec/MAUI/scratch_LE/projects/LE/WRF_5_layers_W1/run
cd /home/slopezr2/scratch/projects/LEKF_OSSE/opcion_1/run

./lekf.x lekf.rc
