#!/bin/bash

#PBS -P bp72
#PBS -q normal
#PBS -j oe
#PBS -l wd
#PBS -l ncpus=25
#PBS -l mem=20GB
#PBS -l walltime=12:00:00
#PBS -l storage=scratch/bp72

export R_LIBS_USER=~/R_lib
export TMPDIR=~/tmp

module load nci-parallel/1.0.0a
module load intel-compiler/2021.10.0
module load openmpi/4.1.0
module load gcc/12.2.0
module load R/4.4.2

Rscript soft_simulation.R
