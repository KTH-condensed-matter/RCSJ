#!/bin/bash

# Submit as sbatch -t 60 -a 1 -n 1 slurm r xxiv

#SBATCH -A SNIC2020-5-92

#SBATCH -n 1

##SBATCH --time=01:00:00

# For Jobarray:
##SBATCH -a 0-1

# Task id is $SLURM_ARRAY_TASK_ID
# Processor core is $SLURM_LOCALID
# Node number is $SLURM_NODEID

dir=$1
shift

dir=$dir-$d$SLURM_ARRAY_TASK_ID

mkdir -p $dir
cp $1 $dir

cd $dir
ln -s ../rcsj . 

srun -n $SLURM_NTASKS -W 4800 $*

cd ..
