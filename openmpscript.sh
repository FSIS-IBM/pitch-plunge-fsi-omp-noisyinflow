#!/bin/bash
#PBS -N "job_1"
#PBS -q workq
#PBS -l select=1:ncpus=24
#PBS -o Output.out
#PBS -e error.err
#PBS -V


cd $PBS_O_WORKDIR

ulimit -s unlimited
export OMP_NUM_THREADS=24
./run_executable.out > output.txt
