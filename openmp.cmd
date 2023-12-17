#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=16
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .

ulimit -s unlimited
export OMP_NUM_THREADS=16
./run_executable.out > output.txt
mv ../job$tpdir $PBS_O_WORKDIR/. 
