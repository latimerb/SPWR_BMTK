#!/bin/bash
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A TG-DBS180005
#SBATCH --job-name=full_build
#SBATCH --output=full_build.out
#SBATCH --time 0-02:30


module purge
#module load python
module load intel
module load openmpi_ib
export PYTHONPATH=$HOME/nrn/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$HOME/nrn/x86_64/lib:$LD_LIBRARY_PATH
export PATH=$HOME/nrn/x86_64/bin:$PATH

rm PN_bg_spikes.h5
rm -rf network/*

echo "Building model at $(date)"

python build_network.py

echo "Done building model at $(date)"


