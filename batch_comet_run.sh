#!/bin/bash

#SBATCH --partition compute
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=24
#SBATCH -A TG-DBS180005
#SBATCH --job-name=full_run
#SBATCH --output=full_run.out
#SBATCH --time 0-03:00

module purge
#module load python
module load intel
module load openmpi_ib
export LD_LIBRARY_PATH=$HOME/nrn/x86_64/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/nrn/lib/python:$PYTHONPATH
export PATH=$HOME/nrn/x86_64/bin:$PATH

rm -rf output


echo "Running model at $(date)"

#mpirun nrniv -mpi -quiet -python3 run_network.py simulation_config.json
ibrun nrniv -mpi -python run_network.py simulation_configECP.json
#python run_network.py simulation_configECP.json

echo "Done running model at $(date)"
