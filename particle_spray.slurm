#!/bin/bash
#SBATCH --job-name=spray     # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1     # number of tasks per node
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G         # memory per cpu-core (4G is default)
#SBATCH --time=03:00:00          # total run time limit (HH:MM:SS)
#SBATCH -p conroy,itc_cluster,shared,sapphire       #partition to submit to
#SBATCH --output=/n/home02/amphillips/p27_nbody/data/particle_sprays/logs/out.%a  # output and error log
#SBATCH --array=0-4 # one sub-job per orbit
#SBATCH --constraint="intel"

# usage: sbatch particle_spray.slurm
# (data/particle_sprays/logs/ has to exist before submitting, or slurm can't write the logs)
# existing spray files are skipped; add --overwrite to the python line to regenerate them.

module purge
module load python/3.12.5-fasrc01
module load gsl/2.8-fasrc01 gcc/14.2.0-fasrc01

mamba activate petar_env

orbits=(gd1 aau pa5 jet c19) #<-- same order as motivation.py
orbit=${orbits[$SLURM_ARRAY_TASK_ID]}
python -u ~/p27_nbody/particle_spray.py --orbit=$orbit --overwrite --progress
