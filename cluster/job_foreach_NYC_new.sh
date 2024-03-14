#!/bin/bash
#SBATCH --account=def-bnosyk   # replace this with your own account
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=4000M      # memory; default unit is megabytes
#SBATCH --time=1-00:00           # time (DD-HH:MM)
#SBATCH --mail-user=xza162@sfu.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)

module load r/3.5.0
export NODESLIST=$(echo $(srun hostname))
R -f ODEsolve_calib_cluster_NYC_new.R