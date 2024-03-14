#!/bin/bash
#SBATCH --account=def-bnosyk   # replace this with your own account
#SBATCH --ntasks=200               # number of processes
#SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#SBATCH --time=5-00:00           # time (DD-HH:MM)
#SBATCH --mail-user=xza162@sfu.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)

module load r/3.5.0
export R_LIBS=~/home/xiaozang
R CMD BATCH --no-save --no-restore ODEsolve_calib_cluster_BAL.R