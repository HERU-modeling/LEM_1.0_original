#### Open MobaXterm: search for application MobaXterm, load existing session

#### load required modules
module load gcc/5.4.0
module load r/3.5.0


#### run program
sbatch job_foreach_NYC.sh
sbatch job_foreach_ATL.sh
sbatch job_foreach_BAL.sh
sbatch job_foreach_MIA.sh
sbatch job_foreach_LA.sh

sbatch job_foreach_NYC_new.sh
sbatch job_foreach_NYC_new2.sh


sbatch job_foreach_NYC_new6.sh

sbatch job_foreach_BAL_new6.sh
sbatch job_foreach_BAL_new6_2.sh
sbatch job_foreach_BAL_new6_3.sh

sbatch job_foreach_MIA_new6.sh
sbatch job_foreach_MIA_new6_2.sh
sbatch job_foreach_MIA_new6_3.sh

sbatch job_foreach_SEA_new6.sh
sbatch job_foreach_SEA_new6_2.sh

sbatch job_foreach_ATL_new6.sh
sbatch job_foreach_ATL_new6_2.sh

sbatch job_foreach_LA_new6.sh

#### see status of program
squeue -u xiaozang


#### see resouce use of program (replace JOBID with actual job ID)
seff JOBID
seff 

#### cancel program
scancel JOBID
scancel 14835661

scancel 14836514


scancel -u xiaozang
