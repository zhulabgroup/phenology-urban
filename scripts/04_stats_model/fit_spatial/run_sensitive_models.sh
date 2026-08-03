#!/bin/bash
#SBATCH --job-name=cities_MCMC_mat_NY_nooutlier
#SBATCH --account=zhukai1                   # Account name
#SBATCH --output=slurm_job%j_task%a.out
#SBATCH --error=slurm_job%j_task%a.err
#SBATCH --array=1
#SBATCH --time=36:00:00
#SBATCH --partition=standard
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jializ@umich.edu        # Where to send mail
#SBATCH --mail-type=BEGIN,END               # Mail notifications at the start and end of the job

# Usage:
# Submit the script to Slurm with:
#    sbatch run_sensitive_models.sh
# When Slurm job array is completed, create a symbolic link to the desirable output with:
#    ln -si "$WORKING_DIR/intermediate/slurm_job123456/jags_mod.rds" \
#          "$WORKING_DIR/intermediate/jags_mod.rds"
#    The `-si` option in `ln` will prompt for confirmation before replacing any existing link or file.


set -uo pipefail

# Define the R script and input/output files
WORKING_DIR=../phenology-urban/
R_SCRIPT_SEN=../phenology-urban/script/04_stats_model/fit_spatial/inla/run_sensitive_models.R
OUTPUT_DIR=${WORKING_DIR}/temp/slurm_logs

# Only for spBayes
module load Rgeospatial/4.4.0-2024-05-11

# Record the start time
start_time=$(date +%s)

# Run the R script
Rscript $R_SCRIPT_SEN


# Record the end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))

# Convert the duration to a human-readable format
hh=$((duration / 3600))
mm=$(((duration % 3600) / 60))
ss=$((duration % 60))
formatted_duration=$(printf "%02d:%02d:%02d" $hh $mm $ss)

# Convert memory to a human-readable format
mem_gb=$((SLURM_MEM_PER_NODE / 1024))

# Log job information
job_info=$OUTPUT_DIR/slurm_job${SLURM_ARRAY_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.log
echo "Job ID: $SLURM_ARRAY_JOB_ID" > $job_info
echo "Job Name: $SLURM_JOB_NAME" >> $job_info
echo "Partition: $SLURM_JOB_PARTITION" >> $job_info
echo "Array Task ID: $SLURM_ARRAY_TASK_ID" >> $job_info
echo "Number of CPUs: $SLURM_CPUS_PER_TASK" >> $job_info
echo "Memory: ${mem_gb}GB" >> $job_info
echo "Time Limit: $SLURM_TIMELIMIT" >> $job_info
echo "Start Time: $(date -d @$start_time '+%Y-%m-%d %H:%M:%S')" >> $job_info
echo "End Time: $(date -d @$end_time '+%Y-%m-%d %H:%M:%S')" >> $job_info
echo "Duration: $formatted_duration" >> $job_info

# Move Slurm output and error files to the output directory
# Move Slurm output and error files to the output directory
mv slurm_job${SLURM_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.out $R_OUTPUT_DIR/
mv slurm_job${SLURM_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.err $R_OUTPUT_DIR/

