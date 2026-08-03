#!/bin/bash
#SBATCH --job-name=download_process_doy_DV
#SBATCH --account=[]                   # Account name
#SBATCH --output=slurm_job%j_task%a.out
#SBATCH --error=slurm_job%j_task%a.err
#SBATCH --array=1
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --mem=15G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=[]        # Where to send mail
#SBATCH --mail-type=BEGIN,END               # Mail notifications at the start and end of the job

# Define the R script and input/output files
WORKING_DIR=../phenology-urban/
R_SCRIPT=../phenology-urban/scripts/download_process_PS.R

# Load necessary modules
module load Rgeospatial/4.2.0-2023-04-25

# Record the start time
start_time=$(date +%s)

# Run the R script
Rscript $R_SCRIPT

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


