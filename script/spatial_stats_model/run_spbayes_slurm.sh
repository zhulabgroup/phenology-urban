#!/bin/bash
#SBATCH --job-name=bayes_sp_urban_phe
#SBATCH --account=zhukai0                   # Account name
#SBATCH --output=slurm_job%j_task%a.out
#SBATCH --error=slurm_job%j_task%a.err
#SBATCH --array=1-4
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=jializ@umich.edu        # Where to send mail
#SBATCH --mail-type=BEGIN,END               # Mail notifications at the start and end of the job

# Usage:
# Submit the script to Slurm with:
#    sbatch run_jags_array_slurm.sh
# When Slurm job array is completed, create a symbolic link to the desirable output with:
#    ln -si "$WORKING_DIR/intermediate/slurm_job123456/jags_mod.rds" \
#          "$WORKING_DIR/intermediate/jags_mod.rds"
#    The `-si` option in `ln` will prompt for confirmation before replacing any existing link or file.


# Define combinations of chains and iterations
PHE_TYPES=("SOS" "EOS")
TAXON_COLS=("canonicalname" "genus")
NUM_PHE=${#PHE_TYPES[@]}
NUM_COL=${#TAXON_COLS[@]}
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
PHE_INDEX=$((INDEX / NUM_COL))
COL_INDEX=$((INDEX % NUM_COL))
PHE_TYPE=${PHE_TYPES[$PHE_INDEX]}
TAXON_COL=${TAXON_COLS[$COL_INDEX]}

echo "Running combination: PHE_TYPE=$PHE_TYPE, TAXON_COL=$TAXON_COL"


# Define the R script and input/output files
R_SCRIPT=../phenology-urban/script/spatial_stats_model/run_spbayes.R
WORKING_DIR=../phenology-urban/
R_INPUT_DIR=../phenology-urban/data/
R_OUTPUT_DIR=../phenology-urban/data/proc/urban/sp_model/fixed_preseason/spbayes/
CITY="NY"

# Export variables to be used in the R script
export R_INPUT_DIR
export R_OUTPUT_DIR
export PHE_TYPE
export CITY
export TAXON_COL

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
job_info=$WORKING_DIR/slurm_job${SLURM_ARRAY_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.log
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
mv slurm_job${SLURM_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.out $R_OUTPUT_DIR/
mv slurm_job${SLURM_JOB_ID}_task${SLURM_ARRAY_TASK_ID}.err $R_OUTPUT_DIR/