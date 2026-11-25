#!/bin/bash
#SBATCH --job-name=download_daymet_prcp
#SBATCH --account=zhukai0                   # Account name
#SBATCH --output=slurm_job%A_task%a.out
#SBATCH --error=slurm_job%A_task%a.err
#SBATCH --array=0-4
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=jializ@umich.edu        # Where to send mail
#SBATCH --mail-type=BEGIN,END               # Mail notifications at the start and end of the job

set -eo pipefail
export LC_ALL=C

TREE_DIR="/home/jializ/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input"
SCRIPT="$SLURM_SUBMIT_DIR/script/process_Daymet_point.py"

shopt -s nullglob
files=("$TREE_DIR"/*.rds)
if (( ${#files[@]} == 0 )); then
	echo "No city RDS files found in $TREE_DIR" >&2
	exit 1
fi
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"

cities=()
for ((index = SLURM_ARRAY_TASK_ID; index < ${#files[@]}; index += 5)); do
	city=$(basename "${files[index]}" .rds)
	cities+=("$city")
done

if (( ${#cities[@]} == 0 )); then
	echo "No cities assigned to task $SLURM_ARRAY_TASK_ID"
	exit 0
fi

module load python3.9-anaconda/2021.11
source "$SLURM_SUBMIT_DIR/../LilyPythonTest/bin/activate"

echo "Processing ${#cities[@]} cities: ${cities[*]}"
python -u "$SCRIPT" --tree-dir "$TREE_DIR" --cities "${cities[@]}" --workers 2


