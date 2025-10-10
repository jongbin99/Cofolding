#!/bin/bash

#SBATCH --job-name=chai_fold_array     # Job name
#SBATCH --output=chai_fold_array_%A_%a.out  # Standard output log
#SBATCH --error=chai_fold_array_%A_%a.err   # Standard error log
##SBATCH --gpus=all                   # Request all GPUs
#SBATCH --cpus-per-task=8            # Number of CPUs (adjust as needed)
#SBATCH --array=236-560%1
#SBATCH --mem=32G                    # Memory allocation (adjust as needed)
#SBATCH --nodelist=epyc-A40          # Run only on the 'epyc-a40' node

#Define directories
FILE_LIST="/mnt/nfs/home/jkim/work/250115_chaifold/fasta_list.txt"
INPUT_DIR="/mnt/nfs/home/jkim/work/250115_chaifold/Mac1_fasta_2"
OUTPUT_DIR="/mnt/nfs/home/jkim/work/250115_chaifold/Chai_output"

#Output directory is created
mkdir -p "$OUTPUT_DIR"

#Python
source /nfs/home/mrachman/anaconda2/bin/activate
conda activate py3.10

# Free GPU memory optimization
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Debug CUDA errors
export CUDA_LAUNCH_BLOCKING=1

# Get the FASTA file corresponding to this array task
FASTA_FILE="$INPUT_DIR/$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")"

# Debugging output
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "FASTA_FILE: $FASTA_FILE"

# Check if the FASTA file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file not found - $FASTA_FILE"
    exit 1
fi

# Define the output file
OUTPUT_FILE="$OUTPUT_DIR/$(basename "${FASTA_FILE%.fasta}.out")"

# Run the folding command
echo "Running chai fold on $FASTA_FILE"
chai fold "$FASTA_FILE" "$OUTPUT_FILE"
