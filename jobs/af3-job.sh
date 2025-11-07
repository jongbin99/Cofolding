#!/bin/bash
#SBATCH --job-name=alphafold3
#SBATCH --output=alphafold3_%j.log
#SBATCH --error=alphafold3_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=8-982%4
#SBATCH --mem=32G
#SBATCH --nodelist=epyc-A40
#SBATCH --partition=gpu.q
#SBATCH --gres=gpu:1

AF_INPUT_DIR="/nfs/home/jkim/work/projects/AmpC_noncofolded/af_input"
AF_OUTPUT_DIR="/nfs/home/jkim/work/projects/AmpC_noncofolded/af_output"
MODEL_PARAMETERS_DIR="/nfs/home/jkim/work/projects/AmpC_noncofolded/af3_weights"
DATABASES_DIR="/local2/af3_public_db"
FILE_LIST="/nfs/home/jkim/work/projects/AmpC_noncofolded/json_list.txt"

TEMP_OUTPUT="/tmp/alphafold_output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# Ensure SLURM_ARRAY_TASK_ID and FILE_LIST exist
if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set."
    exit 1
fi

if [[ ! -f "$FILE_LIST" ]]; then
    echo "Error: File list $FILE_LIST does not exist."
    exit 1
fi

# Get the JSON file corresponding to this array task
json_file="$AF_INPUT_DIR/$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")"

# Ensure the JSON file exists
if [[ ! -f "$json_file" ]]; then
    echo "Error: JSON file $json_file does not exist."
    exit 1
fi

GPU_ID=$((($SLURM_ARRAY_TASK_ID - 1) % 4))
echo "Using GPU $GPU_ID for task $SLURM_ARRAY_TASK_ID"

export CUDA_VISIBLE_DEVICES=$GPU_ID

mkdir -p "$TEMP_OUTPUT"
chmod 777 "$TEMP_OUTPUT"
mkdir -p "$AF_OUTPUT_DIR"

# Run Docker command
echo "=== RUNNING ALPHAFOLD3 on GPU $GPU_ID ==="
newgrp docker <<EOF
docker run \
    --volume "$AF_INPUT_DIR":/root/af_input \
    --volume "$TEMP_OUTPUT":/root/af_output \
    --volume "$MODEL_PARAMETERS_DIR":/root/models \
    --volume "$DATABASES_DIR":/root/public_databases \
    --gpus '"device='$GPU_ID'"' \
    alphafold3 \
    python run_alphafold.py \
    --json_path="/root/af_input/$(basename "$json_file")" \
    --model_dir="/root/models" \
    --output_dir="/root/af_output"
EOF

if [ $? -eq 0 ]; then
    echo "=== COPYING RESULTS TO FINAL LOCATION ==="
    cp -r "$TEMP_OUTPUT"/* "$AF_OUTPUT_DIR"/
    echo "Results copied to $AF_OUTPUT_DIR"
else
    echo "AlphaFold3 run failed on GPU $GPU_ID, check logs above"
fi

rm -rf "$TEMP_OUTPUT"
