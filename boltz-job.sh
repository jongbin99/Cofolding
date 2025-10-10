#!/bin/bash
#SBATCH --job-name=boltz_aff
#SBATCH --output=logs/boltz_%j.out
#SBATCH --error=logs/boltz_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --partition=gpu.q          
#SBATCH --nodelist=epyc-A40        
#SBATCH --gres=gpu:1

mkdir -p logs ../.boltz_cache

export CUDA_DEVICE_ORDER=PCI_BUS_ID
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

source /nfs/soft/anaconda3/bin/activate
conda activate py3.10

YAML="/nfs/home/jkim/work/projects/sacharin.yaml" #change this to your file
BOLTZ_BIN="/nfs/home/jkim/.conda/envs/py3.10/bin/boltz"   

"$BOLTZ_BIN" predict "$YAML" \
  --use_msa_server \
  --cache /nfs/home/jkim/.boltz_cache
