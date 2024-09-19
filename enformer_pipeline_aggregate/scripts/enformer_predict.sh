#!/bin/bash

#PBS -N enformer_predict
#PBS -l nodes=1
#PBS -l walltime=00:30:00
#PBS -A TFXcan
#PBS -q preemptable
#PBS -l filesystems=home:grand
#PBS -o /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline/enformer_predict.out
#PBS -e /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline/enformer_predict.err

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/enformer-tools

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python /beagle3/haky/users/shared_pipelines/enformer_pipeline_aggregate/scripts/enformer_predict.py --parameters /project/haky/users/temi/projects/bpnet-comparison/data/config/enformer_parameters_bpnet_SOX2_SOX2_ESC.json