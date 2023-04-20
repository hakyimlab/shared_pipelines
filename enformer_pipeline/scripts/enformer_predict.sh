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
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools

cd /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python scripts/enformer_predict.py --param_config config_files/run_on_polaris_reference.json