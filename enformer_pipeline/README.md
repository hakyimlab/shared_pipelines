
# Shared ENFORMER repo
Pipeline used to run ENFORMER 

### Authors: Sai and Temi
### Date: Sometime in February 2023

## Usage

### 1. Clone the repo
clone the repo using:
    `git clone https://github.com/hakyimlab/shared_folder.git`. 

### 2. Install the software
- If you don't have conda installed, please install conda. [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

### 3. Create conda environment
- If you are on polaris, you can use the shared software environment: `/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools`

- If you are not on polaris, you can create your own environment using the [`enformer-predict-tools.yaml`](./software/enformer-predict-tools.yaml) file.

    `conda env create -p <path to where you want to put the environment> (e.g. ./enformer-predict-tools) -f [](./software/enformer_pipeline/software/enformer-predict-tools.yaml)`

This environment contains all the software needed to run the pipeline; and you should remember to use it in the configuration file defined below.

### 4. Activate conda environment
`conda activate <<path to the environment>> (e.g. ./enformer-predict-tools)`

### 5. Edit the config.json file
Instructions are below.

### 6. Run the pipeline
There are two ways to do this (for now, #1 is recommended):

1. After editing your `config.json` file , simply call: [`python3 ./scripts/enformer_predict.py`](./scripts/enformer_predict.py) --parameters `{path to config.json file}`. **Make sure that the `provider` option is set to 'highthroughput" in the config file.**

2. Edit the [`enformer_predict.sh`](./scripts/enformer_predict.sh) file to include the path to your config.json file and then call: [`qsub ./enformer_predict.sh`](./scripts/enformer_predict.sh). **Make sure that the `provider` option is set to 'local" in the config file.** 

#### Difference between "highthroughput" and "local" providers

The "highthroughput" provider is used when you want to run the pipeline on a cluster from the login nodes. To make sure your job does not stop when logged out, you can use screen or tmux. The "local" provider is used when you want to run the pipeline as if on your local machine. Here, you have direct access to some gpus, or have requested an interactive session. So, when submitting a pbs script or job, you need to use the "local" provider.

An example of a config.json file is [here](./config_files). You should choose one depending on if you want to predict on the reference genome or on personalized genomes, and on what cluster you are on. Instructions for the config.json file are below. Template files for individuals and regions are [here](./metadata/). 


## Options
### General (and fatal if not provided)
- `project_dir`: (relative or absolute path) Provide a directory where files, logs, and outputs will be created. 
- `sub_dir`: (Bool: true or false) To make the folders clean, predictions will be saved in `{project_dir}/predictions_folder/`
- `model_path`: (relative or absolute path) Path to ENFORMER model.
- `fasta_file`: (Absolute path) Path to the fasta file containing the reference sequences.
- `interval_list_file`: (Relative path to this tab-delimited file) Path to the file where the intervals files are saved.
- `output_dir`: (str) name of the folder where predictions will be saved i.e `{project_dir}/predictions_folder/{prediction_data_name}_{prediction_id}/predictions_{date}/{output_dir}`
- `sequence_source`: (str: "reference" or "personalized", "random") A name for the type of dataset. This is necessary so that the pipeline knows whether to use vcf files or not. 
- `reverse_complement`: (bool) Do you want to make predictions for the reverse complement of the extracted sequence too. Currently, this only works when doing personalized predictions. THIS DOES NOT WORK FOR REFERENCE OR RANDOM ONLY PREDICTIONS YET.
- `prediction_data_name`: (str) A unique id for the predictions that will be used to create folders. Can also use the name of the dataset e.g. "freedman", or "kawakami". If `dataset_type` is "reference", this will be the name of the folder within which predictions will be made. If `dataset_type` is personalized, the unique ids of the individuals will be used.
- `prediction_id`: (str) A unique id
- `predictions_log_dir`: (path to folder): Where should predictions be logged? In the event of the job not completing and you intend to re-run, the file `{project_dir}/predictions_folder/{prediction_data_name}_{prediction_id}/predictions_{date}/{predictions_log_dir}/{individual or id}_log.csv` will be read and used such that if a prediction for a region exists and the `{region}_predictions.h5` file exists, predictions will not be made for that region anylonger.
- `batch_regions`: (int) Predictions for intervals or regions will be split into `batch_regions`.  
- `n_regions`: (int) how many regions should be predicted for at a time. If using parsl, each gpu/parsl app will get at most `n_regions` at a time.
- `tracks_to_save`: indices of ENFORMER tracks to be saved; -1 for all.
- `bins_to_save`: indices of ENFORMER bins to be saved; -1 for all


### Personalized predictions (optional; used only when predicting on individuals)
- `individuals`:(path, list or str) The unique ids of the individuals whose predictions are to be made. If providing a file, the ids should be written row-wise. A list of ids can be supplied and a single id can be supplied as a string.
- `n_individuals`: (int) how many individuals to predict on; -1 for all.
- `batch_individuals`: (int) how many individuals to predict on at a time; -1 for all. Individuals will be split into `batch_individuals`. 
- `vcf_files`: (path) A nested json of the path to the phased vcf file that contains the genotypes/variants of the individuals.
    - `folder`: (path) the directory where the vcf files are located
    - `files`: a nested json of `{chr_n: vcf_file_n}`

### Optional
- `exclude_regions`: (path to file or null): During predictions if there are issues with some regions, those regions are logged here, and will be excluded if the job is re-run. 
- `date`: (str in the form "YYYY-MM-DD") The date these predictions are made. If not supplied, today's date will be used.
- `use_parsl`: (bool: true or false) Should parsl be used to make predicting run faster? If `false` only one GPU is used. 
- `write_log`: (dict of bools) A nested json directive 
    - `logdir`: (str) name of the log directory
    - `logtypes`: Controls what logs should be written. Any of `memory`, `cache`, `error` or `time`, and possible values are `true` or `false`.
- `parsl_parameters`: (dict of bools) If `use_parsl` is true, these parameters are passed to parsl. Dictionary keys are `job_name`, `num_of_full_nodes`, `walltime`, `min_num_blocks`, `queue`, and `max_num_blocks`
    - `job_name`: (str) e.g. "my_predictions"
    - `num_of_full_nodes`: (int) e.g. 10
    - `walltime`: (str) "HH:MM:SS" e.g. "01:10:59"
    - `init_blocks`: (int) e.g. 1
    - `min_num_blocks`: (int) e.g. 0
    - `max_num_blocks`: (int) e.g. 1
    - `queue`: (str) e.g. "preemptable", "full-node"
    - `account`: (str) e.g. "haky", "pi-haky", "TFXcan", "covid-ct"
    - `hpc`: (str) e.g. "beagle3", "polaris", "theta"
    - `provider`: (str) e.g. "highthroughput", "local"
    - `worker_init`: (str) e.g. "conda activate ./enformer-predict-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./enformer-predict-tools/lib"


## To-do
- [ ] > :heavy_exclamation_mark: Allow the submission of a bash script to the queue.
- [X] Check that `os.path.join` paths work correctly.
- [X] User should supply the intervals list within the predictions folder.
- [ ] Provide a module to ensure that inputs are appropriate and available
- [ ] > :heavy_exclamation_mark: Provide checks to ensure that necessary files and folders are available.
- [X] Change `hg38_fasta_file` to a better, more general name like `fasta_file`
- [X] Ensure that the pipeline works with relative paths.
- [X] Change "motif" to "region" in the predictions_log files names
- [X] Pipeline should be able to predict on reverse complements
- [X] Create config templates for different servers (beagle3, polaris, theta e.t.c.)


## Updates

#### Sat Apr 19 2023
- [X] Changed argument passed to `enformer_predict.py` from `--param_config` to `--parameters` to make it more understandable.
- [X] Removed the need to create a temporary config when predicting. Now, the config file is passed as a parameter to the necessary modules.

#### Mon Jun 5 2023
- [X] Added a new function to separate successful from unsuccessful predictions. This is to allow for re-running of unsuccessful predictions.