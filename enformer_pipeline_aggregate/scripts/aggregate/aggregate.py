# import os
# os.chdir('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/scripts/collect_model_data/')

import os
import pandas as pd
import os, sys, json
from datetime import date
import parsl
from parsl.app.app import python_app
import subprocess
import argparse
import shutil

import multiprocessing

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--metadata_file", help="Path to the metadata file", type=str)
parser.add_argument("--agg_types", nargs='+', help='<Required> aggregation type to be used', required=True)
parser.add_argument("--output_directory", help='<Required> folder where aggregation will be saved', required=True)
parser.add_argument("--hpc", required=True)
parser.add_argument("--parsl_executor", required=True)
parser.add_argument("--delete_enformer_outputs", action='store_true')
args = parser.parse_args()

print(f'INFO - Available CPUs are {multiprocessing.cpu_count()}')
if args.delete_enformer_outputs:
    print(f'INFO - will delete enformer predictions after aggregating')
# use parsl if the num of rows of log_data is more than 10000
use_parsl = True

#todays_date = date.today().strftime("%Y-%m-%d")

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
utils_path = os.path.join(script_path, 'utilities')
sys.path.append(utils_path)

#print(sys.path)

try:
    import aggregate_tools
    import parslConfiguration
except ModuleNotFoundError:
    print(f'ERROR - aggregate_tools module not found.')

with open(f'{args.metadata_file}') as f:
    parameters = json.load(f)

    enformer_predictions_path = parameters['enformer_prediction_path']
    #log_file = parameters['log_file']
    sequence_source = parameters['sequence_source']# e.g. "kawakami" or "cistrome"
    todays_date = parameters['run_date']
    base_path = parameters['predictions_folder']
    prediction_logfiles_folder = parameters['prediction_logfiles_folder']
    prediction_id = parameters['prediction_id']
    individuals = parameters['individuals']
    n_individuals = parameters['n_individuals']
    prediction_data_name = parameters['prediction_data_name']

# determine what individuals to predict on and all that
if sequence_source == 'personalized':
    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            if n_individuals == -1:
                individuals = pd.read_table(individuals, header=None)[0].tolist()[0:]
            elif n_individuals > 0:
                individuals = pd.read_table(individuals, header=None)[0].tolist()[0:(n_individuals)]
            print(type(individuals))
        else:
            individuals = [individuals]
            print(individuals)

agg_types = args.agg_types
agg_types = agg_types[0].split(' ')
print(f'INFO - Aggregating these: {agg_types}')

print(f'INFO - Currently on {prediction_data_name}')
save_dir = os.path.abspath(args.output_directory) #os.path.join(base_path, 'aggregated_predictions') 
print(save_dir)
if not os.path.isdir(save_dir):
    os.makedirs(save_dir, exist_ok=True)

if individuals is None:
    ids_names = [prediction_data_name]
elif isinstance(individuals, list):
    ids_names = individuals

predict_utils_one = f'{script_path}/aggregate_tools.py'
exec(open(predict_utils_one).read(), globals(), globals())

if use_parsl == True:
    #bpath = os.path.join(base_path, 'modeling_pipeline')
    print(f'INFO - Using parsl.')

    if args.hpc == 'beagle3':
        parsl_params = {'working_dir':base_path, 'job_name':'aggregate_predictions', 'queue':"preemptable", 'walltime':"01:00:00", 'num_of_full_nodes': 1, 'min_num_blocks':0, 'max_num_blocks':2}
        #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
        #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
        if args.parsl_executor == 'local':
            parsl.load(parslConfiguration.beagle3_localParslConfig(parsl_params))
            #parsl.load(parslConfiguration.beagle3_tpParslConfig(parsl_params))
        elif args.parsl_executor == 'highthroughput':
            parsl.load(parslConfiguration.beagle3_htParslConfig(parsl_params))

    if args.hpc == 'midway3' or args.hpc == 'caslake':
        parsl_params = {'working_dir':base_path, 'job_name':'aggregate_predictions', 'queue':"preemptable", 'walltime':"01:00:00", 'num_of_full_nodes': 1, 'min_num_blocks':0, 'max_num_blocks':2, "worker_init": "source ~/.bashrc; conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib", "init_blocks": 1}
        #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
        #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
        if args.parsl_executor == 'local':
            parsl.load(parslConfiguration.midway3_localParslConfig(parsl_params))
            #parsl.load(parslConfiguration.beagle3_tpParslConfig(parsl_params))
        elif args.parsl_executor == 'highthroughput':
            print(parsl_params)
            parsl.load(parslConfiguration.midway3_htParslConfig(parsl_params))

    if args.hpc == 'polaris':
        parsl_params = {'working_dir':base_path, "job_name": "aggregate_predictions", "num_of_full_nodes": 4, "walltime": "06:00:00", "min_num_blocks": 0, "max_num_blocks": 8, "queue": "preemptable", "init_blocks": 1, "hpc": "polaris", "account": "covid-ct", "provider": "highthroughput", "worker_init": "source ~/.bashrc; conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/lib"}
        
        parsl.load(parslConfiguration.polaris_htParslConfig(parsl_params))

collection_fxn = return_prediction_function(use_parsl)

app_futures = []
for each_id in ids_names:
    if os.path.isfile(os.path.join(prediction_logfiles_folder, f'{each_id}_log.csv')):
        log_data = pd.read_csv(os.path.join(prediction_logfiles_folder, f'{each_id}_log.csv'))
    else:
        continue
    
    #print(log_data.head())
    for each_agg in agg_types:
        log_data = log_data.loc[log_data['sample'] == each_id, ]
        log_data = log_data.drop_duplicates(subset=['region'])

        ind_path = os.path.join(enformer_predictions_path, each_id)
        # check that the folder exists
        if not os.path.exists(ind_path):
            raise Exception(f'WARNING - {each_id} does not exist')
        else:
            app_futures.append(collection_fxn(each_id=each_id, log_data=log_data, predictions_path=ind_path, prediction_id=prediction_id, agg_types=[each_agg], save_dir=save_dir, sequence_source=sequence_source, batch_num=None))

print(f'INFO - Executing {len(app_futures)} app_futures')
if use_parsl == True:
    app_execs = [r.result() for r in app_futures]

print(f'INFO - Aggregation complete for all.')

if args.delete_enformer_outputs:
    print(f'INFO - Deleting enformer predictions')
    if os.path.isdir(enformer_predictions_path):
        shutil.rmtree(enformer_predictions_path)