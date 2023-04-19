# Usage: This script is used to predict on batches using ENFORMER on individuals' regions
# Author: Temi
# Date: Wed 25 Jan 2023

#from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import time
import parsl
from datetime import date
import argparse
import pathlib

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--param_config", help="Path to JSON file of parameters and directives to be used by ENFORMER", type=str)
args = parser.parse_args()

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
batch_utils_path = os.path.join(script_path, 'batch_utils')
sys.path.append(batch_utils_path)
#print(sys.path)

import loggerUtils
import directives

# main 
def main():

    params_path = args.param_config

    with open(f'{params_path}') as f:
        parameters = json.load(f)

        prediction_data_name = parameters['prediction_data_name']
        prediction_id = parameters['prediction_id']
        run_date = parameters['date'] if parameters['date'] is not None else date.today().strftime("%Y-%m-%d")

        if parameters['sub_dir'] == True:
            project_dir = os.path.join(parameters['project_dir'], 'predictions_folder', f'{prediction_data_name}_{prediction_id}', f'predictions_{run_date}')
        elif parameters['sub_dir'] == False:
            project_dir = os.path.join(parameters['project_dir'], f'{prediction_data_name}_{prediction_id}', f'predictions_{run_date}')
        else:
            raise Exception('ERROR - `sub_dir` argument must be a boolean, either true or false')

        interval_list_file = parameters['interval_list_file']
        predictions_log_dir = os.path.join(project_dir, parameters['predictions_log_dir'])
        job_log_dir = os.path.join(project_dir, parameters['write_log']['logdir'])
        n_regions = parameters["n_regions"]
        batch_regions = int(parameters['batch_regions'])
        use_parsl = parameters['use_parsl']
        parsl_parameters = parameters['parsl_parameters']
        sequence_source = parameters['sequence_source']
        exclude_regions = parameters["exclude_regions"]
        reverse_complement = parameters["reverse_complement"]
        create_hdf5_file = parameters["create_hdf5_file"]
    
        metadata_dir = parameters['metadata_dir']
        if not os.path.isdir(metadata_dir):
            os.makedirs(metadata_dir)

        output_dir = os.path.join(project_dir, parameters['output_dir'])
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        if int(n_regions) == -1:
            n_regions = None
        elif int(n_regions) > 0:
            n_regions = (n_regions) if isinstance(n_regions, int) else None

        # personalized parameters 
        individuals = parameters['individuals'] if sequence_source == 'personalized' else None
        vcf_files_dict = parameters['vcf_files'] if sequence_source == 'personalized' else None

        if sequence_source == 'personalized':
             # use only the chromosomes that have been made available in the config file vcf params
            print(f'INFO - Sequence source is {sequence_source}. Using a reference genome + vcf files.')
            chromosomes = list(vcf_files_dict['files'].keys())

            batch_individuals = parameters["batch_individuals"]
            n_individuals = int(parameters['n_individuals'])
        # list of chromosomes (if the sequence source is reference)
        elif sequence_source == 'reference':
            print(f'INFO - Sequence source is {sequence_source}. Using a reference genome.')
            chromosomes = [f'chr{i}' for i in range(1, 23)]
            chromosomes.extend(['chrX'])

        if reverse_complement:
            print(f'INFO - Predicting on reverse complements too')

    # write the params_path to a config.json file in a predefined folder
    tmp_config_data = {'params_path': params_path}
    with open(os.path.join(batch_utils_path, 'tmp_config.json'), mode='w') as cj:
        json.dump(tmp_config_data, cj)

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = project_dir

    if not os.path.isdir(job_log_dir):
        os.makedirs(job_log_dir)

    # set parsl directives
    directives.parsl_directives(use_parsl, parsl_parameters)
    
    # importing this module does not work; best to execute it here
    predict_utils_one = f'{script_path}/batch_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    # decorate the prediction function with or without parsl
    prediction_fxn = return_prediction_function(use_parsl)

    # determine what individuals to predict on and all that
    if sequence_source == 'personalized':
        
        if isinstance(individuals, list):
            id_list = individuals
            pass
        elif isinstance(individuals, type('str')):
            if os.path.isfile(individuals):
                if n_individuals == -1:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:]
                elif n_individuals > 0:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:(n_individuals)]
            else:
                id_list = [individuals]
        print(f'INFO - Found {len(id_list)} individuals to predict on')

    elif sequence_source == 'reference':
        id_list = [prediction_data_name]
        print(f'INFO - Found one reference set named {id_list[0]} to predict on')
    elif sequence_source == 'random':
        id_list = [prediction_data_name]
        print(f'INFO - Prediction will be on a randomly generated set')

    # set log files to be put in a folder and touch the log files per sample
    prediction_logfiles_folder = predictions_log_dir
    if not os.path.isdir(prediction_logfiles_folder):
        os.makedirs(prediction_logfiles_folder)
        
    # list of intervals to be predicted on
    a = pd.read_table(interval_list_file, sep=' ', header=None).dropna(axis=0) #.drop_duplicates(subset=['motif', 'sample', 'status', 'sequence_source'], keep='last')
    list_of_regions = a[0].tolist()[0:(n_regions)] # a list of queries
    print(f'INFO - Found {len(list_of_regions)} regions to be split into batches with at most {batch_regions} regions in each batch.')

    # filter the list of chromosomes to be compatible with the available regions
    chromosomes = list(set([r.split('_')[0] for r in list_of_regions]))
    #print(f'INFO - Chromosomes to predict on are: {chromosomes}')

    # should some regions be excluded?
    if exclude_regions == True:
        # seach for the invalid_regions.csv file
        exclude_file = os.path.join(job_log_dir, 'invalid_queries.csv')
        if os.path.isfile(exclude_file):
            exclude_these_regions = pd.read_csv(exclude_file)['motif'].tolist()
            print(f'INFO - Found regions to be excluded from the input regions.')
            list_of_regions = [l for l in list_of_regions if l not in exclude_these_regions]  
            print(f'INFO - Updated number of regions to predict on is {len(list_of_regions)}')
        else:
            print(f'INFO - No regions to exclude yet. You either did not supply a file, this is the first run, or there are truly no regions to exclude')
            exclude_these_regions = None
    else:
        exclude_file = None
    
    # batch the samples too
    # if you have 1000 individuals, it may be too much
    if len(id_list) > 5:
        if batch_individuals is not None:
            if isinstance(batch_individuals, int):
                sample_batches = list(generate_batch_n_elems(id_list, n = batch_individuals)) # 5 samples in each batch
                print(f'INFO - There are more than 10 individuals. Predictions will be done for every {batch_individuals} individuals.')
            else:
                raise Exception(f'ERROR - argument `batch_individuals` is not a str type. You supplied a {type(batch_individuals).__name__}')
        else:
            print(f'INFO - You have multiple individuals/samples and have not supplied how to batch them. For efficient use of resources, use the `batch_individuals` argument.')
    else:
        sample_batches = [id_list] # put the list in a list
        print(f'INFO - There seem to be just one sample i.e. {sample_batches}. No need to batch.')

    # to make this fast, pass multiple regions to one parsl app
    sample_app_futures = []
    for sample_list in sample_batches:
        for chromosome in chromosomes:
            #print(chromosome)
            chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]
            if sequence_source == 'personalized':
                chr_vcf_file = os.path.join(vcf_files_dict['folder'], vcf_files_dict['files'][chromosome])
            elif sequence_source == 'reference':
                chr_vcf_file = None

            if not chr_list_of_regions:
                print(f'WARNING - {chromosome} motif sites are not available.')
                continue

            # I want many regions to be put in a parsl app
            if len(chr_list_of_regions) > batch_regions:
                region_batches = generate_batch_n_elems(chr_list_of_regions, n=batch_regions) # batch_regions total batches
            else:
                region_batches = [chr_list_of_regions]
            
            count = 0
            for region_list in region_batches:
                #print(len(sample_list))
                #print(f'{len(region_list)} regions in {chromosome} for {len(sample_list)} samples')
                sample_app_futures.append(prediction_fxn(batch_regions=list(region_list), samples=list(sample_list), path_to_vcf = chr_vcf_file, batch_num = count, script_path=script_path, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source))   

                count = count + 1 

    if use_parsl == True:
        print(f'INFO - Executing parsl futures for {len(sample_app_futures)} parsl apps')
        exec_futures = [q.result() for q in sample_app_futures] 
        #print(sample_app_futures)
        print(f'INFO - Finished predictions for all')
    elif use_parsl == False:
        print(f'INFO - Finished predictions for: {sample_app_futures} ...')

    # just so I don't have to deal with having too many resources, I can request a small amount of resource
    check_fxn = return_check_function(use_parsl)
    SUMMARY_FILE = os.path.join(job_log_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.summary')
    summary_exec = []
    for sample in id_list:
        if os.path.isfile(os.path.join(prediction_logfiles_folder, f"{sample}_log.csv")):
            summary_exec.append(check_fxn(sample=sample, predictions_folder=output_dir, log_folder=prediction_logfiles_folder, interval_list_file=interval_list_file, exclude_csv=exclude_file, sequence_source=sequence_source))

    if use_parsl:
        summary_exec = [q.result() for q in summary_exec]
        parsl.clear() # end parsl

    #summary_exec = list(set(summary_exec))
    for i, qr in enumerate(summary_exec):
        loggerUtils.write_logger(log_msg_type=qr['logtype'], logfile=SUMMARY_FILE, message=qr['logmessage'])
    print(f'INFO - Check {SUMMARY_FILE} for a summary of the entire run.')

    # == After predictions are complete, a json file will be written out to help with aggregation
    print(f'INFO - Writing `aggregation_config_{prediction_data_name}_{prediction_id}.json` file to {metadata_dir}')
    agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'prediction_logfiles_folder':prediction_logfiles_folder, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'prediction_id':prediction_id, 'individuals': None if sequence_source in ['reference', 'random'] else individuals, 'n_individuals':n_individuals if sequence_source == 'personalized' else None}

    with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{prediction_id}.json', mode='w')) as wj:
        json.dump(agg_dt, wj)

    # remove temporatry config file
    print(f"INFO - Removing temporary config file at {os.path.join(batch_utils_path, 'tmp_config.json')}")
    os.remove(os.path.join(batch_utils_path, 'config.json'))
              
if (__name__ == '__main__') or (__name__ == 'enformer_predict'):
    #check_input_parameters.check_inputs(args.param_config)

    # job stats will always be written
    import time
    job_start = time.perf_counter()
    main()
    job_end = time.perf_counter()

    job_runtime = job_end - job_start
    print(f'INFO - Completed job in {job_runtime} seconds.')

    



# if create_hdf5_file == True:
#         print(f'[INFO] Creating HDF5 database(s)')
#         finished_predictions = pd.read_csv(logfile_csv)
#         make_db = make_h5_db_parsl(use_parsl = use_parsl)

#         db_parsl = []
#         for each_id in id_list:
#             motifs_list = finished_predictions.loc[finished_predictions['individual'] == each_id, ].motif.values.tolist()
#             motifs_list = list(set(motifs_list))

#             print(f'[INFO] Creating HDF5 database for {each_id} for {len(motifs_list)} predictions.')

#             motifs_list_paths = [f'{output_dir}/{each_id}/{i}_predictions.h5' for i in motifs_list]
#             csv_file = f'{output_dir}/{each_id}_{prediction_id}_predictions.csv'
#             h5_file = f'{output_dir}/{each_id}_{prediction_id}_predictions.hdf5'
#             db_parsl.append(make_db(h5_file = h5_file, csv_file = csv_file, files_list = motifs_list, files_path = motifs_list_paths, dataset = each_id))

        
#         print(db_parsl)
#         if use_parsl == True:
#             exec_parsl = [q.result() for q in tqdm.tqdm(db_parsl, desc=f'[INFO] Executing database futures.')] 
#             print(exec_parsl)



    # #chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]
    # for each_id in id_list:
    #     app_futures = [] # collect futures here
        
    #     # filter where each_id is in the log file
    #     if not logfile is None:
    #         id_logfile = logfile.loc[logfile['individual'] == each_id, : ]
    #     elif logfile is None:
    #         id_logfile = logfile
        
    #     if sequence_source == 'personalized':
    #         if not os.path.exists(f'{output_dir}/{each_id}'):
    #             print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
    #             os.makedirs(f'{output_dir}/{each_id}/haplotype1')
    #             os.makedirs(f'{output_dir}/{each_id}/haplotype2')
    #     elif sequence_source == 'reference':
    #         if not os.path.exists(f'{output_dir}/{each_id}'):
    #             print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
    #             os.makedirs(f'{output_dir}/{each_id}/haplotype0')
    #     elif sequence_source == 'random':
    #         if not os.path.exists(f'{output_dir}/{each_id}'):
    #             print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
    #             os.makedirs(f'{output_dir}/{each_id}/haplotype0')

    #     print(f'[INFO] Collecting appfutures for {each_id}')
    #     for chromosome in chromosomes:
    #         chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]

    #         if not chr_list_of_regions:
    #             continue

    #         chr_vcf_file = os.path.join(vcf_files_dict['folder'], vcf_files_dict['files'][chromosome])
    #         if sequence_source == 'personalized':
    #             def make_cyvcf_object(vcf_file=chr_vcf_file, sample=each_id):
    #                 import cyvcf2
    #                 return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))
    #         elif sequence_source == 'reference':
    #             make_cyvcf_object = None
    #         elif sequence_source == 'random':
    #             make_cyvcf_object = None

    #         #print(f'[INFO] Collecting parsl appfuture batches for {chromosome}: {len(chr_list_of_regions)}')
    #         #for batch_query in batches:

    #         batches = generate_batch(chr_list_of_regions, batch_n=batch_regions)
    #         count = 0
    #         #app_futures = []
    #         for batch_query in batches:
    #             count = count + 1
    #             app_futures.append(prediction_fxn(batch_regions=batch_query, batch_num = count, id=each_id, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=id_logfile, predictions_log_file=logfile_csv))
    



    #for sample_list in sample_batches:
    #     # collect all chromosomes
    #     sample_app_futures = []
    #     for chromosome in chromosomes:
    #         chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]

    #         if sequence_source == 'personalized':
    #             chr_vcf_file = os.path.join(vcf_files_dict['folder'], vcf_files_dict['files'][chromosome])
    #         elif sequence_source == 'reference':
    #             chr_vcf_file = None

    #         if not chr_list_of_regions:
    #             print(f'[WARNING] {chromosome} motif sites are not available.')
    #             continue

    #         region_batches = generate_n_batches(chr_list_of_regions, batch_n=batch_regions) # batch_regions total batches
    #         count = 0
    #         for region_list in region_batches:
    #             sample_app_futures.append(prediction_fxn(batch_regions=region_list, samples=sample_list, path_to_vcf = chr_vcf_file, batch_num = count, script_path=script_path, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source))

    #  if sequence_source == 'reference':
    #     print(f'INFO - Writing `aggregation_config_{prediction_data_name}_{prediction_id}.json` file to {metadata_dir}')

    #     agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'prediction_logfiles_folder':prediction_logfiles_folder, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'prediction_id':prediction_id, 'individuals':None, 'n_individuals':n_individuals}

    #     with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{prediction_id}.json', mode='w')) as wj:
    #         json.dump(agg_dt, wj) 

    # elif sequence_source == 'personalized':
    #     print(f'INFO - Writing `aggregation_config_{prediction_data_name}_{prediction_id}.json` file to {metadata_dir}')

    #     agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'prediction_logfiles_folder':prediction_logfiles_folder, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'prediction_id':prediction_id, 'individuals':individuals, 'n_individuals':n_individuals}

            

    # elif sequence_source == 'random':
    #     print(f'INFO - Writing `aggregation_config_{prediction_data_name}_{prediction_id}.json` file to {metadata_dir}')

    #     agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'prediction_logfiles_folder':prediction_logfiles_folder, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'prediction_id':prediction_id, 'individuals':None, 'n_individuals':n_individuals}

    #     with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{prediction_id}.json', mode='w')) as wj:
    #         json.dump(agg_dt, wj)      