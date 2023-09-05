# Description: This script is used for ENFORMER inference
# Author: Temi
# Date: Wed 25 Jan 2023
# Usage: 

import os, sys, json, re
import pandas as pd 
import time
import parsl
from datetime import date
import argparse

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--parameters", help="Path to JSON file of parameters and directives to be used by ENFORMER", type=str, required=True)
args = parser.parse_args()

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
batch_utils_path = os.path.join(script_path, 'modules')
sys.path.append(batch_utils_path)

import loggerUtils
import directives

# main 
def main():

    params_path = args.parameters

    if not os.path.isabs(params_path):
        params_path = os.path.abspath(params_path)

    # I need predict_utils_two.py downstream and and it is quite tricky to use with parsl
    p_two = os.path.join(script_path, 'modules', 'predictUtils_two.py')

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

    # # write the params_path to a config.json file in a predefined folder
    # tmp_config_data = {'params_path': params_path}
    # tmp_config_file = os.path.join(batch_utils_path, f'tmp_config_{prediction_data_name}_{prediction_id}.json')
    # with open(tmp_config_file, mode='w') as cj:
    #     json.dump(tmp_config_data, cj)

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = project_dir

    if not os.path.isdir(job_log_dir):
        os.makedirs(job_log_dir)

    # set parsl directives
    directives.parsl_directives(use_parsl, parsl_parameters)
    
    # importing this module does not work; best to execute it here
    predict_utils_one = os.path.join(script_path, 'modules', 'predictUtils_one.py')
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
    a = pd.read_table(interval_list_file, sep=' ', header=None).dropna(axis=0) #.drop_duplicates(subset=['region', 'sample', 'status', 'sequence_source'], keep='last')
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
            exclude_these_regions = pd.read_csv(exclude_file)['region'].tolist()
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
                print(f'WARNING - {chromosome} sites are not available.')
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
                sample_app_futures.append(prediction_fxn(batch_regions=list(region_list), samples=list(sample_list), path_to_vcf = chr_vcf_file, batch_num = count, script_path=script_path, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source, tmp_config_path=params_path, p_two=p_two))   

                count = count + 1 

    if use_parsl == True:
        print(f'INFO - Executing parsl futures for {len(sample_app_futures)} parsl apps')
        exec_futures = [q.result() for q in sample_app_futures] 
        #print(sample_app_futures)
        print(f'INFO - Finished predictions for all')
    elif use_parsl == False:
        print(f'INFO - Finished predictions for: {sample_app_futures} ...')

    # # just so I don't have to deal with having too many resources, I can request a small amount of resource
    # check_fxn = return_check_function(use_parsl)
    # SUMMARY_FILE = os.path.join(job_log_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.summary')
    # summary_exec = []
    # for sample in id_list:
    #     if os.path.isfile(os.path.join(prediction_logfiles_folder, f"{sample}_log.csv")):
    #         summary_exec.append(check_fxn(sample=sample, predictions_folder=output_dir, log_folder=prediction_logfiles_folder, interval_list_file=interval_list_file, exclude_csv=exclude_file, sequence_source=sequence_source))

    # if use_parsl:
    #     summary_exec = [q.result() for q in summary_exec]
    #     parsl.clear() # end parsl

    # #summary_exec = list(set(summary_exec))
    # for i, qr in enumerate(summary_exec):
    #     loggerUtils.write_logger(log_msg_type=qr['logtype'], logfile=SUMMARY_FILE, message=qr['logmessage'])

    # # regex the summary file and save the failed ones e.t.c to csv
    # # --- there is a better way to do this but for now, this will do

    # warning_pattern = r"^\[WARNING.*For\s(\w+|\d+).*"
    # success_pattern = r"^\[INFO.*For\s(\w+|\d+).*"
    # with open(SUMMARY_FILE, 'r') as f:
    #     lines = list(set(f.readlines()))
    # # print(line)
    # warning_result = [re.search(warning_pattern, l).group(1) for l in lines if not re.search(warning_pattern, l) is None]
    # success_result = [re.search(success_pattern, l).group(1) for l in lines if not re.search(success_pattern, l) is None]
    # pd.DataFrame(list(set(warning_result))).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.unsuccessful_predictions.csv'), index=False, header=False)
    # pd.DataFrame(list(set(success_result))).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.successful_predictions.csv'), index=False, header=False)

    # collect the successfule predictions
    # successful_predictions = list(set([q['sample'] for q in summary_exec if q['logtype'] == 'INFO']))
    # unsuccessful_predictions = list(set([q['sample'] for q in summary_exec if q['logtype'] == 'WARNING']))
    # pd.DataFrame({'successful_predictions':successful_predictions}).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.successful_predictions.csv'), index=False, header=False)
    # pd.DataFrame({'unsuccessful_predictions':unsuccessful_predictions}).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.unsuccessful_predictions.csv'), index=False, header=False)

    # print(f'INFO - Check {SUMMARY_FILE} for a summary of the entire run.')
    print(f'INFO - Check `{metadata_dir}` for successful and unsuccessful predictions.')

    # == After predictions are complete, a json file will be written out to help with aggregation
    print(f'INFO - Writing `aggregation_config_{prediction_data_name}_{prediction_id}.json` file to {metadata_dir}')
    agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'prediction_logfiles_folder':prediction_logfiles_folder, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'prediction_id':prediction_id, 'individuals': None if sequence_source in ['reference', 'random'] else individuals, 'n_individuals':n_individuals if sequence_source == 'personalized' else None}

    with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{prediction_id}.json', mode='w')) as wj:
        json.dump(agg_dt, wj)

    # remove temporatry config file
    # print(f"INFO - Cleaning up: Removing temporary config file at {tmp_config_file}")
    # os.remove(tmp_config_file)
              
if (__name__ == '__main__') or (__name__ == 'enformer_predict'):
    #check_input_parameters.check_inputs(args.param_config)

    # job stats will always be written
    import time
    job_start = time.perf_counter()
    main()
    job_end = time.perf_counter()

    job_runtime = job_end - job_start
    print(f'INFO - Completed job in {job_runtime} seconds.')





# def slice_bins(locus, bin_size=128, nbins=896):
#     import math
#     locus = locus.split('_')
#     start = int(locus[1])
#     end = int(locus[2])
#     midn = math.ceil((start + end) / 2)
#     print(f'Middle locations is: {midn}')
#     nstart = midn - 57344 # (128*896) / 2
#     #print(f'New start is: {nstart}')
#     nend = midn + (57344 - 1)
#     #print(f'New end is: {nend}')
#     cnt_start = nstart
#     slice_start = 0
#     while cnt_start <= start:
#         cnt_start += bin_size
#         slice_start += 1
#     cnt_start = nstart
#     slice_end = 0
#     while cnt_start <= end:
#         cnt_start += bin_size
#         slice_end += 1
#     return(slice_start, slice_end)

#     # if start % 128 != 0:
#     #     slice_start = math.ceil(start // 128)
#     # else: 
#     #     slice_start = math.floor(start // 128)
#     # return(slice_start)

# a = 'chr2_223164458_223164688'
# sum_bins(a)

# 'chr6_86160197_86176777'             
# 'chr2_223164458_223164688'

