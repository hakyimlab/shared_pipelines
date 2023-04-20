# Usage: This module is used to predict with ENFORMER on genomic sequences
# Author: Temi
# Date: Thurs Feb 2 2023

import json
import os
import tensorflow as tf
import time
import numpy as np
from datetime import date

global module_path, write_log, sequence_source, grow_memory

# read in the config_file
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
#module_path = os.path.abspath(whereis_script)

#if __name__ == '__main__':
if __name__ == 'predictUtils_two':

    with open(f'{whereis_script}/tmp_config.json') as f:
        parameters = json.load(f)
        parameters_file = parameters['params_path']
    
    with open(f'{parameters_file}') as f:
        parameters = json.load(f)

        prediction_data_name = parameters['prediction_data_name']
        prediction_id = parameters['prediction_id']
        run_date = parameters['date'] if parameters['date'] is not None else date.today().strftime("%Y-%m-%d")

        if parameters['sub_dir'] == True:
            project_dir = os.path.join(parameters['project_dir'], 'predictions_folder', f'{prediction_data_name}_{prediction_id}', f'predictions_{run_date}')
        else:
            project_dir = os.path.join(parameters['project_dir'], f'{prediction_data_name}_{prediction_id}', f'predictions_{run_date}')

        sequence_source = parameters['sequence_source']
        write_logdir = parameters['write_log']['logdir']
        write_log = parameters["write_log"]
        fasta_file = parameters['fasta_file']
        model_path = parameters['model_path']
        bins_indices_raw = parameters['bins_to_save']
        tracks_indices_raw = parameters['tracks_to_save']
        reverse_complement = parameters["reverse_complement"]

if any([write_log['logtypes'][ltype] for ltype in ['memory', 'error', 'time', 'cache']]):
    write_log['logdir'] = os.path.join(project_dir, write_logdir)
else:
    write_log['logdir'] = os.path.join(project_dir, write_logdir)

#import checksUtils
import loggerUtils
import sequencesUtils
import predictionUtils
import checksUtils
import collectUtils

global enformer_model
global fasta_extractor
global predictions_expected_shape

#enformer_model = predictionUtils.get_model(model_path)
grow_memory = True
#print(f'GPU Memory before calling batch predict function is {loggerUtils.get_gpu_memory()}')

bins_indices, tracks_indices = collectUtils.parse_bins_and_tracks(bins_indices_raw,tracks_indices_raw)

# Check prediction size for correctness
if bins_indices == None:
    if tracks_indices == None:
        predictions_expected_shape = (896,5313)
    else:
        predictions_expected_shape = (896,len(tracks_indices))
else:
    if tracks_indices == None:
        predictions_expected_shape = (len(bins_indices), 5313)
    else:
        predictions_expected_shape = (len(bins_indices),len(tracks_indices))

#print(bins_indices,tracks_indices)

def enformer_predict_on_batch(batch_regions, samples, logging_dictionary, path_to_vcf, batch_num, output_dir, prediction_logfiles_folder, sequence_source):

    # this could mean
    # - check_queries returned nothing because
        # - nothing should be returned - good ; and it should return none
        # - something is wrong with check_queries ; and I should fix that
    if (not batch_regions) or (batch_regions is None):
        raise Exception(f'[INFO] There are no regions in this batch {batch_num}.')

    # print(f'batch_regions are: {batch_regions}')
    # # print(f'samples are: {samples}')
    # # print(f'path_to_vcf are: {path_to_vcf}')
    # print(f'output_dir are: {output_dir}')
    # print(f'prediction_logfiles_folder are: {prediction_logfiles_folder}')
    # print(f'sequence_source are: {sequence_source}')

    #print(f'GPU Memory at start of batch {batch_num} predict function is {loggerUtils.get_gpu_memory()}')

    if grow_memory == True:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                # Currently, memory growth needs to be the same across GPUs
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
            except RuntimeError as e:
                raise Exception(f'[RUNTIME ERROR] Batch number: {batch_num} of {type(e)} in module {__name__}')

    try:
        #model = predictionUtils.get_model(model_path)
        #model = enformer_model # check global definitions
        enformer_model = predictionUtils.get_model(model_path)
        fasta_extractor = sequencesUtils.get_fastaExtractor(fasta_file)

        #print('Fasta and model successfully loaded')
        logger_output = []
        dlist = {sample: [elem['query'] for elem in logging_dictionary[sample]] for sample in logging_dictionary.keys()}

        for input_region in batch_regions: # input_region is chr1_10_20
            #print(input_region)
            # filter the samples
            v_samples = checksUtils.return_samples_to_predict_on(query=input_region, logging_list_per_sample=dlist)

            #print(f'Creating sequences for {input_region}')
            tic = time.perf_counter()

            samples_enformer_inputs = sequencesUtils.create_input_for_enformer(query_region=input_region, samples=v_samples, path_to_vcf=path_to_vcf, fasta_func=fasta_extractor, hap_type = 'both', resize_for_enformer=True, resize_length=None, write_log=write_log, sequence_source=sequence_source, reverse_complement=reverse_complement)

            toc = time.perf_counter()

            retrieve_time = (toc - tic)/len(v_samples) #len(list(samples_enformer_inputs['sequence'].keys()))

            #print(f'Region {input_region} sequences successfully created')

            # check that all the samples are accounted for
            #print(sorted(list(samples_enformer_inputs['sequence'].keys())))
            if samples_enformer_inputs is None:
                logger_output.append(2)
                continue

            elif samples_enformer_inputs is not None:
                if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                    if sorted(v_samples) != sorted(list(samples_enformer_inputs['sequence'].keys())):
                        missing_samples = [s for s in sorted(v_samples) if s not in sorted(list(samples_enformer_inputs['sequence'].keys()))]
                        # remove missing samples from v_samples
                        if missing_samples:
                            print(f"WARNING - Removing {len(missing_samples)} missing samples from the input samples list.")
                            v_samples = [s for s in v_samples if s not in missing_samples]
                        print(f"[WARNING] Some samples cannot be found. But this job will continue")
                #logging_info_list = [] # collect all logging information here
                for sample in v_samples:
                    if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                        tic = time.perf_counter()

                        unfiltered_sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'][sample])
                    elif samples_enformer_inputs['metadata']['sequence_source'] in ['ref', 'random']:
                        tic = time.perf_counter()

                        unfiltered_sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'])
                    
                    toc = time.perf_counter()
                    predict_time = toc - tic
                    
                    # check that the predictions have the appropriate shapes
                    sample_predictions = {}
                    for hap in unfiltered_sample_predictions.keys():

                        haplo_prediction_cur = np.squeeze(unfiltered_sample_predictions[hap], axis=0)
                        #print(haplo_prediction_cur)
                        sample_predictions[hap] = collectUtils.collect_bins_and_tracks(haplo_prediction_cur, bins_indices, tracks_indices)
                        #print(sample_predictions[hap])

                        if sample_predictions[hap].shape != predictions_expected_shape:
                            #print(sample_predictions[hap])
                            raise Exception(f'ERROR - {sample}\'s {hap} predictions shape is {sample_predictions[hap].shape} and is not equal to expected shape {predictions_expected_shape}.')
                        else:
                            print(f'Sample {sample} {input_region} {hap} predictions are of the correct shape:  {sample_predictions[hap].shape}')
                        
                    # otherwise, you can save the predictions ; prediction will be reshaped to (17, 5313) here
                    sample_logging_info = loggerUtils.save_haplotypes_h5_prediction(haplotype_predictions=sample_predictions, metadata=samples_enformer_inputs['metadata'], output_dir=output_dir, sample=sample)

                    print(f'Sample {sample} {input_region} haplotypes predictions have been saved.')

                    # check logging info/dictionary for the sample and the region
                    logging_type = checksUtils.return_sample_logging_type(sample=sample, query_region=input_region, logging_dictonary=logging_dictionary)
                    #print(logging_type)

                    if logging_type == 'y':
                        if (sample_logging_info is not None) and (len(sample_logging_info) == 4):
                            predictions_log_file = os.path.join(prediction_logfiles_folder, f'{sample}_log.csv')
                            sample_logging_info.extend([predict_time, retrieve_time])
                            logger_output.append(loggerUtils.log_predictions(predictions_log_file=predictions_log_file, what_to_write=sample_logging_info))
                        print(f'Sample {sample} {input_region} haplotypes predictions have been logged.')
                    elif logging_type == 'n':
                        logger_output.append(1)
                        continue
                
            if write_log['logtypes']['memory']:
                mem_use = loggerUtils.get_gpu_memory() # [123, 456]#
                msg_mem_log = f"[MEMORY] (GPU) at the end of batch {batch_num} prediction: free {mem_use[0]} mb, used {mem_use[1]} mb on " 
                if tf.config.list_physical_devices('GPU'):
                    MEMORY_LOG_FILE = os.path.join(write_log['logdir'], "memory_usage.log")
                    loggerUtils.write_logger(log_msg_type = 'memory', logfile = MEMORY_LOG_FILE, message = msg_mem_log)
            else:
                mem_use = loggerUtils.get_gpu_memory() # [123, 456]#
                msg_mem_log = f"[MEMORY] (GPU) at the end of batch {batch_num} prediction: free {mem_use[0]} mb, used {mem_use[1]} mb on " 
                print(msg_mem_log)

            if write_log['logtypes']['cache']:
                msg_cac_log = f'[CACHE] (model) at batch {batch_num}: [{predictionUtils.get_model.cache_info()}]'
                CACHE_LOG_FILE = os.path.join(write_log['logdir'], 'cache_usage.log')
                loggerUtils.write_logger(log_msg_type = 'cache', logfile = CACHE_LOG_FILE, message = msg_cac_log)
            else:
                msg_cac_log = f'[CACHE] (model) at batch {batch_num}: [{predictionUtils.get_model.cache_info()}]'
                print(msg_cac_log)

        return(logger_output)
    
    except (TypeError, AttributeError) as tfe:
        if write_log['logtypes']['error']:
            if tf.config.list_physical_devices('GPU'):
                mem_use = loggerUtils.get_gpu_memory()
                err_mem_log = f"[ERROR] GPU memory error of type {type(tfe).__name__} for batch {batch_num}): free {mem_use[0]} mb, used {mem_use[1]} mb on {loggerUtils.get_gpu_name()}"
                MEMORY_ERROR_FILE = os.path.join(write_log['logdir'], 'error_details.log')
                loggerUtils.write_logger(log_msg_type = 'error', logfile = MEMORY_ERROR_FILE, message = err_mem_log)
        else:
            raise Exception(f'[ERROR] of {type(tfe).__name__} at batch {batch_num}')



