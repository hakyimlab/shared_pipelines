# Usage: This module is used to predict with ENFORMER on individual genomes
# Author: Temi
# Date: Thurs Feb 2 2023


#@python_app
def run_batch_predictions(batch_regions, samples, path_to_vcf, batch_num, script_path, output_dir, prediction_logfiles_folder, sequence_source, tmp_config_path = None, p_two=None): #
    """
    Predict and save on a given batch of regions in the genome

    This function also filters the regions in the batch for those that have been predicted and logged

    Parameters:
        batch_regions: list
            A list of regions with each element (region) in the form `chr_start_end`.
        batch_num: num
            The number of the batch e.g. batch 1, batch 2, e.t.c.
        samples: list
            a list of samples: [a, b, c]
        output_dir: str (path)
            Where should the predictions be saved? Predictions are saved as `{sample}/{haplotype0, haplotype1, haplotype2}/{region}_predictions.h5` 
        path_to_vcf: str (path)
            The path to the vcf file
        prediction_logfiles_folder: str (path)
            When predictions are made, they are logged in this folder and saved as {sample}_log.csv. Also useful for checking the log to prevent re-predicting an existing prediction. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.
        sequence_source: one of 'personalized' or 'reference'

    Returns: num
        A single value of either 0 (if predictions were successful) or 1 (if there was a error).
        Check the call logs or stacks for the source of the error. 
    """

    print("Running batch predictions:", sequence_source)

    import sys, os, faulthandler, time, importlib

    mpath = os.path.join(script_path, 'modules') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from

    # I want to import predictUtils_two but I need it to be dynamic and imported with the path to the config file
    try:
        if p_two is not None:
            spec = importlib.util.spec_from_file_location("predictUtils_two", p_two)
        else:
            raise Exception('p_two is None')

        predictUtils_two = importlib.util.module_from_spec(spec)
        predictUtils_two.tmp_config_path = tmp_config_path
        spec.loader.exec_module(predictUtils_two)

        import checksUtils
        #import predictUtils_two
    except ModuleNotFoundError as merr:
        raise Exception(f'[ERROR] {type(merr).__name__} at run_batch_predictions. Cannot locate either of `checkUtils` or `predictUtils_two` modules.')

    # check_these = itertools.product(samples, [batch_regions])
    # check_results = [checksUtils.check_queries(sample=cq[0], queries=cq[1], output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source) for cq in check_these]

    check_results = {sample: checksUtils.check_queries(sample=sample, queries=batch_regions, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source) for sample in samples}
    
    filtered_check_result = {k: v for k, v in check_results.items() if k is not None}
    print("filtered check results:", check_results)

    # remove empty dict values
    for sample in list(filtered_check_result.keys()):
        if not filtered_check_result[sample]:
            del filtered_check_result[sample]

    # filter out nones
    # filtered_check_result = [r for r in check_results if r is not None]
    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:
        # unlist the list of dictionaries
        # filtered_check_result = [f for f in filtered_check_result for f in f]
        # filtered_check_result = list({v['query']: v for v in filtered_check_result}.values())

        pqueries = [v for k, v in filtered_check_result.items()]
        pqueries = [l for l in pqueries for l in l]
        pqueries = list(set([d['query'] for d in pqueries]))
        # print(f'{len(pqueries)}')
        print("queries", pqueries)

        if pqueries:
            samples = list(filtered_check_result.keys())
            tic = time.perf_counter()

            reg_prediction = predictUtils_two.enformer_predict_on_batch(batch_regions=pqueries, samples=samples, logging_dictionary=filtered_check_result, path_to_vcf = path_to_vcf, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, batch_num=batch_num, sequence_source=sequence_source)
            
            toc = time.perf_counter()
            print(f'[INFO] (time) to predict on batch {batch_num} is {toc - tic}')
            return(reg_prediction) # returns 0 returned by enformer_predict
        else:
            return(1)

def return_prediction_function(use_parsl, fxn=run_batch_predictions):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    ''' 
    if use_parsl == True:
        from parsl.app.app import python_app
        return python_app(fxn)
    elif use_parsl == False:
        return fxn

def generate_n_batches(lst, batch_n, len_lst = None):
    """
    Given a list, this function yields batches of an unspecified size but the number of batches is equal to `batch_n`
    E.g. generate_batch([0, 1, 2, 3, 4, 5, 6], batch_n=2) -> (0, 1, 2, 3), (4, 5, 6)
    
    Parameters:
        lst: list
        batch_n: int
            Number of batches to return
        len_lst: None or num (length of the input list)

    Yields
        `batch_n` batches of the list
    """
    import math
    # how many per batch
    if len_lst is not None:
        n_elems = math.ceil(len_lst/batch_n)
    else:
        n_elems = math.ceil(len(lst)/batch_n)
        
    for i in range(0, len(lst), n_elems):
        yield lst[i:(i + n_elems)]

def generate_batch_n_elems(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def make_h5_db(h5_file, csv_file, files_list, files_path, dataset):
    import h5py
    import pandas as pd

    with h5py.File(f"{h5_file}", "w") as f_dst:
        #h5files = [f for f in os.listdir(f'{project_dir}/') if f.endswith(".h5")]

        dset = f_dst.create_dataset(f'{dataset}_dataset', shape=(len(files_list), 17, 5313), dtype='f4')
        for i, filename in enumerate(files_list):
            with h5py.File(files_path[i]) as f_src:
                dset[i] = f_src[filename]

    pd.DataFrame(files_list, columns=['region']).to_csv(f"{csv_file}")

    return(0)

def make_h5_db_parsl(use_parsl, fxn=make_h5_db):
    '''
    Decorate or not the `make_h5_db` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    if use_parsl == True:
        from parsl.app.app import python_app
        return python_app(fxn)
    elif use_parsl == False:
        return fxn





def check_predictions_and_logs(sample, predictions_folder, log_folder, interval_list_file, exclude_csv, sequence_source):

    import os
    import pandas as pd
    import numpy as np

    if not isinstance(sample, str):
        raise Exception(f'[ERROR] Sample argument should be a str of a valid sample name. You supplied a {type(sample).__name__} type')
    if not os.path.isdir(predictions_folder):
        raise Exception(f'[ERROR] Predictions folder does not exist. You supplied a {predictions_folder}')
    if not os.path.isdir(log_folder):
        raise Exception(f'[ERROR] Predictions log folder does not exist. You supplied a {log_folder}')
    if not isinstance(sequence_source, str):
        raise Exception(f'[ERROR] `sequence_source` argument should be a str of a valid sequence source (either `personalized` or `reference`). You supplied a {sequence_source}')

    if exclude_csv is None:
        exclude_these_regions = None
    elif os.path.isfile(exclude_csv):
        exclude_these_regions = pd.read_csv(exclude_csv)['region'].tolist()
    else:
        exclude_these_regions = None

    if os.path.isfile(interval_list_file):
        queries = pd.read_table(interval_list_file, sep=' ', header=None)[0].tolist()
        if exclude_these_regions is not None:
            queries = [q for q in queries if q not in exclude_these_regions]
    else:
        raise Exception(f'[ERROR] Interval list file {interval_list_file} does not exist.')

    # temporary solutions => remove nans
    queries = [q for q in queries if str(q) != 'nan']
    
    logged_queries = pd.read_csv(os.path.join(log_folder, f'{sample}_log.csv'))['region'].tolist()

    queries_logged = np.array([(query in logged_queries) for query in queries])
    
    if sequence_source == 'personalized': # prediction must be present in two folders
        queries_saved_h1 = [str(f'{predictions_folder}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
        queries_saved_h2 = [str(f'{predictions_folder}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
        queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
    elif sequence_source == 'reference':
        queries_saved = [str(f'{predictions_folder}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
        queries_saved = np.array([os.path.isfile(q) for q in queries_saved])

    if queries_logged.shape[0] != queries_saved.shape[0]:
        raise Exception("Lengths of queries logged and saved conditions are not the same")
    queries_condition = queries_saved * queries_logged # true should not be predicted or logged; False should be
    queries_condition = queries_condition.tolist()

    if all(queries_condition):
        message = f'SUCCESS - For {sample}, all predictions match all logged queries in the interval list files minus excluded regions, if any.'
        return({'logtype': 'info', 'logmessage':message, 'sample':sample})
    else:
        message = f'WARNING - For {sample}, either all predictions don\'t match all logged queries in the interval list files minus excluded regions or vice versa. This can happen if you have supplied a list of intervals but have chosen to predict on a subset. If this is the case, this behavior is normal. If you are unsure, please re-run the enformer prediction pipeline with the same parameters. You may supply a csv file of regions to exclude if available, but this should not matter.'
        return({'logtype': 'warning', 'logmessage':message, 'sample':sample})

def return_check_function(use_parsl, fxn=check_predictions_and_logs):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    if use_parsl == True:
        from parsl.app.app import python_app
        return python_app(fxn)
    elif use_parsl == False:
        return fxn
