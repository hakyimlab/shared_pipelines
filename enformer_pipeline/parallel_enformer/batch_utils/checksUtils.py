# Usage: A module to check that things are fine and all that
# Author: Temi
# Date: Mon Feb 6 2023

def return_samples_to_predict_on(query, logging_list_per_sample):
    ss = [sample for sample in logging_list_per_sample.keys() if query in logging_list_per_sample[sample]]
    return(ss)

def return_sample_logging_type(sample, query_region, logging_dictonary):
    result = [elem['logtype'] for elem in logging_dictonary[sample] if elem['query'] == query_region][0]

    if not result:
        return('n')

    return(result)


def check_queries(sample, queries, output_dir, prediction_logfiles_folder, sequence_source):
    """
    Check whether a given region, for an individual has been predicted and logged.

    Parameters:
        sample: str 
            The name/id of an individual/sample
        queries: an iterable
            A batch of regions each in the genome in the form `chr_start_end`.
        output_dir: str (path)
            The folder where the predictions should have been logged. 
        prediction_logfiles_folder: a directory
            A directory within which the {sample}_log.csv file should have been logged. 
            If the file is not found, the regions are returned and will be logged after prediction.
        sequence_source: str ('personalized' or 'reference')
            where is the sequence sourced from? 
            This argument is needed to know what folders to look into within the output directory.
    
    Returns: dict
        'query': the query region if it has not been logged or predictions don't exist
        'logtype': whether it should be logged if it has not been logged i.e. 'y' or 'n'

    If predictions exist and the query has been logged, this function returns None.
    """
    import pandas as pd
    import os
    import numpy as np

    if isinstance(prediction_logfiles_folder, type(None)):
        output = [{'query':query, 'logtype':'y'} for query in queries]
        return(output)

    else:
        id_logfile = os.path.join(prediction_logfiles_folder, f'{sample}_log.csv')
        if os.path.isfile(id_logfile):
            #print(f'Logfile found for {sample}')
            try:
                id_logfile = pd.read_csv(id_logfile) 
                id_logfile = id_logfile.loc[id_logfile['sample'] == sample, : ]
                # check if the file is saved
                if sequence_source == 'personalized': # prediction must be present in two folders
                    queries_saved_h1 = [str(f'{output_dir}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
                    queries_saved_h2 = [str(f'{output_dir}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
                elif sequence_source == 'reference':
                    queries_saved = [str(f'{output_dir}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved])
                # check if the file is in the logfile
                queries_logged = np.array([(query in id_logfile.region.values) for query in queries])
                # s
                if not queries_logged.shape[0] == queries_saved.shape[0]:
                    raise Exception("Lengths of queries logged and saved conditions are not the same")
                queries_condition = queries_saved * queries_logged # true should not be predicted or logged; False should be
                queries_condition = queries_condition.tolist()
                # by default , log type is yes
                # id = queries.index('chr15_54740749_54740758')
                # print(queries_saved[id])
                # print(queries_logged[id])
                
                output = [{'query': queries[i], 'logtype':'y'} if qc is False else None for i, qc in enumerate(queries_condition)]

                # get indices that are nones
                none_ids = [i for i, q in enumerate(output) if q is None]

                # safely filter out the nones 
                #queries_saved = [q for i, q in enumerate(queries_saved) if i not in none_ids]
                queries_logged = [q for i, q in enumerate(queries_logged) if i not in none_ids]
                output = [q for i, q in enumerate(output) if i not in none_ids]
                #print(output)
                # assuming all predictions have not been logged hence logtype is y

                # change those regions where queries_logged == True to 'n'
                refined_output = [q_details.update({'logtype': 'n'}) if queries_logged[i] == True else q_details.update({'logtype': 'y'}) for i, q_details in enumerate(output)]

                return(output)

            except pd.errors.EmptyDataError:
                id_logfile = None
                output = [{'query':query, 'logtype':'y'} for query in queries]
                return(output)
        else:
            output = [{'query':query, 'logtype':'y'} for query in queries]
            return(output)





