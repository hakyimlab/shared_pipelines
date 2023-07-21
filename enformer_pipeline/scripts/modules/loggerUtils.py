
import logging

# I think logger utils should check for write log directives and cache/save those

#def get_gpu_name():
#    import subprocess
#    cmd = "cat $COBALT_NODEFILE"
#    #cmd = "cat $PBS_NODEFILE"
#    a = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')
#    # cmd = "nvidia-smi -L" #"nvidia-smi --query-gpu=gpu_bus_id --format=csv"
#    # b = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')
#
#    # return(f"{a}-{b}")
#    return(a)
def get_gpu_name():
    import subprocess
    # cmd = "cat $COBALT_NODEFILE"
    # Uncomment the line above if you want to read GPU information from $COBALT_NODEFILE
    cmd = "nvidia-smi -L"  # Use this command to get GPU information using nvidia-smi
    try:
        output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        gpu_info = output.stdout.strip()

        # If the command returns GPU information, return it
        if gpu_info:
            return gpu_info

        # If the command does not return any GPU information,
        # handle the case by returning a default value or None
        return "No NVIDIA GPU found"  # You can modify the default value as needed
    except FileNotFoundError:
        # If 'nvidia-smi' is not found, handle the case gracefully
        return "nvidia-smi not found, no NVIDIA GPU information available"




#def get_gpu_memory():
#    import subprocess
#    command = "nvidia-smi --query-gpu=memory.free,memory.used --format=csv"
#    memory_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[1].split(',')
#    memory_values = [int(x.strip().split(' ')[0]) for i, x in enumerate(memory_info)]
#    return memory_values
import subprocess
import psutil
def get_gpu_memory():
    try:
        # Try running 'nvidia-smi' to get GPU memory on systems with NVIDIA GPUs
        command = "nvidia-smi --query-gpu=memory.free,memory.used --format=csv"
        memory_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[1].split(',')
        memory_values = [int(x.strip().split(' ')[0]) for i, x in enumerate(memory_info)]
    except FileNotFoundError:
        # If 'nvidia-smi' is not found, use the macOS version with Metal
        mem_info = psutil.virtual_memory()
        memory_values = [mem_info.available / (1024 * 1024), (mem_info.total - mem_info.available) / (1024 * 1024)] 

    return memory_values


def setup_logger(logger_name, log_file, level=logging.INFO):
    log_setup = logging.getLogger(logger_name)
    formatter = logging.Formatter('[%(levelname)s: %(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fileHandler = logging.FileHandler(log_file, mode='a')
    fileHandler.setFormatter(formatter)
    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)
    log_setup.setLevel(level)
    log_setup.addHandler(fileHandler)
    log_setup.addHandler(streamHandler)

    return None

def logger(msg, level, logfile):
    if logfile == 'memory'   : log = logging.getLogger('memory_log')
    if logfile == 'cache'   : log = logging.getLogger('cache_log') 
    if logfile == 'run_error'    : log = logging.getLogger('error_log')
    if logfile == 'time'    : log = logging.getLogger('time_log')
    if logfile == 'run_summary' : log = logging.getLogger('summary_log')

    if level == 'info'    : log.info(msg) 
    if level == 'warning' : log.warning(msg)
    if level == 'error'   : log.error(msg)
    if level == 'time'    : log.info(msg)
    #if level == 'summary' : log.info(msg)

    return None

def write_logger(log_msg_type, logfile, message):
    if log_msg_type == 'memory': setup_logger('memory_log', logfile) ; logger(message, 'info', 'memory')
    if log_msg_type == 'cache': setup_logger('cache_log', logfile) ; logger(message, 'info', 'cache')
    if log_msg_type == 'error': setup_logger('error_log', logfile) ; logger(message, 'error', 'run_error')
    if log_msg_type == 'time' : setup_logger('time_log', logfile) ; logger(message, 'info', 'time')
    if log_msg_type == 'info' : setup_logger('summary_log', logfile) ; logger(message, 'info', 'run_summary')
    if log_msg_type == 'warning' : setup_logger('summary_log', logfile) ; logger(message, 'warning', 'run_summary')



def save_haplotypes_h5_prediction(haplotype_predictions, metadata, output_dir, sample):

    import h5py
    import os
    import numpy

    region = metadata['region']
    for key, values in haplotype_predictions.items():
        #print(f'[INFO] This is what is being saved {values.shape}')
        houtput = os.path.join(output_dir, sample, key)
        if not os.path.exists(houtput): os.makedirs(houtput)
        h5save = str(f'{houtput}/{region}_predictions.h5')
        # for i in range(0, values.shape[0]):
            #print(values[i, :].shape)
        with h5py.File(h5save, 'w') as hf:
            hf.create_dataset(region, data=values)

    output = [metadata['region'], sample, 'completed', metadata['sequence_source']]

    return(output)

def log_predictions(predictions_log_file, what_to_write):
    #print('Writing log file')

    import os, csv
    #logfile_csv = f'{predictions_log_dir}/{id}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(predictions_log_file) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(predictions_log_file, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['region', 'sample', 'status', 'sequence_source', 'predict_time', 'retrieve_time']) # 
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

    return(0)

def log_error_sequences(error_folder, what_to_write):
    import os, csv

    open_mode = 'a' if os.path.isfile(error_folder) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(error_folder, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['region', 'reason']) # 
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

    return(0)




