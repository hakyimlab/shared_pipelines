
# Description: Provides the parsl configuration for the different systems
# Author: Temi
# Date: Sometime in early 2023


# 'source /home/temi/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib; echo Running on host `hostname`; echo Running on nodes `cat $PBS_NODEFILE`'

def polaris_htParslConfig(params):
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname
    from parsl.providers import PBSProProvider
    import os
    print(f'INFO - Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    job_name = params['job_name']
    
    # I want to put the cobalt directives 
    sch_options = ['#PBS -l filesystems=home:grand:eagle',
                    f'#PBS -N {job_name}',
                    f'#PBS -k doe'
    ]
    sch_options = '\n'.join(sch_options)
    user_opts = {
        'polaris': {
            # Node setup: activate necessary conda environment and such.
            'worker_init': params['worker_init'],
            'scheduler_options': f'{sch_options}',
            # ALCF allocation to use
            'account': params['account'],
        }
    }
    pbs_htex = Config(
        executors=[
            HighThroughputExecutor(
                label='htex_PBS',
                available_accelerators=4,  # Pin each worker to a different GPU
                max_workers=4,
                cpu_affinity = 'alternating',
                address=address_by_hostname(),
                provider=PBSProProvider(
                    launcher=MpiExecLauncher(
                        bind_cmd="--cpu-bind", overrides="--depth=64 --ppn 1"
                    ),
                    account=user_opts['polaris']['account'],
                    queue=params['queue'], #preemptable',
                    cpus_per_node=64,
                    select_options='ngpus=4',
                    worker_init=user_opts['polaris']['worker_init'],
                    scheduler_options=user_opts['polaris']['scheduler_options'],
                    walltime=params['walltime'],
                    nodes_per_block=params['num_of_full_nodes'],
                    init_blocks=1,
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                ),
            )
        ],
        run_dir=rundir,
        retries=6
    )
    return pbs_htex

def theta_htParslConfig(params):
    import parsl
    from parsl.config import Config
    from parsl.providers import CobaltProvider
    from parsl.launchers import MpiExecLauncher
    from parsl.executors import HighThroughputExecutor
    import os
    print(f'INFO - Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    job_name = params['job_name']
    
    # I want to put the cobalt directives 
    sch_options = ['#COBALT --attrs filesystems=theta-fs0,grand:enable_ssh=1',
                    f'#COBALT --jobname={job_name}'
    ]
    sch_options = '\n'.join(sch_options)
    user_opts = {
        'theta': {
            # Node setup: activate necessary conda environment and such.
            'worker_init': params['worker_init'],
            'scheduler_options': f'{sch_options}',
            # ALCF allocation to use
            'account': params['account'],
        }
    }
    #just in case there is a config loaded already
    cobalt_htex = Config(
        executors=[
            HighThroughputExecutor(
                label='htex-Cobalt',
                max_workers=8, # vs max_workers
                available_accelerators=8,
                worker_debug=True,
                cores_per_worker=24, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=CobaltProvider(
                    queue='full-node',
                    account=user_opts['theta']['account'],
                    launcher=MpiExecLauncher(),
                    walltime=params['walltime'],
                    nodes_per_block=params['num_of_full_nodes'], # number of full-nodes - 3 will launch 3 full nodes at a t   ime for one instance for each `cores_per_worker`
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                    worker_init=user_opts['theta']['worker_init'],
                    cmd_timeout=120,
                    scheduler_options=user_opts['theta']['scheduler_options']
                ),
            )   
        ],
        run_dir=rundir,
        retries=6
    )
    return(cobalt_htex)

def polaris_localParslConfig(params):
    import parsl
    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname
    import os
    print(f'INFO - Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    job_name = params['job_name']
    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=4, # vs max_workers
                available_accelerators=4,
                worker_debug=True,
                cores_per_worker=32, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                address=address_by_hostname(),
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    nodes_per_block=params['num_of_full_nodes'],
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                    launcher=MpiExecLauncher(bind_cmd="--cpu-bind", overrides="--depth=64 --ppn 1"),
                    worker_init=params['worker_init']
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)

def theta_localParslConfig(params):

    import parsl
    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import MpiExecLauncher
    from parsl.addresses import address_by_hostname

    import os

    print(f'INFO - Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    job_name = params['job_name']
    #parsl.clear()

    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=8, # vs max_workers
                available_accelerators=8,
                worker_debug=True,
                cores_per_worker=32, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    nodes_per_block=params['num_of_full_nodes'],
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                    launcher=MpiExecLauncher(),
                    worker_init=params['worker_init']
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)

def beagle3_htParslConfig(params):

    import parsl
    from parsl.config import Config as ParslConfig
    from parsl.providers import SlurmProvider
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SrunLauncher
    import os

    print(f'INFO - Parsl version: {parsl.__version__}')
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')

    scheduler_options = [f"#SBATCH --gres=gpu:4", 
            "#SBATCH --partition=beagle3"
    ]
    scheduler_options = '\n'.join(scheduler_options)

    #SBATCH --constraint=rtx6000		# Only RTX 6000 (can set to v100 or a100)
    #SBATCH --cpus-per-task=1		# Number of threads
    #SBATCH --ntasks-per-node=1		# Number of CPU cores to drive GPU
    # ,"#SBATCH --constraint=rtx6000"

    user_opts = {
        'beagle3': {
            # Node setup: activate necessary conda environment and such.
            'worker_init': params['worker_init'],
            # ALCF allocation to use
            'account': params['account'],
        }
    }

    config = ParslConfig(
        executors=[
            HighThroughputExecutor(
                label="htex_slurm",
                available_accelerators=4,  # Pin each worker to a different GPU
                max_workers=4,
                provider=SlurmProvider(
                    launcher=SrunLauncher(),  # Ensures 1 manger per node, work on all 64 cores
                    account=user_opts['beagle3']['account'],
                    worker_init=user_opts['beagle3']['worker_init'],
                    walltime=params['walltime'],
                    init_blocks=params['init_blocks'],
                    scheduler_options=scheduler_options,
                    nodes_per_block=params['num_of_full_nodes'], 
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],  
                ),
            )
        ],
        run_dir=rundir,
        retries=6
    )
    return config

def beagle3_localParslConfig(params):

    import parsl
    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import SrunLauncher

    import os

    print(f'INFO - Parsl version: {parsl.__version__}')
    # I defined these locations otherwise parsl will use the current directory to output the run informations and log messages
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')
    # job_name = params['job_name']
    #parsl.clear()

    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=4, # vs max_workers
                available_accelerators=4,
                worker_debug=True,
                cores_per_worker=8, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    nodes_per_block=params['num_of_full_nodes'],
                    min_blocks=params['min_num_blocks'],
                    max_blocks=params['max_num_blocks'],
                    launcher=SrunLauncher(),
                    worker_init=params['worker_init']
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)