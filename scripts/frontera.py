from parsl.config import Config
from parsl.channels import LocalChannel
from parsl.providers import SlurmProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.addresses import address_by_hostname


""" This config assumes that it is used to launch parsl tasks from the login nodes
of Frontera at TACC. Each job submitted to the scheduler will request 2 nodes for 10 minutes.
"""
config = Config(
    retries=1,
    executors=[
        HighThroughputExecutor(
            label="frontera_htex",
            prefetch_capacity=1,
            address=address_by_hostname(),
            max_workers=30,          # Set number of workers per node
            provider=SlurmProvider(
                cmd_timeout=120,     # Add extra time for slow scheduler responses
                nodes_per_block=300,
                walltime='04:00:00',
                # partition='development',                                 # Replace with partition name
                partition='normal',

                init_blocks=2,
                min_blocks=1,
                max_blocks=2,

                scheduler_options='''#SBATCH -A FTA-Jha''',   # Enter scheduler_options if needed

                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='source ~/anaconda3/bin/activate; conda activate candle_py3.7',
                launcher=SrunLauncher(),
            ),
        )
    ],
    strategy=None,
)
