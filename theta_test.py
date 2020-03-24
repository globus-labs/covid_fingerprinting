from parsl.config import Config
from parsl.providers import CobaltProvider
from parsl.launchers import AprunLauncher
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname
#from parsl.monitoring.monitoring import MonitoringHub

config = Config(
    executors=[
        HighThroughputExecutor(
            label='theta_local_htex_multinode',

            # this is max_workers per node. Controls node level concurrency
            max_workers=64, # The target process itself if a Multiprocessing application. We do not
            # need to overload the compute node with parsl workers.
            address="10.236.1.195",
            prefetch_capacity=2,
            provider=CobaltProvider(
                # Switch between debug-flat-quad and debug-cache-quad based on queue busy-ness
                queue='debug-flat-quad',
                #queue='default',
                #queue='CVD_Research',
                #account='candle_aesp',
                account='CVD_Research',
                launcher=AprunLauncher(overrides=" -d 64"),
                walltime='01:00:00',
                nodes_per_block=8,
                init_blocks=1,
                min_blocks=1,
                max_blocks=1,
                # Make sure cond
                worker_init='source ~/anaconda3/bin/activate; conda activate candle_py3.7;',
                cmd_timeout=300,
            ),
        )
    ],
    strategy=None,
)

