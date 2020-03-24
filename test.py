# from candle_apps.candle_apps import run_intranode
# from candle_apps.candle_apps import ModelInferer
#from candle_apps.candle_node_local import run_local
#import pandas as pd
import os
import glob
import argparse
import traceback
# from fingerprint.compute_fingerprints import process_files
from fingerprint.fp2 import process_files

import parsl
import pickle

def generate_batch(filename, start=0, batchsize=10, max_batches=10):

    counter = 0
    if max_batches == 0:
        max_batches = 999999999

    x = 'Hello'
    with open(filename) as current:
        yield current.tell()
        counter += 1

        while x and counter < max_batches:
            counter += 1
            for i in range(batchsize):
                x = current.readline()

            yield current.tell()
        return



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--num_batches", default=0,
                        help="Number of batches to load and run. Default=8, if set to 0 the entire file will be used")
    parser.add_argument("-b", "--batch_size", default="4",
                        help="Size of the batch of smiles to send to each node for processing. Default=4, should be 10K")
    parser.add_argument("-s", "--smile_glob", default=".",
                        help="Glob pattern that points to all .smi smile files")
    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    #for smile_dir in glob.glob(args.smile_dir):
    #    print(smile_dir)
    #exit(0)

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 1
    elif args.config == "theta":
        from theta import config
    elif args.config == "theta_test":
        from theta_test import config
    elif args.config == "comet":
        from comet import config

    # Most of the app that hit the timeout will complete if retried.
    # but for this demo, I'm not setting retries.
    # config.retries = 2
    parsl.load(config)

    parsl_runner = parsl.python_app(process_files)

    if args.debug:
        parsl.set_stream_logger()


    os.makedirs(args.outdir, exist_ok=True)

    all_smile_files = glob.glob(args.smile_glob)
    counter = 0
    batch_futures = {}
    chunksize = int(args.batch_size)

    for smile_file in all_smile_files:
        if not smile_file.endswith('smi'):
            print(f"Ignoring {smile_file} not smile file")
            continue

        print("Processing smile_file: {} {}/{}".format(smile_file, counter, len(all_smile_files)))
        counter+=1
        batch_futures[smile_file] = []

        outdir = "{}/{}".format(args.outdir, os.path.basename(smile_file))
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(outdir + '/logs' , exist_ok=True)

        batch_generator = generate_batch(smile_file, start=0,
                                         batchsize=int(args.batch_size),
                                         max_batches=int(args.num_batches))

        i = 0
        for batch_index in batch_generator:        
            print(f"Trying to launch {smile_file} index {i}-{i+chunksize} index:{batch_index}")
            fname = os.path.basename(smile_file)
            csv_file = "{}/{}".format(outdir, 
                                      fname.replace('.smi', f'.chunk-{i}-{i+chunksize}.tsv'))
            log_file = "{}/logs/{}".format(outdir, 
                                           fname.replace('.smi', f'.chunk-{i}-{i+chunksize}.log'))

            if os.path.exists(csv_file):
                # Skip compute entirely if output file already exists
                continue

            # In this case the application expects a single file, we just give it a list with 
            # a single file
            x = parsl_runner(smile_file,
                             csv_file,
                             log_file,
                             batch_index,
                             chunksize,
                             debug=False)
            batch_futures[smile_file].append(x)
            i += chunksize

        # Waiting for all futures
        print("Waiting for all futures from {}".format(smile_file))

        for i in batch_futures[smile_file]:
            try:
                x = i.result()
                print(x)
            except Exception as e:
                print("Exception : {} Traceback : {}".format(e, traceback.format_exc()))
                print(f"Chunk {i} failed")
        print(f"Completed {smile_file}")    

    print("All done!")

