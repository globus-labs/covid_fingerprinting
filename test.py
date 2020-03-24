# from candle_apps.candle_apps import run_intranode
# from candle_apps.candle_apps import ModelInferer
#from candle_apps.candle_node_local import run_local
#import pandas as pd
import os
import glob
import argparse
import traceback
from fingerprint.fp2 import process_files

import parsl
import pickle



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--num_files", default=10000,
                        help="Number of files to load and run. Default=all, if set to 0 the entire file will be used")
    parser.add_argument("-s", "--smile_dir", default=".",
                        help="File path to the smiles csv file")
    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    #for smile_dir in glob.glob(args.smile_dir):
    #    print(smile_dir)
    #exit(0)
    print(f"Loading pkl files from {args.smile_dir}")

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

    all_smile_dirs = glob.glob(args.smile_dir)
    counter = 0
    batch_futures = {}

    for smile_dir in all_smile_dirs:
        print("Processing smile_dir: {} {}/{}".format(smile_dir, counter, len(all_smile_dirs)))
        counter+=1
        batch_futures[smile_dir] = []

        outdir = "{}/{}".format(args.outdir, os.path.basename(smile_dir))
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(outdir + '/logs' , exist_ok=True)

        for smi_file in os.listdir(smile_dir)[:int(args.num_files)]:
            if not smi_file.endswith('.smi'):
                continue
        
            # print("Trying to launch : ", smi_file)
            fname = os.path.basename(smi_file)
            csv_file = "{}/{}".format(outdir, 
                                      fname.replace('.smi', '.csv'))
            log_file = "{}/logs/{}".format(outdir, 
                                           fname.replace('.smi', '.log'))

            if os.path.exists(csv_file):
                # Skip compute entirely if output file already exists
                continue

            smi_file_path = f"{smile_dir}/{smi_file}"

            # In this case the application expects a single file, we just give it a list with 
            # a single file
            x = process_files(smi_file_path,
                              '\t',
                              None,
                              csv_file,
                              log_file,
                              debug=True)
            print(x)
            """
            x = parsl_runner(smi_file_path,
                             '\t',
                             None,
                             csv_file,
                             log_file)

            batch_futures[smile_dir].append(x)
            """
        # Waiting for all futures
        print("Waiting for all futures from {}".format(smile_dir))

        for i in batch_futures[smile_dir]:
            try:
                x = i.result()

            except Exception as e:
                print("Exception : {} Traceback : {}".format(e, traceback.format_exc()))
                print(f"Chunk {i} failed")
        print(f"Completed {smile_dir}")    

    print("All done!")

