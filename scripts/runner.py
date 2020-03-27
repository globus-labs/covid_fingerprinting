# from candle_apps.candle_apps import run_intranode
# from candle_apps.candle_apps import ModelInferer
#from candle_apps.candle_node_local import run_local
#import pandas as pd
import os
import glob
import argparse
import traceback
import parsl
import pickle
import gc
from itertools import islice
from rdkit import Chem
from rdkit.Chem import AllChem
from operator import itemgetter
from targets import target_smiles
import logging

def set_file_logger(filename: str, name: str = 'runner', level: int = logging.DEBUG, format_string = None):
    """Add a stream log handler.

    Args:
        - filename (string): Name of the file to write logs to
        - name (string): Logger name
        - level (logging.LEVEL): Set the logging level.
        - format_string (string): Set the format string

    Returns:
       -  None
    """
    if format_string is None:
        format_string = "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger




@parsl.python_app
def process_old_target(file, targets, N):
    import rdkit
    import base64
    from rdkit import RDLogger
    from rdkit import DataStructs
    from pstats import SortKey
    import pickle
    from operator import itemgetter
    import csv
    

    target_results = {}

    for smile_target in targets:
        best_so_far = [('', 0.0) for index in range(N)]
        bit_target = targets[smile_target]

        if file.endswith('pkl'):
            read = pickle.load( open(file, 'rb') )
        else:
            with open(file, 'r') as f:
                reader = csv.reader(f)
                read = list(map(tuple, reader))

        fingerprint_set = []
        for sm, _, fp in read:
            # Next TRY here as some do not convert
            try:
                bv = DataStructs.ExplicitBitVect(base64.b64decode(fp))
            except:
                bv = None
            fingerprint_set += [(sm, bv)]
            
        # Find scores for non-None fingerprints in fingerprint set
        scores = []
        for (smile, fingerprint) in fingerprint_set:
            try:
                score = DataStructs.TanimotoSimilarity(fingerprint, bit_target)
                scores += [(smile, score)]
            except:
                pass

        sorted_scores = sorted(scores, key=itemgetter(1))
        new_list = sorted_scores[-N:]
        target_results[smile_target] = new_list

    return(target_results)

@parsl.python_app
def process_one_target(file, targets, N, outfile=None):
    import rdkit
    import base64
    from rdkit import RDLogger
    from rdkit import DataStructs
    from pstats import SortKey
    import pickle
    from operator import itemgetter
    import csv
    

    target_results = {}

    if file.endswith('pkl'):
        read = pickle.load( open(file, 'rb') )
    else:
        with open(file, 'r') as f:
            reader = csv.reader(f)
            read = list(map(tuple, reader))
    fingerprint_set = []
    for sm, _, fp in read:
        # Next TRY here as some do not convert
        try:
            bv = DataStructs.ExplicitBitVect(base64.b64decode(fp))
        except:
            bv = None
        fingerprint_set += [(sm, bv)]

    for smile_target in targets:
        best_so_far = [('', 0.0) for index in range(N)]
        bit_target = targets[smile_target]
            
        # Find scores for non-None fingerprints in fingerprint set
        scores = []
        for (smile, fingerprint) in fingerprint_set:
            try:
                score = DataStructs.TanimotoSimilarity(fingerprint, bit_target)
                scores += [(smile, score)]
            except:
                pass

        sorted_scores = sorted(scores, key=itemgetter(1))
        new_list = sorted_scores[-N:]
        target_results[smile_target] = new_list

    if outfile:
        with open(outfile, 'wb') as f:
            pickle.dump(target_results, f)
        return outfile
    else:
        return(target_results)


def launch_slice(all_pickle_files, target_subset, csv_file_handle, index):

    target_futures = []

    count  = 0
    for pickle_file in all_pickle_files:
        count += 1 
        if not pickle_file.endswith('.pkl'):
            logger.warning(f"Ignoring {pickle_file} not smile file")
            continue

        logger.info("Processing pickle_file: {} Target:{}-{}".format(pickle_file, 
                                                                     index,
                                                                     index+100))
        outfile =  "{}/{}.{}-{}.pkl".format(args.outdir, 
                                            os.path.basename(pickle_file),
                                            index, 
                                            index+100)
                                  
        x = process_one_target(pickle_file,
                               target_subset,
                               N = int(args.top_n),
                               outfile=outfile
                           )
        target_futures.append(x)

        # Wait for all futures

        logger.info(f"Waiting for all futures from target_subset {index}-{index+100}")

        all_results = {smile: [] for smile in targets}

    count = 0
    for fu in target_futures:
        logger.debug("Waiting on {}/{} of futures".format(count, len(target_futures)))
        count+=1
        try:
            x = fu.result()
            
            #for smile in x:
            #    all_results[smile].extend(x[smile])
                
        except Exception as e:
            logger.exception(f"Computing failed")
            print("Exception : {} Traceback : {}".format(e, traceback.format_exc()))
    print("{} - {} done".format(i, i+100))

    """
    for smile in all_results:
        logger.info(f"Collating result for target : {smile}")
        sorted_scores = sorted(all_results[smile], key=itemgetter(1), reverse=True)[:int(args.top_n)]
        logger.info("Top {} : {}".format(args.top_n, sorted_scores))
        for match, score,  in sorted_scores:
            print("{},{},{},{}".format(args.name,
                                       '%.6f'%score,
                                       smile,
                                       match), file=csv_file_handle)
        csv_file_handle.flush()
    """
    return [fu.result() for fu in target_futures]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--name", required=True,
                        help="Name to use in csv")
    parser.add_argument("-s", "--smile_glob", default=".",
                        help="Glob pattern that points to all .smi smile files")
    parser.add_argument("-t", "--top_n", default="100",
                        help="Top N items to report, default  = 100")
    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    parser.add_argument("-l", "--logfile", default="runner.log",
                        help="log file path ")
    parser.add_argument("-p", "--pickle", action='store_true', default=False, 
                        help="Output data in pickled format. Default: False")
    args = parser.parse_args()

    #for smile_dir in glob.glob(args.smile_dir):
    #    print(smile_dir)
    #exit(0)

    logger = set_file_logger(args.logfile)

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 1
    elif args.config == "theta":
        from theta import config
    elif args.config == "frontera":
        from frontera import config
    elif args.config == "theta_test":
        from theta_test import config
    elif args.config == "comet":
        from comet import config

    # Most of the app that hit the timeout will complete if retried.
    # but for this demo, I'm not setting retries.
    # config.retries = 2
    parsl.load(config)


    if args.debug:
        parsl.set_stream_logger()


    os.makedirs(args.outdir, exist_ok=True)

    logger.info("Computing fingerprints for all {} targets".format(len(target_smiles)))
    targets = {}
    bit_target = None
    for smile in target_smiles:
        bit_target = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048)
        targets[smile] = bit_target

    logger.info("Targets computed")

    all_pickle_files = glob.glob(args.smile_glob)
    batch_futures = {}
    counter = 0

    csv_file_handle = open(args.logfile.replace('.log', '.csv'), 'w')
        
    """
    all_outfiles = []
    for i in range(0,len(targets),100):
    
        target_subset = {key:targets[key]  for key in islice(targets.keys(), i, i+100)}
        target_futures = []
        outfiles = launch_slice(all_pickle_files, target_subset, csv_file_handle, i)    
        all_outfiles.extend(outfiles)
        print("Garbage collect : ", gc.collect())
        #break
        
    """

    outfiles = launch_slice(all_pickle_files, targets, csv_file_handle, 1000)    
    all_outfiles.extend(outfiles)
    print("Garbage collect : ", gc.collect())
    #break
                                               
    print("All done!")

