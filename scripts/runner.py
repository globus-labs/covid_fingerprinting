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
from rdkit import Chem
from rdkit.Chem import AllChem

from targets import target_smiles

@parsl.python_app
def process_one_target(file, targets, N):
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-s", "--smile_glob", default=".",
                        help="Glob pattern that points to all .smi smile files")
    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    parser.add_argument("-p", "--pickle", action='store_true', default=False, 
                        help="Output data in pickled format. Default: False")
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


    if args.debug:
        parsl.set_stream_logger()


    os.makedirs(args.outdir, exist_ok=True)

    print("Computing fingerprints for all {} targets".format(len(target_smiles)))
    targets = {}
    bit_target = None
    for smile in target_smiles:
        bit_target = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048)
        targets[smile] = bit_target

    print("Targets computed")
    print(targets)

    all_pickle_files = glob.glob(args.smile_glob)
    batch_futures = {}
    counter = 0

    target_smile = 'O=C(Nc1cccc(c1)S(=O)(=O)N1CCCCC1)CN1Cc2c(C1)cccc2'

    for pickle_file in all_pickle_files:
        if not pickle_file.endswith('.pkl'):
            print(f"Ignoring {pickle_file} not smile file")
            continue

        print("Processing pickle_file: {} {}/{}".format(pickle_file, counter, len(all_pickle_files)))


        #outdir = "{}/{}".format(args.outdir, os.path.basename(pickle_file))
        #os.makedirs(outdir, exist_ok=True)
        #os.makedirs(outdir + '/logs' , exist_ok=True)


        x = process_one_target(pickle_file,
                               targets,
                               N = 100)
        batch_futures[pickle_file] = x
        counter += 1

    # Waiting for all futures
    print("Waiting for all futures from {}".format(pickle_file))

    all_results = {smile: [] for smile in targets}
    for pkl_file in batch_futures:
        i = batch_futures[pkl_file]
        try:
            x = i.result()
            for smile in x:
                all_results[smile].extend(x[smile])
            print(x)
        except Exception as e:
            print("Exception : {} Traceback : {}".format(e, traceback.format_exc()))
            print(f"Computing on {pkl_file} failed")
    print(f"Completed {pickle_file}")    

    print("All done!")

