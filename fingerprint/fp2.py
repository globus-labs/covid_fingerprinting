from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import RDLogger
from pathlib import Path
import glob, sys
import argparse
RDLogger.DisableLog('rdApp.*')
import pickle

import logging

def set_file_logger(filename: str, name: str = 'candle', level: int = logging.DEBUG, format_string = None):
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

def process_files(smile_file, out_file, log_file, index_start, batchsize, pickle_out, debug=False):
    import time
    import shutil
    import os
    start = time.time()
    import logging

    from rdkit.Chem import AllChem

    count = 0
    bad_count = 0
    logger = set_file_logger(log_file, level=logging.DEBUG if debug else logging.INFO)
    logger.info(f"Running fingerprint on {smile_file} from index_start:{index_start}")

    tmp_out_file = '/dev/shm/{}'.format(os.path.basename(out_file))

    logger.info(f"Writing output temporarily to {tmp_out_file}")
    smiles = []
    with open(smile_file) as current:        
        current.seek(index_start)            
        smiles = [current.readline() for i in range(batchsize)]

    with open(tmp_out_file, 'w') as out_handle:
        if pickle_out:
            outputs = []
        for line in smiles:
            smile, *remainder = line.split()
            print('XXXX', line, smile, remainder)
            try:
                logger.debug(f"Processing smile {smile}")
                fprint = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048).ToBase64()
            except:
                logger.exception("Caught exception")
                fprint = None
                bad_count += 1
            id = remainder[0] if remainder else ''
            if pickle_out:
                outputs += [(smile, id, fprint)]
            else:
                print('{},{},{}'.format(smile, id, fprint), file=out_handle)
            count += 1

    if pickle_out:
        pickle.dump( outputs, tmp_out_file )
    shutil.move(tmp_out_file, out_file)

    logger.info("Bad smile count {}".format(bad_count))
    logger.info("Completed {} smiles from {} in {:8.3f}s".format(count,
                                                                 smile_file,
                                                                 time.time() - start))
    logger.handlers.pop()
    return out_file


def main(argv):
    parser = argparse.ArgumentParser(description='Program to process Moyer maps')
    parser.add_argument('-d', '--debug', help='Debug level', required=False, type=int, choices=[0, 1, 2])
    parser.add_argument('-P', '--profile', help='Profile code', required=False, default=False, action='store_true')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    args = parser.parse_args()

    # Default debug level is 1: set -d 0 for no info at all, -d 2 for verbose info
    if args.debug != None:
        global debug
        debug = args.debug

    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    if args.input == 'Enamine':
        files = glob.glob('SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only.smi')
        sep   = '\t'
        process_files(files, sep, None)
    elif args.input == 'pubchem':
        files = ['SMILEs/pubchem_canonical.csv']
        sep   = ','
        process_files(files, sep, 'pubchem')

    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

if __name__ == '__main__':
    process_files('/projects/CVD_Research/foster/SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_16_canonical_only.smi', 
                  '/projects/candle_aesp/yadu/covid_fingerprinting/test.csv', 
                  '/projects/candle_aesp/yadu/covid_fingerprinting/test.log', 
                  251,
                  5,
                  debug=True)
    #main(sys.argv[1:])
