from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import RDLogger
from pathlib import Path
import glob, sys
import argparse
RDLogger.DisableLog('rdApp.*')

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

def process_files(smile_file, csv_file, log_file, debug=False):
    import time
    import shutil
    import os
    start = time.time()
    import logging

    from rdkit.Chem import AllChem

    count = 0
    logger = set_file_logger(log_file, level=logging.DEBUG if debug else logging.INFO)
    logger.info(f"Running fingerprint on {smile_file}")

    tmp_csv_file = '/dev/shm/{}'.format(os.path.basename(csv_file))

    logger.info(f"Writing output temporarily to {tmp_csv_file}")
    with open(tmp_csv_file, 'w') as csv_handle:
        with open(smile_file) as f:
            for line in f:
                clean_line = line.strip()
                smile, *remainder = line.split()
                try:
                    logger.debug(f"Processing smile {smile}")
                    fprint = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048).ToBase64()
                except:
                    logger.exception("Caught exception")
                    fprint = None
                print('{}, {}'.format(smile, fprint), file=csv_handle)
                count += 1

    shutil.move(tmp_csv_file, csv_file)

    logger.info("Completed {} smiles from {} in {:8.3f}s".format(count,
                                                                 smile_file,
                                                                 time.time() - start))
    logger.handlers.pop()


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
    process_files('/projects/candle_aesp/yadu/covid_fingerprinting/2019q3-4_Enamine_REAL_01.smi', 
                  '/projects/candle_aesp/yadu/covid_fingerprinting/test.csv', 
                  '/projects/candle_aesp/yadu/covid_fingerprinting/test.log', 
                  debug=True)
    #main(sys.argv[1:])
