from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import RDLogger
import pandas as pd
from pathlib import Path
import glob, sys
import argparse
RDLogger.DisableLog('rdApp.*')


def process_files(files, sep, headers):
    for file in files:
        if headers == None:
            smiles = pd.read_csv(file, sep=sep, header=None)
            smiles.columns = ['canonical_smile', 'id']
        else:
            smiles = pd.read_csv(file, sep=sep)
        fps = []
        for i, row in smiles.iterrows():
            try:
                fps += [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(row['canonical_smile']), 2, nBits=2048).ToBase64()]
            except:
                fps += ['None']
        smiles['ECFP4'] = fps
        smiles.to_csv(str(Path(file).with_suffix(''))+'_fp.csv', index=False)


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
   main(sys.argv[1:])
