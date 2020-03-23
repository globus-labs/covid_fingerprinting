from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import RDLogger
import pandas as pd
from pathlib import Path
import glob, sys
import argparse
import pickle
RDLogger.DisableLog('rdApp.*')
import cProfile, pstats, io
from pstats import SortKey

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
                fps += [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(row['canonical_smile']), 2, nBits=2048)]
            except:
                fps += ['None']
        smiles['ECFP4'] = fps
        pickle.dump( smiles, open( str(Path(file).with_suffix(''))+'_ecfp4.pkl', "wb" ) )


def main(argv):
    parser = argparse.ArgumentParser(description='Program to generate fingerprints for SMILEs')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-P', '--profile', help='Profile code', required=False, default=False, action='store_true')
    args = parser.parse_args()

    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    if args.input == 'Enamine':
        #files = glob.glob('SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only.smi')
        files = glob.glob('SMILEs/Enamine_Real_SMILEs/*.smi')
        process_files(files, '\t', None)
    else:
        print('Unknown input')

    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

if __name__ == '__main__':
   main(sys.argv[1:])
