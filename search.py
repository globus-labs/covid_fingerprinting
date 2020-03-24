import rdkit
import random
import base64
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import RDLogger
from rdkit import DataStructs
import numpy as np
#from rdkit import DataStructs.cDataStructs.ExplicitBitVect
import pandas as pd
from pathlib import Path
import glob, sys
import argparse
import cProfile, pstats, io
from pstats import SortKey
import pickle
from operator import itemgetter
import matplotlib.pyplot as plt
import csv
import datetime


RDLogger.DisableLog('rdApp.*')


def search_files(files, fingerprints):
    with open('score.csv', 'w') as csvfile:
        score_writer = csv.writer(csvfile, delimiter=' ')
        for file in files:
            print('Processing file %s at %s'%(file, str(datetime.datetime.now())))
            smiles = pickle.load( open(file, 'rb') )
    
            # Precompute bitvectors
            fps = []
            for i, row in smiles.iterrows():
                try:
                    fps += [DataStructs.ExplicitBitVect(base64.b64decode(row['fingerprint']))]
                except:
                    fps += [None]
                    print('None')
                if i%100000==0 and i>0:
                    print('.', end='', flush=True)
                if i%1000000==0 and i>0:
                    if i%10000000==0:
                        print('M', end='', flush=True)
                    else:
                        print('m', end='', flush=True)

            smiles['fp'] = fps
            print('\n  Precomputed at %s '%(str(datetime.datetime.now())), end='')
    
            plt.rcParams["figure.figsize"] = (12,10)
            plt.figure()
            plt.title(file, fontsize=12)
            # For each of our target SMILE strings
            for (insmile, fp2) in fingerprints:
                print('\n  %s '%insmile, sep='')
                scores = []
                # For each row in the file we are comparing against
                for i, row in smiles.iterrows():
                    try:
                        score = DataStructs.TanimotoSimilarity(row['fp'], fp2)
                        scores += [(row['canonical_smile'], score)]
                    except:
                        pass
                    if i%100000==0 and i>0:
                        print('.', end='', flush=True)
                    if i%1000000==0 and i>0:
                        if i%10000000==0:
                            print('M', end='', flush=True)
                        else:
                            print('m', end='', flush=True)
                print(' %s'%(str(datetime.datetime.now())))
                # Add line to graph
                sorted_scores = sorted(scores, key=itemgetter(1))
                scores_only = [x[1] for x in sorted_scores]
                plt.step(np.arange(len(sorted_scores)), np.array(scores_only), label=insmile, linewidth=0.5)
                
                # Select top 200
                lastN = sorted_scores[-200:]
                lastN.reverse()
                for (smile, score) in lastN:
                    score_writer.writerow([file, '%.6f'%score, insmile, smile])
            plt.legend(fontsize=6)
            plt.savefig('fig%d.pdf'%(random.randint(1,100000)))
    
smiles1 = ['O=C(Nc1cccc(c1)S(=O)(=O)N1CCCCC1)CN1Cc2c(C1)cccc2'] # First element of pubchem, I think

smilesAll =\
       ['COC1CC2CC1CC2N(C)C1=NC[C@H]2OC(C(C)C)O[C@@H]2CN1', # First element of 15M
        'S=C1SC(=Cc2ccccn2)C(=O)N1c1ccccc1', # First element of enamine/01
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(/C=C/C4=CC=CC=C4)=O',
        'CC(C)(C)OC(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=C(F)C=C3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC(C)C)NC(/C=C/C3=CC=CC=C3)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CCCC)NC(/C=C/C3=CC=CC=C3)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC#C)NC(/C=C/C3=CC=CC=C3)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](C3CC3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3CCCCC3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3CC3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3CCC3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3CCCC3)NC(/C=C/C4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCCCC)=O)=O)C[C@@H]1CCNC1=O)[C@H](CC2=CC=CC=C2)NC(/C=C/C3=CC=CC=C3)=O',
        'O=C(N[C@H](C(C(NC(C)(C)C)=O)=O)C[C@@H]1CCNC1=O)[C@H](CC2=CC=CC=C2)NC(/C=C/C3=CC=CC=C3)=O',
        'N[C@H](C(C(NCC(C)C)=O)=O)C[C@@H]1CCNC1=O)[C@H](CC(C)C)NC(/C=C/C2=CC=CC=C2)=O',
        'O=C(N[C@H](C(C(NCC(OC)=O)=O)=O)C[C@@H]1CCNC1=O)[C@H](CC2=CC=CC=C2)NC(/C=C/C3=CC=CC=C3)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(CCC4=CC=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(C(C=C4)=CC5=C4OCCO5)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(C4=C5C(OCCO5)=CC=C4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(C4=CC(C=CC=C5)=C5O4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(C4=CC(C=CC=C5)=C5S4)=O',
        'O=C(N[C@H](C(C(NCC1=CC=CC=C1)=O)=O)C[C@@H]2CCNC2=O)[C@H](CC3=CC=CC=C3)NC(C4=CN(C=C(Br)C=C5)C5=N4)=O'
    ]

def compute_fingerprints_for_smiles(the_smiles):
    fingerprints = []
    for smile in the_smiles:
        try:
            #fingerprints += [(smile, AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048).ToBase64())]
            fingerprints += [(smile, AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048))]
        except:
            fingerprints += [(smile, 'None')]
    return(fingerprints)


def main(argv):
    parser = argparse.ArgumentParser(description='Program to process Moyer maps')
    parser.add_argument('-d', '--debug', help='Debug level', required=False, type=int, choices=[0, 1, 2])
    parser.add_argument('-P', '--profile', help='Profile code', required=False, default=False, action='store_true')
    #parser.add_argument('-i', '--input', help='Input file', required=True)
    args = parser.parse_args()

    # Default debug level is 1: set -d 0 for no info at all, -d 2 for verbose info
    if args.debug != None:
        global debug
        debug = args.debug

    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    # files = ['SMILEs/pubchem_canonical.csv']
    #files = glob.glob('SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only.smi')

    #files = ['SMILEs/Enamine_REAL_diversity_set_15.5M_fp.pkl']
    #smiles = smilesAll

    #files = ['test_fp.pkl', 'test2_fp.pkl', 'test3_fp.pkl' ]
    #smiles = smilesAll

    files = ['SMILEs/Enamine_REAL_diversity_set_15.5M_fp.pkl', 'SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only_fp.pkl']
    smiles = smilesAll
    
    fingerprints = compute_fingerprints_for_smiles(smiles)
    search_files(files, fingerprints)

    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

if __name__ == '__main__':
   main(sys.argv[1:])
