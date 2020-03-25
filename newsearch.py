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

debug = 0

RDLogger.DisableLog('rdApp.*')


def search_fingerprint_set(fingerprint_set, target, N):
    scores = []

    # Find scores for non-None fingerprints in fingerprint set
    for (smile, fingerprint) in fingerprint_set:
        try:
            score = DataStructs.TanimotoSimilarity(fingerprint, target)
            scores += [(smile, score)]
        except:
            pass

    sorted_scores = sorted(scores, key=itemgetter(1))

    return(sorted_scores[-N:])

# Two arguments are sorted lists, in increasing order; so is output, with largest values
def merge_lists(list1, list2, N):
    # Say best so far is: 12 14 16 18
    #     new list is   :  1 10 13 15
    # Updated list      : 14 15 16 18
    index1 = N-1
    index2 = N-1
    merged = []
    for i in range(N):
        element1 = list1[index1]
        element2 = list2[index2]
        if element1[1] < element2[1]:
            merged += [element2]
            index2 -= 1
        else:
            merged += [element1]
            index1 -= 1
    merged.reverse()     
    return merged

def get_fingerprint(smile):
    return( AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048) )

def print_list(list):
    out = ''
    for (key, score) in list:
        out += '    %0.6f %s\n'%(score,key)
    return(out)

def process_one_target(files, target, N):
    try:
        bit_target = get_fingerprint(target)
    except:
        return None
    best_so_far = [('', 0.0) for index in range(N)] 
    for file in files:
        with open(file, newline='') as f:
            reader = csv.reader(f)
            read = [(sm, fp) for [sm, id, fp] in reader]
            fingerprint_set = []
            for sm, fp in read:
                try:
                    bv = DataStructs.ExplicitBitVect(base64.b64decode(fp))
                except:
                    if debug > 0:
                        print('    fingerprint fail:', sm)
                    bv = None
                fingerprint_set += [(sm, bv)]
            #fingerprint_set = [(sm, DataStructs.ExplicitBitVect(base64.b64decode(fp))) for [sm, id, fp] in reader]
            new_list        = search_fingerprint_set(fingerprint_set, bit_target, N)
            if debug > 1:
                print('Merge:')
                print('  Pre:\n%s'%print_list(best_so_far))
                print('  New:\n%s'%print_list(new_list))
            best_so_far     = merge_lists(best_so_far, new_list, N)
            if debug > 1:
                print('  Now:\n%s'%print_list(best_so_far))
    return(best_so_far)


smiles1 = ['O=C(Nc1cccc(c1)S(=O)(=O)N1CCCCC1)CN1Cc2c(C1)cccc2'] # First element of pubchem, I think

smilesAll =\
       ['COC1CC2CC1CC2N(C)C1=NC[C@H]2OC(C(C)C)O[C@@H]2CN1', # First element of 15M
        'S=C1SC(=Cc2ccccn2)C(=O)N1c1ccccc1', # First element of enamine/01
        'O[C@H]1CCCCN(C1)S(=O)(=O)c1cccc(c1)NC(=O)C', #element in Zinc15
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

def plot_figure(label, target, results, N, count):
    plt.rcParams["figure.figsize"] = (8,6)
    plt.figure()
    plt.title('Top %d matches in %s for\n%s'%(N, label, target), fontsize=8)
    plt.ylim(0,1)
    scores_only = [x[1] for x in results]
    plt.step(np.arange(len(scores_only)), np.array(scores_only), linewidth=1)
    alphas = ''.join(c for c in target if c.isalpha() or c.isdigit() or c=='=' or c=='@')
    figure_name = 'fig_%s__%s__%d.pdf'%(label, alphas[0:30], count) # random.randint(1,100))
    print('  creating %s'%figure_name)
    plt.savefig(figure_name)
    plt.close()


def process_targets_on_files(label, files, smiles, N, figures):
    num_smiles = len(smiles)
    if files == []:
        print('No files for',label)
        exit(1)
    with open(label+'_output.csv', 'w') as csvfile:
        result_writer = csv.writer(csvfile, delimiter=' ')
        count = 0
        for target in smiles:
            print('%s (%d of %d) : %s'%(label, count, num_smiles, target))
            results = process_one_target(files, target, N)
            if results == None:
                print('  Bad target')
            else:
                if figures:
                    plot_figure(label, target, results, N, count)
                # Put higher results first in output file
                results.reverse()
                # Note that we don't record, or write, the file from which the match came
                print('  top 10:')
                for i in range(10):
                    print('    %0.6f : %s'%(results[i][1], results[i][0]))
                for (match, score) in results:
                    result_writer.writerow([label, '%.6f'%score, target, match])
            count += 1

def main(argv):
    parser = argparse.ArgumentParser(description='Program to process Moyer maps')
    parser.add_argument('-d', '--debug', help='Debug level', required=False, type=int, choices=[0, 1, 2])
    parser.add_argument('-P', '--profile', help='Profile code', required=False, default=False, action='store_true')
    parser.add_argument('-n', '--number', help='Number of candidates to select', required=False, type=int, default=1000)
    parser.add_argument('-f', '--figures', help='Whether to create figures', required=False, default=False, action='store_true')
    #parser.add_argument('-i', '--input', help='Input file', required=True)
    args = parser.parse_args()

    # Default debug level is 1: set -d 0 for no info at all, -d 2 for verbose info
    if args.debug != None:
        global debug
        debug = args.debug

    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    pubchem_label  = 'pubchem'
    pubchem_files  = glob.glob('outputs/pubchem_canonical.smi/*.csv')
    pubchem_smiles = smilesAll

    enamine_label  = 'enamine'
    enamine_files  = glob.glob('outputs/2019q3-4_Enamine*/*.smi')
    enamine_smiles = smilesAll

    zinc_15_label  = 'zinc_15'
    zinc_15_files  = glob.glob('outputs/zinc15_unique_smile_name*/*.csv')
    zinc_15_smiles = smilesAll

    #label  = 'test'
    #files  = ['f1.csv', 'f2.csv', 'f3.csv']
    #smiles = smilesAll

    process_targets_on_files(pubchem_label, pubchem_files, pubchem_smiles, args.number, args.figures)
    #process_targets_on_files(enamine_label, enamine_files, enamine_smiles, args.number, args.figures)
    #process_targets_on_files(zinc_15_label, zinc_15_files, zinc_15_smiles, args.number, args.figures)

    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

if __name__ == '__main__':
   main(sys.argv[1:])
