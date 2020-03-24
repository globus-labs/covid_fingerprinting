import pandas as pd
import pickle

def pickle_file(from_file, sep, to_file):
    #smiles = pd.read_csv(from_file, sep=sep, header=None)
    smiles = pd.read_csv(from_file, sep=sep)
    smiles.columns = ['canonical_smile', 'id', 'fingerprint']
    pickle.dump( smiles, open( to_file, "wb" ) )

#pickle_file('SMILEs/Enamine_REAL_diversity_set_15.5M_fp.csv', '\t', 'SMILEs/Enamine_REAL_diversity_set_15.5M_fp.pkl')
pickle_file('test_fp.csv', '\t', 'test_fp.pkl')
#pickle_file('SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only_fp.csv', ',',
#            'SMILEs/Enamine_Real_SMILEs/2019q3-4_Enamine_REAL_01_canonical_only_fp.pkl')
