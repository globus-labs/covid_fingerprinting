# target_smiles = ['O=C(Nc1cccc(c1)S(=O)(=O)N1CCCCC1)CN1Cc2c(C1)cccc2']
# target_smiles = ["O=C(NCCO)c1cccc(O)c1"]
target_smiles = None
with open('top.7.5k.ml.PLPro_pocket23_dock.csv') as f:
    smiles = f.readlines()[:1000]
    target_smiles = [smile.strip().split(',')[-1] for smile in smiles]
    
#print(target_smiles)
    
