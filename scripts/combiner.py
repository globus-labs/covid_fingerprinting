import os
import glob
import pickle
import argparse
from operator import itemgetter

def combine(pkl_files, name, output=None, top_n_matches=100):

    all_results = {}
    for pkl_file in pkl_files:
        # print(f"Loading {pkl_file}")
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
            for target in data:
                # print(target, data[target][0])
                scores = all_results.get(target, [])
                scores.extend(data[target])
                sorted_scores = sorted(scores, key=itemgetter(1), reverse=True)[:top_n_matches]
                all_results[target] = sorted_scores

    
    with open(output, 'w') as csv_file_handle:
        for smile in all_results:
            sorted_scores = sorted(all_results[smile], key=itemgetter(1), reverse=True)
            for match, score, identifier in sorted_scores:
                print("{},{},{},{},{}".format(name,
                                              '%.6f'%score,
                                              smile,
                                              match,
                                              identifier), file=csv_file_handle)    
    
    print(f"Writing csv file to {output} done")


if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("--top_n_matches", default="100",
                        help="Top N matches to the target to calculate")
    parser.add_argument("-i", "--input_dir", default=".",
	                help="Inputs directory")
    parser.add_argument("-t", "--target_dir", help="Target dir")
    args = parser.parse_args()

    #outdir="top_100_similar_1000_targets"
    # outdir = "/home1/02551/yadunand/ScreenPilot/covid_fingerprinting/scripts/drugbank/all_targets_similar_10_matches"
    outdir = "test_out"
    os.makedirs(outdir, exist_ok=True)


    print(args.input_dir)
    # for source in ['enamine_diversity', 'ZINC15', 'SureChEMBL', 'pubchem', 'GDB13', 'Enamine_Real']:    
    for source in os.listdir(args.input_dir):
        #print("Source : ", source)
        targets = [f.strip('.csv') for f in glob.glob(args.target_dir)]
        #print("Target : ", targets)
        # for target in ['3CLPro_pocket1', 'ADRP_pocket1', 'CoV_pocket1', 'PLPro_pocket3']:
        #targets = [t.strip('ml.').rsplit('_', 1)[0] for t in os.listdir(args.target_dir)]
        for target in targets:
            base_target = os.path.basename(target)
            # print("Searching for : ", f"{args.input_dir}/{source}/*{base_target}*.pkl")
            files = glob.glob(f"{args.input_dir}/{source}/*{base_target}*.pkl")
            # print("Source files : ", files)
            output=f"{outdir}/{source}/{base_target}_50_targets_top_100_similar.csv"

            print(f"Starting {source} for target:{target} at ${output}") 
            
            os.makedirs(os.path.dirname(output), exist_ok=True)
            combine(files, source, 
                    output=f"{outdir}/{source}_{target}_50_targets_top_100_similar.csv",
                    top_n_matches=int(args.top_n_matches))
            #break
        #break
