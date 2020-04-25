#!/bin/bash

set -x
TOP_N=50
LOGSUFFIX="V3.April_9_all_targets_top_50_similarity"
echo $LOGSUFFIX

#TARGET_GLOB="./targets/*csv"
# TARGET_GLOB="./March_30/ml.ADRP-ADPR_pocket1_dock.csv.top10k.csv"
# TARGET_GLOB="adrp_aps_leads.csv"
# TARGET_GLOB="./DrugBank_sliced/*.csv"
TARGET_GLOB="./V3.April_9/*top50"

for i in $(ls /scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints)
do 
    source=$i
    source_path="/scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints/$source/*pkl"
    
    python3 runner.py -s "$source_path" \
	-l $source_$LOGSUFFIX.log \
	-n $source \
	-o $PWD/$source \
	--target_glob="$TARGET_GLOB" \
	--top_n_targets=$TOP_N \
	--top_n_matches=100 \
	-c frontera_small

    ls -thor $PWD/$source | mail -s "[Frontera] $source done" yadudoc1729@gmail.com
done

mail -s "[Frontera] All done" yadudoc1729@gmail.com

exit 0
python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/DrugBank/*pkl" \
    -l DrugBank_$LOGSUFFIX.log \
    -n DrugBank \
    -o $PWD/DrugBank \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/DrugBank | mail -s "DrugBank done" yadudoc1729@gmail.com


python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/Enamine_REAL_diversity_set_15.5M.smi/*pkl" \
    -l enamine_diversity_$LOGSUFFIX.log \
    -n enamine_diversity \
    -o $PWD/enamine_diversity \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/enamine_diversity | mail -s "Enamine_diversity done" yadudoc1729@gmail.com

python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/pubchem*/*pkl" \
    -l pubchem_$LOGSUFFIX.log \
    -n pubchem \
    -o $PWD/pubchem \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/pubchem | mail -s "Pubchem done" yadudoc1729@gmail.com


python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/SureChEMBL-64/*/*.pkl" \
    -l SureChEMBL_$LOGSUFFIX.log \
    -n SureChEMBL \
    -o $PWD/SureChEMBL \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/SureChEMBL | mail -s "SureChEMBL done" yadudoc1729@gmail.com

python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/GDB-13-64/*/*pkl" \
    -l GDB13_$LOGSUFFIX.log \
    -n GDB13 \
    -o $PWD/GDB13 \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/GDB13 | mail -s "GDB13 done" yadudoc1729@gmail.com

python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Enamine_Real_Fingerprints/*/*" \
    -l Enamine_Real_$LOGSUFFIX.log \
    -n Enamine_Real \
    -o $PWD/Enamine_Real \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/Enamine_Real | mail -s "Enamine_Real done" yadudoc1729@gmail.com

python3 runner.py -s "/home1/02551/yadunand/ScreenPilot/Fingerprints/ZINC15/*/*pkl" \
    -l ZINC15_$LOGSUFFIX.log \
    -n ZINC15 \
    -o $PWD/ZINC15 \
    --target_glob="$TARGET_GLOB" \
    --top_n_targets=$TOP_N \
    -c frontera_small

ls -thor $PWD/pubchem | mail -s "ZINC15 done" yadudoc1729@gmail.com

python3 combiner.py -t DrugBank_sliced/ --top_n_matches=10

ls -thor drugbank/*/ | mail -s "Combiner done" yadudoc1729@gmail.com

exit 0


for target_path in /home1/02551/yadunand/ScreenPilot/covid_fingerprinting/scripts/targets/*
#for target_path in /home1/02551/yadunand/ScreenPilot/covid_fingerprinting/scripts/targets/top.7.5k.ml.PLPro_pocket3_dock.csv
do
    echo $target_path
    target=$(basename $target_path)
    s1=${target%"_dock.csv"}
    PREFIX=${s1#"top.7.5k.ml."}

    LOGSUFFIX="${PREFIX}_100_targets_top_100_similarity"
    #echo "python3 run_only100.py -s /scratch1/02551/yadunand/ScreenPilot/Fingerprints/SureChEMBL-64/SureChEMBL_all_sorted_canonical.smi/SureChEMBL_all_sorted_canonical.chunk-0-1000000.pkl -l SureChEMBL_0-1M.$LOGSUFFIX.log -n SureChEMBL -c frontera -o SureChEMBL --target=$target_path"


    #python3 run_only100.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/SureChEMBL-64/*/*.pkl" -l SureChEMBL_$LOGSUFFIX.log -n SureChEMBL -c frontera_small -o SureChEMBL --target=$target_path
    #python3 run_only100.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/pubchem_canonical.smi/*pkl" -l pubchem_$LOGSUFFIX.log -n pubchem -o PubChem -c frontera_small --target=$target_path

    #python3 run_only100.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/GDB-13-64/*/*pkl"      -l GDB13_$LOGSUFFIX.log -n GDB13 -c frontera -o GDB13 --target=$target_path
    #python3 run_only100.py -s "/scratch1/02551/yadunand/ScreenPilot/Enamine_Real_Fingerprints/*/*" -l Enamine_Real_$LOGSUFFIX.log -n Enamine_Real -o Enamine_Real -c frontera --target=$target_path
    
    #python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/SureChEMBL-64/*/*.pkl" -l SureChEMBL_$LOGSUFFIX.log -n SureChEMBL -c frontera -o SureChEMBL --
    #python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/GDB-13-64/*/*pkl"      -l GDB13_$LOGSUFFIX.log -n GDB13 -c frontera
    #python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Fingerprints/pubchem_canonical.smi/*pkl" -l pubchem_$LOGSUFFIX.log -n pubchem -c frontera
    #python3 runner.py -s "/scratch1/02551/yadunand/ScreenPilot/Enamine_Real_Fingerprints/*/*" -l Enamine_Real_$LOGSUFFIX.log -n Enamine_Real -c frontera

done



