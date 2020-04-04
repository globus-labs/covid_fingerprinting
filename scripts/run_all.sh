
TOP_N=100
LOGSUFFIX="March_30_${TOP_N}_targets_top_100_similarity"
echo $LOGSUFFIX

#TARGET_GLOB="./targets/*csv"
# TARGET_GLOB="./March_30/ml.ADRP-ADPR_pocket1_dock.csv.top10k.csv"
# TARGET_GLOB="adrp_aps_leads.csv"
TARGET_GLOB="./March_30/*.csv"

cat <<EOF
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
    -c frontera

ls -thor $PWD/GDB13 | mail -s "GDB13 done" yadudoc1729@gmail.com
EOF

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



