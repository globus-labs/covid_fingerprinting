#!/bin/bash

#set -x
TOP_N=1000
LOGSUFFIX="hong_MCU_top_100_similarity"
echo $LOGSUFFIX

# TARGET_GLOB="./V5.April_24/agg_split_*"
# TARGET_GLOB="./V5.April_24/agg_split_*"
#TARGET_GLOB="./top.7.5k.ml.PLPro_pocket23_dock.csv"
TARGET_GLOB="./Revised_Screen_Hit_List_Nov20.csv"
# for i in $(ls /scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints)
for i in MCU
do 
    source=$i
    dataset_path="/scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints/$source/*csv"

    if [[ -f "$PWD/$source/done" ]]
    then
	echo "$source is completed"
	continue
    fi
	
     python3 runner.py \
	-l $source_$LOGSUFFIX.log \
	--dataset_glob="$dataset_path" \
	-o $PWD/$source \
	--target_glob="$TARGET_GLOB" \
	--top_n_targets=$TOP_N \
	--top_n_matches=100 \
	-c frontera_small

    status="$?"
    if [[ "$status" == "0" ]]
    then
        touch "$PWD/$source/done"
    fi
    ls -thor $PWD/$source | mail -s "[Frontera] $source done with $status" yadudoc1729@gmail.com

done

#mail -s "[Frontera] All done" yadudoc1729@gmail.com

exit 0



