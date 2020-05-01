#!/bin/bash

#set -x
TOP_N=1000
LOGSUFFIX="V5.April_24_all_targets_top_100_similarity"
echo $LOGSUFFIX

TARGET_GLOB="./V5.April_24/agg_split_*"

for i in $(ls /scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints)
do 
    source=$i
    source_path="/scratch1/02551/yadunand/ScreenPilot/Fingerprints/fingerprints/$source/*pkl"

    if [[ -f "$PWD/$source/done" ]]
    then
	echo "$source is completed"
	continue
    fi
	
     python3 runner.py -s "$source_path" \
	-l $source_$LOGSUFFIX.log \
	-n $source \
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



