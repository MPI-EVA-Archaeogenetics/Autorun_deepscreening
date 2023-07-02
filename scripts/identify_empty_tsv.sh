#!/bin/bash

input=/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_inputs/Bacterial_Viral_Prescreening/No_Pathogen_Capture

outputFile=/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/debugging/List_empty_batches.tsv

allTSV=($(find ${input} -name "*_eager_input.tsv"))

for tsv in ${allTSV[@]}; do
    lineCount=`grep -v -c 'Sample' ${tsv}`
    if [[ ${lineCount} -le 0 ]]
    then
        #Check if already present in ../debugging/List_empty_batches.tsv
        countDebug=`grep -c "${tsv}" ${outputFile}` 
        #echo ${countDebug}
        if [[ ${countDebug} -eq 0 ]]
        then
            echo "${tsv} is empty, will be added to ${outputFile}"
            echo $tsv >> ${outputFile}
        fi
    fi
done