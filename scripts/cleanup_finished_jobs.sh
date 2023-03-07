#!/bin/bash

##Script to clean up finished jobs due to issue with copying of the maltextract results

finishedJobs=($(find /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_outputs/Bacterial_Viral_Prescreening/No_Pathogen_Capture/*/multiqc/ -name "*multiqc_report.*"))
nextflowPath='/home/srv_pestis_screening/.conda/envs/nextflow/bin/'

for i in ${finishedJobs[@]}; do
    isInFile=0
    base=$(echo ${i} | awk 'BEGIN{FS=OFS="/";} {$NF=$(NF-1)=""; print $0}' | sed 's/\/$//')
    isInFile=$(grep -Fc "${base}" /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/debugging/finished_cleaned_runs.txt)
    if [[ ${isInFile} -eq 0 ]]; then
        DIR=${base}/maltextract/
        if [ -d "${DIR}" ]
            then
            if [ "$(ls -A ${DIR}/results)" ]; then
                echo "Take action $DIR/results is not Empty"
                cd ${base}
                ${nextflowPath}/nextflow clean -f -k -q
                echo ${base} >> /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/debugging/finished_cleaned_runs.txt
            else
                echo "$DIR is Empty"
	        fi
        else
	        echo "Directory $DIR not found."
        fi
    fi
done