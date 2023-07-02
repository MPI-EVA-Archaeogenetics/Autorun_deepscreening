#!/bin/bash

analysis=$1
type=$2
year=$3

re='^[0-9]+$'
if [[ ${year} =~ ${re} ]]
then
    allRuns=(`ls /mnt/archgen/Autorun/Results/Bacterial_Viral_Prescreening/ | grep "^${year}"`)
else
    allRuns=(`ls /mnt/archgen/Autorun/Results/Bacterial_Viral_Prescreening/`)
fi


for i in ${allRuns[@]}; do
    echo "Checking ${i}"
    tsv=`ls /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_inputs/Bacterial_Viral_Prescreening/No_Pathogen_Capture/${i}*`
    if [ -z ${tsv} ]
    then
        Rscript /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/prepare_eager_deepScreening_tsv.R \
        -s ${i} \
        -a ${analysis} \
        -o /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_inputs \
        -t ${type} \
        -d \
        /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/.credentialsCluster
    else
        echo "${tsv} already exists, checking if ${analysis} has been rerun in Autorun"
        dateTsv=`date +%F -r ${tsv}`
        numberModifiedFiles=`find /mnt/archgen/Autorun/Results/Bacterial_Viral_Prescreening/${i} -type f -newermt ${dateTsv} | wc -l`
        if [ ${numberModifiedFiles} -eq 0 ]
        then
            echo "${tsv} is up to date, my job here is done!"
        else
            Rscript /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/prepare_eager_deepScreening_tsv.R \
            -s ${i} \
            -a ${analysis} \
            -o /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_inputs \
            -t ${type} \
            -d \
            /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/.credentialsCluster
            echo "${tsv} is updated"
        fi
    fi
done
