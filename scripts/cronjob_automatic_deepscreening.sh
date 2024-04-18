#!/bin/bash

# Run by cron. 
# Will do the following:
# - Create tsv for any of new runs
# - Identify batches without any libraries in the tsv and add their names to List_empty_batches.tsv
# - Run script to initiate nf-core/eager in array mode containing all unrun tsv

##2023-02-2024: it will only automatically run jobs from 2023 until we catched up, in future all runs with new data

#Change to scripts directories
cd /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts

#Create tsv 
export R_LIBS_USER=/home/srv_pestis_screening/r_libs_screening/4.0/

/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/create_tsv_for_all_runs.sh Bacterial_Viral_Prescreening No_Pathogen_Capture 24 > /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/.eager_input_log 2>&1

#Identify batches without any libraries in the tsv
bash /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/identify_empty_tsv.sh

#Run nf-core/eager in array mode
/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/run_eager.sh -a -y 24
