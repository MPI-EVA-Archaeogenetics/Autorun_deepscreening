#!/bin/bash

. /etc/profile

# Help message
Help()
{
   # Display Help
   echo "Starts nf-core/eager for deep screening."
   echo
   echo "Syntax: run_eager.sh [-a|h|y]"
   echo "options:"
   echo "a     Runs in array mode."
   echo "h     Print this Help."
   echo "y     Creates queue for specific year."
   echo
}


## Defaults to create arrays and temp_file that will contain the queue
array=''
temp_file=''
year=''

# Get the options
while getopts ":hay:" option; do
   case $option in
        h) # display Help
            Help
            exit;;
        a) # Activate array option
            array='TRUE';;
        y) # Activate array option
            year=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

if [[ ${array} == 'TRUE' ]]; then
    echo "Array mode on"
    if [[ -z ${year} ]]; then
        temp_file="/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/queue/$(date +'%y%m%d_%H%M')_Autorun_deepscreening_queue.txt"
    else
        temp_file="/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/queue/$(date +'%y%m%d_%H%M')_${year}_Autorun_deepscreening_queue.txt"
    fi
    echo ${temp_file}
fi

##Path to nextflow in cluster
nextflowPath='/home/srv_pestis_screening/.conda/envs/nextflow/bin/'
## Path to input directory
mainInDir='/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_inputs'
## Path to output directory
mainOutDir='/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_outputs'
##Fix eager version to 2.4.6 to use the correct malt version
eager_version='2.4.7'
##Path to config file for autorun deep screening
autorun_config='/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/conf/Autorun_deepscreening.config'

## Set base profiles for EVA cluster.
nextflow_profiles="eva,archgen,autorun,local_paths"

##Create a loop to pick up all unrun eager_inputs
##loop through the different types of analysis
for analysis_type in "Bacterial_Viral_Prescreening"; do
    ##Set up the profiles for the specific analysis type
    analysis_profiles="${nextflow_profiles},${analysis_type}"
    ##Loop through No_Pathogen_Capture & Pathogen_Capture, 
    ## right now Pathogen_Capture won't be included
    for capture_type in "No_Pathogen_Capture"; do
        ##Loop thorough all the TSV and set up the outdir
        tsv_inputs=$(find ${mainInDir}/${analysis_type}/${capture_type} -name "${year}*.tsv")
        for eager_input in ${tsv_inputs}; do
            line_count=$(wc -l ${eager_input} | awk '{print $1}')
            if [[ $line_count == 1 ]]; then
                echo "Empty ${eager_input}, nf-core/eager not executed"
                echo ${eager_input} >> ../debugging/List_empty_batches_nodeepscreening_runned.tsv
            else
                ##Extract sequencing batch id from the eager_input.tsv
                sequencing_batch_id=$(basename $eager_input _eager_input.tsv)
                ##Set the output directory for each of the sequencing batch id
                eager_output_dir="${mainOutDir}/${analysis_type}/${capture_type}/${sequencing_batch_id}"
                ##Ideally run name should start with p followed by the sequencing_batch_id
                ##But right not not easy to resume, so setting name to -resume
                if [[ -d ${eager_output_dir}/work  ]]; then
                    run_name="-resume"
                else
                    run_name="-name p${sequencing_batch_id}"
                fi
                mkdir -p ${eager_output_dir}
                ##Start eager run if either the multiqc_report does not exist, or if TSV newer than the report 
                ##(useful if new data due to wrong index combinations in demultiplexing).
                ##This is check with -nt in the if statement
                multiqc_file=$(find ${eager_output_dir}/multiqc -name "*multiqc_report.html")
                if [[ ${eager_input} -nt ${multiqc_file} ]]; then
                    if [[ ${array} == 'TRUE' ]]; then
                        ## For array submissions, the commands to be run will be added one by one to the temp_file
                        ## Then once all jobs have been added, submit that to qsub with each line being its own job.
                        ## Use `continue` to avoid running eager interactivetly for arrayed jobs.
                        ##Spawner to be run within the parent of the output directory 
                        echo "cd ${eager_output_dir} ; ${nextflowPath}/nextflow run nf-core/eager \
                        -r ${eager_version} \
                        -c ${autorun_config} \
                        -profile ${analysis_profiles} \
                        --input ${eager_input} \
                        --outdir ${eager_output_dir} \
                        -w ${eager_output_dir}/work \
                        ${run_name}" | tr -s " " >> ${temp_file}
                        continue ## Skip running eager interactively if arrays are requested
                    fi

                    ## NON-ARRAY RUNS
                    ## Change to input directory to run from, to keep one cwd per run.
                    cd ${eager_output_dir}
                    ## Debugging info.
                    echo "Running eager on ${eager_input}:"
                    echo "${nextflowPath}/nextflow run nf-core/eager \
                    -r ${eager_version} \
                    -c ${autorun_config} \
                    -profile ${analysis_profiles} \
                    --input ${eager_input} \
                    --outdir ${eager_output_dir} \
                    -w ${eager_output_dir}/work \
                    ${run_name}"
            

                    ## Actually run eager now.
                    ## Monitor run in nf tower. Only works if TOWER_ACCESS_TOKEN is set.
                    ## Runs show in the Autorun_Eager workspace on tower.nf
                    ${nextflowPath}/nextflow run nf-core/eager \
                    -r ${eager_version} \
                    -c ${autorun_config} \
                    -profile ${analysis_profiles} \
                    --input ${eager_input} \
                    --outdir ${eager_output_dir} \
                    -w ${eager_output_dir}/work \
                    ${run_name}
            
                    cd ${mainInDir} ## Then back to root dir
                fi
            fi
        done
    done
done

## If array is requested submit the created array file to qsub below
if [[ ${array} == 'TRUE' ]]; then
    mkdir -p /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/array_Logs/$(basename ${temp_file}) ## Create new directory for the logs for more traversable structure
    jn=$(wc -l ${temp_file} | cut -f 1 -d " ") ## number of jobs equals number of lines
    export NXF_OPTS='-Xms4G -Xmx4G' ## Set 4GB limit to Nextflow VM
    export JAVA_OPTS='-Xms8G -Xmx8G' ## Set 8GB limit to Java VM
    ## -V Pass environment to job (includes nxf/java opts)
    ## -S /bin/bash Use bash
    ## -l v_hmem=40G ## 30GB memory limit (8 for java + the rest for garbage collector)
    ## -pe smp 2 ## Use two cores. one for nextflow, one for garbage collector
    ## -n AE_spawner ## Name the job
    ## -cwd Run in currect run directory (ran commands include a cd anyway, but to find the files at least)
    ## -j y ## join stderr and stdout into one output log file
    ## -b y ## Provided command is a binary already (i.e. executable)
    ## -o /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/array_Logs/ ## Keep all log files in one directory.
    ## -tc 2 ## Number of concurrent spawner jobs (10)
    ## -t 1-${jn} ## The number of array jobs (from 1 to $jn)
    echo "qsub -V -S /bin/bash -l h_vmem=40G -pe smp 2 -N AD_spawner_$(basename ${temp_file}) -cwd -j y -b y -o /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/array_Logs/$(basename ${temp_file}) -tc 2 -t 1-${jn} /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/submit_as_array.sh ${temp_file}"
    /opt/sge/bin/lx-amd64/qsub -V -S /bin/bash -l h_vmem=40G -pe smp 2 -N AD_spawner_$(basename ${temp_file}) -cwd -j y -b y -o /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/array_Logs/$(basename ${temp_file}) -tc 2 -t 1-${jn} /mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/scripts/submit_as_array.sh ${temp_file}
fi








