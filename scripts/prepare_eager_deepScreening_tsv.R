#!/usr/bin/env Rscript

## Need sidora.core to pull sequencing metainfo
if (!require('sidora.core')) {
  if(!require('remotes')) install.packages('remotes')
  remotes::install_github('sidora-tools/sidora.core', quiet=T)
} else {library(sidora.core)}

## Need pandora2eager to generate TSV table
if (!require('pandora2eager')) {
  if(!require('remotes')) install.packages('remotes')
  remotes::install_github('sidora-tools/pandora2eager', quiet=T)
} else {library(pandora2eager)}

require(purrr)
require(tidyverse, warn.conflicts = F)
require(optparse)
require(readr)
require(stringr)

## Validate analysis type option input
validate_analysis_type <- function(option, opt_str, value, parser) {
  valid_entries <- c("Bacterial_Viral_Prescreening")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[prepare_eager_tsv.R] error: Invalid analysis type: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_sequencing_type <- function(option, opt_str, value, parser) {
  valid_entries <- c("No_Pathogen_Capture", "Pathogen_Capture", "All")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[prepare_eager_tsv.R] error: Invalid analysis type: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
parser <- add_option(parser, c("-s", "--sequencing_batch_id"), type = 'character', 
                     action = "store", dest = "sequencing_batch_id", 
                     help = "The Pandora sequencing batch ID for which Deep Screening will be runned. A TSV file will be prepared
			for each library in this run, containing all relevant processed fastq files resulting from the Autorun prescreening pipelines",
                     default = NA)
parser <- add_option(parser, c("-a", "--analysis_type"), type = 'character',
                     action = "callback", dest = "analysis_type",
                     callback = validate_analysis_type, default=NA,
                     help = "The analysis type to compile the data from. Should be one of: 'Bacterial_Viral_Prescreening'.")
parser <- add_option(parser, c("-o", "--outDir"), type = 'character',
                     action = "store", dest = "outdir",
                     help= "The desired output directory. Within this directory, one subdirectory will be 
			created per analysis type, within that one subdirectory per sequencing batch ID. By default, it is the current directory.",
                     default = ".")
parser <- add_option(parser, c("-t", "--sequencing_type"), type = 'character',
                     action = "callback", dest = "sequencing_type",
                     callback = validate_sequencing_type,
                     help= "The sequencing type one wants to deep screen. By default, it will output all the libraries in the tsv. Should be one of: 'No_Pathogen_Capture', 'Pathogen_Capture', 'All'",
                     default = "Allx")

arguments <- parse_args(parser, positional_arguments = 1)
opts <- arguments$options

cred_file <- arguments$args
sequencing_batch_id <- opts$sequencing_batch_id
analysis_type <- opts$analysis_type
sequencing_type <- opts$sequencing_type

if (is.na(sequencing_batch_id)) {
  stop(call.=F, "\n[prepare_eager_deepScreening_tsv.R] ERROR: No sequencing batch ID type provided with -s. Please see --help for more information.\n")
}

if (opts$outdir == ".") {
  write(paste0("[prepare_eager_deepScreening_tsv.R] No output directory specified, TSV will be stored in the currect directory under: ",analysis_type,"/",sequencing_type,"/",sequencing_batch_id), stdout())
}

output_dir <- paste0(opts$outdir,"/",analysis_type,"/",sequencing_type)
output_file <- paste0(output_dir,"/",sequencing_batch_id,"_eager_input.tsv")


#Check if output_dir exists 
if (!dir.exists(output_dir)) {
  write(paste0("[prepare_eager_deepScreening_tsv.R] Creating output directory '",output_dir,"'"), stdout())
  dir.create(output_dir, showWarnings = F, recursive = T) ## Create output directory and subdirs if they do not exist.
}

if(file.exists(output_file)){
  stop(call. = F, paste0("\n[prepare_eager_deepScreening_tsv.R] ",output_file," already exists, my job here is done!\nIf you want to remake the TSV please remove the already existing one"))
}
con <- sidora.core::get_pandora_connection(cred_file)

complete_pandora_table <- join_pandora_tables(
  get_df_list(
    c(make_complete_table_list(
      c("TAB_Site", "TAB_Analysis")
    )), con = con,
    cache = F
  )
) %>% 
  convert_all_ids_to_values(., con = con)

sequencingAll <- complete_pandora_table %>%
  filter(sequencing.Run_Id == sequencing_batch_id) %>%
  filter(sample.Ethically_culturally_sensitive == FALSE) %>%
  filter(!is.na(raw_data.Full_Raw_Data_Id))

autorunFinished <- sequencingAll %>%
  filter(analysis.Analysis_Id == analysis_type) %>%
  filter(analysis.Title == "Fastq mapped reads") 

if (length(unique(sequencingAll$raw_data.Full_Raw_Data_Id)) == length(unique(autorunFinished$raw_data.Full_Raw_Data_Id))) {
  write("Autorun has finished", stdout())
  results <- autorunFinished %>%
    select(individual.Full_Individual_Id,individual.Organism,library.Full_Library_Id,library.Protocol,analysis.Result,sequencing.Sequencing_Id,sequencing.Full_Sequencing_Id,sequencing.Single_Stranded) %>%
    mutate(
      Colour_Chemistry=4,
      SeqType="SE",
      Lane=row_number(),
      Strandedness=case_when(
        sequencing.Single_Stranded == 'yes' ~ "single",
        sequencing.Single_Stranded == 'no' ~ "double",
        is.na(sequencing.Single_Stranded) ~ "Unknown"), ## So far, any NAs are for libraries that were never sequenced, but just in case.
      inferred_udg=map_chr(library.Protocol, function(.){pandora2eager::infer_library_specs(.)[2]}),
      ## If UDG treatment cannot be assigned, but library is ssDNA, then assume none (since no trimming anyway)
      ## Also set to none for UDG treatment that could not be determined, the samples should still run through pathogen
      UDG_Treatment=case_when(
        Strandedness == 'single' & inferred_udg == 'Unknown' ~ "none",
        inferred_udg == 'Unknown' ~ "none",
        TRUE ~ inferred_udg
      ),
      R2=NA,
      BAM=NA,
      ## Add `_ss` to sample name for ssDNA libraries. Avoids file name collisions and allows easier merging of genotypes for end users.
      Sample_Name = case_when(
        sequencing.Single_Stranded == 'yes' ~ paste0(individual.Full_Individual_Id, "_ss"),
        TRUE ~ individual.Full_Individual_Id
      ),
      ## Also add the suffix to the Sample_ID part of the Library_ID. This ensures that in the MultiQC report, the ssDNA libraries will be sorted after the ssDNA sample.
      Library_ID = case_when(
        sequencing.Single_Stranded == 'yes' ~ paste0(Sample_Name, ".", stringr::str_split_fixed(library.Full_Library_Id, "\\.", 2)[,2]),
        TRUE ~ library.Full_Library_Id
      )
    ) %>%
    select(
      "Sample_Name",
      "Library_ID",
      "Lane",
      "Colour_Chemistry",
      "SeqType",
      "Organism"=individual.Organism,
      "Strandedness",
      "UDG_Treatment",
      "R1"=analysis.Result,
      "R2",
      "BAM"
    )

  if(sequencing_type == "No_Pathogen_Capture"){
    pathogenCaptureCodes <- read.csv("Probe_Set_Pathogens.csv",stringsAsFactors = F)
    results %>%
      filter(!str_detect(R1, paste(pathogenCaptureCodes$Short.Name, collapse = '|'))) %>%
      write_tsv(file= output_file)
  } 
  else if(sequencing_type == "Pathogen_Capture"){
    pathogenCaptureCodes <- read.csv("Probe_Set_Pathogens.csv",stringsAsFactors = F)
    results %>%
      filter(str_detect(R1, paste(pathogenCaptureCodes$Short.Name, collapse = '|'))) %>%
      write_tsv(file= output_file)
  } 
  else if(sequencing_type == "All"){
    write_tsv(results, file= output_file)
  }
  write(paste0("TSV created in: ",output_file), stdout())
  write("[prepare_eager_deepScreening_tsv.R] WARNING: if UDG treatment could not be determine, it was assumed that no UDG treatment was performed at all", stdout())
} else {
  write("[prepare_eager_deepScreening_tsv.R] WARNING: Prescreening has not finished yet!! No TSV produced.", stdout())
}

