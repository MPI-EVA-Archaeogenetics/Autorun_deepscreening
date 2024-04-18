#Script to crawl Autorun_deepscreening nf-core/eager results
#Authors: Aida Andrades Valtueña

## Load necessary libraries ##
library(tidyverse)
library(sidora.core)

if (!require('sidora.core')) {
  if(!require('remotes')) install.packages('remotes')
  remotes::install_github('sidora-tools/sidora.core', quiet=T)
} else {library(sidora.core)}

require(purrr)
require(tidyverse, warn.conflicts = F)
require(optparse)
require(readr)
require(stringr)

#Validation of arguments
validate_file <- function(option, opt_str, value, parser) {
  ifelse(!is.na(value), 
         return(value), 
         stop(call.=F, "\n[pandora2AutorunDeepscreeningResults.R] error: No file with libraries provided!\n Please provide a file with library IDs using the -f option"))
}
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

## Functions
extractLibraries <- function(nameTable, folder, matchExpression, outdir, pathDeepScreening, librariesdeepscreening, column) {
  list_tables <- list()
  for (sequencing in unique(librariesdeepscreening$sequencing.Run_Id)) {
    table <- paste(pathDeepScreening,sequencing,"maltextract","results",folder,nameTable, sep = "/")
    if (file.exists(table)) {
      list_tables <- append(list_tables, table)
    }
  }
  
  print(list_tables)
  table_combined <- list_tables %>% 
    lapply(read_tsv) %>% 
    reduce(full_join, by=column) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
    select(column,matches(matchExpression))
  
  if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = F, recursive = T) ## Create output directory and subdirs if they do not exist.
}

  output_table <- paste(outdir, nameTable, sep = "/")
  write_tsv(table_combined, output_table)
  
}

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
##Make -f mandatory with valid_entries
parser <- add_option(parser, c("-f", "--file_libraries_id"), type = 'character', 
                     action = "callback", dest = "file_libraries_id",
                     callback = validate_file,
                     help = "CSV file containing the libraries IDs to gather results for",
                     default = NA
                     )
parser <- add_option(parser, c("-a", "--analysis_type"), type = 'character',
                     action = "callback", dest = "analysis_type",
                     callback = validate_analysis_type, default=NA,
                     help = "The analysis type to compile the data from. Should be one of: 'Bacterial_Viral_Prescreening'.")
parser <- add_option(parser, c("-t", "--sequencing_type"), type = 'character',
                     action = "callback", dest = "sequencing_type",
                     callback = validate_sequencing_type,
                     help= "The sequencing type one wants to deep screen. By default, it will output all the libraries in the tsv. Should be one of: 'No_Pathogen_Capture', 'Pathogen_Capture', 'All'",
                     default = "All")
parser <- add_option(parser, c("-o", "--outDir"), type = 'character',
                     action = "store", dest = "outdir",
                     help= "The desired output directory. By default, it is the current directory.",
                     default = "."
                     )

arguments <- parse_args(parser, positional_arguments = 1)
opts <- arguments$options

cred_file <- arguments$args
con <- sidora.core::get_pandora_connection(cred_file)
analysis <- opts$analysis_type
type <- opts$sequencing_type

#Set path to the deep screening results
results_path_deepscreening <- paste("/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_outputs",analysis,type,sep = "/")

#Read file containing libraries IDs to extract deep screening results
libraries <- read.csv(opts$file_libraries_id, header = T)

#Creating output directory if it doesn't exist
output_dir <- paste(opts$outdir,"maltextract", sep = "/")
if (!dir.exists(output_dir)) {
  write(paste0("[pandora2AutorunDeepscreeningResults.R] Creating output directory '",output_dir,"'"), stdout())
  dir.create(output_dir, showWarnings = F, recursive = T) ## Create output directory and subdirs if they do not exist.
}

complete_pandora_table <- join_pandora_tables(
  get_df_list(
    c(make_complete_table_list(
      c("TAB_Site", "TAB_Raw_Data")
    )), con = con
  )
)

libraries_deepscreening <- complete_pandora_table %>%
  filter(library.Full_Library_Id %in% libraries$Library_Id) %>%
  filter(!is.na(sequencing.Full_Sequencing_Id))%>%
  select(library.Full_Library_Id, individual.Full_Individual_Id, library.Full_Library_Id, sequencing.Full_Sequencing_Id, sequencing.Run_Id)

write_tsv(libraries_deepscreening, paste(output_dir, "Libraries_info.tsv", sep = "/"))
#libraries_deepscreening

#Set a match expression to search in the different tables
match_expression <- paste(unique(libraries_deepscreening$individual.Full_Individual_Id), collapse = "|")


extractLibraries("heatmap_overview_Wevid.tsv", "", match_expression, output_dir, results_path_deepscreening, libraries_deepscreening, "node")
extractLibraries("RunSummary.txt", "ancient",match_expression, paste(output_dir, "ancient", sep = "/"), results_path_deepscreening, libraries_deepscreening, "Node")
extractLibraries("TotalCount.txt", "ancient",match_expression, paste(output_dir, "ancient", sep = "/"), results_path_deepscreening, libraries_deepscreening, "Node")
extractLibraries("RunSummary.txt", "default",match_expression, paste(output_dir, "default", sep = "/"), results_path_deepscreening, libraries_deepscreening, "Node")
extractLibraries("TotalCount.txt", "default",match_expression, paste(output_dir, "default", sep = "/"), results_path_deepscreening, libraries_deepscreening, "Node")

# TODO: Copy rma6 files to malt folder
# TODO: Copy pdf profiles to maltextract/pdf_candidates folder 

#TODO complete table dumps
#TODO: Species name dumps