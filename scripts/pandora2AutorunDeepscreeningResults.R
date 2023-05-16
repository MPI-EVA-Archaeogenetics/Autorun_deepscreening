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
validate_sites <- function(option, opt_str, value, parser) {
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

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
##Make -f mandatory with valid_entries
parser <- add_option(parser, c("-f", "--file_sites_id"), type = 'character', 
                     action = "callback", dest = "file_libraries_id",
                     callback = validate_sites,
                     help = "CSV file containing all the site IDs",
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

results_path_deepscreening <- paste("/mnt/archgen/pathogen_resources/screening/Autorun_deepscreening/eager_outputs",analysis,type,sep = "/")

sites <- read.csv(opts$file_libraries_id, header = T)

complete_pandora_table <- join_pandora_tables(
  get_df_list(
    c(make_complete_table_list(
      c("TAB_Site", "TAB_Raw_Data")
    )), con = con
  )
)

sites_deepscreening <- complete_pandora_table %>%
  filter(site.Full_Site_Id %in% sites$Site_ID) %>%
  filter(!is.na(sequencing.Full_Sequencing_Id))%>%
  select(site.Full_Site_Id, individual.Full_Individual_Id, library.Full_Library_Id, sequencing.Full_Sequencing_Id, sequencing.Run_Id)

list_heatmap <- list()
for (sequencing in unique(sites_deepscreening$sequencing.Run_Id)) {
  heatmap <- paste(results_path_deepscreening,sequencing,"maltextract","results","heatmap_overview_Wevid.tsv", sep = "/")
  list_heatmap <- append(list_heatmap, heatmap)
}

matchExpression <- paste(unique(sites_deepscreening$individual.Full_Individual_Id), collapse = "|")

print(list_heatmap)
heatmap_combined <- list_heatmap %>% 
  lapply(read_tsv) %>% 
  reduce(full_join, by="node") %>% 
  mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
  select(node,matches(matchExpression))

head(heatmap_combined)



