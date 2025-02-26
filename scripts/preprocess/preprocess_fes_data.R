# Data is on disk and downloaded
# from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1132.
# It contains data that has been processed with Agilent Feature
# Extraction software.

library(dplyr)
setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')

# 1. Confirm no missing files during download (should be 246 files).
# 2. Create a separate file for each sample that includes the 
#    processed signal data from both the green and red channels.
# 3. Remove unwanted columns.

# 1
print_missing_files <- function() {
  # Print files that are missing after download.
  all_246_files <- unique(readLines(ALL_SAMP_NAMES_FILE))
  downloaded_files <- list.files(AGIL_DATA_FES_DIR)
  missing_files <- list()
  for (file_name in all_246_files) {
    samp_name <- strsplit(file_name, "\\.txt")[[1]][1]
    n_matches <- grep(samp_name, downloaded_files)
    if (length(n_matches) == 0) {
      missing_files <- append(missing_files, samp_name)
    }
  }
  num_missing = length(missing_files)
  message(num_missing, " missing files: ", missing_files, ".")
}
print_missing_files()

################

# 2
make_empty_df <- function() {
  # Make empty df that will be used as template.
  return(
    data.frame(
      GeneName = character(),
      gProcessedSignal = numeric(0),
      rProcessedSignal = numeric(0),
      gMedianSignal = numeric(0),
      rMedianSignal = numeric(0),
      gIsSaturated = logical(0),
      rIsSaturated = logical(0),
      gIsFeatNonUnifOL = logical(0),
      rIsFeatNonUnifOL = logical(0),
      gIsBGNonUnifOL = logical(0),
      rIsBGNonUnifOL = logical(0),
      gIsFeatPopnOL = logical(0),
      rIsFeatPopnOL = logical(0),
      gIsBGPopnOL = logical(0),
      rIsBGPopnOL = logical(0),
      gIsWellAboveBG = logical(0),
      rIsWellAboveBG = logical(0),
      PositionX = numeric(0),
      PositionY = numeric(0)
    )
  )
}

filter_fes_data1A <- function(agil_probe_ids) {
  # Load each file into R and extract the processed signal data.
  # Don't do this for genes not in the Agilent Gene List. 
  # Save text file to agilent_data_fes_filt1 dir.
  #
  # agil_probe_ids: Probe IDs to keep because documented by agilent.
  fes_files <- list.files(AGIL_DATA_FES_DIR, full.names=TRUE)
  n_fes_files <- length(fes_files)
  for (n in 1:n_fes_files) {
    message("File ", n, " / ", n_fes_files)
    fes_file <- fes_files[n]
    fes_df <- read.delim(fes_file, header=FALSE)
    rows <- list()
    for (idx in seq(START_ROW, nrow(fes_df))) {
      row_data <- fes_df[idx, ]
      if (idx == START_ROW) {
        # Confirm data frame is aligned as expected
        stopifnot(row_data$V1 == "FEATURES")
        stopifnot(row_data$V7 == "ProbeName")
        stopifnot(row_data$V14 == "gProcessedSignal")
        stopifnot(row_data$V15 == "rProcessedSignal")
        stopifnot(row_data$V18 == "gMedianSignal")
        stopifnot(row_data$V19 == "rMedianSignal")
        stopifnot(row_data$V24 == "gIsSaturated")
        stopifnot(row_data$V25 == "rIsSaturated")
        stopifnot(row_data$V26 == "gIsFeatNonUnifOL")
        stopifnot(row_data$V27 == "rIsFeatNonUnifOL")
        stopifnot(row_data$V28 == "gIsBGNonUnifOL")
        stopifnot(row_data$V29 == "rIsBGNonUnifOL")
        stopifnot(row_data$V30 == "gIsFeatPopnOL")
        stopifnot(row_data$V31 == "rIsFeatPopnOL")
        stopifnot(row_data$V32 == "gIsBGPopnOL")
        stopifnot(row_data$V33 == "rIsBGPopnOL")
        stopifnot(row_data$V39 == "gIsWellAboveBG")
        stopifnot(row_data$V40 == "rIsWellAboveBG")
        stopifnot(row_data$V9 == "PositionX")
        stopifnot(row_data$V10 == "PositionY")
      } else {
        # Potential gene row
        probe_name <- row_data$V7
        if (probe_name %in% agil_probe_ids | probe_name == NEG_CONTROL) {
          # Add data to data frame
          new_row <- list(
            GeneName = row_data$V7,
            gProcessedSignal = row_data$V14,
            rProcessedSignal = row_data$V15,
            gMedianSignal = row_data$V18,
            rMedianSignal = row_data$V19,
            gIsSaturated = row_data$V24,
            rIsSaturated = row_data$V25,
            gIsFeatNonUnifOL = row_data$V26,
            rIsFeatNonUnifOL = row_data$V27,
            gIsBGNonUnifOL = row_data$V28,
            rIsBGNonUnifOL = row_data$V29,
            gIsFeatPopnOL = row_data$V30,
            rIsFeatPopnOL = row_data$V31,
            gIsBGPopnOL = row_data$V32,
            rIsBGPopnOL = row_data$V33,
            gIsWellAboveBG = row_data$V39,
            rIsWellAboveBG = row_data$V40,
            PositionX = row_data$V9,
            PositionY = row_data$V10
          )
          rows <- append(rows, list(new_row))
        }
      }
    }
    # Convert to data frame now
    df <- do.call(rbind, rows)
    save_name <- basename(fes_file)
    save_path <- file.path(AGIL_DATA_FES_FILT1A_DIR, save_name)
    write.table(df, file = save_path, sep = "\t", row.names = FALSE)
    n <- n + 1 
  }
}

# Takes about 32 hours to run. Will be 246 files saved.
agil_gene_df <- read.delim(AGIL_GENE_LIST_FILE)
agil_probe_ids <- agil_gene_df$xProbeID
message(length(agil_probe_ids), " probe ids in gene list.")
dir.create(AGIL_DATA_FES_FILT1A_DIR, showWarnings = TRUE)
timestamp()
filter_fes_data1A(agil_probe_ids)

################

# 3
filter_fes_data1B <- function () {
  # Remove unwanted columns because filter_fes_data1A saved all
  # columns that are potentially needed.
  filt1A_files <- list.files(AGIL_DATA_FES_FILT1A_DIR, full.names=TRUE)
  n_filt1A_files <- length(filt1A_files)
  n <- 1
  for (filt1A_file in filt1A_files) {
    message("File ", n, " / ", n_filt1A_files)
    filt1A_df <- read.delim(filt1A_file)
    cols_to_keep <- c("GeneName", "gProcessedSignal", "rProcessedSignal",
                      "gIsSaturated", "rIsSaturated", "gIsFeatNonUnifOL",
                      "rIsFeatNonUnifOL", "gIsBGNonUnifOL", "rIsBGNonUnifOL",
                      "gIsFeatPopnOL", "rIsFeatPopnOL", "gIsBGPopnOL",
                      "rIsBGPopnOL", "gIsWellAboveBG", "rIsWellAboveBG")
    filt1A_df <- filt1A_df %>% select(all_of(cols_to_keep))
    save_name <- basename(filt1A_file)
    save_path <- file.path(AGIL_DATA_FES_FILT1B_DIR, save_name)
    write.table(filt1A_df, file = save_path, sep = "\t", row.names = FALSE)
    n <- n + 1
  }
}

# Takes about 2 minute to run.
dir.create(AGIL_DATA_FES_FILT1B_DIR, showWarnings = TRUE)
filter_fes_data1B()
