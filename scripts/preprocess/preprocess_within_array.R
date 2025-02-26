# gProcessedSignal, and rProcessedSignal data at this point are background
# subtracted and dye normalized by Agilent Feature Extraction software.

library(purrr)
setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')

# Still need to:
#   1. Extract both gProcessed channels for single patient into single df.
#   2a. Set signal of probes that are outliers == 1 to NA
#   2b. Set signal of probes that are IsSaturated == 1 to NA
#   2b. Subtract the negative control signal from other probe signals.
#   3. Split each file into signal file and well above file.

# 1
filter_fes_data2A <- function() {
  # Only keep both green channel for each patient. This will
  # ensure that regardless of microarray design, one green
  # channel represents one condition, and the other green
  # channel represents the other condition.
  sample_pairs_df <- read.csv(SAMP_PAIRS_FILE)
  unique_individuals <- unique(sample_pairs_df$Individual)
  num_inds <- length(unique_individuals)
  stopifnot(length(unique_individuals) == nrow(sample_pairs_df) / 2)
  filt1B_files <- list.files(AGIL_DATA_FES_FILT1B_DIR, full.names = TRUE)
  n <- 1
  for (unique_ind in unique_individuals) {
    message("Individual ", n, " / ", num_inds)
    samps <- sample_pairs_df %>%
      filter(Individual == unique_ind) %>% select(Sample)
    samps <- samps[, 1]
    samp_files <- filt1B_files[sapply(filt1B_files, function(file)
      {any(sapply(samps, function(sample) grepl(sample, file)))})]
    stopifnot(length(samp_files) == 2)
    file1_df <- read.delim(samp_files[1])
    file2_df <- read.delim(samp_files[2])
    # Use inverse because GeneName col
    file1_df <- file1_df[, !grepl("^r", colnames(file1_df))]
    file2_df <- file2_df[, !grepl("^r", colnames(file2_df))]
    # Save files
    file1_df_save_path <- file.path(AGIL_DATA_FES_FILT2A_DIR,
                                    basename(samp_files[1]))
    file2_df_save_path <- file.path(AGIL_DATA_FES_FILT2A_DIR,
                                    basename(samp_files[2]))
    write.table(file1_df, file = file1_df_save_path, sep = "\t",
                row.names = TRUE)
    write.table(file2_df, file = file2_df_save_path, sep = "\t",
                row.names = TRUE)
    n <- n + 1
  }
}

dir.create(AGIL_DATA_FES_FILT2A_DIR, showWarnings = TRUE)
filter_fes_data2A()

################

# 2a
replace_outliers <- function(df) {
  # Use the outlier columns (IsFeatNonUnifOL, IsBGNonUnifOL, IsFeatPopnOL,
  # IsBGPopnOL) to flag probes that have any outlier value equal to 1 by
  # setting corresponding signal values to NA.
  # IsFeatNonUnifOL
  df$gProcessedSignal[df$gIsFeatNonUnifOL == 1] <- NA
  df$gIsWellAboveBG[df$gIsFeatNonUnifOL == 1] <- NA
  # IsBGNonUnifOL
  df$gProcessedSignal[df$gIsBGNonUnifOL == 1] <- NA
  df$gIsWellAboveBG[df$gIsBGNonUnifOL == 1] <- NA
  # IsFeatPopnOL
  df$gProcessedSignal[df$gIsFeatPopnOL == 1] <- NA
  df$gIsWellAboveBG[df$gIsFeatPopnOL == 1] <- NA
  # IsBGPopnOL
  df$gProcessedSignal[df$gIsBGPopnOL == 1] <- NA
  df$gIsWellAboveBG[df$gIsBGPopnOL == 1] <- NA
  # Remove outlier columns
  df <- df[, (names(df) %in% c("GeneName",
                               "gProcessedSignal",
                               "gIsWellAboveBG",
                               "gIsSaturated"))]
  return(df)
}

# 2b
replace_saturated <- function(df) {
  # Use the IsSaturated column to flag probes that have IsSaturated
  # equal to 1 by setting corresponding signal values to NA.
  df$gProcessedSignal[df$gIsSaturated == 1] <- NA
  df$gIsWellAboveBG[df$gIsSaturated == 1] <- NA
  df <- df[, !(names(df) %in% c("gIsSaturated"))]
  return(df)
}

g_nc_medians <- numeric()
get_neg_control_medians <- function(df) {
  # About 600 negative control probe repeats, so get median for each channel.
  neg_controls <- df[df$GeneName == NEG_CONTROL, ]
  g_median <- median(neg_controls$gProcessedSignal, na.rm = TRUE)
  g_nc_medians <<- c(g_nc_medians, g_median)
  return(g_median)
}

# 2c
subtract_neg_control <- function(bg_df) {
  # Subtract the negative control medians from other probes. Do each
  # channel separately. Should be left with 41K rows.
  g_neg_control <- get_neg_control_medians(bg_df)
  # Remove negative control rows
  bg_df <- bg_df[bg_df$GeneName != NEG_CONTROL, ]
  # Subtract negative control values from other probes
  bg_df$gProcessedSignal <- pmax(0, bg_df$gProcessedSignal - g_neg_control)
  new_n_rows <- nrow(bg_df)
  if (new_n_rows == NUM_AGIL_GENES) {
    message("Success removing neg controls. ", new_n_rows, " rows")
  } else {
    message("Failed removing neg controls. ", new_n_rows, " rows")
  }
  return(bg_df)
}

# 2 full
# Function used to run filtering 2B data steps above.
filter_fes_data2B <- function() {
  filt2A_files <- list.files(AGIL_DATA_FES_FILT2A_DIR, full.names = TRUE)
  walk(filt2A_files, ~{
    message("Processing file: ", .x)
    filt2A_df <- read.delim(.x)
    filt2A_df <- replace_outliers(filt2A_df)
    filt2A_df <- replace_saturated(filt2A_df)
    filt2A_df <- subtract_neg_control(filt2A_df)
    # Make row names the gene names
    rownames(filt2A_df) <- filt2A_df$GeneName
    filt2A_df <- filt2A_df %>% select(-GeneName)
    save_path <- file.path(AGIL_DATA_FES_FILT2B_DIR, basename(.x))
    write.table(filt2A_df, file = save_path, sep = "\t", row.names = TRUE)
  })
}

# Takes about 2 minutes to run.
dir.create(AGIL_DATA_FES_FILT2B_DIR, showWarnings = TRUE)
filter_fes_data2B()
message('green neg control median: ', median(g_nc_medians))

################

# 3
filter_fes_data2C <- function () {
  # Each file at this point contains only gene names as row names, a
  # ProcessedSignal column, IsWellAboveBG, and IsFeatNonUnifOL column.
  # Split each file into two files with only ProcessedSignal, IsWellAboveBG.
  filt2B_files <- list.files(AGIL_DATA_FES_FILT2B_DIR, full.names=TRUE)
  n_filt2B_files <- length(filt2B_files)
  n <- 1
  for (filt2B_file in filt2B_files) {
    message("File ", n, " / ", n_filt2B_files)
    filt2B_df <- read.delim(filt2B_file)
    signal_df <- filt2B_df[, "gProcessedSignal", drop = FALSE]
    well_above_df <- filt2B_df[, "gIsWellAboveBG", drop = FALSE]
    save_name <- basename(filt2B_file)
    sig_save_path <- file.path(AGIL_DATA_FES_FILT2C_SIG_DIR, save_name)
    well_above_save_path <- file.path(AGIL_DATA_FES_FILT2C_WA_DIR, save_name)
    write.table(signal_df, file = sig_save_path, sep = "\t")
    write.table(well_above_df, file = well_above_save_path, sep = "\t")
    n <- n + 1
  }
}

# Takes about 1 minute to run.
dir.create(AGIL_DATA_FES_FILT2C_DIR, showWarnings = TRUE)
dir.create(AGIL_DATA_FES_FILT2C_SIG_DIR, showWarnings = TRUE)
dir.create(AGIL_DATA_FES_FILT2C_WA_DIR, showWarnings = TRUE)
filter_fes_data2C()
