# Utility functions used through project.

setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')


split_extension <- function(file_name) {
  # Remove extension from file name
  #
  # file_name: The file name with extension
  #
  # Returns the name, ext
  parts <- strsplit(file_name, "\\.")
  name <- parts[[1]][1]
  ext <- parts[[1]][2]
  return(list(name, ext))
}


sample_id_from_file_path <- function(file_path) {
  # Get the sample unique sample id from file path
  #
  # file_path: File path or file name
  #
  # Returns the sample id (###)
  file_name <- basename(file_path)
  parts = strsplit(file_name, "_")
  full_id <- parts[[1]][2]
  sample_id <- substr(full_id, nchar(full_id) - 2, nchar(full_id))
  return(paste0("samp_", sample_id))
}


add_response_var <- function(df, all_classes = FALSE) {
  # Add the response variable to the data frame first row.
  # If response variable already present, just return df
  #
  # df: The data frame without response variable
  #
  # Returns the data frame with response variable in fist row
  if (rownames(df)[1] == "disease") {
    message("Already added response variable to this data frame!")
    return(df)
  }
  response_vars_df <- read.csv(RESPONSE_VARS_FILE, row.names=1)
  response_row <- factor(character(0))
  for (samp_name in colnames(df)) {
    disease_state <- response_vars_df[samp_name, ]
    if (!all_classes) {
      if (disease_state != "normal") {
        disease_state <- "nsclc"
      }
    }
    response_row <- c(response_row, disease_state)
  }
  df <- rbind(response_row, df)
  rownames(df)[1] <- "disease"
  return(df)
}


remove_response_var <- function(df, in_row = TRUE) {
  # Remove response variable from data frame
  # Assumes response variable in first column or row.
  # Also used to remove gender variable.
  #
  # df: The data frame with response variable
  # in_row: Whether the response variable is in first row
  #
  # Returns the data frame without the response variable
  if (colnames(df)[1] != "disease" & rownames(df)[1] != "disease") {
    if (colnames(df)[1] != "gender" & rownames(df)[1] != "gender") {
      message("Already removed response variable from this data frame!")
      return(df)
    }
  }
  if (in_row) {
    df <- df[-1, ]
  } else {
    df <- df[, -1]
  }
  df[] <- lapply(df, as.numeric)
  return(df)
}


split_by_class <- function(df) {
  # Split data frame into 2 by response variable.
  # Response variable needs to be in first row.
  #
  # df: The data frame with all samples
  #
  # Returns the data frame with normal samples
  # and data frame with nsclc samples
  response_classes <- as.character(df[1, ])
  normal_df <- df[, response_classes == "normal"]
  nsclc_df <- df[, response_classes != "normal"]
  return(list(normal_df = normal_df, nsclc_df = nsclc_df))
}


load_processed_df <- function(log2 = TRUE, add_resp = TRUE) {
  # Load the processed data frame which is the filtered data frame
  # with variables removed and 246 samples.
  # 
  # log2: Whether to log2 transform the data
  # add_resp: Whether to add response var
  #
  # Returns the data frame
  processed_df <- read.delim(FINAL_FILT_DF_FILE)
  processed_df[] <- lapply(processed_df, as.numeric)
  if (log2) {
    processed_df <- log2_df(processed_df)
  }
  if (add_resp) {
    processed_df <- add_response_var(processed_df)
  }
  return(processed_df)
}


combine_files_into_df <- function(dir_path, add_resp) {
  # Combine files from dir_path into a single data frame.
  # Each file is a column.
  #
  # dir_path: Directory where files are
  # add_resp: Whether to add response variable to data frame
  #
  # Returns: The data frame
  start_files <- list.files(dir_path, full.names=TRUE)
  n_start_files <- length(start_files)
  n <- 1
  combined_df <- NULL
  for (start_file in start_files) {
    message("Adding file: ", n, " / ", n_start_files)
    samp_name <- basename(start_file)
    samp_name <- split_extension(samp_name)[[1]]
    start_df <- read.csv(start_file, row.names=1, sep="")
    if (is.null(combined_df)) {
      # Add row names and column
      combined_df <- start_df
    } else {
      # Did not have row names column
      if (grepl(AGIL_DATA_FES_FILT2C_SIG_DIR, start_file)) {
        combined_df <- cbind(combined_df, start_df$gProcessedSignal)
      } else if (grepl(AGIL_DATA_FES_FILT2C_WA_DIR, start_file)) {
        combined_df <- cbind(combined_df, start_df$gIsWellAboveBG)
      }
    }
    colnames(combined_df)[ncol(combined_df)] <- samp_name
    n <- n + 1
  }
  if (add_resp) {
    combined_df <- add_response_var(combined_df)
  }
  return(combined_df)
}


get_stats_for_df <- function(df) {
  # Get basic stats for a data frame. Does not matter if variables are rows
  # or columns because stats are taken over full data frame.
  #
  # df: The data frame
  #
  # Returns: The stats for df
  vec <- unlist(df)
  describe_vec <- describe(vec, IQR = TRUE, skew = TRUE, range = TRUE)
  print(describe_vec)
  q1 <- quantile(vec, 0.25)
  q3 <- quantile(vec, 0.75)
  print(paste("q1: ", q1))
  print(paste("q3: ", q3))
  n_high_outliers <- sum(vec > q3 + (1.5 * describe_vec$IQR))
  print(paste("Num high outliers: ", n_high_outliers))
  n_low_outliers <- sum(vec < q1 - (1.5 * describe_vec$IQR))
  print(paste("Num low outliers: ", n_low_outliers))
}


log2_df <- function(df) {
  # Log the data frame. Set values less than
  # 1 to 1 so that they become 0 after log2.
  #
  # df: The data frame to log transform
  #
  # Returns: The log2 transformed data
  df[df < 1] <- 1
  df <- log2(df)
  return(df)
}


unlog2_df <- function(df) {
  # Unlog the data frame.
  # 
  # df: The data frame
  #
  # Returns: The unlogged data frame
  df <- (2 ^ df)
  return(df)
}


impute_df_NAs <- function(df) {
  # Replace NA values with median of row.
  # Start from last row because if too many
  # NA, then remove row.
  #
  # df: Data frame with samples as columns, genes as rows
  #
  # Returns the data frame with imputed NA
  num_rows <- nrow(df)
  num_NA_replaced <- 0
  num_rows_removed <- 0
  for (idx in num_rows:1) {
    idx_row <- df[idx, ]
    idx_row <- as.numeric(idx_row)
    NA_idxs <- which(is.na(idx_row), arr.ind = TRUE)
    if (length(NA_idxs > 0)) {
      if (length(NA_idxs) > ncol(df) * 0.25) {
        df <- df[-idx, ]
        num_rows_removed <- num_rows_removed + 1
      }
      else {
        row_median <- median(idx_row, na.rm = TRUE)
        for (NA_idx in NA_idxs) {
          # NA_idx will have row at first idx and col at second
          df[idx, NA_idx] <- row_median
          num_NA_replaced <- num_NA_replaced + 1 
        }
      }
    }
  }
  message(num_NA_replaced, " replaced NA vals.")
  message(num_rows_removed, " rows removed.")
  stopifnot(!anyNA(df))
  return(df)
}
