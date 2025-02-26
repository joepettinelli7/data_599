# Split into training data and test data.
# Keep pairs of normal and NSCLC in the same data set.

setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')

# Data was logged before saving during filtering
processed_df <- load_processed_df(log2 = FALSE, add_resp = TRUE)

# Get splits
sample_pairs_df <- read.csv(SAMP_PAIRS_FILE)
unique_individuals <- unique(sample_pairs_df$Individual)
training_size = round(length(unique_individuals) * 0.80)
set.seed(SEED)
training_individuals <- sample(x = unique_individuals,
                               size = training_size,
                               replace = FALSE)
stopifnot(length(training_individuals) == training_size)
training_indices <- which(sample_pairs_df$Individual %in% training_individuals)
stopifnot(length(training_indices) == 2 * length(training_individuals))
training_samples <- sample_pairs_df$Sample[training_indices]

training_df <- processed_df[, colnames(processed_df) %in% training_samples]
test_df <- processed_df[, !(colnames(processed_df) %in% training_samples)]
stopifnot(length(intersect(colnames(training_df), colnames(test_df))) == 0)

# Save training df and test df that is log2 transformed
write.table(training_df, file = TRAINING_DF_FILE, sep = "\t", row.names = TRUE)
write.table(test_df, file = TEST_DF_FILE, sep = "\t", row.names = TRUE)

# Use the custom function to get stats for training and test set
no_response_training_df <- remove_response_var(training_df, in_row = TRUE)
no_response_test_df <- remove_response_var(test_df, in_row = TRUE)
get_stats_for_df(no_response_training_df)
get_stats_for_df(no_response_test_df)
rm(no_response_training_df)
rm(no_response_test_df)

# Get training data samples by class
training_samp_names <- names(training_df)
normal_idxs <- which(training_df[1, ] == "normal")
nsclc_idxs <- which(training_df[1, ] == "nsclc")
training_norm_samples <- training_samp_names[normal_idxs]
training_nsclc_samples <- training_samp_names[nsclc_idxs]
stopifnot(length(training_norm_samples) + length(training_nsclc_samples) 
          == ncol(training_df))

get_samp_ids <- function(full_names) {
  samp_ids <- character(0)
  for (full_name in full_names) {
    samp_id <- sample_id_from_file_path(full_name)
    samp_ids <- append(samp_ids, samp_id)
  }
  return(samp_ids)
}

training_norm_ids <- get_samp_ids(training_norm_samples)
training_nsclc_ids <- get_samp_ids(training_nsclc_samples)
training_norm_ids
training_nsclc_ids

# Get test set samples by class
test_samp_names <- names(test_df)
normal_idxs <- which(test_df[1, ] == "normal")
nsclc_idxs <- which(test_df[1, ] == "nsclc")
test_norm_samples <- test_samp_names[normal_idxs]
test_nsclc_samples <- test_samp_names[nsclc_idxs]
stopifnot(length(test_norm_samples) + length(test_nsclc_samples) 
          == ncol(test_df))
test_normal_ids <- get_samp_ids(test_norm_samples)
test_nsclc_ids <- get_samp_ids(test_nsclc_samples)
test_normal_ids
test_nsclc_ids

rm(processed_df)
rm(sample_pairs_df)
rm(test_df)
rm(training_df)
