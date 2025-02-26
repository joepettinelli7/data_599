# Probe sets / genes need to be filtered out that represent biological noise.

setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')
source('scripts/plotting.R')

# 1. Use the well_above_df to filter out probes by fraction IsWellAboveBG.
#    Remove any probes that are not IsWellAboveBG in at least 25% of samples
#    in at least one class (normal or nsclc).
# 2a. Get noise threshold with trimmed means.
# 2b. Use normalized_df to filter out probes with low amplitude.

# Load normalized signal data. Undo log2 transform from normalization.
# Make a single data frame with the IsWellAboveBG data. Add class labels to
# IsWellAboveBG data frame. Takes about 2 minutes to run.
# Should have 41,001 rows in each with response variable.
normalized_df <- read.delim(NORMALIZED_DF_FILE)
normalized_df <- unlog2_df(normalized_df)
well_above_df <- combine_files_into_df(AGIL_DATA_FES_FILT2C_WA_DIR, FALSE)
well_above_df <- impute_df_NAs(well_above_df)
well_above_df <- add_response_var(well_above_df)

# 1
get_probes_to_keep_single <- function(class_df) {
  # Get probes to keep for single class by fraction IsWellAboveBG.
  # Potentially twice as fast as finding probes to remove.
  # 
  # class_df: The data frame with samples from single class
  #
  # Returns: The probe names to keep
  num_samps <- ncol(class_df)
  stopifnot(num_samps == ncol(well_above_df) / 2)
  # Convert class_df to numeric
  class_df <- remove_response_var(class_df)
  # Calculate the number of well above BG for each probe (row)
  well_above_counts <- rowSums(class_df)
  keep_thresh <- num_samps * 0.25
  probes_to_keep <- rownames(class_df)[well_above_counts > keep_thresh]
  return(probes_to_keep)
}

# Split data frame into two by class.
class_dfs <- split_by_class(well_above_df)
probes_to_keep_norm <- get_probes_to_keep_single(class_dfs$normal_df)
probes_to_keep_nsclc <- get_probes_to_keep_single(class_dfs$nsclc_df)
probes_to_keep <- union(probes_to_keep_norm, probes_to_keep_nsclc)
frac_probes_to_rem <- setdiff(row.names(normalized_df), probes_to_keep)
n_frac_probes_to_rem <- length(frac_probes_to_rem)
message("Filter by fraction IsWellAboveBG found ", n_frac_probes_to_rem, 
        " genes.")
rm(well_above_df)
rm(class_dfs)

################

# 2a
get_trimmed_means <- function() {
  # Calculate trimmed means of normalized df before removing any genes
  # Get 4% trimmed mean for each sample (each microarray).
  #
  # Returns: The 246 trimmed means
  n_cols <- ncol(normalized_df)
  trimmed_means <- apply(normalized_df, 2, function(x) {
    lower_thresh <- quantile(x, 0.02)
    upper_thresh <- quantile(x, 0.98)
    inside_thresh_signals <- x[x > lower_thresh & x < upper_thresh]
    mean(inside_thresh_signals)
  })
  return(as.numeric(trimmed_means))
}

# Get mean of all trimmed means. This would be estimate of biological noise.
# Assert range of trimmed means is small.
trimmed_means <- get_trimmed_means()
stopifnot(length(trimmed_means) == ncol(normalized_df))
trimmed_means_range <- max(trimmed_means) - min(trimmed_means)
message("Trimmed means range is ", round(trimmed_means_range, 3), ".")
percent_range <- (trimmed_means_range / max(trimmed_means)) * 100
message("Trimmed means range is ", round(percent_range, 3), "% of max.")
trimmed_means_mean <- mean(trimmed_means)
message("Trimmed means mean is ", round(trimmed_means_mean, 3), ".")
barplot(trimmed_means, names.arg = NULL, main = "Trimmed Means",
        xlab = "Samples", ylab = "Log2 Expression" )
noise_thresh <- ceiling(trimmed_means_mean / 100) * 100
message("Noise thresh: ", noise_thresh)
message("Max value: ", max(normalized_df))

# Keep only probes_to_keep probes in the the normalized_df.
normalized_df <- normalized_df[!(rownames(normalized_df)
                                 %in% frac_probes_to_rem),]
stopifnot(nrow(normalized_df) == length(probes_to_keep))

# 2b
get_low_amp_genes <- function() {
  # Find low amplitude genes to remove now using noise_thresh from above.
  # Find genes with amplitude (over all samples) less than biological noise.
  #
  # Returns: The gene names
  amplitudes <- apply(normalized_df, 1, function(x) max(x) - min(x))
  low_amp_genes <- rownames(normalized_df)[amplitudes < noise_thresh]
  return(low_amp_genes)
}

low_amp_genes_to_rem <- get_low_amp_genes()
n_genes_to_rem <- length(low_amp_genes_to_rem)
message("Filter by low amplitude found ", n_genes_to_rem, " genes.")

# Count total number of removed genes
unique_remove_genes <- union(low_amp_genes_to_rem, frac_probes_to_rem)
message("Filtering methods found ", length(unique_remove_genes), " total.")

# Remove low_amp_genes_to_rem.
normalized_df_filt <- normalized_df[!(rownames(normalized_df)
                                       %in% low_amp_genes_to_rem),]
rm(normalized_df)
num_neg <- sum(normalized_df_filt < 0)
message("There are ", num_neg, " negative values. Should be 0.")
if (num_neg > 0) {
  paste("at", which(normalized_df_filt < 0, arr.ind = TRUE))
}

message("This leaves ", nrow(normalized_df_filt), " genes.")
normalized_df_filt <- log2_df(normalized_df_filt)
get_stats_for_df(normalized_df_filt)
vec <- unlist(normalized_df_filt)
vec <- as.numeric(vec)
make_hist(vec)
make_boxplot(normalized_df_filt, ylim = c(0, 20))
# Data is saved logged
write.table(normalized_df_filt, file =  FINAL_FILT_DF_FILE, sep = "\t")
rm(normalized_df_filt)
