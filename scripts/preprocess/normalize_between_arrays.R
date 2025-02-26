# Need to normalize between arrays to make them comparable for
# filtering in the next step. Need to load all files into single data frame.

library(limma)
library(psych)
setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')
source('scripts/plotting.R')

# 1
# Combine all files into single data frame.
# Load signal data frame and log2 transform
signal_df <- combine_files_into_df(AGIL_DATA_FES_FILT2C_SIG_DIR, FALSE)
signal_df <- log2_df(signal_df)
signal_df <- impute_df_NAs(signal_df)
get_stats_for_df(signal_df)
vec <- unlist(signal_df)
vec <- as.numeric(vec)
make_hist(vec)
make_boxplot(signal_df, ylim=c(0, 20))
rm(vec)

################

# 2
# Now normalize signal data frame. Use quantile because cyclic loess
# introduced many negative values.
normalized_matrix <- normalizeBetweenArrays(signal_df, method = "quantile")
normalized_df <- as.data.frame(normalized_matrix)
rm(normalized_matrix)
rm(signal_df)
stopifnot(min(normalized_df) >= 0)

################

# 3
# Histogram, boxplot
get_stats_for_df(normalized_df)
vec <- unlist(normalized_df)
vec <- as.numeric(vec)
make_hist(vec)
make_boxplot(normalized_df, ylim = c(0, 20))
write.table(normalized_df, file =  NORMALIZED_DF_FILE, sep = "\t")
rm(vec)
