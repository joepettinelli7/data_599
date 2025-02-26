# This script is used to do recursive feature selection
# with random forests and SVMs

setwd(file.path(path.expand("~/Desktop"), "data_599"))
library(parallel)
library(doParallel)
source('scripts/constants.R')
source('scripts/modeling/modeling_funcs.R')
source('scripts/modeling/train_e1071_svm.R')
source('scripts/df_funcs.R')
source('scripts/plotting.R')
library(dplyr)  # library dplyr after plyr

NUM_RF_TREES <- 2000
NUM_TUNE_RESAMPS <- 5

# Load training df with response variable. Already log2 transformed.
# Convert to rows have samples, columns have variables.
training_df <- read.delim(TRAINING_DF_FILE)
training_df <- as.data.frame(t(training_df))
train_disease_vec <- as.factor(training_df$disease)
stopifnot(length(levels(train_disease_vec)) == 2)
training_df <- remove_response_var(training_df, in_row=FALSE)
stopifnot(!any(duplicated(colnames(training_df))))
agil_gene_df <- read.delim(AGIL_GENE_LIST_FILE)

get_set_sizes <- function(fss) {
  # Min size is 2 because bug with rfe when subset size is 1.
  # Make in ascending order.
  first_val <- fss[1]
  if (first_val == 2) {
    return(fss)
  } else if (first_val > 32) {
    new_val <- first_val %/% 2
  } else {
    new_val <- first_val - 1
  }
  fss <- c(new_val, fss)
  return(get_set_sizes(fss))
}

# Get sizes vec
# Full set of variables automatically used.
start_size <- ncol(training_df) %/% 2
sizes_vec <- get_set_sizes(c(start_size))
sizes_vec

# Get elements for each external resampling iteration (Each parallel RFE).
set.seed(SEED)
resamps <- createResample(y = train_disease_vec, times = NUM_RFEs)
# Get seeds for each resampling iteration
all_seeds <- vector(mode = "list", length = NUM_RFEs + 1)
for(idx in 1:NUM_RFEs){
  all_seeds[[idx]] <- sample.int(n = 1000, size = (length(sizes_vec) + 1))
}
all_seeds[[NUM_RFEs + 1]] <- sample.int(n = 1000, 1)

# recursive feature elimination parameters for both RF and SVM
base_rfe_ctrl <- rfeControl(rerank = TRUE,
                            method ="boot",
                            saveDetails = TRUE,
                            number = NUM_RFEs,
                            verbose = TRUE,
                            returnResamp = "all",
                            index = resamps,
                            seeds = all_seeds,
                            allowParallel = TRUE
)

# Detect number of cores, not threads.
num_cores <- detectCores(logical = FALSE)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

########## Feature selection for RF.
base_rfe_ctrl$functions <- rf_rfe_funcs
timestamp()
rf_rfe <- rfe(x = training_df,
              y = train_disease_vec,
              rfeControl = base_rfe_ctrl,
              metric = "ROC",
              maximize = TRUE,
              ntree = NUM_RF_TREES,
              sizes = sizes_vec
)
save(rf_rfe, file = RF_RFE_FILE)
rf_rfe
rf_man_size <- 13
plot_5_stats(rf_rfe$results, rf_rfe$optsize, final_size = rf_man_size)
# Update with chosen subset size
rf_final <- update(rf_rfe, training_df, train_disease_vec, rf_man_size)
rf_opt_rows <- agil_gene_df %>% dplyr::filter(ProbeID %in% rf_final$bestVar)
rf_opt_rows
rf_rfe$optsize<- rf_man_size
save(rf_final, file = RF_FINAL_FILE)

##########
# Train control for linear and RBF SVM
tr_control = trainControl(method = "boot",
                          number = NUM_TUNE_RESAMPS,
                          summaryFunction = fivestats_summary,
                          selectionFunction = "best",
                          classProbs = TRUE,
                          verboseIter = FALSE)

########## Feature selection for linear SVM.
training_matrix <- as.matrix(training_df)
base_rfe_ctrl$functions <- svm_rfe_funcs_linear
lin_svm_tune_grid <- expand.grid(C = 2^c(-6:-2))
lin_svm_tune_grid
timestamp()
lin_svm_rfe <- rfe(x = training_matrix,
                   y = train_disease_vec,
                   rfeControl = base_rfe_ctrl,
                   metric = "ROC",
                   maximize = TRUE,
                   sizes = sizes_vec,
                   tuneGrid = lin_svm_tune_grid,
                   method = "svmLinear",
                   preProc = c("center", "scale"),
                   trControl = tr_control
)
save(lin_svm_rfe, file = SVM_RFE_LIN_FILE)
lin_svm_man_size <- 2

# Get most chosen C from inner resamples
lin_C_counts <- lin_svm_rfe$variables %>%
    dplyr::filter(Variables == lin_svm_rfe$optsize) %>%
    dplyr::count(C) %>%
    dplyr::arrange(desc(n))
lin_C_counts
lin_svm_rfe$fit$bestTune$C <- lin_C_counts[1, 1]

plot_5_stats(lin_svm_rfe$results, lin_svm_rfe$optsize,
             final_size = lin_svm_man_size)
lin_svm_opt_rows <- agil_gene_df %>% dplyr::filter(ProbeID %in% 
                                                       lin_svm_rfe$optVariables)
lin_svm_opt_rows
lin_svm_final <- update(lin_svm_rfe, training_matrix, train_disease_vec,
                        lin_svm_man_size)
lin_svm_rfe$optsize <- lin_svm_man_size
save(lin_svm_final, file = SVM_FINAL_LIN_FILE)

########## Feature selection for RBF SVM.
base_rfe_ctrl$functions <- svm_rfe_funcs_rbf
gamma_vals <- as.numeric(sigest(training_matrix, frac = 1.0, scaled = TRUE))
rbf_svm_tune_grid <- expand.grid(cost = 2^c(-4:-1), gamma = gamma_vals)
rbf_svm_tune_grid
timestamp()
rbf_svm_rfe <- rfe(x = training_matrix,
                   y = train_disease_vec,
                   rfeControl = base_rfe_ctrl,
                   metric = "ROC",
                   maximize = TRUE,
                   sizes = sizes_vec,
                   tuneGrid = rbf_svm_tune_grid,
                   method = e1071_rbf_svm,
                   preProc = c("center", "scale"),
                   trControl = tr_control
)
save(rbf_svm_rfe, file = SVM_RFE_RBF_FILE)
rbf_svm_man_size <- 29

# Get most chosen cost (C) and gamma (sigma)
rbf_cost_counts <- rbf_svm_rfe$variables %>%
    dplyr::filter(Variables == lin_svm_rfe$optsize) %>%
    dplyr::count(cost) %>%
    dplyr::arrange(desc(n))
rbf_cost_counts
rbf_gamma_counts <- rbf_svm_rfe$variables %>%
    dplyr::filter(Variables == lin_svm_rfe$optsize) %>%
    dplyr::count(gamma) %>%
    dplyr::arrange(desc(n))
rbf_gamma_counts

plot_5_stats(rbf_svm_rfe$results, rbf_svm_rfe$optsize,
             final_size = rbf_svm_man_size)
# Update with chosen subset size
rbf_svm_rfe$fit$bestTune$cost <- rbf_cost_counts[1, 1]
rbf_svm_rfe$fit$bestTune$gamma <- rbf_gamma_counts[1, 1]
rbf_svm_final <- update(rbf_svm_rfe, training_matrix, train_disease_vec,
                        rbf_svm_man_size)
rbf_svm_opt_rows <- agil_gene_df %>% dplyr::filter(ProbeID %in% 
                                                       rbf_svm_final$bestVar)
rbf_svm_opt_rows
rbf_svm_rfe$optsize <- rbf_svm_man_size
save(rbf_svm_final, file = SVM_FINAL_RBF_FILE)

stopCluster(cl)
