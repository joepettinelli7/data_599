# Test the final model on the test data set.

setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')
source('scripts/modeling/modeling_funcs.R')
source('scripts/plotting.R')

# Load the test data
test_df <- read.delim(TEST_DF_FILE)
test_df <- as.data.frame(t(test_df))
test_disease_vec <- as.factor(test_df$disease)
stopifnot(length(levels(test_disease_vec)) == 2)
test_df <- remove_response_var(test_df, in_row=FALSE)

# Load final models
load(RF_FINAL_FILE)
load(SVM_FINAL_LIN_FILE)
load(SVM_FINAL_RBF_FILE)
set.seed(SEED)

# RF-RFE eval
rf_test_df <- test_df[, rf_final$bestVar]
rf_preds <- rf_pred(rf_final$fit, rf_test_df)
rf_preds$obs <- test_disease_vec
fivestats_summary(rf_preds, plot_roc = TRUE)

# Linear SVM-RFE eval
lin_svm_test_df <- test_df[, lin_svm_final$bestVar]
lin_svm_preds <- svm_pred(lin_svm_final$fit, lin_svm_test_df)
lin_svm_preds$obs <- test_disease_vec
fivestats_summary(lin_svm_preds, plot_roc = TRUE)

# RBF SVM-RFE eval
rbf_svm_test_df <- test_df[, rbf_svm_final$bestVar]
rbf_svm_preds <- svm_pred(rbf_svm_final$fit, rbf_svm_test_df)
rbf_svm_preds$obs <- test_disease_vec
fivestats_summary(rbf_svm_preds, plot_roc = TRUE)

# Get confusion matrices
rf_conf_matrix <- confusionMatrix(rf_preds$pred, test_disease_vec)
rf_conf_matrix$table
lin_svm_conf_matrix <- confusionMatrix(lin_svm_preds$pred, test_disease_vec)
lin_svm_conf_matrix$table
rbf_svm_conf_matrix <- confusionMatrix(rbf_svm_preds$pred, test_disease_vec)
rbf_svm_conf_matrix$table
