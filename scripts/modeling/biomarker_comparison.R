# Compare the biomarkers from the three feature selection experiments.
# Use the original rfe object, not updated, to avoid internal validation.

setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')
source('scripts/modeling/modeling_funcs.R')
source('scripts/plotting.R')

# Plot summaries against each other at chosen subset size.
# optsize should be the manually chosen subset size.
rf_opt <- rf_rfe$results[rf_rfe$results$Variables == rf_rfe$optsize, ]

lin_svm_opt <- lin_svm_rfe$results[lin_svm_rfe$results$Variables == 
                                       lin_svm_rfe$optsize, ]

rbf_svm_opt <- rbf_svm_rfe$results[rbf_svm_rfe$results$Variables == 
                                       rbf_svm_rfe$optsize, ]

plot_summary_comparison(rf_opt, lin_svm_opt, rbf_svm_opt)

# Plot variable importance. Must use updated model with correct variables.
# RF-RFE
load(RF_FINAL_FILE)
ranked_rf_vars <- rf_rank(rf_final$fit)
sorted_rf_vars <- dplyr::arrange(ranked_rf_vars, desc(rank_metric))
plot_var_imp(sorted_rf_vars,
             nrow(sorted_rf_vars),
             'Mean Decrease Accuracy',
             "slategray1")

# Linear SVM-RFE
load(SVM_FINAL_LIN_FILE)
ranked_lin_svm_vars <- svm_rank_linear(lin_svm_final$fit)
sorted_lin_svm_vars <- dplyr::arrange(ranked_lin_svm_vars, desc(rank_metric))
plot_var_imp(sorted_lin_svm_vars,
             nrow(sorted_lin_svm_vars),
             'Squared Weight',
             rgb(249/255, 226/255, 213/255))

# RBF SVM-RFE
load(SVM_FINAL_RBF_FILE)
ranked_rbf_svm_vars <- svm_rank_rbf(rbf_svm_final$fit)
sorted_rbf_svm_vars <- dplyr::arrange(ranked_rbf_svm_vars, desc(rank_metric))
plot_var_imp(sorted_rbf_svm_vars,
             nrow(sorted_rbf_svm_vars),
             'Change in Cost Function',
             rgb(249/255, 226/255, 213/255))
