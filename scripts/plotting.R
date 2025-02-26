# This file is used for plotting data (boxplot, pca, etc)

library(ggplot2)
library(tidyr)
setwd(file.path(path.expand("~/Desktop"), "data_599"))

make_hist <- function(vec) {
  hist(vec, main = "", ylab = "Count", xlab = "Log_2 Gene Expression")
}

make_boxplot <- function(df, ylim, response_vec = NULL) {
  # Make a box plot of df. Separate box plot by class.
  # df: Data frame should have samples as columns, variables
  #     as the rows. With or without response variable.
  # ylim: Numeric vector of length two for y axis range
  # response_vec: The response variable vector. If not NULL, skip
  # straight to plotting because df already sorted by class
  #
  # Returns:
  df <- add_response_var(df)
  class_dfs <- split_by_class(df)
  normal_df <- class_dfs$normal_df
  nsclc_df <- class_dfs$nsclc_df
  normal_df_no_res <- remove_response_var(normal_df)
  nsclc_df_no_res <- remove_response_var(nsclc_df)
  resp_sorted_df <- cbind(normal_df_no_res, nsclc_df_no_res)
  resp_sorted_df <- add_response_var(resp_sorted_df)
  response_vec <- unlist(resp_sorted_df[1, ])
  names(response_vec) <- NULL
  resp_sorted_df <- remove_response_var(resp_sorted_df)
  colors <- ifelse(response_vec == unique(response_vec)[1], "black", "salmon")
  boxplot(resp_sorted_df, 
          col = colors,
          outline = FALSE,
          ylim = ylim,
          main = "",
          xaxt = "n",
          xlab = "Samples",
          ylab = "Log_2 Expression")
  legend("topright", legend = unique(response_vec), fill = c("black", "salmon"))
}


plot_5_stats <- function(df, select_size, final_size = NULL) {
  # Plot the five summary stats for recursive feature elimination
  #
  # df: The rfe result
  # select_size: The size selected by rfe()
  # final_size: Size chosen manually
  data_long <- pivot_longer(df, cols = c(ROC, Accuracy, Sens, Spec), 
                            names_to = "Metric", values_to = "Value")
  data_long$Metric[data_long$Metric == "ROC"] <- "ROC-AUC"
  data_long$Variables <- factor(data_long$Variables,
                                levels = unique(data_long$Variables))
  ggplot(data_long, aes(x = Variables, y = Value, color = Metric,
                        group = Metric)) +
    geom_line() +
    geom_point() +
    labs(title = "Subset Size Performance",
         x = "Number of Variables",
         y = "Performance",
         color = "Summary Metric") +
    scale_color_manual(values = c("ROC-AUC" = "darkgreen", "Accuracy" = "blue",
                                  "Sens" = "red", "Spec" = "darkgray")) +
    theme_minimal() +
    theme(legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 16),
          axis.text.x = element_text(angle = -90, hjust = 1, size = 16),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16)) +
    scale_y_continuous(
      expand = c(0, 0), breaks = round(seq(from = min(data_long$Value), 
                                     to = max(data_long$Value),
                                     length.out = 15), 2)) +
    geom_vline(xintercept = as.factor(select_size), 
               color = "black", linetype = "dashed") +
    geom_vline(xintercept = as.factor(final_size), 
               color = "green", linetype = "dashed")
}


plot_summary_comparison <- function(rf_df, lin_svm_df, rbf_svm_df) {
  # Plot final fit model performance against each other
  #
  # rf_df: random forests data frame
  # lin_svm_df: linear kernel SVMs data frame
  # rbf_svm_df: RBF kernel SVMs data frame
  rf_df <- rf_df[, c("ROC", "Accuracy", "Sens", "Spec")]
  lin_svm_df <- lin_svm_df[, c("ROC", "Accuracy", "Sens", "Spec")]
  rbf_svm_df <- rbf_svm_df[, c("ROC", "Accuracy", "Sens", "Spec")]
  rf_df$Class <- "RF"
  lin_svm_df$Class <- "Linear SVM"
  rbf_svm_df$Class <- "RBF SVM"
  combined_df <- rbind(rf_df, lin_svm_df, rbf_svm_df)
  data_long <- pivot_longer(combined_df,
                            cols = c("ROC", "Accuracy", "Sens", "Spec"),
                            names_to = "Metric", values_to = "Value")
  data_long$Metric[data_long$Metric == "ROC"] <- "ROC-AUC"
  ggplot(data_long, aes(x = Metric, y = Value, color = Class)) +
    geom_point(size = 3) +
    labs(title = "Biomarker Comparison", x = "Metric", y = "Performance",
         color = "Model Class") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
}


plot_var_imp <- function(df, n, xaxis, color) {
  # Plot variable importance of final models
  #
  # df: The variable importance data frame
  # n: The number of variables to plot
  # xaxis: The name of x-axis. MDA for random forest and
  #         squared weights for lin SVM, change cost function for rbf svm. 
  # color: The color of the plot bars.
  
  
  
  df$bar_color <- ifelse(seq_along(df$rank_metric) %% 2 == 0, color, "white")
  top_df <- head(df, n)
  ggplot(top_df, aes(x = rank_metric, y = reorder(var, rank_metric),
                     fill = bar_color)) +
    geom_col(color = "black") +
    labs(title = "Variable Importance",
         x = xaxis,
         y = "Gene Transcript") +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18)) + 
    scale_x_continuous(labels = scales::comma) +
    scale_fill_identity()
}


plot_roc <- function(roc_df) {
  # Plot the ROC-AUC and misclassifications
  #
  # roc_df: The ROC df
  ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line() +
    labs(title = "ROC Curve", x = "FPR", y = "TPR") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16)) +
    geom_abline(intercept = 0, slope = 1, alpha = .3, linetype = "dashed") +
    annotate("point", x = roc_df$fpr[500], y = roc_df$tpr[500], color = "red")
}


plot_misclass <- function(df, thresholds, pos_class, neg_class) {
  # Plot misclassification vs threshold
  #
  # df: The data
  # thresholds: The thresholds 0-1
  # pos_class: The positive class
  # neg_class: The negative class
  misclass <- function(df, threshold) {
    fp <- sum(df[, pos_class] >= threshold & df$obs == neg_class)
    fn <- sum(df[, pos_class] < threshold & df$obs == pos_class)
    return(fp + fn)
  }
  df_for_plot <- data.frame(threshold = thresholds,
                            num_mis <- sapply(thresholds,
                                              function(x) misclass(df, x)))
  ggplot(df_for_plot, aes(x = threshold, y = num_mis)) + geom_line() +
    labs(title="Misclassifcations", x="Threshold", y="Number Misclassified") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16)) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed")
}
