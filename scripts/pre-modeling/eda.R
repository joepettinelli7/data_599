# Perform EDA on filtered data.

library(psych)
library(ggplot2)
library(plotly)
setwd(file.path(path.expand("~/Desktop"), "data_599"))
source('scripts/constants.R')
source('scripts/df_funcs.R')
source('scripts/plotting.R')

# Data was logged before saving during filtering
processed_df <- load_processed_df(log2 = FALSE, add_resp = TRUE)

# Get stats
class_dfs <- split_by_class(processed_df)
normal_df <- class_dfs$normal_df
nsclc_df <- class_dfs$nsclc_df
normal_df_no_res <- remove_response_var(normal_df)
nsclc_df_no_res <- remove_response_var(nsclc_df)
# Get stats by class
get_stats_for_df(normal_df_no_res)
get_stats_for_df(nsclc_df_no_res)

# Do PCA to visualize the classes.
# Transpose data frame to have genes as columns and samples as rows
resp_sorted_df <- cbind(normal_df_no_res, nsclc_df_no_res)
resp_sorted_df <- add_response_var(resp_sorted_df)
resp_sorted_df <- as.data.frame(t(resp_sorted_df))
response_vec <- as.factor(resp_sorted_df$disease)
resp_sorted_df <- remove_response_var(resp_sorted_df, in_row = FALSE)
pca_result <- prcomp(resp_sorted_df)
summary(pca_result)
pca_df <- as.data.frame(pca_result$x[, 1:3])
# Add class labels
response_vec <- as.factor(gsub("normal", "Normal", response_vec))
response_vec <- as.factor(gsub("nsclc", "NSCLC", response_vec))
pca_df$Class <- response_vec
ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) + 
  geom_point(alpha = 0.7, size = 2) + theme_minimal() + 
  labs(title = "PCA - First Two Principal Components by Diagnosis Class",
       x = "PC1", y = "PC2") +
  scale_color_manual(values = c("Normal" = "black", "NSCLC" = "red")) + 
  theme(title = element_text(size = 18),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))

plot_ly(
  pca_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Class,
  colors = c("Normal" = "black", "NSCLC" = "red"),
  type = 'scatter3d',
  mode = 'markers'
) %>%
  layout(
    scene = list(
      xaxis = list(title = 'PC1'),
      yaxis = list(title = 'PC2'),
      zaxis = list(title = 'PC3')
    ),
    title = "PCA - First Three Principal Components by Diagnosis Class"
)

# Get most influential variables
pca_loadings = pca_result$rotation
top_genes_PC1 <- pca_loadings[order(abs(pca_loadings[, "PC1"]),
                                    decreasing = TRUE), "PC1"]
head(top_genes_PC1, 5)

top_genes_PC2 <- pca_loadings[order(abs(pca_loadings[, "PC2"]),
                                    decreasing = TRUE), "PC2"]
head(top_genes_PC2, 5)

top_genes_PC3 <- pca_loadings[order(abs(pca_loadings[, "PC3"]),
                                     decreasing = TRUE), "PC3"]
head(top_genes_PC3, 5)

rm(class_dfs)
rm(normal_df)
rm(nsclc_df)
rm(normal_df_no_res)
rm(nsclc_df_no_res)
rm(pca_df)
rm(pca_loadings)
rm(pca_result)
rm(processed_df)
rm(resp_sorted_df)

# Do PCA will all NSCLC classes separate and no normal class
processed_df <- load_processed_df(log2 = FALSE, add_resp = FALSE)
all_class_df <- add_response_var(processed_df, all_classes = TRUE)
# all_class_df <- all_class_df[, all_class_df[1, ] != "normal"]
all_class_df <- as.data.frame(t(all_class_df))
response_vec <- as.factor(all_class_df$disease)
response_vec <- as.factor(gsub("normal", "Normal", response_vec))
all_class_df <- remove_response_var(all_class_df, in_row = FALSE)
pca_result <- prcomp(all_class_df)
summary(pca_result)
pca_df <- as.data.frame(pca_result$x[, 1:3])
pca_df$Class <- response_vec
ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) + 
  geom_point(alpha = 0.7, size = 2) + theme_minimal() + 
  labs(title = "PCA - First Two Principal Components by NSCLC Class",
       x = "PC1", y = "PC2") +
  scale_color_manual(values = c("Normal" = "black", "NSCLC" = "red",
                                "SCC" = "blue", "AC" = "gold",
                                "LCC" = "green")) + 
  theme(title = element_text(size = 18),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))

plot_ly(
  pca_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Class,
  colors = c("Normal" = "black", "NSCLC" = "red",
             "SCC" = "blue", "AC" = "gold",
             "LCC" = "green"),
  type = 'scatter3d',
  mode = 'markers'
) %>%
  layout(
    scene = list(
      xaxis = list(title = 'PC1'),
      yaxis = list(title = 'PC2'),
      zaxis = list(title = 'PC3')
    ),
    title = "PCA - First Three Principal Components by NSCLC Class"
)

rm(all_class_df)
rm(pca_df)
rm(pca_result)
rm(processed_df)

# Do PCA using gender instead of disease
add_gender_var <- function(df) {
  # Add the gender variable to the data frame first row.
  # If disease or gender variable already present, just return df
  #
  # df: The data frame without gender variable
  #
  # Returns the data frame with response variable in fist row
  if (rownames(df)[1] == "disease" | rownames(df)[1] == "gender") {
    message("Already have binary variable in this data frame!")
    return(df)
  }
  gender_vars_df <- read.csv(GENDER_VARS_FILE, row.names=1)
  gender_row <- factor(character(0))
  for (samp_name in colnames(df)) {
    gender <- gender_vars_df[samp_name, ]
    gender_row <- c(gender_row, gender)
  }
  df <- rbind(gender_row, df)
  rownames(df)[1] <- "gender"
  return(df)
}

processed_df <- load_processed_df(log2 = FALSE, add_resp = FALSE)
processed_df <- add_gender_var(processed_df)
processed_df <- as.data.frame(t(processed_df))
gender_vec <- as.factor(processed_df$gender)
processed_df <- remove_response_var(processed_df, in_row = FALSE)
pca_result <- prcomp(processed_df)
summary(pca_result)
pca_df <- as.data.frame(pca_result$x[, 1:3])
# Add class labels
pca_df$Class <- gender_vec
ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) + 
  geom_point(alpha = 0.7, size = 2) + theme_minimal() + 
  labs(title = "PCA - First Two Principal Components by Gender",
       x = "PC1", y = "PC2") +
  scale_color_manual(values = c("M" = "black", "F" = "red")) + 
  theme(title = element_text(size = 18),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))

rm(pca_df)
rm(pca_result)
