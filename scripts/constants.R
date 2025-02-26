# Paths from current working directory
setwd(file.path(path.expand("~/Desktop"), "data_599"))

SEED <- 100

# Directories
AGIL_DATA_FES_DIR <- file.path("data", "agilent_data_fes")
AGIL_DATA_FES_FILT1A_DIR <- file.path("data", "agilent_data_fes_filt1A")
AGIL_DATA_FES_FILT1B_DIR <- file.path("data", "agilent_data_fes_filt1B")
AGIL_DATA_FES_FILT2A_DIR <- file.path("data", "agilent_data_fes_filt2A")
AGIL_DATA_FES_FILT2B_DIR <- file.path("data", "agilent_data_fes_filt2B")
AGIL_DATA_FES_FILT2C_DIR <- file.path("data", "agilent_data_fes_filt2C")
AGIL_DATA_FES_FILT2C_SIG_DIR <- file.path(AGIL_DATA_FES_FILT2C_DIR, "signal")
AGIL_DATA_FES_FILT2C_WA_DIR <- file.path(AGIL_DATA_FES_FILT2C_DIR, "well_ab")

# Files
ALL_SAMP_NAMES_FILE <- file.path("data", "all_agilent_sample_names.txt")
AGIL_GENE_LIST_FILE <- file.path("data", "014850_D_GeneList_20240521.txt")
COMMON_PROBES_FILE <- file.path("data", "common_probes.RData")
RESPONSE_VARS_FILE <- file.path("data", "response_variables.txt")
SAMP_PAIRS_FILE <- file.path("data", "sample_pairs.txt")
GENDER_VARS_FILE <- file.path("data", "gender_variables.txt")
NORMALIZED_DF_FILE <- file.path("data", "normalized_df.txt")
FINAL_FILT_DF_FILE <- file.path("data", "processed_df.txt")
TRAINING_DF_FILE <- file.path("data", "training_df.txt")
TEST_DF_FILE <- file.path("data", "test_df.txt")
RF_RFE_FILE <- file.path("data", "rf_rfe.RData")
SVM_RFE_LIN_FILE <- file.path("data", "svm_rfe_linear.RData")
SVM_RFE_RBF_FILE <- file.path("data", "svm_rfe_rbf.RData")
RF_FINAL_FILE <- file.path("data", "rf_final.RData")
SVM_FINAL_LIN_FILE <- file.path("data", "svm_final_linear.RData")
SVM_FINAL_RBF_FILE <- file.path("data", "svm_final_rbf.RData")

# Agilent data
START_ROW <- 10
NEG_CONTROL <- "(-)3xSLv1"
NUM_AGIL_GENES <- 41000

# Feature selection
NUM_RFEs <- 100
ONE_SE_THRESH <- 20
