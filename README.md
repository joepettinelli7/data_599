# data_599

## What does this project do?
- This project contains the code that was used for my MS thesis. The thesis_paper.pdf document contains the written part and results.
  The main objective of the project was to identify a transcriptomic biomarker for predicting non-small cell lung cancer using random forests and
  support vector machines as the learning algorithms.

## Why is this project useful?
- This project is useful because it contains a full explanation and blueprint for identifying a transcriptomic biomarker.

## How can users get started?
1. Clone this repository to desktop.
2. Download data from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1132 (about 17GB and 246 .txt files).
3. Put the files in Desktop/data_599/data/agilent_data_fes/ folder.

## General workflow:
Run scripts in the following order (would take about 7 days to complete full workflow using iMac with 6-Core Intel i5 processor and 8GB DDR4):

preprocess/
1. preprocess_fes_data.R
2. preprocess_within_array.R
3. normalize_between_arrays.R
4. filter_vars.R

pre-modeling/
1. eda.R
2. partition_data.R

modeling/
1. feature_selection.R
2. biomarker_comparison.R
3. model_evaluation.R
