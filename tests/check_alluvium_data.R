#!/usr/bin/env Rscript
# Test script to examine alluvium plot data

library(readxl)
library(tidyverse)

# Read the data
file_path <- "/Users/lixiang/Downloads/PlotCase/物种组成冲积图/NCOMMS-24-29043C_Source Data to Main Figures.xlsx"

cat("Reading data from:\n", file_path, "\n\n")

df <- read_excel(file_path, sheet = "Fig.1a", skip = 1)

cat("Data dimensions:", nrow(df), "rows x", ncol(df), "columns\n\n")

cat("Column names:\n")
print(names(df))

cat("\n\nFirst few rows:\n")
print(head(df))

cat("\n\nData structure:\n")
str(df)

cat("\n\nColumn types:\n")
sapply(df, class)
