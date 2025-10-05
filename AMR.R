# Loading Libraries
library(readxl)        
library(ggplot2)      
library(dplyr)         
library(tidyr)         


# Setting Up Working Directories and Files
setwd("C:/Users/User/Desktop/Project7/")
genes_data <- read_excel("Gene_Nucleotide_Composition.xlsx")

head(genes_data)
str(genes_data)

# 3. Compare GC Content Distributions

ggplot(genes_data, aes(x = Type, y = `GC%`, fill = Type)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(
    title = "GC Content Comparison Between Gene Types",
    x = "Gene Category",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "orange"))


#Compare Nucleotide Bias

# Gather nucleotides into a long format

nuc_data <- genes_data %>%
  select(Gene_ID, Type, `A%`, `T%`, `G%`, `C%`) %>%
  pivot_longer(cols = c(`A%`, `T%`, `G%`, `C%`),
               names_to = "Nucleotide", values_to = "Percentage")

ggplot(nuc_data, aes(x = Nucleotide, y = Percentage, fill = Type)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  labs(
    title = "Nucleotide Composition in Resistant vs Housekeeping Genes",
    x = "Nucleotide Base",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "orange"))

