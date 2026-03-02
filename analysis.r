# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

library(tidyverse)

# Load data here ----------------------
# Load each file with a meaningful variable name.
expr <- read_csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")
meta <- read_csv("data/GSE60450_filtered_metadata.csv")


# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(expr)
# cols = 14
# rows = 23735

## Metadata
dim(meta)
# cols = 4
# rows = 12

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?
expr_long <- expr %>% 
  pivot_longer(cols = contains("GSM"),
               names_to = "sample",
               values_to = "counts")

meta <- meta %>% 
  mutate(sample = ...1) %>% 
  select(-...1) %>% 
  relocate(sample, .before = "characteristics")

combined <- meta %>% 
  left_join(expr_long, join_by(sample)) %>% 
  mutate(gene_id = ...1) %>% 
  select(-...1) %>% 
  relocate(gene_id, .before = "gene_symbol")

#check if joined correctly
unique(combined$sample)
unique(expr_long$sample)
               

# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2
combined <- combined %>% 
  mutate(immunotype = case_when(immunophenotype == "luminal cell population" ~ "lum",
                                immunophenotype == "basal cell population" ~ "bas"),
         devestage = case_when(`developmental stage` == "18.5 day pregnancy" ~ "preg",
                               `developmental stage` == "2 day lactation" ~ "lac",
                               .default = "virg")) %>% 
  mutate(celltype = str_c(immunotype, devestage)) %>% 
  mutate(logcounts = log2(counts))

# chose 4 genes to visualize
# random 4
gene_vis <- c("Gnai3", "Pbsn", "Cdc45", "H19")

p <- combined %>%
  filter(gene_symbol %in% gene_vis) %>% 
  ggplot(aes(celltype, logcounts, colour = celltype)) +
  geom_boxplot() +
  facet_wrap(~ gene_symbol)

p


## Save the plot
### Show code for saving the plot with ggsave() or a similar function

ggsave("plots/celltype_vs_logcounts.png", p, width = 8, height = 8)
