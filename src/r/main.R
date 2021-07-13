# main.R

# backup

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# The idea of this main script is to manage all other scrips from here.

# !!!!!!!!!!!!!!!!!!!!!!!! review !!!!!!!!!!!!!!!
# However, in this data pooling project no script can run without manual
# intervention. Scripts which needs manual interventions are prefixed with
# keyword "s*_manual". Hence in effect this main script functions as a place
# holder where in the details of other sripts and how they are connected are
# documented.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions
# 3. Summary of the purpose of each scripts and their ordinal positions
# 4. Details of each scripts and the main R-objects/results created within




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/tables")){
  dir.create(path = "results/tables", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/figures")){
  dir.create(path = "results/figures", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/data")){
  dir.create(path = "results/data", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/documents")){
  dir.create(path = "results/documents", recursive = TRUE, mode = "0700")
}


out_tables <- "results/tables/"
out_figures <- "results/figures/"
out_data <- "results/data/"
out_documents <- "results/documents/"


# # Output directory
# if(!dir.exists("results/main")){
#   dir.create(path = "results/main", recursive = TRUE, mode = "0700")
# }
# outdir <- "results/main/"

# # Tmp directory
# if(!dir.exists("results/tmp")){
#   dir.create(path = "results/tmp", recursive = TRUE, mode = "0700")
# }
# tmpdir <- "results/tmp/"

#
# ==============================================================================




# 2. Library: Loading of depended R packages and local functions
# ==============================================================================

# Public
library(tidyverse)
library(ggstance) # position_dodgev() # vertical dodging
library(MCPcounter) # MCPcounter.estimate: produce cell type estimates
# library(hablar) # retype: Transforms all elements into simple classes
library(rmarkdown) # render: To render the Rmarkdown file into specific output format
# library(openxlsx) #

library(writexl) # write_xlsx: wrirte xlsx files in multiple sheets
# library(readxl) # read_excel: Read xls and xlsx files
library(altmeta) # metahet: Meta-Analysis Heterogeneity Measures
# library(meta) # metabin: Meta-analysis of binary outcome data
# library(MetaGxBreast)
library(ggpubr) # ggarrange: Arrange nultiple ggplots
                # annotate_figure: Annotate figures including: i) ggplots, ii) arranged ggplots from ggarrange(), grid.arrange() and plot_grid().
library(pdftools) # qpdf::pdf_subset: Split, Combine and Compress PDF Files
library(genefu) # rescale() based on quantiles
library(lmtest) # lrtest: likilihood ratio test
library(survival) # pspline(), coxph()
library(splines) # ns()
# library(coxme)
# library(lme4)
library(survminer) # ggsurvplot()
# NOTE: surminer depends on  ggpubr package which has ggarrange() function
# Try using ggarrange to arrange ggaurvplot()
library(grid) # unit()
library(survcomp) # combine.test()
library(pROC) # auc

# library(ggsci) # scale_color_jco() Ref : https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#:~:text=In%20this%20case%2C%20we'll,that%20are%20color%2Dblind%20friendly.



# Private
source("src/r/functions.R")

#
# ==============================================================================




# 3. Summary of the pupose of each scripts and their ordinal positions
# ==============================================================================

# Note:
# The keyword "manual" in script filename implies that some part of the script
# requires manual cleaning/formatting/curation etc.

# main.R:
#       This script !!!!!!!
#       Main script from which other scripts are managed.
# functions.R:
#       Supporting functions for main and other scripts.
#       Sourced under "Library" section.
# s1_pooled_dataset_cleaning.R:
#       Filter pooled dataset based on the criteria specified in the text.
# tilsig.R
# s2_module_score_and_celltype_estimation.R:
#       Clean gene-modules of biological processes and celltype.
#       Validate the cleaned gene-modules using TCGA/METABRIC(MeatGxBreast) datasets.
#       Estimate module and celltype scores.
# s3_exploring_prognosis.R:
#       Prognosis and heterogenity testing.
# s4_exploring_interaction.R:
#       Interaction and heterogenity testing.
# s5_figures_tables_data.R:
#       Generate journal specific figures, tables and dataset.
# manuscript.Rmd:
#       Rmarkdown version of final manuscript.

#
# ==============================================================================




# 4. Details of each scripts and the main Robjects/results created within
# ==============================================================================


# s1_manual_pooled_dataset_cleaning.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("r/s1_manual_pooled_dataset_cleaning.R")

# Robject created: 1) clin_neoadj, 2) expr_neoadj
# Relevant Robjects: 1) clin_neoadj, 2) expr_neoadj




# s2_module_score_and_celltype_estimation.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# source("r/s2_module_score_and_celltype_estimation.R")

# Robject created: 1) module_consolidated, 2) module_list, 3) module_list_subset,
#                  4) module_stat
# Robject updated: 1) clin_neoadj (updated with module-score and celltype estimates)
# Relevant Robjects:  1) clin_neoadj, 2) expr_neoadj,
#                     3) module_consolidated, 4) module_list,
#                     5) module_list_subset, 6) module_stat





# s3_exploring_prognosis.R
# >>>>>>>>>>>>>>>>>>>>>>>>

# source("r/s3_exploring_prognosis.R")

# Robject created: 1) plot_data_prognosis
# Relevant Robjects:  1) clin_neoadj, 2) expr_neoadj,
#                     3) module_consolidated, 4) module_list,
#                     5) module_list_subset, 6) module_stat,
#                     7) plot_data_prognosis





# s4_exploring_interaction.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>

# source("r/s4_exploring_interaction.R")

# Robject created: 1) plot_data_interaction
# Relevant Robjects:  1) clin_neoadj, 2) expr_neoadj,
#                     3) module_consolidated, 4) module_list,
#                     5) module_list_subset, 6) module_stat,
#                     7) plot_data_prognosis, 8) plot_data_interaction



# s5_figures_tables_data.R
# >>>>>>>>>>>>>>>>>>>>>>>>

# source("r/s5_figures_tables_data.R")

# Output: Prognosis and interaction forest plots.
# Relevant Robjects:  1) clin_neoadj, 2) expr_neoadj,
#                     3) module_consolidated, 4) module_list,
#                     5) module_list_subset, 6) module_stat,
#                     7) plot_data_prognosis, 8) plot_data_interaction




# manuscript.Rmd
# >>>>>>>>>>>>>>

# # Rmarkdown version of the manuscript after review.
# render("src/r/manuscript.Rmd", output_dir = out_documents)
# # Ref: https://community.rstudio.com/t/is-it-possible-to-save-the-html-output-in-a-directory-which-is-not-the-one-where-the-rmd-file-resides/3588/5


#
# ==============================================================================


