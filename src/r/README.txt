# ==============================================================================
# This README details the following
# 1. Summary of the purpose of each scripts and their ordinal positions
# 2. Details of each scripts and the main R-objects/results created within
# ==============================================================================


# 1. Summary of the pupose of each scripts and their ordinal positions
# ==============================================================================

# Note:
# The keyword "manual" in script filename implies that some part of the script
# requires manual cleaning/formatting/curation etc.


# @@@@ Updated !!!
# main.R:
#       This script !!!!!!!
#       Main script from which other scripts are managed.
# functions.R:
#       Supporting functions for main and other scripts.
#       Sourced under "Library" section.
# s1_manual_pooled_dataset_cleaning.R:
#       Filter pooled dataset based on the criteria specified in the text.
# s2_manual_finher_dataset_cleaning.R:
#       Clean finher dataset based on the criteria specified in the text.
# s3_pub_gene_module_processing.R # gruossoTILloc_generalTIL_control.R
# s4_denovo_gene_module_generation.R

# !!!!!!!!!!!!!!
#       s3_TIL_localization_modules.R
#       s4_TIL_modules.R
#       s5_General_modules_representing_TIL_associated_BP.R
#       s6_Celltype_estimation.R
#       s7_Control_modules.R
# !!!!!!!!!!!!

# s3_finher_tilsig.R: (! RENAME finher_tilsig.R)
#       De-novo TILsig generation.
# s4_module_score_and_celltype_estimation.R: (! RENAME s2_module_score_and_celltype_estimation.R)
#       Clean gene-modules of biological processes and celltype.
#       Validate the cleaned gene-modules using TCGA/METABRIC(MeatGxBreast) datasets.
#       Jaccard-index, Co-correlation heatmaps.
#       Justify redundant signature use or merge signatures.
#       Estimate module and celltype scores.
# s3_exploring_prognosis.R:
#       Prognosis and heterogenity testing.
# s4_exploring_interaction.R:
#       Interaction and heterogenity testing.
# s5_figures_tables_data.R:
#       Generate journal specific figures, tables and dataset.
# manuscript.Rmd:
#       Rmarkdown version of final manuscript.


# @@@@ Backup old !!!
# main.R:
#       This script !!!!!!!
#       Main script from which other scripts are managed.
# functions.R:
#       Supporting functions for main and other scripts.
#       Sourced under "Library" section.
# s1_manual_pooled_dataset_cleaning.R:
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




# 2. Details of each scripts and the main Robjects/results created within
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

