
# This README details the summary and purpose of each script and their ordinal positions
# ==============================================================================

Initialize.R: Initialize the folder structure for the project and load R packages.

functions.R: Support functions for other scripts. Sourced under the "Library" section.

README.txt: This file.

s1_pooled_dataset_cleaning.R: Filter pooled GEO dataset based on the criteria specified in the manuscript.

s2_finher_cleaning.R: Prepare FinHER gene expression data and clinical data for analysis.

s3_pub_gene_module_processing.R: Public gene module processing, gene agreement, subsetting, correlation analysis, and pooling.

s4_finher_tilsig.R: De novo TILsig derivation from FinHER TIL counts.

s5_module_score_and_celltype_estimation.R: Compute module scores and celltype estimates.

s6.1_finher_tn.R: Prognosis and MTA*Immune*Survival analysis in FinHER TNBC.

s6.2_supp_firth_finher_tn.R: Checking if Firth penalization reduces the wide confidence interval of estimates derived from FinHER TNBC.

s7.1_finher_her2.R: Prognosis and MTA*Immune*Survival analysis in FinHER HER2+BC.

s7.2_supp_finher_her2_ph_fixed.R: PH fixing for PH failed models using a step function of time approach.

s7.3_supp_firth_finher_her2.R: Checking if Firth penalization reduces the wide confidence interval of estimates derived from FinHER HER2+BC.

s8_finher_plots_mta_paper.R: Generated analysis plots from the FinHER dataset.

s9.1_geo_exploring_interaction _finneo.R: MTA*Immune*pCR analysis in GEO dataset using FinGEO module subsets.

s9.2_geo_exploring_interaction _neo.R: MTA*Immune*pCR analysis in GEO dataset using GEO module subsets.

s9.3_supp_firth_geo_exploring_interaction _finneo.R: Checking if Firth penalization reduces the wide confidence interval of interaction estimates derived from GEO.

s10.1_geo_exploring_prognosis_finneo.R: Prognosis analysis in GEO dataset using FinGEO module subsets.

s10.2_geo_exploring_prognosis_neo.R: Prognosis analysis in GEO dataset using GEO module subsets.

s10.3_supp_firth_geo_exploring_prognosis_finneo.R: Checking if Firth penalization reduces the wide confidence interval of prognosis estimates derived from GEO.

S11.1_geo_plots_finneo.R: Generate analysis plots from the GEO dataset using FinGEO module subsets.

S11.2_geo_plots_neo.R: Generate analysis plots from the GEO dataset using GEO module subsets.

S12_manuscript_tables.R: Generate tables for manuscripts.

#
# ==============================================================================
