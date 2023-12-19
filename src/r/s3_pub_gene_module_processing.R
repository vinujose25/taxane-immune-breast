# s3_pub_gene_module_processing.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Process and validate
# 1) TIL localization based GO-Biological-Process gene modules (Grusso et al)
# 2) General gene modules representing all BPs associated with TIL localization
# 3) Control gene modules; proliferation and visual_perception BP
# 4) MCPcounter cellmarker

# Processing includes
# a) Gene overlap and correlation analysis between similar gene modules (original modules)
# b) and pooling of highly correlated gene modules using TCGA/Metabric dataset.
# The idea is to eliminate dataset or transcriptome-snaphost effect from
# - gene modules representing similar BPs.

# Validation includes
# a) High correlation between full (original) module and finher/neoadj module subset in
# -TCGA/Metabric datasets.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data (geo, finher, expr_tcga, expr_metabric)
# 2. Load, clean, subset gene-modules.
# 3. Original module gene agreement and correlation analysis and pooling if needed
# 4. Validation of processed gene modues (Full to subset correlation)
# 5. Save Robjects




# 1. Load data (geo, finher, expr_tcga, expr_metabric)
# ==============================================================================

# clinical (to update module sore and celltype estimates)
load("results/data/clin_neoadj.RData")
load("results/data/clin_finher.RData")

# expression (to subset original gene modules and validate in expr_tcga/expr_metabric)
load("results/data/expr_neoadj.RData")
load("results/data/expr_finher.RData")

# tcga/metabric from metagxbreast R-package
# Ref: inhouse project "expr_tcga-expr_metabric-metagxbreast" (https://osf.io/jb8za/)
load("data/tcga.RData")
load("data/metabric.RData")

dim(expr_neoadj) # 9184 1500
dim(expr_finher) # 3350 301

dim(tcga$expr) # [1] 19405  1074
dim(metabric$expr) # [1] 24924  2115


# Convert expression matrix to samples x genes (for module score algorithm)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
expr_neoadj <- expr_neoadj %>% t_tibble(names_x_desc = "Sample_id")
expr_finher <- expr_finher %>% t_tibble(names_x_desc = "Sample_id")
expr_tcga <- tcga$expr %>% t_tibble(names_x_desc = "Sample_id")
expr_metabric <- metabric$expr %>% t_tibble(names_x_desc = "Sample_id")



# Discard genes with atleast one NA expression values in any samples
# (NAs will creat problems in module score algorithm)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

expr_neoadj %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE
# 9185
# expr_neoadj: No genes with NAs, No genes to discard

expr_finher %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE
# 3351
# expr_finher: No genes with NAs, No genes to discard

expr_tcga %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE
# 19406
# expr_tcga: No genes with NAs, No genes to discard

expr_metabric %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE  TRUE
# 24918     7
# expr_metabric: 7 genes with NAs, 7 genes to discard
idx <- expr_metabric %>% purrr::map_lgl(~!any(is.na(.x))) %>% which()
expr_metabric <- expr_metabric[,idx]

#
# ==============================================================================




# 2. Load, clean, subset gene-modules.
# ==============================================================================


# gene modules (vectorised)
module_consolidated <- read_tsv(
  file = "data/gene_modules/Gene_modules_consolidated_clean_ncbi_hugo.tsv"
)


# clean gene modules
module_consolidated <- module_consolidated %>%
  dplyr::mutate(

    # Set non-id/missing values to NAs

    Ncbi_gene_id = if_else(
      Ncbi_gene_id == "-",
      NA_character_,
      Ncbi_gene_id
    ),

    Hugo_gene_symbol = if_else(
      Hugo_gene_symbol == "-",
      NA_character_,
      Hugo_gene_symbol
    )

    # Direction = if_else(Coefficient < 0 , -1, 1)
    # Error: not all > 0 is positive; eg Sotiriou 2006
    # This code is moved to in-house code preparing gene modules.
  )


# vector to list
module_list <- purrr::map(
  module_consolidated$Module_id %>% unique(),
  function(nme, module_consolidated){

    module_consolidated %>%
      dplyr::filter(Module_id == nme & (!is.na(Ncbi_gene_id))) %>%
      dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE) %>%
      dplyr::select(Ncbi_gene_id, Hugo_gene_symbol, Direction) %>%
      dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id)) # to match expr gene-ids

  },
  module_consolidated
)
names(module_list) <- module_consolidated$Module_id %>% unique()

# Discard tissue signature as negative control,
# as tissue signature may contains connective tissue based BPs.

# Select and order gene modules
nme <-  c(
  # TIL localization
  "Gruosso_2019_Immune.cdsig1",
  "Gruosso_2019_Interferon.edsig2",
  "Gruosso_2019_Cholesterol.edsig5",
  "Gruosso_2019_Fibrosis.cdsig3",
  # General immune
  "Hamy_2016_Immune",
  "Teschendorff_2007_Immune",
  "Yang_2018_Immune",
  # General interferon
  "Desmedt_2008_Immune",
  "Farmer_2009_Interferon.mx1",
  "Hamy_2016_Interferon",
  "Nirmal_2018_Interferon",
  # General fibrosis
  "Bergamaschi_2008_Ecm1",
  "Hamy_2016_Ecm",
  "Naba_2014_Ecmcore",
  "Triulzi_2013_Ecm",
  # General cholestrol
  "Ehmsen_2019_Cholesterol",
  "Kimbung_2016_Cholesterol",
  "Kuzu_2016_Cholesterol",
  "Simigdala_2016_Cholesterol",
  "Sorrentino_2014_Cholesterol.mevalonate",
  # MCP-counter
  "Becht_2016_T.Cells",
  "Becht_2016_Cd8.T.Cells",
  "Becht_2016_Cytotoxic.Lymphocytes",
  "Becht_2016_B.Lineage",
  "Becht_2016_Nk.Cells",
  "Becht_2016_Monocytic.Lineage",
  "Becht_2016_Myeloid.Dendritic.Cells",
  "Becht_2016_Neutrophils",
  "Becht_2016_Endothelial.Cells",
  "Becht_2016_Fibroblasts",
  # Positive control: proliferation
  "Dai_2005_Proliferation",
  "Desmedt_2008_Proliferation",
  "Farmer_2009_Proliferation.tpx2",
  "Lundberg_2017_Proliferation",
  "Nirmal_2018_Proliferation",
  "Yang_2018_Proliferation",
  "Sotiriou_2006_Proliferation",
  # Negative control:visual perception / tissue / Behaviour signatures
  # "GO_0007601_Visual.perception",
  "Dezso_2008_Retina",
  "Dezso_2008_Tonsil",
  "Dezso_2008_Small.Intestine", "Dezso_2008_Thymus",
  "Dezso_2008_Mammary.Gland", "Dezso_2008_Brain",
  "Dezso_2008_Prostate", "Dezso_2008_Ovary",
  "Dezso_2008_Placenta", "Dezso_2008_Kidney",
  "Dezso_2008_Colon", "Dezso_2008_Fetal.Liver",
  "Dezso_2008_Fetal.Kidney", "Dezso_2008_Thyroid",
  "Dezso_2008_Lung", "Dezso_2008_Trachea",
  "Dezso_2008_Spinal.Cord", "Dezso_2008_Spleen",
  "Dezso_2008_Liver", "Dezso_2008_Uterus",
  "Dezso_2008_Skeletal.muscle", "Dezso_2008_Adrenal.Gland",
  "Dezso_2008_Testis", "Dezso_2008_Bone.Marrow",
  "Dezso_2008_Fetal.Brain",
  "Dezso_2008_Pbls", "Dezso_2008_Skin",
  "Dezso_2008_Salivary.Gland", "Dezso_2008_Pancreas",
  "Dezso_2008_Fetal.Thymus", "Dezso_2008_Heart",

  "MSigDB_Abnormal.aggressive.impulsive.or.violent.behavior",
  "MSigDB_Abnormal.consumption.behavior",
  "MSigDB_Abnormal.cry",
  "MSigDB_Abnormal.drinking.behavior",
  "MSigDB_Abnormal.eating.behavior",
  "MSigDB_Abnormal.emotion.affect.behavior",
  "MSigDB_Abnormal.fear.anxiety.related.behavior",
  "MSigDB_Abnormal.social.behavior",
  "MSigDB_Addictive.behavior",
  "MSigDB_Aggressive.behavior",
  "MSigDB_Autistic.behavior",
  "MSigDB_Inappropriate.behavior",
  "MSigDB_Inappropriate.crying",
  "MSigDB_Inappropriate.laughter",
  "MSigDB_Inappropriate.sexual.behavior",
  "MSigDB_Irritability",
  "MSigDB_Low.self.esteem",
  "MSigDB_Mental.deterioration",
  "MSigDB_Mood.swings",
  "MSigDB_Obsessive.compulsive.behavior",
  "MSigDB_Personality.changes",
  "MSigDB_Personality.disorder",
  "MSigDB_Repetitive.compulsive.behavior",
  "MSigDB_Restrictive.behavior",
  "MSigDB_Self.injurious.behavior"
  )

module_list <- module_list[nme]

# rename gene modules
nme <-  c(
  # TIL localization
  "Gruosso2019_Immune",
  "Gruosso2019_Interferon",
  "Gruosso2019_Cholesterol",
  "Gruosso2019_Fibrosis",
  # General immune
  "Hamy2016_Immune",
  "Teschendorff2007_Immune",
  "Yang2018_Immune",
  # General interferon
  "Desmedt2008_STAT1",
  "Farmer2009_MX1",
  "Hamy2016_Interferon",
  "Nirmal2018_Interferon",
  # General fibrosis
  "Bergamaschi2008_Ecm1",
  "Hamy2016_Ecm",
  "Naba2014_Ecmcore",
  "Triulzi2013_Ecm",
  # General cholestrol
  "Ehmsen2019_Chol",
  "Kimbung2016_Chol",
  "Kuzu2016_Chol",
  "Simigdala2016_Chol",
  "Sorrentino2014_Chol",
  # MCP-counter
  "Tcell",
  "CD8.Tcell",
  "Cyto.Lymphocyte",
  "B.Lineage",
  "NK.Cells",
  "Monocytic.Lineage",
  "Myeloid.Dendritic",
  "Neutrophils",
  "Endothelial",
  "Fibroblasts",
  # Positive control: proliferation
  "Dai2005_Prolif",
  "Desmedt2008_AURKA",
  "Farmer2009_TPX2",
  "Lundberg2017_Prolif",
  "Nirmal2018_Prolif",
  "Yang2018_Prolif",
  "Sotiriou2006_Prolif",
  # Negative control: visual perception / tissue signatures
  # "Visual.perception",
  "Retina",
  "Tonsil",
  "Small.Intestine", "Thymus",
  "Mammary.Gland", "Brain",
  "Prostate", "Ovary",
  "Placenta", "Kidney",
  "Colon", "Fetal.Liver",
  "Fetal.Kidney", "Thyroid",
  "Lung", "Trachea",
  "Spinal.Cord", "Spleen",
  "Liver", "Uterus",
  "Skeletal.muscle", "Adrenal.Gland",
  "Testis", "Bone.Marrow",
  "Fetal.Brain",
  "Pbls", "Skin",
  "Salivary.Gland", "Pancreas",
  "Fetal.Thymus", "Heart",

  "Abnormal.aggressive",
  "Abnormal.consumption",
  "Abnormal.cry",
  "Abnormal.drinking",
  "Abnormal.eating",
  "Abnormal.emotion.affect",
  "Abnormal.fear.anxiety",
  "Abnormal.social",
  "Addictive",
  "Aggressive",
  "Autistic",
  "Inappropriate.behavior",
  "Inappropriate.crying",
  "Inappropriate.laughter",
  "Inappropriate.sexual",
  "Irritability",
  "Low.self.esteem",
  "Mental.deterioration",
  "Mood.swings",
  "Obsessive.compulsive",
  "Personality.changes",
  "Personality.disorder",
  "Repetitive.compulsive",
  "Restrictive",
  "Self.injurious"
  )

# renaming module list
names(module_list) <- nme


# # append tilsig versions
# # >>>>>>>>>>>>>>>>>>>>>>
#
# module_list <- c(tilsig, module_list)
#
#
#
# # Neoadj module subset
# module_list_neoadj <- purrr::map(
#   module_list,
#   function(sig, genes){
#     sig %>%
#       dplyr::filter(Ncbi_gene_id %in% all_of(genes))
#     # Gene-ids are already formatted to ncbi_*
#   },
#   genes = names(expr_neoadj)[-1]
# )
#
#
# # finher module subset
# module_list_finher <- purrr::map(
#   module_list,
#   function(sig, genes){
#     sig %>%
#       dplyr::filter(Ncbi_gene_id %in% all_of(genes))
#     # Gene-ids are already formatted to ncbi_*
#   },
#   genes = names(expr_finher)[-1]
# )

#
# ==============================================================================




# 3. Original module gene agreement and correlation analysis, pooling (if needed),
#   and subsetting according to finher and neoadj genes
# ==============================================================================


# tilsig_bp geneset agreement plot

module_overlap <- module_gene_overlap(lst = module_list)

module_annot <-  module_overlap %>%
  dplyr::select(Module_name, Module_name2) %>%
  dplyr::mutate(
    Module_group = purrr::map_chr(
      Module_name,
      ~case_when(
        .x %in% c("Gruosso2019_Immune",
                  "Gruosso2019_Interferon",
                  "Gruosso2019_Cholesterol",
                  "Gruosso2019_Fibrosis") ~ "TIL-localization",
        .x %in% c("Hamy2016_Immune",
                  "Teschendorff2007_Immune",
                  "Yang2018_Immune") ~ "General-immune",
        .x %in% c( "Desmedt2008_STAT1",
                   "Farmer2009_MX1",
                   "Hamy2016_Interferon",
                   "Nirmal2018_Interferon") ~ "Genereal-interferon",
        .x %in% c("Bergamaschi2008_Ecm1",
                  "Hamy2016_Ecm",
                  "Naba2014_Ecmcore",
                  "Triulzi2013_Ecm") ~ "General-ECM",
        .x %in% c("Ehmsen2019_Chol",
                  "Kimbung2016_Chol",
                  "Kuzu2016_Chol",
                  "Simigdala2016_Chol",
                  "Sorrentino2014_Chol") ~ "General-cholesterol",
        .x %in% c("Tcell",
                  "CD8.Tcell",
                  "Cyto.Lymphocyte",
                  "B.Lineage",
                  "NK.Cells",
                  "Monocytic.Lineage",
                  "Myeloid.Dendritic",
                  "Neutrophils",
                  "Endothelial",
                  "Fibroblasts") ~ "MCP-counter",
        .x %in% c("Dai2005_Prolif",
                  "Desmedt2008_AURKA",
                  "Farmer2009_TPX2",
                  "Lundberg2017_Prolif",
                  "Nirmal2018_Prolif",
                  "Yang2018_Prolif",
                  "Sotiriou2006_Prolif") ~ "Positive-control",
        .x %in% c( "Visual.perception",
                   "Retina",
                   "Tonsil",
                   "Small.Intestine", "Thymus",
                   "Mammary.Gland", "Brain",
                   "Prostate", "Ovary",
                   "Placenta", "Kidney",
                   "Colon", "Fetal.Liver",
                   "Fetal.Kidney", "Thyroid",
                   "Lung", "Trachea",
                   "Spinal.Cord", "Spleen",
                   "Liver", "Uterus",
                   "Skeletal.muscle", "Adrenal.Gland",
                   "Testis", "Bone.Marrow",
                   "Fetal.Brain",
                   "Pbls", "Skin",
                   "Salivary.Gland", "Pancreas",
                   "Fetal.Thymus", "Heart",

                   "Abnormal.aggressive",
                   "Abnormal.consumption",
                   "Abnormal.cry",
                   "Abnormal.drinking",
                   "Abnormal.eating",
                   "Abnormal.emotion.affect",
                   "Abnormal.fear.anxiety",
                   "Abnormal.social",
                   "Addictive",
                   "Aggressive",
                   "Autistic",
                   "Inappropriate.behavior",
                   "Inappropriate.crying",
                   "Inappropriate.laughter",
                   "Inappropriate.sexual",
                   "Irritability",
                   "Low.self.esteem",
                   "Mental.deterioration",
                   "Mood.swings",
                   "Obsessive.compulsive",
                   "Personality.changes",
                   "Personality.disorder",
                   "Repetitive.compulsive",
                   "Restrictive",
                   "Self.injurious"

                   ) ~ "Negative-control",
        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(
        .x == "TIL-localization" ~ "tomato",
        .x == "General-immune" ~ "forestgreen",
        .x == "Genereal-interferon" ~ "purple",
        .x == "General-ECM" ~ "orange",
        .x == "General-cholesterol" ~ "deeppink",
        .x == "MCP-counter" ~ "steelblue",
        .x == "Positive-control" ~ "brown",
        .x == "Negative-control" ~ "gray20",
        TRUE ~ "Error"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  module_group_hline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 2.5,
  ticksize.y = 2.5,
  angle = 90,
  axis.text.size = 5.5,
  axis.text.x.vjust = .5,
  plotarea_text_size = 1

)

pdf(file = str_c(out_figures,"Module_pub_gene_overlap_proportion2.pdf"),
    width = 7.5, height = 9)
print(p_agreement) # +  guides(color = "none", fill = "none"))
dev.off()



# Module-score: Metabric

score <- get_module_score_2(
  x = expr_metabric,
  module_list = module_list,
  by = "Ncbi_gene_id"
)


p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  module_group_hline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 2.5,
  ticksize.y = 2.5,
  angle = 90,
  axis.text.size = 5.5,
  axis.text.x.vjust = .5,
  plotarea_text_size = 1
)


pdf(file = str_c(out_figures,"Module_pub_correlation_metabric2.pdf"),
    width = 7.5, height = 9)
print(p_correlation) # +  guides(color = "none", fill = "none"))
dev.off()


# Module-score: TCGA

score <- get_module_score_2(
  x = expr_tcga,
  module_list = module_list,
  by = "Ncbi_gene_id"
)


p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  module_group_hline = c(0.5, 4.5, 7.5, 11.5, 15.5, 20.5, 30.5, 37.5, 93.5),
  line_col = "gray20",
  ticksize.x = 2.5,
  ticksize.y = 2.5,
  angle = 90,
  axis.text.size = 5.5,
  axis.text.x.vjust = .5,
  plotarea_text_size = 1
)


pdf(file = str_c(out_figures,"Module_pub_correlation_tcga2.pdf"),
    width = 7.5, height = 9)
print(p_correlation )#+  guides(color = "none", fill = "none"))
dev.off()




# Negative control module selection
# #################################

# TCGA

score <- get_module_score_2(
  x = expr_tcga,
  module_list = module_list,
  by = "Ncbi_gene_id"
)

cor_tcga <- score %>%
  dplyr::select(-1) %>%
  cor(method = "pearson")
diag(cor_tcga) <- 0

cor_tcga <- cor_tcga[1:37,] %>% as_tibble(rownames = "Row")

nme <- purrr::map_lgl(cor_tcga[,-1],function(x){
  all(abs(round(x, digits = 1))<=0.2)
})
nme[nme]

# Trachea Adrenal.Gland
# TRUE          TRUE



# Metabric

score <- get_module_score_2(
  x = expr_metabric,
  module_list = module_list,
  by = "Ncbi_gene_id"
)

cor_metabric <- score %>%
  dplyr::select(-1) %>%
  cor(method = "pearson")
diag(cor_metabric) <- 0

cor_metabric <- cor_metabric[1:37,] %>% as_tibble(rownames = "Row")

nme <- purrr::map_lgl(cor_metabric[,-1],function(x){
  all(abs(round(x, digits = 1))<=0.2)
})
nme[nme]

# Brain            Fetal.Liver                Thyroid        Skeletal.muscle
# TRUE                   TRUE                   TRUE                   TRUE
# Fetal.Brain             Skin Inappropriate.behavior        Low.self.esteem
# TRUE                   TRUE                   TRUE                   TRUE

# No common tissue signature which is poorly correlated (absolute correltion <= 0.2) with
# TIL-localization, General-immune/interferon/fibrosis/cholesterol, and
# MCP-counter gene modules.
# Literature review showed taht among the poorly correlated modules, Brain tissue
# signature seems most distinct from breast expression profiles.
# Hence brain tissue signature is considered as negative control.

# Negative control: Brain module from Dezso et al.2008

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Updated negative control
# Pool all poorly correlated non-breast tissue signatures
# Negative control:
# Brain, Fetal.Liver, Thyroid, Skeletal.muscle, Fetal.Brain, Skin, Trachea, and Adrenal.Gland
# tissue modules from Dezso et al.2008.





# Positive control module selection
# #################################

# All proliferation modules are correlated.
# Proliferation correaltion range: 0.7 to 1 in tcga, 0.9 to 1 in metabric

# Gene overlap proportion ranges between 0 to 0.9, reflecting the different
# dataset (transcriptomic snapshots) used to derive the gene modules.

# Hence all proliferation modules are pooled together.


# Positive control: Pooled proliferation module


# General-Immune/Interferon/Cholesterol/Fibrosis module selection
# ###############################################################

# All general modules are correlated within their BPs.
# Cholesterol correlation range: 0.5 to 0.9 in tcga, 0.7 to 0.9 in metabric
# Ecm correlation range: 0.3 to 1.0 in tcga, 0.2 to 0.9 in metabric
# Interferon correlation range: 0.6 to 1.0 in tcga, 0.5 to 1.0 in metabric
# Immune correlation range: 0.6 to 1 in tcga, 0.2 to 1.0 in metabric

# Gene overlap proportion ranges between
# 0.2 to 0.8 for cholesteros modules,
# 0.0 to 0.6 for ecm modules,
# 0.1 to 0.8 for interferon modules, and
# 0.0 to 0.4 for immune modules, reflecting the different dataset
# (transcriptomic snapshots) used to derive the gene modules.

# Hence all general modules are pooled within their respective BPs.


# Pooled general immune modules.
# Pooled general interferon modules.
# Pooled general ecm modules.
# Pooled general cholesterol modules.


# TIL-localization (Gruosso etal.) module selection
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 70% of interferon module genes are present in immune module, which is expected
# as both immune and interferon BPs are high with CD8 Tcell infiltration.
# No common genes between other modules.

# Mtrabric correlation:

# Fibrosis is weakly anti-correlated with immune/interferon/cholesterol modules as
# Fibrosis is high in Immune desert/margin restricted TIL.
# Cholesterol is weakly positively correlated with immune/intereron as cholesterol
# is high in Stroma-restricted TIL.
# Immune and interferon is highly correlated as both are high in highly inflamed tumors.

# TCGA correlation:

# Fibrosis is weakly anti-correlated with cholesterol and not correlated with immune/interferon
# A strong anti-correlation between Fibrosis and immune/interferon is expected.
# The observed no correlation could be due to platfomr effect (RNAseq).
# Cholesterol is not correlated with immune/intereron as cholesterol.
# A weak positive correlation between cholesterol and immune/intereron is expected,
# as cholesterol is aslo high in stroma restricted TIL.
# Again, the observed no correlation could be due to platfomr effect (RNAseq).
# Immune and interferon is highly correlated as both are high in fully inflamed tumors.

# Althogh immune and interferon modules are highly correlated, they are not
# pooled as did with general modules due to the uniqness of module derivation by grusso et al.


# Pooled modules
# >>>>>>>>>>>>>>

names(module_list)

module_list_merged <- c(
  # Gruosso modules intact
  module_list[c("Gruosso2019_Immune",
                "Gruosso2019_Interferon",
                "Gruosso2019_Cholesterol",
                "Gruosso2019_Fibrosis")],

  list(
    General_Immune = bind_rows(
      module_list[c("Hamy2016_Immune",
                    "Teschendorff2007_Immune",
                    "Yang2018_Immune")]),

    General_Interferon = bind_rows(
      module_list[c("Desmedt2008_STAT1",
                    "Farmer2009_MX1",
                    "Hamy2016_Interferon",
                    "Nirmal2018_Interferon")]),

    General_ECM = bind_rows(
      module_list[c("Bergamaschi2008_Ecm1",
                    "Hamy2016_Ecm",
                    "Naba2014_Ecmcore",
                    "Triulzi2013_Ecm")]),

    General_Cholesterol = bind_rows(
      module_list[c("Ehmsen2019_Chol",
                    "Kimbung2016_Chol",
                    "Kuzu2016_Chol",
                    "Simigdala2016_Chol",
                    "Sorrentino2014_Chol")]),

    General_Proliferation = bind_rows(
      module_list[c("Dai2005_Prolif",
                    "Desmedt2008_AURKA",
                    "Farmer2009_TPX2",
                    "Lundberg2017_Prolif",
                    "Nirmal2018_Prolif",
                    "Yang2018_Prolif",
                    "Sotiriou2006_Prolif")]),

    Tissue_Non.breast = bind_rows(
    module_list[c(
      "Brain",
      "Fetal.Liver",
      "Thyroid",
      "Skeletal.muscle",
      "Fetal.Brain",
      "Skin",
      "Trachea",
      "Adrenal.Gland"
    )]),

    Human_Behaviour = bind_rows(
      module_list[c(
        # "Abnormal.aggressive",
        # "Abnormal.consumption",
        # "Abnormal.cry",
        # "Abnormal.drinking",
        # "Abnormal.eating",
        # "Abnormal.emotion.affect",
        # "Abnormal.fear.anxiety",
        # "Abnormal.social",
        # "Addictive",
        # "Aggressive",
        # "Autistic",
        "Inappropriate.behavior",
        # "Inappropriate.crying",
        # "Inappropriate.laughter",
        # "Inappropriate.sexual",
        # "Irritability",
        "Low.self.esteem"
        # "Mental.deterioration",
        # "Mood.swings",
        # "Obsessive.compulsive",
        # "Personality.changes",
        # "Personality.disorder",
        # "Repetitive.compulsive",
        # "Restrictive",
        # "Self.injurious"
      )]
      )
  ),

  # module_list[c("Brain")],

  # MCP counter modules are included for completness
  # Celltype estimation is done seperatly
  module_list[c("Tcell", "CD8.Tcell", "Cyto.Lymphocyte", "B.Lineage",
                "NK.Cells", "Monocytic.Lineage", "Myeloid.Dendritic",
                "Neutrophils", "Endothelial", "Fibroblasts")]
)


module_list_merged <- purrr::map(module_list_merged,
                                 ~.x %>% dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE))



# Module subset
# >>>>>>>>>>>>>


# Neoadj module subset
module_list_merged_neoadj <- purrr::map(
  module_list_merged,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_neoadj)[-1]
)


# finher module subset
module_list_merged_finher <- purrr::map(
  module_list_merged,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_finher)[-1]
)


# Module size
tibble(
  Module_name = names(module_list_merged),
  Original = purrr::map(module_list_merged,nrow) %>% unlist(),
  Finher = purrr::map(module_list_merged_finher,nrow) %>% unlist(),
  Neoadj = purrr::map(module_list_merged_neoadj,nrow) %>% unlist()
) %>%
  as.data.frame()
#                Module_name Original Finher Neoadj
# 1       Gruosso2019_Immune      716    212    493
# 2   Gruosso2019_Interferon      111     36     80
# 3  Gruosso2019_Cholesterol       49     11     41
# 4     Gruosso2019_Fibrosis      316     95    231
# 5           General_Immune       43     19     26
# 6       General_Interferon      151     55    104
# 7              General_ECM      146     82    105
# 8      General_Cholesterol       27     13     20
# 9    General_Proliferation      710    172    490
# 10       Tissue_Non.breast      259     12    103
# 11         Human_Behaviour       63     13     46
# 12                   Tcell       16      1     10
# 13               CD8.Tcell        1      0      0
# 14         Cyto.Lymphocyte        9      0      3
# 15               B.Lineage        9      2      8
# 16                NK.Cells        9      0      4
# 17       Monocytic.Lineage        7      2      5
# 18       Myeloid.Dendritic        6      0      3
# 19             Neutrophils       15      3     10
# 20             Endothelial       33      0     22
# 21             Fibroblasts        8      7      6


#
# ==============================================================================




# 4. Validation of processed gene modues (Full to subset correlation)
# ==============================================================================


# Neoadj module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# validation on expr_tcga
validation_neoadj_expr_tcga <- validate_gene_modules(
  module_full_list = module_list_merged,
  module_subset_list = module_list_merged_neoadj,
  validation_data = expr_tcga)


# validation on expr_metabric
validation_neoadj_expr_metabric <- validate_gene_modules(
  module_full_list = module_list_merged,
  module_subset_list = module_list_merged_neoadj,
  validation_data = expr_metabric)

# validation stat
validation_neoadj_stat <- bind_cols(
  validation_neoadj_expr_tcga %>%
    dplyr::select(Module_id, Module_full_size, Module_subset_size,
                  Module_subset_size_percent),
  validation_neoadj_expr_tcga %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_tcga_",.x)),
  validation_neoadj_expr_metabric %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_metabric_",.x))
) %>%
  dplyr::mutate(
    Valid_modules = ((expr_tcga_Pearson %>% round(digits = 1)) >= 0.8 &
                       (expr_metabric_Pearson %>% round(digits = 1)) >= 0.8),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )

# Relevant valid modules
validation_neoadj_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::mutate(out = str_c(Module_id, "(N:", Module_subset_size,",",
                            round(Module_subset_size_percent, digits = 1),"%,P:",
                            round(expr_tcga_Pearson, digits = 2),")")) %>%
  dplyr::select(out) %>%
  tibble::deframe()

# [1] "Gruosso2019_Immune(N:493,68.9%,P:1)"        "Gruosso2019_Interferon(N:80,72.1%,P:0.99)"
# [3] "Gruosso2019_Cholesterol(N:41,83.7%,P:0.99)" "Gruosso2019_Fibrosis(N:231,73.1%,P:1)"

# [5] "General_Immune(N:26,60.5%,P:1)"             "General_Interferon(N:104,68.9%,P:0.98)"
# [7] "General_ECM(N:105,71.9%,P:1)"               "General_Cholesterol(N:20,74.1%,P:0.87)"

# [9] "General_Proliferation(N:490,69%,P:1)"       "Tissue_Non.breast(N:103,39.8%,P:0.86)"
# [11] "Human_Behaviour(N:46,73%,P:0.96)"           "Tcell(N:10,62.5%,P:1)"

# [13] "Cyto.Lymphocyte(N:3,33.3%,P:0.94)"          "B.Lineage(N:8,88.9%,P:1)"
# [15] "NK.Cells(N:4,44.4%,P:0.96)"                 "Monocytic.Lineage(N:5,71.4%,P:0.99)"
# [17] "Myeloid.Dendritic(N:3,50%,P:0.99)"          "Neutrophils(N:10,66.7%,P:0.85)"
# [19] "Endothelial(N:22,66.7%,P:0.99)"             "Fibroblasts(N:6,75%,P:0.99)"


# FinHER module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# validation on expr_tcga
validation_finher_expr_tcga <- validate_gene_modules(
  module_full_list = module_list_merged,
  module_subset_list = module_list_merged_finher,
  validation_data = expr_tcga)

# validation in expr_metabric
validation_finher_expr_metabric <- validate_gene_modules(
  module_full_list = module_list_merged,
  module_subset_list = module_list_merged_finher,
  validation_data = expr_metabric)

# validation stat
validation_finher_stat <- bind_cols(
  validation_finher_expr_tcga %>%
    dplyr::select(Module_id, Module_full_size, Module_subset_size,
                  Module_subset_size_percent),
  validation_finher_expr_tcga %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_tcga_",.x)),
  validation_finher_expr_metabric %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_metabric_",.x))
) %>%
  dplyr::mutate(
    Valid_modules = ((expr_tcga_Pearson %>% round(digits = 1)) >= 0.8 &
                       (expr_metabric_Pearson %>% round(digits = 1)) >= 0.8),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )


# valid modules
validation_finher_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::mutate(out = str_c(Module_id, "(N:", Module_subset_size,",",
                            round(Module_subset_size_percent, digits = 1),"%,P:",
                            round(expr_tcga_Pearson, digits = 2),")")) %>%
  dplyr::select(out) %>%
  tibble::deframe()
# [1] "Gruosso2019_Immune(N:212,29.6%,P:0.94)"     "Gruosso2019_Interferon(N:36,32.4%,P:0.97)"
# [3] "Gruosso2019_Cholesterol(N:11,22.4%,P:0.91)" "Gruosso2019_Fibrosis(N:95,30.1%,P:0.95)"

# [5] "General_Immune(N:19,44.2%,P:0.99)"          "General_Interferon(N:55,36.4%,P:0.83)"
# [7] "General_ECM(N:82,56.2%,P:0.98)"             "General_Cholesterol(N:13,48.1%,P:0.91)"

# [9] "General_Proliferation(N:172,24.2%,P:0.94)"  "Tcell(N:1,6.2%,P:0.93)"
# [11] "Monocytic.Lineage(N:2,28.6%,P:0.84)"        "Fibroblasts(N:7,87.5%,P:0.99)"


# The negative control Tissue_Non.breast and Human_Behavour modules did not pass
# validation process.
# However, it is decided to keep the both modules as negative control, as Tissue_Non.breast and Human_Behavour
# signatures were not supposed to present in breast normal/tumor tissue.
# Further the number of genes present in FinHER dataset (12/13) is very small compared
# to Neadj dataset (n=103/46).



# # cleaning
# rm(validation_neoadj_expr_tcga, validation_neoadj_expr_metabric)
# rm(validation_finher_expr_tcga, validation_finher_expr_metabric)
# rm(expr_tcga, expr_metabric)


#
# ==============================================================================




# 5. Save Robjects
# ==============================================================================


# Robject created:
save(module_list_merged, file = str_c(out_data, "module_list_merged.RData"))
save(module_list_merged_neoadj, file = str_c(out_data, "module_list_merged_neoadj.RData"))
save(module_list_merged_finher, file = str_c(out_data, "module_list_merged_finher.RData"))

save(validation_neoadj_stat, file = str_c(out_data, "validation_neoadj_stat.RData"))
save(validation_finher_stat, file = str_c(out_data, "validation_finher_stat.RData"))


#
# ==============================================================================


