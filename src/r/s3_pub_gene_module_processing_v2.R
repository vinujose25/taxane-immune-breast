# s3_pub_gene_module_processing.R

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Identical to s3_pub_gene_module_processing.R but with
# strict validation cutoffs.



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
  "Desmedt_2008_Immune",
  # General interferon
  "Farmer_2009_Interferon.mx1",
  "Hamy_2016_Interferon",
  "Nirmal_2018_Interferon",
  # General fibrosis
  # "Bergamaschi_2008_Ecm1", # the original paper contain ECM signature
  # representing four subtypes, it is not sure whic one is associated with TIL
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
  "Becht_2016_Fibroblasts" #,
  # # Positive control: proliferation
  # "Dai_2005_Proliferation",
  # "Desmedt_2008_Proliferation",
  # "Farmer_2009_Proliferation.tpx2",
  # "Lundberg_2017_Proliferation",
  # "Nirmal_2018_Proliferation",
  # "Yang_2018_Proliferation",
  # "Sotiriou_2006_Proliferation",
  # # Negative control:visual perception / tissue / Behaviour signatures
  # # "GO_0007601_Visual.perception",
  # "Dezso_2008_Retina",
  # "Dezso_2008_Tonsil",
  # "Dezso_2008_Small.Intestine", "Dezso_2008_Thymus",
  # "Dezso_2008_Mammary.Gland", "Dezso_2008_Brain",
  # "Dezso_2008_Prostate", "Dezso_2008_Ovary",
  # "Dezso_2008_Placenta", "Dezso_2008_Kidney",
  # "Dezso_2008_Colon", "Dezso_2008_Fetal.Liver",
  # "Dezso_2008_Fetal.Kidney", "Dezso_2008_Thyroid",
  # "Dezso_2008_Lung", "Dezso_2008_Trachea",
  # "Dezso_2008_Spinal.Cord", "Dezso_2008_Spleen",
  # "Dezso_2008_Liver", "Dezso_2008_Uterus",
  # "Dezso_2008_Skeletal.muscle", "Dezso_2008_Adrenal.Gland",
  # "Dezso_2008_Testis", "Dezso_2008_Bone.Marrow",
  # "Dezso_2008_Fetal.Brain",
  # "Dezso_2008_Pbls", "Dezso_2008_Skin",
  # "Dezso_2008_Salivary.Gland", "Dezso_2008_Pancreas",
  # "Dezso_2008_Fetal.Thymus", "Dezso_2008_Heart",
  #
  # "MSigDB_Abnormal.aggressive.impulsive.or.violent.behavior",
  # "MSigDB_Abnormal.consumption.behavior",
  # "MSigDB_Abnormal.cry",
  # "MSigDB_Abnormal.drinking.behavior",
  # "MSigDB_Abnormal.eating.behavior",
  # "MSigDB_Abnormal.emotion.affect.behavior",
  # "MSigDB_Abnormal.fear.anxiety.related.behavior",
  # "MSigDB_Abnormal.social.behavior",
  # "MSigDB_Addictive.behavior",
  # "MSigDB_Aggressive.behavior",
  # "MSigDB_Autistic.behavior",
  # "MSigDB_Inappropriate.behavior",
  # "MSigDB_Inappropriate.crying",
  # "MSigDB_Inappropriate.laughter",
  # "MSigDB_Inappropriate.sexual.behavior",
  # "MSigDB_Irritability",
  # "MSigDB_Low.self.esteem",
  # "MSigDB_Mental.deterioration",
  # "MSigDB_Mood.swings",
  # "MSigDB_Obsessive.compulsive.behavior",
  # "MSigDB_Personality.changes",
  # "MSigDB_Personality.disorder",
  # "MSigDB_Repetitive.compulsive.behavior",
  # "MSigDB_Restrictive.behavior",
  # "MSigDB_Self.injurious.behavior"
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
  "Desmedt2008_STAT1",
  # General interferon
  "Farmer2009_MX1",
  "Hamy2016_Interferon",
  "Nirmal2018_Interferon",
  # General fibrosis
  # "Bergamaschi2008_Ecm1",
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
  "Fibroblasts"#,
  # # Positive control: proliferation
  # "Dai2005_Prolif",
  # "Desmedt2008_AURKA",
  # "Farmer2009_TPX2",
  # "Lundberg2017_Prolif",
  # "Nirmal2018_Prolif",
  # "Yang2018_Prolif",
  # "Sotiriou2006_Prolif",
  # # Negative control: visual perception / tissue signatures
  # # "Visual.perception",
  # "Retina",
  # "Tonsil",
  # "Small.Intestine", "Thymus",
  # "Mammary.Gland", "Brain",
  # "Prostate", "Ovary",
  # "Placenta", "Kidney",
  # "Colon", "Fetal.Liver",
  # "Fetal.Kidney", "Thyroid",
  # "Lung", "Trachea",
  # "Spinal.Cord", "Spleen",
  # "Liver", "Uterus",
  # "Skeletal.muscle", "Adrenal.Gland",
  # "Testis", "Bone.Marrow",
  # "Fetal.Brain",
  # "Pbls", "Skin",
  # "Salivary.Gland", "Pancreas",
  # "Fetal.Thymus", "Heart",
  #
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
  # "Inappropriate.behavior",
  # "Inappropriate.crying",
  # "Inappropriate.laughter",
  # "Inappropriate.sexual",
  # "Irritability",
  # "Low.self.esteem",
  # "Mental.deterioration",
  # "Mood.swings",
  # "Obsessive.compulsive",
  # "Personality.changes",
  # "Personality.disorder",
  # "Repetitive.compulsive",
  # "Restrictive",
  # "Self.injurious"
  )

# renaming module list
names(module_list) <- nme


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
                  "Gruosso2019_Fibrosis") ~ "TIL-Localization",
        .x %in% c("Hamy2016_Immune",
                  "Teschendorff2007_Immune",
                  "Yang2018_Immune",
                  "Desmedt2008_STAT1") ~ "Immune", #"General-Immune",
        .x %in% c("Farmer2009_MX1",
                   "Hamy2016_Interferon",
                   "Nirmal2018_Interferon") ~ "Interferon", #"General-Interferon",
        .x %in% c(#"Bergamaschi2008_Ecm1",
                  "Hamy2016_Ecm",
                  "Naba2014_Ecmcore",
                  "Triulzi2013_Ecm") ~ "Fibrosis", #"General-Fibrosis",
        .x %in% c("Ehmsen2019_Chol",
                  "Kimbung2016_Chol",
                  "Kuzu2016_Chol",
                  "Simigdala2016_Chol",
                  "Sorrentino2014_Chol") ~"Cholesterol",  #"General-Cholesterol",
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
        # .x %in% c("Dai2005_Prolif",
        #           "Desmedt2008_AURKA",
        #           "Farmer2009_TPX2",
        #           "Lundberg2017_Prolif",
        #           "Nirmal2018_Prolif",
        #           "Yang2018_Prolif",
        #           "Sotiriou2006_Prolif") ~ "Positive-Control",
        # .x %in% c( "Visual.perception",
        #            "Retina",
        #            "Tonsil",
        #            "Small.Intestine", "Thymus",
        #            "Mammary.Gland", "Brain",
        #            "Prostate", "Ovary",
        #            "Placenta", "Kidney",
        #            "Colon", "Fetal.Liver",
        #            "Fetal.Kidney", "Thyroid",
        #            "Lung", "Trachea",
        #            "Spinal.Cord", "Spleen",
        #            "Liver", "Uterus",
        #            "Skeletal.muscle", "Adrenal.Gland",
        #            "Testis", "Bone.Marrow",
        #            "Fetal.Brain",
        #            "Pbls", "Skin",
        #            "Salivary.Gland", "Pancreas",
        #            "Fetal.Thymus", "Heart",
        #
        #            "Abnormal.aggressive",
        #            "Abnormal.consumption",
        #            "Abnormal.cry",
        #            "Abnormal.drinking",
        #            "Abnormal.eating",
        #            "Abnormal.emotion.affect",
        #            "Abnormal.fear.anxiety",
        #            "Abnormal.social",
        #            "Addictive",
        #            "Aggressive",
        #            "Autistic",
        #            "Inappropriate.behavior",
        #            "Inappropriate.crying",
        #            "Inappropriate.laughter",
        #            "Inappropriate.sexual",
        #            "Irritability",
        #            "Low.self.esteem",
        #            "Mental.deterioration",
        #            "Mood.swings",
        #            "Obsessive.compulsive",
        #            "Personality.changes",
        #            "Personality.disorder",
        #            "Repetitive.compulsive",
        #            "Restrictive",
        #            "Self.injurious"
        #
        #            ) ~ "Negative-Control",
        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(

        # color palatte reference
        # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9

        .x == "TIL-Localization" ~ "#e41a1c", # red "tomato",
        .x == "Immune" ~ "#377eb8", # blue "forestgreen",
        .x == "Interferon" ~ "#4daf4a", # green "purple",
        .x == "Fibrosis" ~ "#984ea3",  #"purple",
        .x == "Cholesterol" ~ "#ff7f00", #"orange",
        .x == "MCP-counter" ~ "#a65628", # brown "steelblue",

        # .x == "General-Immune" ~ "forestgreen",
        # .x == "General-Interferon" ~ "purple",
        # .x == "General-ECM" ~ "orange",
        # .x == "General-Cholesterol" ~ "deeppink",
        # .x == "Positive-Control" ~ "brown",
        # .x == "Negative-Control" ~ "gray20",

        TRUE ~ "Error"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  module_group_hline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 4.5,
  ticksize.y = 3.5,
  angle = 90,
  axis.text.size = 7,
  axis.text.x.vjust = .5,
  plotarea_text_size =  1.75

)

pdf(file = str_c(out_figures,"Module_pub_gene_overlap_proportion_v2.pdf"),
    width = 5.5, height = 6)
print(
  p_agreement +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(face = "bold"))
)
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
  module_group_vline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  module_group_hline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 4.5,
  ticksize.y = 3.5,
  angle = 90,
  axis.text.size = 7,
  axis.text.x.vjust = .5,
  plotarea_text_size = 1.75
)


pdf(file = str_c(out_figures,"Module_pub_correlation_metabric_v2.pdf"),
    width = 5.5, height = 6)
print(
  p_correlation +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(face = "bold"))
  )
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
  module_group_vline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  module_group_hline = c(0.5, 4.5, 8.5, 11.5, 14.5, 19.5, 29.5),
  line_col = "gray20",
  ticksize.x = 4.5,
  ticksize.y = 3.5,
  angle = 90,
  axis.text.size = 7,
  axis.text.x.vjust = .5,
  plotarea_text_size = 1.75
)


pdf(file = str_c(out_figures,"Module_pub_correlation_tcga_v2.pdf"),
    width = 5.5, height = 6)
print(
  p_correlation +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(face = "bold"))
)
dev.off()



# Module subset
# >>>>>>>>>>>>>


# Finher+ Neoadj common module subset
module_list_finneo <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = intersect(names(expr_finher)[-1], names(expr_neoadj)[-1])
)


# Neoadj module subset
module_list_neo <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_neoadj)[-1]
)



# Module size
module_size_filter <- tibble(
  Module_name = names(module_list),
  Original = purrr::map(module_list,nrow) %>% unlist(),

  Finneo = purrr::map(module_list_finneo,nrow) %>% unlist(),
  Finneo_size_select = (Finneo>=3),

  Neo = purrr::map(module_list_neo,nrow) %>% unlist(),
  Neo_size_select = (Neo>=3)

) %>%
  as.data.frame()

module_size_filter %>%
  dplyr::filter(Finneo_size_select | Neo_size_select)
#                Module_name Original Finneo Finneo_size_select Neo Neo_size_select
# 1       Gruosso2019_Immune      716    144               TRUE 493            TRUE
# 2   Gruosso2019_Interferon      111     25               TRUE  80            TRUE
# 3  Gruosso2019_Cholesterol       49     10               TRUE  41            TRUE
# 4     Gruosso2019_Fibrosis      316     74               TRUE 231            TRUE
# 5          Hamy2016_Immune       26      6               TRUE  13            TRUE
# 6  Teschendorff2007_Immune        7      2              FALSE   4            TRUE
# 7          Yang2018_Immune       17      6               TRUE  13            TRUE
# 8        Desmedt2008_STAT1       95     19               TRUE  73            TRUE
# 9           Farmer2009_MX1       40     18               TRUE  30            TRUE
# 10     Hamy2016_Interferon       11      6               TRUE   8            TRUE
# 11   Nirmal2018_Interferon       66     19               TRUE  39            TRUE
# 12            Hamy2016_Ecm       21      9               TRUE  11            TRUE
# 13        Naba2014_Ecmcore       55     28               TRUE  40            TRUE
# 14         Triulzi2013_Ecm       57     25               TRUE  41            TRUE
# 15         Ehmsen2019_Chol        8      5               TRUE   7            TRUE
# 16        Kimbung2016_Chol       20      9               TRUE  16            TRUE
# 17           Kuzu2016_Chol        7      3               TRUE   5            TRUE
# 18      Simigdala2016_Chol        8      4               TRUE   7            TRUE
# 19     Sorrentino2014_Chol        9      5               TRUE   7            TRUE
# 20                   Tcell       16      1              FALSE  10            TRUE
# 21         Cyto.Lymphocyte        9      0              FALSE   3            TRUE
# 22               B.Lineage        9      1              FALSE   8            TRUE
# 23                NK.Cells        9      0              FALSE   4            TRUE
# 24       Monocytic.Lineage        7      1              FALSE   5            TRUE
# 25       Myeloid.Dendritic        6      0              FALSE   3            TRUE
# 26             Neutrophils       15      3               TRUE  10            TRUE
# 27             Endothelial       33      0              FALSE  22            TRUE
# 28             Fibroblasts        8      5               TRUE   6            TRUE

#
# ==============================================================================






# 4. Validation of processed gene modues (Full to subset correlation)
# ==============================================================================



# Finneo module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


module_size_selected <- module_size_filter %>%
  dplyr::filter(Finneo_size_select) # !!!!!!!!!!!!

module_size_selected <- module_size_selected$Module_name

# validation on expr_tcga
validation_finneo_expr_tcga <- validate_gene_modules(
  module_full_list = module_list[module_size_selected],
  module_subset_list = module_list_finneo[module_size_selected],
  validation_data = expr_tcga)

# validation in expr_metabric
validation_finneo_expr_metabric <- validate_gene_modules(
  module_full_list = module_list[module_size_selected],
  module_subset_list = module_list_finneo[module_size_selected],
  validation_data = expr_metabric)

# validation stat
validation_finneo_stat <- bind_cols(
  validation_finneo_expr_tcga %>%
    dplyr::select(Module_id, Module_full_size, Module_subset_size,
                  Module_subset_size_percent),
  validation_finneo_expr_tcga %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_tcga_",.x)),
  validation_finneo_expr_metabric %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_metabric_",.x))
) %>%
  dplyr::mutate(
    Valid_modules = ((expr_tcga_Pearson %>% round(digits = 2)) > 0.95 &
                       (expr_metabric_Pearson %>% round(digits = 2)) > 0.95),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )


# size selected valid modules in finher dataset
validation_finneo_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::mutate(out = str_c(Module_id, "(N:", Module_subset_size,",",
                            round(Module_subset_size_percent, digits = 1),"%,P:",
                            round(expr_tcga_Pearson, digits = 2),")")) %>%
  dplyr::select(out) %>%
  tibble::deframe()
# [1] "Gruosso2019_Interferon(N:25,22.5%,P:0.96)" "Hamy2016_Immune(N:6,23.1%,P:0.97)"
# [3] "Yang2018_Immune(N:6,35.3%,P:0.98)"         "Farmer2009_MX1(N:18,45%,P:0.99)"
# [5] "Hamy2016_Interferon(N:6,54.5%,P:0.99)"     "Nirmal2018_Interferon(N:19,28.8%,P:0.99)"
# [7] "Hamy2016_Ecm(N:9,42.9%,P:0.99)"            "Naba2014_Ecmcore(N:28,50.9%,P:0.98)"
# [9] "Triulzi2013_Ecm(N:25,43.9%,P:0.99)"        "Sorrentino2014_Chol(N:5,55.6%,P:0.96)"
# [11] "Fibroblasts(N:5,62.5%,P:0.98)"

module_list_sizecor_selected_finneo <- c(

  module_list_finneo[
    validation_finneo_stat %>%
      dplyr::filter(Valid_modules == TRUE) %>%
      dplyr::select(Module_id) %>%
      tibble::deframe()
  ],

  # merged
  list(
    Pooled_Immune = bind_rows(
      module_list_finneo[c("Hamy2016_Immune",
                           "Yang2018_Immune")]),

    Pooled_Interferon = bind_rows(
      module_list_finneo[c("Gruosso2019_Interferon",
                           "Farmer2009_MX1",
                           "Hamy2016_Interferon",
                           "Nirmal2018_Interferon")]),

    Pooled_Fibrosis = bind_rows(
      module_list_finneo[c("Hamy2016_Ecm",
                           "Naba2014_Ecmcore",
                           "Triulzi2013_Ecm")]),

    Pooled_Cholesterol = module_list_finneo$Sorrentino2014_Chol
  )
)


# eliminate redundant genes from merged modules
module_list_sizecor_selected_finneo <- purrr::map(
  module_list_sizecor_selected_finneo,
  ~.x %>% dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE))




# Neoadj module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

module_size_selected <- module_size_filter %>%
  dplyr::filter(Neo_size_select) # !!!!!!!!!!!!

module_size_selected <- module_size_selected$Module_name



# validation on expr_tcga
validation_neo_expr_tcga <- validate_gene_modules(
  module_full_list = module_list[module_size_selected],
  module_subset_list = module_list_neo[module_size_selected],
  validation_data = expr_tcga)


# validation on expr_metabric
validation_neo_expr_metabric <- validate_gene_modules(
  module_full_list = module_list[module_size_selected],
  module_subset_list = module_list_neo[module_size_selected],
  validation_data = expr_metabric)

# validation stat
validation_neo_stat <- bind_cols(
  validation_neo_expr_tcga %>%
    dplyr::select(Module_id, Module_full_size, Module_subset_size,
                  Module_subset_size_percent),
  validation_neo_expr_tcga %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_tcga_",.x)),
  validation_neo_expr_metabric %>%
    dplyr::select(Pearson, Spearman) %>%
    dplyr::rename_all(~str_c("expr_metabric_",.x))
) %>%
  dplyr::mutate(
    Valid_modules = ((expr_tcga_Pearson %>% round(digits = 2)) > 0.95 &
                       (expr_metabric_Pearson %>% round(digits = 2)) > 0.95),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )


# size selected valid modules in neoadj dataset
validation_neo_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::mutate(out = str_c(Module_id, "(N:", Module_subset_size,",",
                            round(Module_subset_size_percent, digits = 1),"%,P:",
                            round(expr_tcga_Pearson, digits = 2),")")) %>%
  dplyr::select(out) %>%
  tibble::deframe()

# [1] "Gruosso2019_Immune(N:493,68.9%,P:1)"        "Gruosso2019_Interferon(N:80,72.1%,P:0.99)"
# [3] "Gruosso2019_Cholesterol(N:41,83.7%,P:0.99)" "Gruosso2019_Fibrosis(N:231,73.1%,P:1)"
# [5] "Hamy2016_Immune(N:13,50%,P:0.99)"           "Teschendorff2007_Immune(N:4,57.1%,P:1)"
# [7] "Yang2018_Immune(N:13,76.5%,P:1)"            "Desmedt2008_STAT1(N:73,76.8%,P:0.98)"
# [9] "Farmer2009_MX1(N:30,75%,P:1)"               "Hamy2016_Interferon(N:8,72.7%,P:0.99)"
# [11] "Nirmal2018_Interferon(N:39,59.1%,P:1)"      "Hamy2016_Ecm(N:11,52.4%,P:0.99)"
# [13] "Naba2014_Ecmcore(N:40,72.7%,P:0.99)"        "Triulzi2013_Ecm(N:41,71.9%,P:1)"
# [15] "Ehmsen2019_Chol(N:7,87.5%,P:0.99)"          "Simigdala2016_Chol(N:7,87.5%,P:0.98)"
# [17] "Sorrentino2014_Chol(N:7,77.8%,P:0.97)"      "Tcell(N:10,62.5%,P:1)"
# [19] "B.Lineage(N:8,88.9%,P:1)"                   "Monocytic.Lineage(N:5,71.4%,P:0.99)"
# [21] "Myeloid.Dendritic(N:3,50%,P:0.99)"          "Endothelial(N:22,66.7%,P:0.99)"
# [23] "Fibroblasts(N:6,75%,P:0.99)"


module_list_sizecor_selected_neo <- c(
  module_list_neo[
    validation_neo_stat %>%
      dplyr::filter(Valid_modules == TRUE) %>%
      dplyr::select(Module_id) %>%
      tibble::deframe()
  ],

  # merged
  list(
    Pooled_Immune = bind_rows(
      module_list_neo[c("Gruosso2019_Immune",
                        "Hamy2016_Immune",
                        "Teschendorff2007_Immune",
                        "Yang2018_Immune",
                        "Desmedt2008_STAT1")]),

    Pooled_Interferon = bind_rows(
      module_list_neo[c("Gruosso2019_Interferon",
                        "Farmer2009_MX1",
                        "Hamy2016_Interferon",
                        "Nirmal2018_Interferon")]),

    Pooled_Fibrosis = bind_rows(
      module_list_neo[c("Gruosso2019_Fibrosis",
                        "Hamy2016_Ecm",
                        "Naba2014_Ecmcore",
                        "Triulzi2013_Ecm")]),

    Pooled_Cholesterol = bind_rows(
      module_list_neo[c("Gruosso2019_Cholesterol",
                        "Ehmsen2019_Chol",
                        "Simigdala2016_Chol",
                        "Sorrentino2014_Chol")])
  )
)


# eliminate redundant genes from merged modules
module_list_sizecor_selected_neo <- purrr::map(
  module_list_sizecor_selected_neo,
  ~.x %>% dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE))







# # cleaning
# rm(validation_neoadj_expr_tcga, validation_neoadj_expr_metabric)
# rm(validation_finher_expr_tcga, validation_finher_expr_metabric)
# rm(expr_tcga, expr_metabric)


#
# ==============================================================================




# 5. Save Robjects
# ==============================================================================


# # Robject created:
# save(module_list_merged, file = str_c(out_data, "module_list_merged.RData"))
# save(module_list_merged_neoadj, file = str_c(out_data, "module_list_merged_neoadj.RData"))
# save(module_list_merged_finher, file = str_c(out_data, "module_list_merged_finher.RData"))
#
# save(validation_neoadj_stat, file = str_c(out_data, "validation_neoadj_stat.RData"))
# save(validation_finher_stat, file = str_c(out_data, "validation_finher_stat.RData"))


# Robject created:
save(module_list, file = str_c(out_data, "module_list.RData"))

save(module_list_finneo, file = str_c(out_data, "module_list_finneo.RData"))
save(module_list_neo, file = str_c(out_data, "module_list_neo.RData"))


save(module_list_sizecor_selected_finneo, file = str_c(out_data, "module_list_sizecor_selected_finneo.RData"))
save(module_list_sizecor_selected_neo, file = str_c(out_data, "module_list_sizecor_selected_neo.RData"))

save(validation_neo_stat, file = str_c(out_data, "validation_neo_stat.RData"))
save(validation_finneo_stat, file = str_c(out_data, "validation_finneo_stat.RData"))


#
# ==============================================================================



# 6. Write out
# ==============================================================================


# Module list: original, finneo and neo subsets

xx <- module_list
idx = names(xx) %>% str_detect("_", negate = T)
names(xx)[idx] <- str_c("MCPcounter_", names(xx)[idx])

lst <- list(
  Original_modules= purrr::map2(
    names(xx), xx,
    ~(.y %>%
        dplyr::mutate(Module_name = .x,
                      Ncbi_gene_id = str_replace(Ncbi_gene_id,"ncbi_","")) %>%
        dplyr::select("Module_name", "Ncbi_gene_id",
                      "Hugo_gene_symbol","Direction")
    )) %>%
    bind_rows()
)


xx <- module_list_sizecor_selected_finneo
idx = names(xx) %>% str_detect("_", negate = T)
names(xx)[idx] <- str_c("MCPcounter_", names(xx)[idx])

lst <- c(lst,
         list(
           FinGEO_module_subsets = purrr::map2(
             names(xx), xx,
             ~(.y %>%
                 dplyr::mutate(Module_name = .x,
                               Ncbi_gene_id = str_replace(Ncbi_gene_id,"ncbi_","")) %>%
                 dplyr::select("Module_name", "Ncbi_gene_id",
                               "Hugo_gene_symbol","Direction")
             )) %>%
             bind_rows()
         )
)

xx <- module_list_sizecor_selected_neo
idx = names(xx) %>% str_detect("_", negate = T)
names(xx)[idx] <- str_c("MCPcounter_", names(xx)[idx])

lst <- c(lst,
         list(
           GEO_module_subsets = purrr::map2(
             names(xx), xx,
             ~(.y %>%
                 dplyr::mutate(Module_name = .x,
                               Ncbi_gene_id = str_replace(Ncbi_gene_id,"ncbi_","")) %>%
                 dplyr::select("Module_name", "Ncbi_gene_id",
                               "Hugo_gene_symbol","Direction")
             )) %>%
             bind_rows()
         )
)

write_xlsx(x = lst,
           path = "results/tables/Module_list.xlsx"
)




# Module validation stat

lst <- list(
  FinGEO_module_subsets = validation_finneo_stat,
  GEO_module_subsets = validation_neo_stat
)

write_xlsx(x = lst,
           path = "results/tables/Module_validation_statistics.xlsx"
)


save(validation_neo_stat, file = str_c(out_data, "validation_neo_stat.RData"))
save(validation_finneo_stat, file = str_c(out_data, "validation_finneo_stat.RData"))



#
# ==============================================================================
