# s2_module_score_and_celltype_estimation.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Clean gene-modules of biological processes and celltype (TIL-localization/
# proliferation/MCPcounter gene modules).
# Validate cleaned gene modules in expr_tcga/expr_metabric (MetaGxBreast).
# Estimate module and celltype scores.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data (geo, finher, expr_tcga, expr_metabric)
# 2. Load, clean, subset gene-modules.
# 3. Module gene agreement testing and module pooling
# 4. Gene-module validation in expr_tcga and expr_metabric.
# 4. Compute module-scores and update clin_neoadj.
# 5. Estimate celltype scores and update clin_neoadj.
# 6. Additional formating of clinincal data to aid in analysis
# 7. Save Robjects



# 1. Load and prepare data (geo, finher, expr_tcga, expr_metabric)
# ==============================================================================

# clinical (to update module sore and celltype estimates)
load("results/data/clin_neoadj.RData")
load("results/data/clin_finher.RData")

# expression (to subset original gene modules and validate in expr_tcga/expr_metabric)
load("results/data/expr_neoadj.RData")
load("results/data/expr_finher.RData")

# # tcga/metabric from metagxbreast R-package
# # Ref: inhouse project "expr_tcga-expr_metabric-metagxbreast" (https://osf.io/jb8za/)
# load("data/tcga.RData")
# load("data/metabric.RData")

dim(expr_neoadj) # 9184 1500
dim(expr_finher) # 3350 301

# dim(tcga$expr) # [1] 19405  1074
# dim(metabric$expr) # [1] 24924  2115


# Convert expression matrix to samples x genes (for module score algorithm)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
expr_neoadj <- expr_neoadj %>% t_tibble(names_x_desc = "Sample_id")
expr_finher <- expr_finher %>% t_tibble(names_x_desc = "Sample_id")



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


#
# ==============================================================================




# 2. Load and, clean gene-modules.
# ==============================================================================

# tilsig versions
load("results/data/tilsig_clean.RData")
load("results/data/tilsig_bp_merged.RData")

tilsig_list <- c(
  list(TILsig = tilsig_clean$ALL %>%
         dplyr::select(Ncbi_gene_id1, Direction) %>%
         dplyr::rename(Ncbi_gene_id = "Ncbi_gene_id1")),
  # tilsig_bp_merged %>% purrr::map(function(x){x %>% dplyr::mutate(Direction = 1)})
  tilsig_bp_merged
)
names(tilsig_list) <- str_c("De.novo_", names(tilsig_list))



# Pub gene modules
load("results/data/module_list_merged_finher.RData")
load("results/data/module_list_merged_neoadj.RData")

load("results/data/validation_finher_stat.RData")
load("results/data/validation_neoadj_stat.RData")


# valid modules
nme <- validation_finher_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
# explicit inclution of non-valid Brain as negative control
module_list_merged_finher <- module_list_merged_finher[c(nme,
                                                         "Tissue_Non.breast",
                                                         "Human_Behaviour" )]



# valid modules
nme <- validation_neoadj_stat %>%
  dplyr::filter(Valid_modules == TRUE) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
module_list_merged_neoadj <- module_list_merged_neoadj[nme]



module_list_finher <- c(tilsig_list,# no filtering as it is derived from finher dataset
                        module_list_merged_finher)

module_list_neoadj <- c(tilsig_list %>% # filtering tilsig based on neoadj dataset genes
                          purrr::map(function(sig, genes){
                            sig %>%
                              dplyr::filter(Ncbi_gene_id %in% all_of(genes))
                          },
                          genes = names(expr_neoadj)[-1]
                          ),
                        module_list_merged_neoadj)


# order modules

names(module_list_finher)
# [1] "De.novo_TILsig"          "De.novo_Immune"          "De.novo_ECM"             "Gruosso2019_Immune"
# [5] "Gruosso2019_Interferon"  "Gruosso2019_Cholesterol" "Gruosso2019_Fibrosis"    "General_Immune"
# [9] "General_Interferon"      "General_ECM"             "General_Cholesterol"     "General_Proliferation"
# [13] "Tcell"                   "Monocytic.Lineage"       "Fibroblasts"             "Tissue_Non.breast"
# [17] "Human_Behaviour"

module_list_finher <- module_list_finher[c(
  "De.novo_TILsig", "De.novo_Immune", "De.novo_ECM",

  "Gruosso2019_Immune", "Gruosso2019_Interferon",
  "Gruosso2019_Cholesterol", "Gruosso2019_Fibrosis",

  "General_Immune", "General_Interferon", "General_ECM", "General_Cholesterol",

  "Tcell", "Monocytic.Lineage", "Fibroblasts",

  "General_Proliferation", "Tissue_Non.breast", "Human_Behaviour"
)
]



names(module_list_neoadj)
# [1] "De.novo_TILsig"          "De.novo_Immune"          "De.novo_ECM"             "Gruosso2019_Immune"
# [5] "Gruosso2019_Interferon"  "Gruosso2019_Cholesterol" "Gruosso2019_Fibrosis"    "General_Immune"
# [9] "General_Interferon"      "General_ECM"             "General_Cholesterol"     "General_Proliferation"
# [13] "Tissue_Non.breast"       "Human_Behaviour"         "Tcell"                   "Cyto.Lymphocyte"
# [17] "B.Lineage"               "NK.Cells"                "Monocytic.Lineage"       "Myeloid.Dendritic"
# [21] "Neutrophils"             "Endothelial"             "Fibroblasts"

module_list_neoadj <- module_list_neoadj[c(
  "De.novo_TILsig", "De.novo_Immune", "De.novo_ECM",

  "Gruosso2019_Immune", "Gruosso2019_Interferon",
  "Gruosso2019_Cholesterol", "Gruosso2019_Fibrosis",

  "General_Immune", "General_Interferon", "General_ECM", "General_Cholesterol",

  "Tcell", "Cyto.Lymphocyte", "B.Lineage", "NK.Cells", "Monocytic.Lineage",
  "Myeloid.Dendritic", "Neutrophils", "Endothelial", "Fibroblasts",

  "General_Proliferation", "Tissue_Non.breast", "Human_Behaviour"
)]


#
# ==============================================================================



# 3. Compute module-scores
# ==============================================================================


# Module-score

score_finher <- get_module_score_2(
  x = expr_finher,
  module_list = module_list_finher,
  by = "Ncbi_gene_id"
)

score_neoadj <- get_module_score_2(
  x = expr_neoadj,
  module_list = module_list_neoadj,
  by = "Ncbi_gene_id"
)



#
# ==============================================================================



# 4. Estimate celltype scores and update score_finher/score_neoadj
# ==============================================================================

# MCPcounter signatures from github
# probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character")
# genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

# MCPcounter input
# expression: matrix or data.frame with features in rows and samples in columns


# FinHER
# >>>>>>

x <- t_tibble(expr_finher, names_x_desc = "Ncbi_gene_id")
xx <- x[, -1] %>% as.matrix()
rownames(xx) <- x$Ncbi_gene_id %>% str_replace("ncbi_", "")

score <- MCPcounter.estimate(
  expression = xx,
  featuresType = "ENTREZ_ID"
) %>%
  as_tibble(rownames = "Cell_type") %>%
  dplyr::mutate(
    Cell_type = str_c(#"MCPcounter_", # validmodule dataframe has different names
                      Cell_type %>% str_to_title() %>% str_replace_all(" ", "."))
  ) %>%
  t_tibble(names_x_desc = "CEL_filename") # id name identical to clin_finher

# Note that only 5 celltypes were estimated out of 10 MCPcounter celltypes.
# Celltypes are missing as the marker genes for the celltypes were
# missing in finher dataset.

# Rename celltypes to aid filtering based on validation dataframe
names(score) <-  purrr::map_chr(
  names(score),
  ~(
    case_when(
    .x == "T.Cells" ~ "Tcell",
    .x == "Cd8.T.Cells" ~ "CD8.Tcell",
    .x == "Cytotoxic.Lymphocytes" ~ "Cyto.Lymphocyte",
    .x == "Nk.Cells" ~ "NK.Cells",
    .x == "Myeloid.Dendritic.Cells" ~ "Myeloid.Dendritic",
    .x == "Endothelial.Cells" ~ "Endothelial",
    TRUE ~ .x # other names not changed
    )
  )
)



# Filtering invalid MCPcounter estimates due missing marker genes
nme <- validation_finher_stat %>%
  dplyr::filter(Valid_modules) %>%
  # "Becht_2016" and "MCPcounter" represents same modules.
  # Modules with "Becht_2016" keyword is computed as an average.
  # Modules with "MCPcounter" keyword is computed as using original algorithm.
  # dplyr::mutate(Module_id = Module_id %>% str_replace("Becht_2016", "MCPcounter")) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()


nme <- intersect(nme, names(score))

score <- score %>%
  dplyr::select(1, all_of(nme)) %>%
  dplyr::rename_with(~(str_c("MCPcounter_", .x)), .cols = -1)


# score_finher updation
score_finher <- score_finher %>%
  dplyr::left_join(
    score %>% dplyr::rename(Sample_id = "CEL_filename"),
    by = "Sample_id")



# Neoadj
# >>>>>>

x <- t_tibble(expr_neoadj, names_x_desc = "Ncbi_gene_id")
xx <- x[, -1] %>% as.matrix()
rownames(xx) <- x$Ncbi_gene_id %>% str_replace("ncbi_", "")

score <- MCPcounter.estimate(
  expression = xx,
  featuresType = "ENTREZ_ID"
) %>%
  as_tibble(rownames = "Cell_type") %>%
  dplyr::mutate(
    Cell_type = str_c(#"MCPcounter_", # validmodule dataframe has different names
                      Cell_type %>% str_to_title() %>% str_replace_all(" ", "."))
    ) %>%
  t_tibble(names_x_desc = "Sample_geo_accession")  # id name identical to clin_neoadj

# Note that only 9 celltypes were estimated out of 10 MCPcounter celltypes.
# Cd8.T.Cells is missing as the single marker gene for this celltype is
# missing in pooled neoadjuvant dataset.


# Rename celltypes to aid filtering based on validation dataframe
names(score) <-  purrr::map_chr(
  names(score),
  ~(
    case_when(
      .x == "T.Cells" ~ "Tcell",
      .x == "Cd8.T.Cells" ~ "CD8.Tcell",
      .x == "Cytotoxic.Lymphocytes" ~ "Cyto.Lymphocyte",
      .x == "Nk.Cells" ~ "NK.Cells",
      .x == "Myeloid.Dendritic.Cells" ~ "Myeloid.Dendritic",
      .x == "Endothelial.Cells" ~ "Endothelial",
      TRUE ~ .x # other names not changed
    )
  )
)




# Filtering invalid MCPcounter estimates due missing marker genes
nme <- validation_neoadj_stat %>%
  dplyr::filter(Valid_modules) %>%
  # "Becht_2016" and "MCPcounter" represents same modules.
  # Modules with "Becht_2016" keyword is computed as an average.
  # Modules with "MCPcounter" keyword is computed as using original algorithm.
  # dplyr::mutate(Module_id = Module_id %>% str_replace("Becht_2016", "MCPcounter")) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()

nme <- intersect(nme, names(score))

score <- score %>%
  dplyr::select(1, all_of(nme)) %>%
  dplyr::rename_with(~(str_c("MCPcounter_", .x)), .cols = -1)


# score_neoadj updation
score_neoadj <- score_neoadj %>%
  dplyr::left_join(
    score %>% dplyr::rename(Sample_id = "Sample_geo_accession"),
    by = "Sample_id")

#
# ==============================================================================



# 5. Module agreement(gene), correlation in finher and neoadj
# ==============================================================================

glimpse(score_finher)
glimpse(score_neoadj)

# FinHER

module_overlap <- module_gene_overlap(lst = module_list_finher)

module_annot <-  module_overlap %>%
  dplyr::select(Module_name, Module_name2) %>%
  dplyr::mutate(
    Module_group = purrr::map_chr(
      Module_name,
      ~case_when(
        .x %in% c(
          "De.novo_TILsig",
          "De.novo_Immune",
          "De.novo_ECM"
        ) ~ "De-novo",
        .x %in% c(
          "Gruosso2019_Immune",
          "Gruosso2019_Interferon",
          "Gruosso2019_Cholesterol",
          "Gruosso2019_Fibrosis"
        ) ~ "TIL-localization",
        .x %in% c(
          "General_Immune",
          "General_Interferon",
          "General_ECM",
          "General_Cholesterol"
        ) ~ "General",
        .x %in% c(
          "Tcell",
          "CD8.Tcell",
          "Cyto.Lymphocyte",
          "B.Lineage",
          "NK.Cells",
          "Monocytic.Lineage",
          "Myeloid.Dendritic",
          "Neutrophils",
          "Endothelial",
          "Fibroblasts"
        ) ~ "MCP-counter",
        .x %in% c(
          "General_Proliferation",
          "Tissue_Non.breast",
          "Human_Behaviour"
        ) ~ "Control",
        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(
        .x == "De-novo" ~ "tomato",
        .x == "TIL-localization" ~ "forestgreen",
        .x == "General" ~ "purple",
        .x == "MCP-counter" ~ "orange",
        .x == "Control" ~ "gray20",
        TRUE ~ "Error"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 3.5, 7.5, 11.5, 14.5, 17.5),
  module_group_hline = c(0.5, 3.5, 7.5, 11.5, 14.5, 17.5),
  module_group_legend_rows = 3,
  line_col = "gray20",
  ticksize.x = 7,
  ticksize.y = 5.5,
  angle = 45,
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2

)



# Module-score

names(score_finher)
# [1] "Sample_id"                 "De.novo_TILsig"               "De.novo_Immune"
# [4] "De.novo_ECM"                  "Gruosso2019_Immune"           "Gruosso2019_Interferon"
# [7] "Gruosso2019_Cholesterol"      "Gruosso2019_Fibrosis"         "General_Immune"
# [10] "General_Interferon"           "General_ECM"                  "General_Cholesterol"
# [13] "Tcell"                        "Monocytic.Lineage"            "Fibroblasts"
# [16] "General_Proliferation"        "Tissue_Non.breast"            "Human_Behaviour"
# [19] "MCPcounter_Tcell"             "MCPcounter_Monocytic.Lineage" "MCPcounter_Fibroblasts"

# Correlation between sigscore and MCPcounter estimates
cor(score_finher[,c("Tcell","Monocytic.Lineage", "Fibroblasts")],
    score_finher[,c("MCPcounter_Tcell", "MCPcounter_Monocytic.Lineage",
                    "MCPcounter_Fibroblasts")]) %>%
  diag()
#  1 1 1 # perfect correlation between sigscore and MCPcounter estimates

# Celltype sigscore is replaced by MCPcounter estimates
score <- score_finher[ , c(1:12,19:21,16:18)] %>%
  dplyr::rename_with(~(str_replace(.x, "MCPcounter_", "")))

p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 3.5, 7.5, 11.5, 14.5, 17.5),
  module_group_hline = c(0.5, 3.5, 7.5, 11.5, 14.5, 17.5),
  module_group_legend_rows = 3,
  line_col = "gray20",
  ticksize.x = 7,
  ticksize.y = 5.5,
  angle = 45,
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2
)


# plot printin
pdf(file = str_c(out_figures,"Module_list_finher_overlap_and_cocorrelation.pdf"),
    width = 7.5, height = 9)
print(
  ggarrange( p_agreement, p_correlation,
             ncol = 1,
             nrow = 2,
             labels = c("A", "B"),
             align = "hv",
             widths = 1,
             # heights = c(.61,.39),
             legend = "right",
             common.legend = F
  ),
  newpage = F
)
dev.off()




# Neoadj

module_overlap <- module_gene_overlap(lst = module_list_neoadj)

module_annot <-  module_overlap %>%
  dplyr::select(Module_name, Module_name2) %>%
  dplyr::mutate(
    Module_group = purrr::map_chr(
      Module_name,
      ~case_when(
        .x %in% c(
          "De.novo_TILsig",
          "De.novo_Immune",
          "De.novo_ECM"
        ) ~ "De-novo",
        .x %in% c(
          "Gruosso2019_Immune",
          "Gruosso2019_Interferon",
          "Gruosso2019_Cholesterol",
          "Gruosso2019_Fibrosis"
        ) ~ "TIL-localization",
        .x %in% c(
          "General_Immune",
          "General_Interferon",
          "General_ECM",
          "General_Cholesterol"
        ) ~ "General",
        .x %in% c(
          "Tcell",
          "CD8.Tcell",
          "Cyto.Lymphocyte",
          "B.Lineage",
          "NK.Cells",
          "Monocytic.Lineage",
          "Myeloid.Dendritic",
          "Neutrophils",
          "Endothelial",
          "Fibroblasts"
        ) ~ "MCP-counter",
        .x %in% c(
          "General_Proliferation",
          "Tissue_Non.breast",
          "Human_Behaviour"
        ) ~ "Control",
        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(
        .x == "De-novo" ~ "tomato",
        .x == "TIL-localization" ~ "forestgreen",
        .x == "General" ~ "purple",
        .x == "MCP-counter" ~ "orange",
        .x == "Control" ~ "gray20",
        TRUE ~ "Error"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 3.5, 7.5, 11.5, 20.5, 23.5),
  module_group_hline = c(0.5, 3.5, 7.5, 11.5, 20.5, 23.5),
  module_group_legend_rows = 3,
  line_col = "gray20",
  ticksize.x = 6.5,
  ticksize.y = 4.5,
  angle = 45,
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2

)



# Module-score

names(score_neoadj)
# [1] "Sample_id"                    "De.novo_TILsig"               "De.novo_Immune"
# [4] "De.novo_ECM"                  "Gruosso2019_Immune"           "Gruosso2019_Interferon"
# [7] "Gruosso2019_Cholesterol"      "Gruosso2019_Fibrosis"         "General_Immune"
# [10] "General_Interferon"           "General_ECM"                  "General_Cholesterol"
# [13] "Tcell"                        "Cyto.Lymphocyte"              "B.Lineage"
# [16] "NK.Cells"                     "Monocytic.Lineage"            "Myeloid.Dendritic"
# [19] "Neutrophils"                  "Endothelial"                  "Fibroblasts"
# [22] "General_Proliferation"        "Tissue_Non.breast"            "Human_Behaviour"
# [25] "MCPcounter_Tcell"             "MCPcounter_Cyto.Lymphocyte"   "MCPcounter_B.Lineage"
# [28] "MCPcounter_NK.Cells"          "MCPcounter_Monocytic.Lineage" "MCPcounter_Myeloid.Dendritic"
# [31] "MCPcounter_Neutrophils"       "MCPcounter_Endothelial"       "MCPcounter_Fibroblasts"


# Correlation between sigscore and MCPcounter estimates
cor(score_neoadj[,c("Tcell", "Cyto.Lymphocyte", "B.Lineage",
  "NK.Cells", "Monocytic.Lineage", "Myeloid.Dendritic",
  "Neutrophils", "Endothelial", "Fibroblasts")],
  score_neoadj[,c("MCPcounter_Tcell", "MCPcounter_Cyto.Lymphocyte", "MCPcounter_B.Lineage",
      "MCPcounter_NK.Cells", "MCPcounter_Monocytic.Lineage", "MCPcounter_Myeloid.Dendritic",
      "MCPcounter_Neutrophils", "MCPcounter_Endothelial", "MCPcounter_Fibroblasts")]) %>%
  diag()
#  1 1 1 1 1 1 1 1 1 # perfect correlation between sigscore and MCPcounter estimates


# Celltype sigscore is replaced by MCPcounter estimates
score <- score_neoadj[ , c(1:12,25:33,22:24)] %>%
  dplyr::rename_with(~(str_replace(.x, "MCPcounter_", "")))


p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 3.5, 7.5, 11.5, 20.5, 23.5),
  module_group_hline = c(0.5, 3.5, 7.5, 11.5, 20.5, 23.5),
  module_group_legend_rows = 3,
  line_col = "gray20",
  ticksize.x = 6.5,
  ticksize.y = 4.5,
  angle = 45,
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2
)




# plot printin
pdf(file = str_c(out_figures,"Module_list_neoadj_overlap_and_cocorrelation.pdf"),
    width = 7.5, height = 9)
print(
  ggarrange( p_agreement, p_correlation,
             ncol = 1,
             nrow = 2,
             labels = c("A", "B"),
             align = "hv",
             widths = 1,
             # heights = c(.61,.39),
             legend = "right",
             common.legend = F
  ),
  newpage = F
)
dev.off()

#
# ==============================================================================



# 6. Additional formating of clinincal data to aid in analysis
# ==============================================================================

# clin_neoadj <- clin_neoadj %>%
#   dplyr::mutate(
#     # Renaming by preserving original variables
#     TILsig_scaled = TILsig %>% genefu::rescale(q = 0.05),
#     TILsig_APP_Fc = APP_Fc %>% genefu::rescale(q = 0.05),
#     TILsig_Immune = Immune %>% genefu::rescale(q = 0.05),
#     TILsig_IFNg = IFN_gamma %>% genefu::rescale(q = 0.05),
#     TILsig_ECM = ECM %>% genefu::rescale(q = 0.05),
#     TILsig_Adhesion = Adhesion %>% genefu::rescale(q = 0.05),
#     Immune1 = Gruosso_2019_Immune.cdsig1 %>% genefu::rescale(q = 0.05),
#     Immune2 = Hamy_2016_Immune %>% genefu::rescale(q = 0.05),
#     # Immune3 = Desmedt_2008_Immune %>% genefu::rescale(q = 0.05),
#     Immune3 = Yang_2018_Immune %>% genefu::rescale(q = 0.05),
#     Interferon1 = Gruosso_2019_Interferon.edsig2 %>% genefu::rescale(q = 0.05),
#     Interferon2 = Hamy_2016_Interferon %>% genefu::rescale(q = 0.05),
#     Interferon3 = Nirmal_2018_Interferon %>% genefu::rescale(q = 0.05),
#     Cholesterol1 = Gruosso_2019_Cholesterol.edsig5 %>% genefu::rescale(q = 0.05),
#     Cholesterol2 = Sorrentino_2014_Cholesterol.mevalonate %>% genefu::rescale(q = 0.05),
#     Cholesterol3 = Simigdala_2016_Cholesterol %>% genefu::rescale(q = 0.05),
#     Fibrosis1 = Gruosso_2019_Fibrosis.cdsig3 %>% genefu::rescale(q = 0.05),
#     Fibrosis2 = Hamy_2016_Ecm %>% genefu::rescale(q = 0.05),
#     Fibrosis3 = Triulzi_2013_Ecm %>% genefu::rescale(q = 0.05),
#     Proliferation1 = Desmedt_2008_Proliferation %>% genefu::rescale(q = 0.05),
#     Proliferation2 = Yang_2018_Proliferation %>% genefu::rescale(q = 0.05),
#     Proliferation3 = Nirmal_2018_Proliferation %>% genefu::rescale(q = 0.05),
#     Tcell = MCPcounter_T.Cells %>% genefu::rescale(q = 0.05),
#     CLymphocyte = MCPcounter_Cytotoxic.Lymphocytes %>% genefu::rescale(q = 0.05),
#     Bcell = MCPcounter_B.Lineage %>% genefu::rescale(q = 0.05),
#     NKcell = MCPcounter_Nk.Cells %>% genefu::rescale(q = 0.05),
#     Monocyte = MCPcounter_Monocytic.Lineage %>% genefu::rescale(q = 0.05),
#     MDendritic = MCPcounter_Myeloid.Dendritic.Cells %>% genefu::rescale(q = 0.05),
#     # Endothelial = Becht_2016_Endothelial.Cells, # reliable but irrelevant !!!
#     Fibroblast = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05)
#   )
#
#
# clin_finher <- clin_finher %>%
#   dplyr::mutate(
#     # Renaming by preserving original variables
#     TILsig_scaled = TILsig %>% genefu::rescale(q = 0.05),
#     TILsig_APP_Fc = APP_Fc %>% genefu::rescale(q = 0.05),
#     TILsig_Immune = Immune %>% genefu::rescale(q = 0.05),
#     TILsig_IFNg = IFN_gamma %>% genefu::rescale(q = 0.05),
#     TILsig_Innate = Innate %>% genefu::rescale(q = 0.05),
#     TILsig_ECM = ECM %>% genefu::rescale(q = 0.05),
#     TILsig_Adhesion = Adhesion %>% genefu::rescale(q = 0.05),
#     Immune1 = Gruosso_2019_Immune.cdsig1 %>% genefu::rescale(q = 0.05),
#     Immune2 = Hamy_2016_Immune %>% genefu::rescale(q = 0.05),
#     # Immune3 = Desmedt_2008_Immune %>% genefu::rescale(q = 0.05),
#     Immune3 = Yang_2018_Immune %>% genefu::rescale(q = 0.05),
#     Interferon1 = Gruosso_2019_Interferon.edsig2 %>% genefu::rescale(q = 0.05),
#     Interferon2 = Hamy_2016_Interferon %>% genefu::rescale(q = 0.05),
#     Interferon3 = Nirmal_2018_Interferon %>% genefu::rescale(q = 0.05),
#     # Cholesterol1 = Gruosso_2019_Cholesterol.edsig5, # unreliable !!!
#     Cholesterol2 = Sorrentino_2014_Cholesterol.mevalonate %>% genefu::rescale(q = 0.05),
#     # Cholesterol3 = Simigdala_2016_Cholesterol, # unreliable !!!
#     Fibrosis1 = Gruosso_2019_Fibrosis.cdsig3 %>% genefu::rescale(q = 0.05),
#     Fibrosis2 = Hamy_2016_Ecm %>% genefu::rescale(q = 0.05),
#     Fibrosis3 = Triulzi_2013_Ecm %>% genefu::rescale(q = 0.05),
#     Proliferation1 = Desmedt_2008_Proliferation %>% genefu::rescale(q = 0.05),
#     Proliferation2 = Yang_2018_Proliferation %>% genefu::rescale(q = 0.05),
#     Proliferation3 = Nirmal_2018_Proliferation %>% genefu::rescale(q = 0.05),
#     Tcell = MCPcounter_T.Cells %>% genefu::rescale(q = 0.05),
#     # CLymphocyte = MCPcounter_Cytotoxic.Lymphocytes, # unreliable !!!
#     # Bcell = MCPcounter_B.Lineage, # unreliable !!!
#     # NKcell = MCPcounter_Nk.Cells, # unreliable !!!
#     # Monocyte = MCPcounter_Monocytic.Lineage, # unreliable !!!
#     # MDendritic = MCPcounter_Myeloid.Dendritic.Cells, # unreliable !!!
#     # Endothelial = Becht_2016_Endothelial.Cells, # unreliable and irrelevant !!!
#     Fibroblast = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05)
#   )

clin_neoadj <- clin_neoadj %>%
  dplyr::select(1:98) %>%
  dplyr::left_join(
    score_neoadj %>%
      dplyr::rename(Sample_geo_accession = "Sample_id"),
    by = "Sample_geo_accession")

clin_finher <- clin_finher %>%
  dplyr::select(1:78) %>%
  dplyr::left_join(
    score_finher %>%
      dplyr::rename(CEL_filename = "Sample_id"),
    by = "CEL_filename")


# Rescaling neoadj score
clin_neoadj <- clin_neoadj %>%
dplyr::mutate(

  # Renaming by preserving original variables
  scaled_Denovo_TILsig = De.novo_TILsig %>% genefu::rescale(q = 0.05),
  scaled_Denovo_Immune = De.novo_Immune %>% genefu::rescale(q = 0.05),
  scaled_Denovo_ECM = De.novo_ECM %>% genefu::rescale(q = 0.05),

  scaled_Gruosso2019_Immune = Gruosso2019_Immune %>% genefu::rescale(q = 0.05),
  scaled_Gruosso2019_Interferon = Gruosso2019_Interferon %>% genefu::rescale(q = 0.05),
  scaled_Gruosso2019_Cholesterol = Gruosso2019_Cholesterol %>% genefu::rescale(q = 0.05),
  scaled_Gruosso2019_Fibrosis = Gruosso2019_Fibrosis %>% genefu::rescale(q = 0.05),

  scaled_General_Immune = General_Immune %>% genefu::rescale(q = 0.05),
  scaled_General_Interferon = General_Interferon %>% genefu::rescale(q = 0.05),
  scaled_General_ECM = General_ECM %>% genefu::rescale(q = 0.05),
  scaled_General_Cholesterol = General_Cholesterol %>% genefu::rescale(q = 0.05),

  # !!! No rescaling on MCPcounter estimates
  # # Only MCPcounter estimates of celltypes are considered
  # MCPcounter_Tcell = MCPcounter_Tcell %>% genefu::rescale(q = 0.05),
  # MCPcounter_Cyto.Lymphocyte = MCPcounter_Cyto.Lymphocyte %>% genefu::rescale(q = 0.05),
  # MCPcounter_B.Lineage = MCPcounter_B.Lineage %>% genefu::rescale(q = 0.05),
  # MCPcounter_NK.Cells = MCPcounter_NK.Cells %>% genefu::rescale(q = 0.05),
  # MCPcounter_Monocytic.Lineage = MCPcounter_Monocytic.Lineage %>% genefu::rescale(q = 0.05),
  # MCPcounter_Myeloid.Dendritic = MCPcounter_Myeloid.Dendritic %>% genefu::rescale(q = 0.05),
  # MCPcounter_Neutrophils = MCPcounter_Neutrophils %>% genefu::rescale(q = 0.05),
  # MCPcounter_Endothelial = MCPcounter_Endothelial %>% genefu::rescale(q = 0.05),
  # MCPcounter_Fibroblasts = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05),

  scaled_Control_Proliferation = General_Proliferation %>% genefu::rescale(q = 0.05),
  scaled_Control_Non.breast.tissue = Tissue_Non.breast %>% genefu::rescale(q = 0.05),
  scaled_Control_Human_Behaviour = Human_Behaviour %>% genefu::rescale(q = 0.05)
)



# Rescaling finher score
clin_finher <- clin_finher %>%
  dplyr::mutate(

    # Renaming by preserving original variables
    scaled_Denovo_TILsig = De.novo_TILsig %>% genefu::rescale(q = 0.05),
    scaled_Denovo_Immune = De.novo_Immune %>% genefu::rescale(q = 0.05),
    scaled_Denovo_ECM = De.novo_ECM %>% genefu::rescale(q = 0.05),

    scaled_Gruosso2019_Immune = Gruosso2019_Immune %>% genefu::rescale(q = 0.05),
    scaled_Gruosso2019_Interferon = Gruosso2019_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Gruosso2019_Cholesterol = Gruosso2019_Cholesterol %>% genefu::rescale(q = 0.05),
    scaled_Gruosso2019_Fibrosis = Gruosso2019_Fibrosis %>% genefu::rescale(q = 0.05),

    scaled_General_Immune = General_Immune %>% genefu::rescale(q = 0.05),
    scaled_General_Interferon = General_Interferon %>% genefu::rescale(q = 0.05),
    scaled_General_ECM = General_ECM %>% genefu::rescale(q = 0.05),
    scaled_General_Cholesterol = General_Cholesterol %>% genefu::rescale(q = 0.05),

    # !!! No rescaling on MCPcounter estimates
    # # Only MCPcounter estimates of celltypes are considered
    # scaled_MCPcounter_Tcell = MCPcounter_Tcell %>% genefu::rescale(q = 0.05),
    # scaled_MCPcounter_Monocytic.Lineage = MCPcounter_Monocytic.Lineage %>% genefu::rescale(q = 0.05),
    # scaled_MCPcounter_Fibroblasts = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05),

    scaled_Control_Proliferation = General_Proliferation %>% genefu::rescale(q = 0.05),
    scaled_Control_Non.breast.tissue = Tissue_Non.breast %>% genefu::rescale(q = 0.05),
    scaled_Control_Human_Behaviour = Human_Behaviour %>% genefu::rescale(q = 0.05)

  )


#
# ==============================================================================



# 7. Save Robjects
# ==============================================================================

# # Robject created:
# save(module_consolidated, file = str_c(out_data, "module_consolidated.RData"))
# save(module_list, file = str_c(out_data, "module_list.RData"))
# # save(module_list_subset, file = str_c(out_data, "module_list_subset.RData")) # obsolete
# # Print consolidated finher and neoadj module subset list in figures_tables_data.R script
# save(module_list_neoadj, file = str_c(out_data, "module_list_neoadj.RData"))
# save(module_list_finher, file = str_c(out_data, "module_list_finher.RData"))

# # save(module_stat, file = str_c(out_data, "module_stat.RData")) # obsolete
# save(validation_neoadj_stat, file = str_c(out_data, "validation_neoadj_stat.RData"))
# save(validation_finher_stat, file = str_c(out_data, "validation_finher_stat.RData"))

# Robject updated:
save(clin_neoadj, file = str_c(out_data, "clin_neoadj.RData"))
save(clin_finher, file = str_c(out_data, "clin_finher.RData"))

#
# ====================================================  ==========================
