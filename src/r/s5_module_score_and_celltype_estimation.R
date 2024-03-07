# s5_module_score_and_celltype_estimation_v2.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Estimate module and celltype scores.



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load data (geo, finher)
# 2. Load validated gene-modules.
# 3. Compute module-scores.
# 4. Estimate celltype scores and update score_finher/score_neoadj
# 5. Module agreement(gene), correlation in finher and neoadj
# 6. Additional formatting of clinical data to aid in analysis
# 7. Save Robjects



# 1. Load and prepare data (geo and finher)
# ==============================================================================

# clinical (to update module sore and celltype estimates)
load("results/data/clin_neoadj.RData")
load("results/data/clin_finher.RData")

# expression (to subset original gene modules and validate in expr_tcga/expr_metabric)
load("results/data/expr_neoadj.RData")
load("results/data/expr_finher.RData")

# To select validated MCPcounter sigs
load("results/data/validation_finneo_stat.RData")
load("results/data/validation_neo_stat.RData")

dim(expr_neoadj) # 9184 1500
dim(expr_finher) # 3350 301


# Convert expression matrix to samples x genes (for module score algorithm)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
expr_neoadj <- expr_neoadj %>% t_tibble(names_x_desc = "Sample_id")
expr_finher <- expr_finher %>% t_tibble(names_x_desc = "Sample_id")



# Discard genes with at least one NA expression values in any samples
# (NAs will create problems in module score algorithm)
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




# 2. Load validated gene-modules.
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
load("results/data/module_list_sizecor_selected_finneo.RData") # merged module prefix: "Pooled_"
load("results/data/module_list_sizecor_selected_neo.RData") # merged module prefix: "Pooled_"


# Update module_list_sizecor_selected_* with tilsig
module_list_sizecor_selected_finneo <- c(tilsig_list %>% # filtering tilsig based on neoadj dataset genes
                                           # to make tilsig identical in both finher and neoadj dataset
                                           purrr::map(function(sig, genes){
                                             sig %>%
                                               dplyr::filter(Ncbi_gene_id %in% all_of(genes))
                                           },
                                           genes = names(expr_neoadj)[-1]
                                           ),
                                         module_list_sizecor_selected_finneo)

module_list_sizecor_selected_neo <- c(tilsig_list %>% # filtering tilsig based on neoadj dataset genes
                                        # to make tilsig identical in both finher and neoadj dataset
                                        purrr::map(function(sig, genes){
                                          sig %>%
                                            dplyr::filter(Ncbi_gene_id %in% all_of(genes))
                                        },
                                        genes = names(expr_neoadj)[-1]
                                        ),
                                      module_list_sizecor_selected_neo)




# Order modules
names(module_list_sizecor_selected_finneo)
# [1] "De.novo_TILsig"         "De.novo_Immune"         "De.novo_ECM"
# [4] "Gruosso2019_Interferon" "Hamy2016_Immune"        "Yang2018_Immune"
# [7] "Farmer2009_MX1"         "Hamy2016_Interferon"    "Nirmal2018_Interferon"
# [10] "Hamy2016_Ecm"           "Naba2014_Ecmcore"       "Triulzi2013_Ecm"
# [13] "Sorrentino2014_Chol"    "Fibroblasts"            "Pooled_Immune"
# [16] "Pooled_Interferon"     "Pooled_Fibrosis"            "Pooled_Cholesterol"

module_list_sizecor_selected_finneo <- module_list_sizecor_selected_finneo[
  c(
    "De.novo_TILsig", "De.novo_Immune", "De.novo_ECM",

    "Hamy2016_Immune", "Yang2018_Immune",

    "Gruosso2019_Interferon", "Farmer2009_MX1", "Hamy2016_Interferon", "Nirmal2018_Interferon",

    "Hamy2016_Ecm", "Naba2014_Ecmcore", "Triulzi2013_Ecm",

    "Sorrentino2014_Chol",

    "Pooled_Immune", "Pooled_Interferon", "Pooled_Fibrosis", "Pooled_Cholesterol",

    "Fibroblasts"
  )
]



names(module_list_sizecor_selected_neo)
# [1] "De.novo_TILsig"          "De.novo_Immune"          "De.novo_ECM"
# [4] "Gruosso2019_Immune"      "Gruosso2019_Interferon"  "Gruosso2019_Cholesterol"
# [7] "Gruosso2019_Fibrosis"    "Hamy2016_Immune"         "Teschendorff2007_Immune"
# [10] "Yang2018_Immune"         "Desmedt2008_STAT1"       "Farmer2009_MX1"
# [13] "Hamy2016_Interferon"     "Nirmal2018_Interferon"   "Hamy2016_Ecm"
# [16] "Naba2014_Ecmcore"        "Triulzi2013_Ecm"         "Ehmsen2019_Chol"
# [19] "Simigdala2016_Chol"      "Sorrentino2014_Chol"     "Tcell"
# [22] "B.Lineage"               "Monocytic.Lineage"       "Myeloid.Dendritic"
# [25] "Endothelial"             "Fibroblasts"             "Pooled_Immune"
# [28] "Pooled_Interferon"       "Pooled_Fibrosis"         "Pooled_Cholesterol"

module_list_sizecor_selected_neo <- module_list_sizecor_selected_neo[c(
  "De.novo_TILsig", "De.novo_Immune", "De.novo_ECM",

  "Gruosso2019_Immune", "Hamy2016_Immune", "Teschendorff2007_Immune", "Yang2018_Immune", "Desmedt2008_STAT1",

  "Gruosso2019_Interferon", "Farmer2009_MX1", "Hamy2016_Interferon", "Nirmal2018_Interferon",

  "Gruosso2019_Cholesterol", "Ehmsen2019_Chol", "Simigdala2016_Chol", "Sorrentino2014_Chol",

  "Gruosso2019_Fibrosis", "Hamy2016_Ecm", "Naba2014_Ecmcore", "Triulzi2013_Ecm",

  "Pooled_Immune", "Pooled_Interferon", "Pooled_Fibrosis", "Pooled_Cholesterol",

  "Tcell", "B.Lineage", "Monocytic.Lineage",
  "Myeloid.Dendritic", "Endothelial", "Fibroblasts"
)]


#
# ==============================================================================



# 3. Compute module-scores
# ==============================================================================


# Module-score

# A: Main comparison
score_finher_finneo <- get_module_score_2(
  x = expr_finher,
  module_list = module_list_sizecor_selected_finneo,
  by = "Ncbi_gene_id"
)

# B: Main comparison
score_neoadj_finneo <- get_module_score_2(
  x = expr_neoadj,
  module_list = module_list_sizecor_selected_finneo,
  by = "Ncbi_gene_id"
)

# C: Supplementary
score_neoadj_neo <- get_module_score_2(
  x = expr_neoadj,
  module_list = module_list_sizecor_selected_neo,
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
nme <- validation_finneo_stat %>%
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
score_finher_finneo <- score_finher_finneo %>%
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
nme <- validation_neo_stat %>%
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
score_neoadj_finneo <- score_neoadj_finneo %>%
  dplyr::left_join(
    score %>% dplyr::rename(Sample_id = "Sample_geo_accession"),
    by = "Sample_id")
# discarding MCPcounter estimates not valid in finher;
# finneo considers common valid modules in both finher and neoadjuvant
score_neoadj_finneo = score_neoadj_finneo[names(score_finher_finneo)]

score_neoadj_neo <- score_neoadj_neo %>%
  dplyr::left_join(
    score %>% dplyr::rename(Sample_id = "Sample_geo_accession"),
    by = "Sample_id")

#
# ==============================================================================





# 5. Module agreement(gene), correlation in finher and neoadj
# ==============================================================================

# glimpse(score_finher_finneo)
# glimpse(score_neoadj_finneo)
# glimpse(score_neoadj_neo)



# Finneo module agreement
# >>>>>>>>>>>>>>>>>>>>>>>

module_overlap <- module_gene_overlap(lst = module_list_sizecor_selected_finneo)

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
          "Gruosso2019_Immune", "Hamy2016_Immune",
          "Teschendorff2007_Immune", "Yang2018_Immune",
          "Desmedt2008_STAT1"
        ) ~ "Immune", #"General-Immune",

        .x %in% c(
          "Gruosso2019_Interferon",
          "Farmer2009_MX1", "Hamy2016_Interferon",
          "Nirmal2018_Interferon"
        ) ~ "Interferon", #"General-Interferon",

        .x %in% c(
          "Gruosso2019_Cholesterol", "Ehmsen2019_Chol",
          "Simigdala2016_Chol", "Sorrentino2014_Chol"
        ) ~ "Cholesterol", # "General-Cholesterol",

        .x %in% c(
          "Gruosso2019_Fibrosis",
          "Hamy2016_Ecm", "Naba2014_Ecmcore",
          "Triulzi2013_Ecm"
        ) ~ "Fibrosis", # "General-Fibrosis",

        .x %in% c(
          "Pooled_Immune",
          "Pooled_Interferon",
          "Pooled_Fibrosis",
          "Pooled_Cholesterol"
        ) ~ "Pooled",

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

        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(

        .x == "De-novo" ~ "#e41a1c", # red "tomato",
        .x == "Immune" ~ "#377eb8", # blue "forestgreen",
        .x == "Interferon" ~ "#4daf4a", # green "purple",
        .x == "Fibrosis" ~ "#984ea3",  #"purple",
        .x == "Cholesterol" ~ "#ff7f00", #"orange",
        .x == "Pooled" ~ "gray20",
        .x == "MCP-counter" ~ "#a65628", # brown "steelblue",

         TRUE ~ "Error"
      ))
  )



p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_hline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 8,
  ticksize.y = 5.5,
  angle = 90,
  axis.text.size = 7,
  axis.text.x.vjust = .5,
  plotarea_text_size =  1.75

)


# plot printin
pdf(file = str_c(out_figures,"Module_list_sizecor_selected_finneo_overlap.pdf"),
    width = 5.5, height = 6)
print(
  p_agreement +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(face = "bold"))
      )
dev.off()





# Finneo module-score in finher and neoadj
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# FinHER

names(score_finher_finneo)
# [1] "Sample_id"              "De.novo_TILsig"         "De.novo_Immune"
# [4] "De.novo_ECM"            "Hamy2016_Immune"        "Yang2018_Immune"
# [7] "Gruosso2019_Interferon" "Farmer2009_MX1"         "Hamy2016_Interferon"
# [10] "Nirmal2018_Interferon"  "Hamy2016_Ecm"           "Naba2014_Ecmcore"
# [13] "Triulzi2013_Ecm"        "Sorrentino2014_Chol"    "Pooled_Immune"
# [16] "Pooled_Interferon"      "Pooled_ECM"             "Pooled_Cholesterol"
# [19] "Fibroblasts"            "MCPcounter_Fibroblasts"

# Correlation between sigscore and MCPcounter estimates
cor(score_finher_finneo[,c("Fibroblasts")],
    score_finher_finneo[,c("MCPcounter_Fibroblasts")]) %>%
  diag()
#  0.9914172 # perfect correlation between sigscore and MCPcounter estimates

# Celltype sigscore is replaced by MCPcounter estimates
score <- score_finher_finneo[ , c(1:18,20)] %>%
  dplyr::rename_with(~(str_replace(.x, "MCPcounter_", "")))

p_correlation_finher <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_hline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 6.5,
  ticksize.y = 5.5,
  angle = 45,
  axis.text.size = 7,
  axis.text.x.vjust = 1,
  plotarea_text_size = 1.75
)


# Neoadj

names(score_neoadj_finneo)
# [1] "Sample_id"              "De.novo_TILsig"         "De.novo_Immune"
# [4] "De.novo_ECM"            "Hamy2016_Immune"        "Yang2018_Immune"
# [7] "Gruosso2019_Interferon" "Farmer2009_MX1"         "Hamy2016_Interferon"
# [10] "Nirmal2018_Interferon"  "Hamy2016_Ecm"           "Naba2014_Ecmcore"
# [13] "Triulzi2013_Ecm"        "Sorrentino2014_Chol"    "Pooled_Immune"
# [16] "Pooled_Interferon"      "Pooled_ECM"             "Pooled_Cholesterol"
# [19] "Fibroblasts"            "MCPcounter_Fibroblasts"

# Correlation between sigscore and MCPcounter estimates
cor(score_neoadj_finneo[,c("Fibroblasts")],
    score_neoadj_finneo[,c("MCPcounter_Fibroblasts")]) %>%
  diag()
#  0.9730348 # perfect correlation between sigscore and MCPcounter estimates

# Celltype sigscore is replaced by MCPcounter estimates
score <- score_neoadj_finneo[ , c(1:18,20)] %>%
  dplyr::rename_with(~(str_replace(.x, "MCPcounter_", "")))

p_correlation_neoadj <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_hline = c(0.5, 3.5, 5.5, 9.5, 12.5, 13.5,17.5,18.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 6.5,
  ticksize.y = 5.5,
  angle = 45,
  axis.text.size = 7,
  axis.text.x.vjust = 1,
  plotarea_text_size = 1.75
)



# plot printing
pdf(file = str_c(out_figures,"Module_list_sizecor_selected_finneo_correlation.pdf"),
    width = 6.5, height = 8.5)
print(
  ggarrange(
    p_correlation_finher +
      theme(plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(face = "bold")), # A

    p_correlation_neoadj +
      theme(plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(face = "bold")), # B

    ncol = 1,
    nrow = 2,
    labels = c("A", "B"),
    align = "hv",
    widths = 1,
    legend = "right",
    common.legend = F
  ),
  newpage = F
)
dev.off()




# Neo module agreement and correlation in neoadj
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

module_overlap <- module_gene_overlap(lst = module_list_sizecor_selected_neo)

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
          "Gruosso2019_Immune", "Hamy2016_Immune",
          "Teschendorff2007_Immune", "Yang2018_Immune",
          "Desmedt2008_STAT1"
        ) ~ "Immune", # "General-Immune",

        .x %in% c(
          "Gruosso2019_Interferon",
          "Farmer2009_MX1", "Hamy2016_Interferon",
          "Nirmal2018_Interferon"
        ) ~ "Interferon", # "General-Interferon",

        .x %in% c(
          "Gruosso2019_Cholesterol", "Ehmsen2019_Chol",
          "Simigdala2016_Chol", "Sorrentino2014_Chol"
        ) ~ "Cholesterol", # "General-Cholesterol",

        .x %in% c(
          "Gruosso2019_Fibrosis",
          "Hamy2016_Ecm", "Naba2014_Ecmcore", "Triulzi2013_Ecm"
        ) ~ "Fibrosis", # "General-Fibrosis",

        .x %in% c(
          "Pooled_Immune",
          "Pooled_Interferon",
          "Pooled_Fibrosis",
          "Pooled_Cholesterol"
        ) ~ "Pooled",

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

        TRUE ~ "Error"
      )
    ),
    Module_group_color = purrr::map_chr(
      Module_group,
      ~case_when(

        .x == "De-novo" ~ "#e41a1c", # red "tomato",
        .x == "Immune" ~ "#377eb8", # blue "forestgreen",
        .x == "Interferon" ~ "#4daf4a", # green "purple",
        .x == "Fibrosis" ~ "#984ea3",  #"purple",
        .x == "Cholesterol" ~ "#ff7f00", #"orange",
        .x == "Pooled" ~ "gray20",
        .x == "MCP-counter" ~ "#a65628", # brown "steelblue",

        TRUE ~ "Error"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_color,
  module_group_vline = c(0.5, 3.5, 8.5, 12.5, 16.5, 20.5,24.5, 30.5),
  module_group_hline = c(0.5, 3.5, 8.5, 12.5, 16.5, 20.5,24.5, 30.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 4.5,
  ticksize.y = 3.5,
  angle = 45,
  axis.text.size = 7,
  axis.text.x.vjust = 1,
  plotarea_text_size = 1.75
)



# Module-score

names(score_neoadj_neo)
# [1] "Sample_id"                    "De.novo_TILsig"
# [3] "De.novo_Immune"               "De.novo_ECM"
# [5] "Gruosso2019_Immune"           "Hamy2016_Immune"
# [7] "Teschendorff2007_Immune"      "Yang2018_Immune"
# [9] "Desmedt2008_STAT1"            "Gruosso2019_Interferon"
# [11] "Farmer2009_MX1"               "Hamy2016_Interferon"
# [13] "Nirmal2018_Interferon"        "Gruosso2019_Cholesterol"
# [15] "Ehmsen2019_Chol"              "Simigdala2016_Chol"
# [17] "Sorrentino2014_Chol"          "Gruosso2019_Fibrosis"
# [19] "Hamy2016_Ecm"                 "Naba2014_Ecmcore"
# [21] "Triulzi2013_Ecm"              "Pooled_Immune"
# [23] "Pooled_Interferon"            "Pooled_Fibrosis"
# [25] "Pooled_Cholesterol"           "Tcell"
# [27] "B.Lineage"                    "Monocytic.Lineage"
# [29] "Myeloid.Dendritic"            "Endothelial"
# [31] "Fibroblasts"                  "MCPcounter_Tcell"
# [33] "MCPcounter_B.Lineage"         "MCPcounter_Monocytic.Lineage"
# [35] "MCPcounter_Myeloid.Dendritic" "MCPcounter_Endothelial"
# [37] "MCPcounter_Fibroblasts"

# Correlation between sigscore and MCPcounter estimates
cor(score_neoadj_neo[,c("Tcell","B.Lineage",
  "Monocytic.Lineage", "Myeloid.Dendritic",
  "Endothelial", "Fibroblasts")],
  score_neoadj_neo[,c("MCPcounter_Tcell", "MCPcounter_B.Lineage",
       "MCPcounter_Monocytic.Lineage", "MCPcounter_Myeloid.Dendritic",
      "MCPcounter_Endothelial", "MCPcounter_Fibroblasts")]) %>%
  diag()
#  1 1 1 1 1 1 1 1 1 # perfect correlation between sigscore and MCPcounter estimates


# Celltype sigscore is replaced by MCPcounter estimates
score <- score_neoadj_neo[ ,c(1:25,32:37)] %>%
  dplyr::rename_with(~(str_replace(.x, "MCPcounter_", "")))


p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_color, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 3.5, 8.5, 12.5, 16.5, 20.5,24.5, 30.5),
  module_group_hline = c(0.5, 3.5, 8.5, 12.5, 16.5, 20.5,24.5, 30.5),
  module_group_legend_rows = 4,
  line_col = "gray20",
  ticksize.x = 4.5,
  ticksize.y = 3.5,
  angle = 45,
  axis.text.size = 7,
  axis.text.x.vjust = 1,
  plotarea_text_size = 1.75
)




# plot printin
pdf(file = str_c(out_figures,"Module_list_sizecor_selected_neo_overlap_correlation.pdf"),
    width = 7.5, height = 9)

print(
  ggarrange(
    p_agreement +
      theme(plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(face = "bold")), # A
    p_correlation  +
      theme(plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(face = "bold")), # B
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



# 6. Additional formatting of clinical data to aid in analysis
# ==============================================================================

clin_finher_finneo <- clin_finher %>%
  dplyr::left_join(
    score_finher_finneo %>%
      dplyr::rename(CEL_filename = "Sample_id"),
    by = "CEL_filename")


clin_neoadj_finneo <- clin_neoadj %>%
  dplyr::left_join(
    score_neoadj_finneo %>%
      dplyr::rename(Sample_geo_accession = "Sample_id"),
    by = "Sample_geo_accession")


clin_neoadj_neo <- clin_neoadj %>%
  dplyr::left_join(
    score_neoadj_neo %>%
      dplyr::rename(Sample_geo_accession = "Sample_id"),
    by = "Sample_geo_accession")



# Rescaling finher score
clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::mutate(

    # Renaming by preserving original variables
    scaled_Denovo_TILsig = De.novo_TILsig %>% genefu::rescale(q = 0.05),
    scaled_Denovo_Immune = De.novo_Immune %>% genefu::rescale(q = 0.05),
    scaled_Denovo_ECM = De.novo_ECM %>% genefu::rescale(q = 0.05),

    scaled_Hamy2016_Immune = Hamy2016_Immune %>% genefu::rescale(q = 0.05),
    scaled_Yang2018_Immune = Yang2018_Immune %>% genefu::rescale(q = 0.05),
    scaled_Gruosso2019_Interferon = Gruosso2019_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Farmer2009_MX1 = Farmer2009_MX1 %>% genefu::rescale(q = 0.05),
    scaled_Hamy2016_Interferon = Hamy2016_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Nirmal2018_Interferon = Nirmal2018_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Hamy2016_Ecm = Hamy2016_Ecm %>% genefu::rescale(q = 0.05),
    scaled_Naba2014_Ecmcore = Naba2014_Ecmcore %>% genefu::rescale(q = 0.05),
    scaled_Triulzi2013_Ecm = Triulzi2013_Ecm %>% genefu::rescale(q = 0.05),
    scaled_Sorrentino2014_Chol = Sorrentino2014_Chol %>% genefu::rescale(q = 0.05),

    scaled_Pooled_Immune = Pooled_Immune %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Interferon = Pooled_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Fibrosis = Pooled_Fibrosis %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Cholesterol = Pooled_Cholesterol %>% genefu::rescale(q = 0.05)

    # Only MCP counter estimates were considered
    # scaled_Fibroblasts = Fibroblasts %>% genefu::rescale(q = 0.05),

    # !!! No rescaling on MCPcounter estimates
    # # Original MCPcounter estimates of celltypes are considered

  )



# Rescaling neoadj score
clin_neoadj_finneo <- clin_neoadj_finneo %>%
  dplyr::mutate(

    # Renaming by preserving original variables
    scaled_Denovo_TILsig = De.novo_TILsig %>% genefu::rescale(q = 0.05),
    scaled_Denovo_Immune = De.novo_Immune %>% genefu::rescale(q = 0.05),
    scaled_Denovo_ECM = De.novo_ECM %>% genefu::rescale(q = 0.05),


    scaled_Hamy2016_Immune = Hamy2016_Immune %>% genefu::rescale(q = 0.05),
    scaled_Yang2018_Immune = Yang2018_Immune %>% genefu::rescale(q = 0.05),

    scaled_Gruosso2019_Interferon = Gruosso2019_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Farmer2009_MX1 = Farmer2009_MX1 %>% genefu::rescale(q = 0.05),
    scaled_Hamy2016_Interferon = Hamy2016_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Nirmal2018_Interferon = Nirmal2018_Interferon %>% genefu::rescale(q = 0.05),

    scaled_Hamy2016_Ecm = Hamy2016_Ecm %>% genefu::rescale(q = 0.05),
    scaled_Naba2014_Ecmcore = Naba2014_Ecmcore %>% genefu::rescale(q = 0.05),
    scaled_Triulzi2013_Ecm = Triulzi2013_Ecm %>% genefu::rescale(q = 0.05),

    scaled_Sorrentino2014_Chol = Sorrentino2014_Chol %>% genefu::rescale(q = 0.05),

    scaled_Pooled_Immune = Pooled_Immune %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Interferon = Pooled_Interferon %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Fibrosis = Pooled_Fibrosis %>% genefu::rescale(q = 0.05),
    scaled_Pooled_Cholesterol = Pooled_Cholesterol %>% genefu::rescale(q = 0.05)

    # Only MCP counter estimates were considered
    # scaled_Fibroblasts = Fibroblasts %>% genefu::rescale(q = 0.05),
    #
    # !!! No rescaling on MCPcounter estimates
    # # Original MCPcounter estimates of celltypes are considered

  )

clin_neoadj_neo <- clin_neoadj_neo %>%
dplyr::mutate(

  # Renaming by preserving original variables
  scaled_Denovo_TILsig = De.novo_TILsig %>% genefu::rescale(q = 0.05),
  scaled_Denovo_Immune = De.novo_Immune %>% genefu::rescale(q = 0.05),
  scaled_Denovo_ECM = De.novo_ECM %>% genefu::rescale(q = 0.05),


  scaled_Gruosso2019_Immune = Gruosso2019_Immune %>% genefu::rescale(q = 0.05),
  scaled_Hamy2016_Immune = Hamy2016_Immune %>% genefu::rescale(q = 0.05),
  scaled_Teschendorff2007_Immune = Teschendorff2007_Immune %>% genefu::rescale(q = 0.05),
  scaled_Yang2018_Immune = Yang2018_Immune %>% genefu::rescale(q = 0.05),
  scaled_Desmedt2008_STAT1 = Desmedt2008_STAT1 %>% genefu::rescale(q = 0.05),

  scaled_Gruosso2019_Interferon = Gruosso2019_Interferon %>% genefu::rescale(q = 0.05),
  scaled_Farmer2009_MX1 = Farmer2009_MX1 %>% genefu::rescale(q = 0.05),
  scaled_Hamy2016_Interferon = Hamy2016_Interferon %>% genefu::rescale(q = 0.05),
  scaled_Nirmal2018_Interferon = Nirmal2018_Interferon %>% genefu::rescale(q = 0.05),

  scaled_Gruosso2019_Cholesterol = Gruosso2019_Cholesterol %>% genefu::rescale(q = 0.05),
  scaled_Ehmsen2019_Chol = Ehmsen2019_Chol %>% genefu::rescale(q = 0.05),
  scaled_Simigdala2016_Chol = Simigdala2016_Chol %>% genefu::rescale(q = 0.05),
  scaled_Sorrentino2014_Chol = Sorrentino2014_Chol %>% genefu::rescale(q = 0.05),

  scaled_Gruosso2019_Fibrosis = Gruosso2019_Fibrosis %>% genefu::rescale(q = 0.05),
  scaled_Hamy2016_Ecm = Hamy2016_Ecm %>% genefu::rescale(q = 0.05),
  scaled_Naba2014_Ecmcore = Naba2014_Ecmcore %>% genefu::rescale(q = 0.05),
  scaled_Triulzi2013_Ecm = Triulzi2013_Ecm %>% genefu::rescale(q = 0.05),

  scaled_Pooled_Immune = Pooled_Immune %>% genefu::rescale(q = 0.05),
  scaled_Pooled_Interferon = Pooled_Interferon %>% genefu::rescale(q = 0.05),
  scaled_Pooled_Fibrosis = Pooled_Fibrosis %>% genefu::rescale(q = 0.05),
  scaled_Pooled_Cholesterol = Pooled_Cholesterol %>% genefu::rescale(q = 0.05)

  # Only MCP counter estimates were considered
  # scaled_Fibroblasts = Fibroblasts %>% genefu::rescale(q = 0.05),
  #
  # !!! No rescaling on MCPcounter estimates
  # Original MCPcounter estimates of celltypes are considered

)

#
# ==============================================================================



# 7. Save Robjects
# ==============================================================================

save(clin_finher_finneo, file = str_c(out_data, "clin_finher_finneo.RData"))
save(clin_neoadj_finneo, file = str_c(out_data, "clin_neoadj_finneo.RData"))
save(clin_neoadj_neo, file = str_c(out_data, "clin_neoadj_neo.RData"))

#
# ==============================================================================


# Clean memory
# ==============================================================================

rm(clin_finher_finneo, clin_neoadj_finneo, clin_neoadj_neo,
   score_finher_finneo, score_neoadj_finneo, score_neoadj_neo,
   clin_finher, clin_neoadj,
   expr_finher, expr_neoadj,
   module_list_sizecor_selected_finneo, module_list_sizecor_selected_neo,
   validation_finneo_stat, validation_neo_stat,
   tilsig_list, tilsig_clean, tilsig_bp_merged)

#
# ==============================================================================
