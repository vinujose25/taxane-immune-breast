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
# 3. Gene-module validation in expr_tcga and expr_metabric.
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



# 2. Load and, clean gene-modules.
# ==============================================================================

# tilsig versions
load("results/data/tilsig_clean.RData")
load("results/data/tilsig_bp_merged.RData")

tilsig <- c(
  list(TILsig = tilsig_clean$ALL %>%
         dplyr::select(Ncbi_gene_id1, Direction) %>%
         dplyr::rename(Ncbi_gene_id = "Ncbi_gene_id1")),
  tilsig_bp_merged %>% purrr::map(function(x){x %>% dplyr::mutate(Direction = 1)})
)



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
    ),

    Direction = if_else(Coefficient < 0 , -1, 1)
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



# append tilsig versions
# >>>>>>>>>>>>>>>>>>>>>>

module_list <- c(tilsig, module_list)



# Neoadj module subset
module_list_neoadj <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_neoadj)[-1]
)


# finher module subset
module_list_finher <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_finher)[-1]
)


#
# ==============================================================================



# 3. Gene-module validation in expr_tcga-expr_metabric.
# ==============================================================================


# Neoadj module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# validation on expr_tcga
validation_neoadj_expr_tcga <- validate_gene_modules(
  module_full_list = module_list,
  module_subset_list = module_list_neoadj,
  validation_data = expr_tcga)


# validation on expr_metabric
validation_neoadj_expr_metabric <- validate_gene_modules(
  module_full_list = module_list,
  module_subset_list = module_list_neoadj,
  validation_data = expr_metabric)


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
    Valid_modules = (expr_tcga_Pearson > 0.9 & expr_metabric_Pearson > 0.9),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )



# FinHER module-subset validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# validation on expr_tcga
validation_finher_expr_tcga <- validate_gene_modules(
  module_full_list = module_list,
  module_subset_list = module_list_finher,
  validation_data = expr_tcga)

# validation in expr_metabric
validation_finher_expr_metabric <- validate_gene_modules(
  module_full_list = module_list,
  module_subset_list = module_list_finher,
  validation_data = expr_metabric)


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
    Valid_modules = (expr_tcga_Pearson > 0.9 & expr_metabric_Pearson > 0.9),
    Valid_modules = if_else(is.na(Valid_modules), FALSE, Valid_modules)
  )

# cleaning
rm(validation_neoadj_expr_tcga, validation_neoadj_expr_metabric)
rm(validation_finher_expr_tcga, validation_finher_expr_metabric)
rm(expr_tcga, expr_metabric)

#
# ==============================================================================



# 4. Compute module-scores and update clin_neoadj/clin_finher.
# ==============================================================================


# Neoadj module scores
# >>>>>>>>>>>>>>>>>>>>

# Valid modules
nme <- validation_neoadj_stat %>%
  dplyr::filter(Valid_modules) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
# 39 valid modules out of 43


# Module-score
score <- get_module_score_2(
  x = expr_neoadj,
  module_list = module_list_neoadj[nme],
  by = "Ncbi_gene_id"
) %>%
  dplyr::rename(Sample_geo_accession = "Sample_id")


# clin updation
clin_neoadj <- clin_neoadj %>%
  dplyr::left_join(score, by = "Sample_geo_accession")



# Finher module scores
# >>>>>>>>>>>>>>>>>>>>

# Valid modules
nme <- validation_finher_stat %>%
  dplyr::filter(Valid_modules) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
# 23 valid modules out of 36


# Module-score
score <- get_module_score_2(
  x = expr_finher,
  module_list = module_list_finher[nme],
  by = "Ncbi_gene_id"
) %>%
  dplyr::rename(CEL_filename = "Sample_id")


# clin updation
clin_finher <- clin_finher %>%
  dplyr::left_join(score, by = "CEL_filename")


#
# ==============================================================================



# 5. Estimate celltype scores and update clin_neoadj.
# ==============================================================================

# MCPcounter signatures from github
# probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character")
# genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

# MCPcounter input
# expression: matrix or data.frame with features in rows and samples in columns


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
    Cell_type = str_c("MCPcounter_",
                      Cell_type %>% str_to_title() %>% str_replace_all(" ", "."))
    ) %>%
  t_tibble(names_x_desc = "Sample_geo_accession")  # id name identical to clin_neoadj

# Note that only 9 celltypes were estimated out of 10 MCPcounter celltypes.
# Cd8.T.Cells is missing as the single marker gene for this celltype is
# missing in pooled neoadjuvant dataset.


# Filtering invalid MCPcounter estimates due missing marker genes
nme <- validation_neoadj_stat %>%
  dplyr::filter(Valid_modules) %>%
  # "Becht_2016" and "MCPcounter" represents same modules.
  # Modules with "Becht_2016" keyword is computed as an average.
  # Modules with "MCPcounter" keyword is computed as using original algorithm.
  dplyr::mutate(Module_id = Module_id %>% str_replace("Becht_2016", "MCPcounter")) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
nme <- intersect(nme, names(score))
score <- score %>%
  dplyr::select(1, all_of(nme))


# clin updation
clin_neoadj <- clin_neoadj %>%
  dplyr::left_join(score, by = "Sample_geo_accession")



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
    Cell_type = str_c("MCPcounter_",
                      Cell_type %>% str_to_title() %>% str_replace_all(" ", "."))
  ) %>%
  t_tibble(names_x_desc = "CEL_filename") # id name identical to clin_finher

# Note that only 5 celltypes were estimated out of 10 MCPcounter celltypes.
# Celltypes are missing as the marker genes for the celltypes were
# missing in finher dataset.


# Filtering invalid MCPcounter estimates due missing marker genes
nme <- validation_finher_stat %>%
  dplyr::filter(Valid_modules) %>%
  # "Becht_2016" and "MCPcounter" represents same modules.
  # Modules with "Becht_2016" keyword is computed as an average.
  # Modules with "MCPcounter" keyword is computed as using original algorithm.
  dplyr::mutate(Module_id = Module_id %>% str_replace("Becht_2016", "MCPcounter")) %>%
  dplyr::select(Module_id) %>%
  tibble::deframe()
nme <- intersect(nme, names(score))
score <- score %>%
  dplyr::select(1, all_of(nme))


# clin updation
clin_finher <- clin_finher %>%
  dplyr::left_join(score, by = "CEL_filename")


#
# ==============================================================================



# 6. Additional formating of clinincal data to aid in analysis
# ==============================================================================

clin_neoadj <- clin_neoadj %>%
  dplyr::mutate(
    # Renaming by preserving original variables
    TILsig_scaled = TILsig %>% genefu::rescale(q = 0.05),
    TILsig_APP_Fc = APP_Fc %>% genefu::rescale(q = 0.05),
    TILsig_Immune = Immune %>% genefu::rescale(q = 0.05),
    TILsig_IFNg = IFN_gamma %>% genefu::rescale(q = 0.05),
    TILsig_ECM = ECM %>% genefu::rescale(q = 0.05),
    TILsig_Adhesion = Adhesion %>% genefu::rescale(q = 0.05),
    Immune1 = Gruosso_2019_Immune.cdsig1 %>% genefu::rescale(q = 0.05),
    Immune2 = Hamy_2016_Immune %>% genefu::rescale(q = 0.05),
    # Immune3 = Desmedt_2008_Immune %>% genefu::rescale(q = 0.05),
    Immune3 = Yang_2018_Immune %>% genefu::rescale(q = 0.05),
    Interferon1 = Gruosso_2019_Interferon.edsig2 %>% genefu::rescale(q = 0.05),
    Interferon2 = Hamy_2016_Interferon %>% genefu::rescale(q = 0.05),
    Interferon3 = Nirmal_2018_Interferon %>% genefu::rescale(q = 0.05),
    Cholesterol1 = Gruosso_2019_Cholesterol.edsig5 %>% genefu::rescale(q = 0.05),
    Cholesterol2 = Sorrentino_2014_Cholesterol.mevalonate %>% genefu::rescale(q = 0.05),
    Cholesterol3 = Simigdala_2016_Cholesterol %>% genefu::rescale(q = 0.05),
    Fibrosis1 = Gruosso_2019_Fibrosis.cdsig3 %>% genefu::rescale(q = 0.05),
    Fibrosis2 = Hamy_2016_Ecm %>% genefu::rescale(q = 0.05),
    Fibrosis3 = Triulzi_2013_Ecm %>% genefu::rescale(q = 0.05),
    Proliferation1 = Desmedt_2008_Proliferation %>% genefu::rescale(q = 0.05),
    Proliferation2 = Yang_2018_Proliferation %>% genefu::rescale(q = 0.05),
    Proliferation3 = Nirmal_2018_Proliferation %>% genefu::rescale(q = 0.05),
    Tcell = MCPcounter_T.Cells %>% genefu::rescale(q = 0.05),
    CLymphocyte = MCPcounter_Cytotoxic.Lymphocytes %>% genefu::rescale(q = 0.05),
    Bcell = MCPcounter_B.Lineage %>% genefu::rescale(q = 0.05),
    NKcell = MCPcounter_Nk.Cells %>% genefu::rescale(q = 0.05),
    Monocyte = MCPcounter_Monocytic.Lineage %>% genefu::rescale(q = 0.05),
    MDendritic = MCPcounter_Myeloid.Dendritic.Cells %>% genefu::rescale(q = 0.05),
    # Endothelial = Becht_2016_Endothelial.Cells, # reliable but irrelevant !!!
    Fibroblast = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05)
  )


clin_finher <- clin_finher %>%
  dplyr::mutate(
    # Renaming by preserving original variables
    TILsig_scaled = TILsig %>% genefu::rescale(q = 0.05),
    TILsig_APP_Fc = APP_Fc %>% genefu::rescale(q = 0.05),
    TILsig_Immune = Immune %>% genefu::rescale(q = 0.05),
    TILsig_IFNg = IFN_gamma %>% genefu::rescale(q = 0.05),
    TILsig_Innate = Innate %>% genefu::rescale(q = 0.05),
    TILsig_ECM = ECM %>% genefu::rescale(q = 0.05),
    TILsig_Adhesion = Adhesion %>% genefu::rescale(q = 0.05),
    Immune1 = Gruosso_2019_Immune.cdsig1 %>% genefu::rescale(q = 0.05),
    Immune2 = Hamy_2016_Immune %>% genefu::rescale(q = 0.05),
    # Immune3 = Desmedt_2008_Immune %>% genefu::rescale(q = 0.05),
    Immune3 = Yang_2018_Immune %>% genefu::rescale(q = 0.05),
    Interferon1 = Gruosso_2019_Interferon.edsig2 %>% genefu::rescale(q = 0.05),
    Interferon2 = Hamy_2016_Interferon %>% genefu::rescale(q = 0.05),
    Interferon3 = Nirmal_2018_Interferon %>% genefu::rescale(q = 0.05),
    # Cholesterol1 = Gruosso_2019_Cholesterol.edsig5, # unreliable !!!
    Cholesterol2 = Sorrentino_2014_Cholesterol.mevalonate %>% genefu::rescale(q = 0.05),
    # Cholesterol3 = Simigdala_2016_Cholesterol, # unreliable !!!
    Fibrosis1 = Gruosso_2019_Fibrosis.cdsig3 %>% genefu::rescale(q = 0.05),
    Fibrosis2 = Hamy_2016_Ecm %>% genefu::rescale(q = 0.05),
    Fibrosis3 = Triulzi_2013_Ecm %>% genefu::rescale(q = 0.05),
    Proliferation1 = Desmedt_2008_Proliferation %>% genefu::rescale(q = 0.05),
    Proliferation2 = Yang_2018_Proliferation %>% genefu::rescale(q = 0.05),
    Proliferation3 = Nirmal_2018_Proliferation %>% genefu::rescale(q = 0.05),
    Tcell = MCPcounter_T.Cells %>% genefu::rescale(q = 0.05),
    # CLymphocyte = MCPcounter_Cytotoxic.Lymphocytes, # unreliable !!!
    # Bcell = MCPcounter_B.Lineage, # unreliable !!!
    # NKcell = MCPcounter_Nk.Cells, # unreliable !!!
    # Monocyte = MCPcounter_Monocytic.Lineage, # unreliable !!!
    # MDendritic = MCPcounter_Myeloid.Dendritic.Cells, # unreliable !!!
    # Endothelial = Becht_2016_Endothelial.Cells, # unreliable and irrelevant !!!
    Fibroblast = MCPcounter_Fibroblasts %>% genefu::rescale(q = 0.05)
  )


#
# ==============================================================================



# 7. Save Robjects
# ==============================================================================

# Robject created:
save(module_consolidated, file = str_c(out_data, "module_consolidated.RData"))
save(module_list, file = str_c(out_data, "module_list.RData"))
# save(module_list_subset, file = str_c(out_data, "module_list_subset.RData")) # obsolete
# Print consolidated finher and neoadj module subset list in figures_tables_data.R script
save(module_list_neoadj, file = str_c(out_data, "module_list_neoadj.RData"))
save(module_list_finher, file = str_c(out_data, "module_list_finher.RData"))

# save(module_stat, file = str_c(out_data, "module_stat.RData")) # obsolete
save(validation_neoadj_stat, file = str_c(out_data, "validation_neoadj_stat.RData"))
save(validation_finher_stat, file = str_c(out_data, "validation_finher_stat.RData"))

# Robject updated:
save(clin_neoadj, file = str_c(out_data, "clin_neoadj.RData"))
save(clin_finher, file = str_c(out_data, "clin_finher.RData"))

#
# ====================================================  ==========================
