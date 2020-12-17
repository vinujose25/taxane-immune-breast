# s2_module_score_and_celltype_estimation.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Clean gene-modules of biological processes and celltype (TIL-localization/
# proliferation/MCPcounter gene modules).
# Validate cleaned gene modules in TCGA/Metabric (MetaGxBreast).
# Estimate module and celltype scores.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data (geo, tcga, metabric)
# 2. Load, clean and summarize gene-modules.
# 3. Gene-module validation in tcga-metabric.
# 4. Compute module-scores and update clin_neoadj.
# 5. Estimate celltype scores and update clin_neoadj.
# 6. Additional formating of clinincal data to aid in analysis
# 7. Save Robjects



# 1. Load data (geo, tcga, metabric)
# ==============================================================================

# geo data
# >>>>>>>>

load("results/data/expr_neoadj.RData")
load("results/data/clin_neoadj.RData")


# tcga-metabric
# >>>>>>>>>>>>>

load("~/projects/on/tcga-metabric-metagxbreast/results/data/tcga.RData")
# Note:
# Use tcga expression and clinical data independetly.
# Ids between expression and clinical data are not matching.
load("~/projects/on/tcga-metabric-metagxbreast/results/data/metabric.RData")

# Convert expression matrix to samples x genes
# (Required for module score algorithm)
tcga$expr <- tcga$expr %>% t_tibble(names_x_desc = "Sample_id")
metabric$expr <- metabric$expr %>% t_tibble(names_x_desc = "Sample_id")


# tcga-metabric: Discard genes with atleast one NA expression values in any samples
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tcga$expr %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE
# 19406
# tcga: No genes with NAs, No genes to discard

metabric$expr %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE  TRUE
# 24918     7
# metabric: 7 genes with NAs, 7 genes to discard
idx <- metabric$expr %>% purrr::map_lgl(~!any(is.na(.x))) %>% which()
metabric$expr <- metabric$expr[,idx]

#
# ==============================================================================



# 2. Load, clean and summarize gene-modules.
# ==============================================================================

# gene modules
module_consolidated <- read_tsv(file = "data/gene_modules/Gene_modules_consolidated_clean_ncbi_hugo.tsv")


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

# # Consolidated and cleaned gene modules
# write_tsv(x = module_consolidated,
#           path = str_c(out_data, "module_consolidated.tsv"))


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



module_list_subset <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = expr_neoadj$Ncbi_gene_id# %>% str_replace("ncbi_", "")
)



module_stat <- tibble(
  Module_id = names(module_list),
  Original_size = purrr::map_int(module_list, nrow),
  Subset_size = purrr::map_int(module_list_subset, nrow)
) %>%
  dplyr::mutate(
    Subset_size_percent = (Subset_size/Original_size) * 100
  )


# Discard modules with subset size ==0
identical(module_stat$Module_id, names(module_list)) # TRUE
module_list <- module_list[module_stat$Subset_size > 0]
module_list_subset <- module_list_subset[module_stat$Subset_size > 0]

#
# ==============================================================================



# 3. Gene-module validation in tcga-metabric.
# ==============================================================================

# tcga corr
# >>>>>>>>>

score_full <- get_module_score(x = tcga$expr, module_list = module_list, by = "Ncbi_gene_id")
score_subset <- get_module_score(x = tcga$expr, module_list = module_list_subset, by = "Ncbi_gene_id")

xx <- cor(score_full[,-1], score_subset[,-1], method = "pearson") %>% diag()
xx <- tibble(
  Module_id = names(xx),
  Cor_tcga = xx
)

# update module_stat with tcga correlation
module_stat <- module_stat %>%
  dplyr::left_join(xx, by = "Module_id")



# metabric corr
# >>>>>>>>>>>>>

score_full <- get_module_score(x = metabric$expr, module_list = module_list, by = "Ncbi_gene_id")
score_subset <- get_module_score(x = metabric$expr, module_list = module_list_subset, by = "Ncbi_gene_id")

xx <- cor(score_full[,-1], score_subset[,-1], method = "pearson") %>% diag()
xx <- tibble(
  Module_id = names(xx),
  Cor_metabric = xx
)

# update module_stat with metabric correlation
module_stat <- module_stat %>%
  dplyr::left_join(xx, by = "Module_id")


# Valid gene-modules
# >>>>>>>>>>>>>>>>>>

module_stat <- module_stat %>%
  dplyr::mutate(Is_reliable = if_else(
    (Cor_tcga > 0.9 & Cor_metabric > 0.9),
    "yes", "no"))

#
# ==============================================================================



# 4. Compute module-scores and update clin_neoadj.
# ==============================================================================

# Module filtering
nme <- module_stat %>%
  dplyr::filter(!is.na(Is_reliable) &
                  (Is_reliable == "yes") &
                  str_detect(string = Module_id, pattern = "Becht_2016", negate = TRUE)) %>%
  # Becht_2016 is signatures from MCPcounter.
  # Need original algorithm for celltype estimation.
  dplyr::select(Module_id) %>%
  tibble::deframe()

# Module-score
x <- t_tibble(x = expr_neoadj, names_x_desc = "Sample_geo_accession")
score <- get_module_score(x = x, module_list = module_list_subset[nme], by = "Ncbi_gene_id")

# clin updation
clin_neoadj <- clin_neoadj %>%
  dplyr::left_join(score, by = "Sample_geo_accession")

#
# ==============================================================================



# 5. Estimate celltype scores and update clin_neoadj.
# ==============================================================================

# MCPcounter signatures from github
# probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character")
# genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

x <- expr_neoadj[, -1] %>% as.matrix()
rownames(x) <- expr_neoadj$Ncbi_gene_id %>% str_replace("ncbi_", "")

score <- MCPcounter.estimate(
  expression = x,
  featuresType = "ENTREZ_ID"
) %>%
  t()
# Note that estimate for Cd8.T.Cells is missing as
# the single marker gene for this celltype is missing in pooled dataset.


# Making celltype nams identical to that in module_list
colnames(score) <- str_c(
  "Becht_2016_",
  (colnames(score) %>%
     str_to_title() %>%
     str_replace_all(" ", "."))
)


# Celltype estimate filtering
nme <- module_stat %>%
  dplyr::filter(!is.na(Is_reliable) &
                  (Is_reliable == "yes") &
                  str_detect(string = Module_id, pattern = "Becht_2016")) %>%
  # Becht_2016 is signatures from MCPcounter.
  dplyr::select(Module_id) %>%
  tibble::deframe() # Discarded Neutrophils signature


score <- score %>%
  as_tibble(rownames = "Sample_geo_accession") %>%
  dplyr::select(Sample_geo_accession, all_of(nme))

# clin updation
clin_neoadj <- clin_neoadj %>%
  dplyr::left_join(score, by = "Sample_geo_accession")


#
# ==============================================================================



# 6. Additional formating of clinincal data to aid in analysis
# ==============================================================================

clin_neoadj <- clin_neoadj %>%
  dplyr::mutate(
    # Renaming by preserving original variables
    Immune1 = Gruosso_2019_Immune.cdsig1,
    Immune2 = Hamy_2016_Immune,
    Immune3 = Desmedt_2008_Immune,
    Interferon1 = Gruosso_2019_Interferon.edsig2,
    Interferon2 = Hamy_2016_Interferon,
    Interferon3 = Nirmal_2018_Interferon,
    Cholesterol1 = Gruosso_2019_Cholesterol.edsig5,
    Cholesterol2 = Sorrentino_2014_Cholesterol.mevalonate,
    Cholesterol3 = Simigdala_2016_Cholesterol,
    Fibrosis1 = Gruosso_2019_Fibrosis.cdsig3,
    Fibrosis2 = Hamy_2016_Ecm,
    Fibrosis3 = Triulzi_2013_Ecm,
    Proliferation1 = Desmedt_2008_Proliferation,
    Proliferation2 = Yang_2018_Proliferation,
    Proliferation3 = Nirmal_2018_Proliferation,
    Tcell = Becht_2016_T.Cells,
    CLymphocyte = Becht_2016_Cytotoxic.Lymphocytes,
    Bcell = Becht_2016_B.Lineage,
    NKcell = Becht_2016_Nk.Cells,
    Monocyte = Becht_2016_Monocytic.Lineage,
    MDendritic = Becht_2016_Myeloid.Dendritic.Cells,
    # Endothelial = Becht_2016_Endothelial.Cells, # not reliable
    Fibroblast = Becht_2016_Fibroblasts
  )


#
# ==============================================================================




# 7. Save Robjects
# ==============================================================================


# Robject created:
save(module_consolidated, file = str_c(out_data, "module_consolidated.RData"))
save(module_list, file = str_c(out_data, "module_list.RData"))
save(module_list_subset, file = str_c(out_data, "module_list_subset.RData"))
save(module_stat, file = str_c(out_data, "module_stat.RData"))

# Robject updated:
save(clin_neoadj, file = str_c(out_data, "clin_neoadj.RData"))


#
# ==============================================================================
