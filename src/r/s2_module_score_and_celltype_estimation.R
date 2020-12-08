# s2_biological_process_and_celltype_estimation.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Clean gene-modules of biological processes and celltype (TIL-localization/
# proliferation/MCPcounter gene modules).
# Validate cleaned gene modules in TCGA/Metabric (MetaGxBreast).
# Estimate biological process and celltype scores.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load geo-data.
# 2. Load and prepare TCGA data for gene-module validation.
# 3. Load, clean and summarize gene-modules.
# 4. Gene-module validation in TCGA.
# 5. Module-score computation of geo data.
# 6. Cell type estimation of geo data.
# 7. Update clin_neoadj with module score and celltype estmiates.



# 1. Load geo-data.
# ==============================================================================

load("results/data/expr_neoadj.RData")
load("results/data/clin_neoadj.RData")

#
# ==============================================================================




# 2. Load and prepare TCGA/Metabric data for gene-module validation.
# ==============================================================================


# In expression set object, the "assayData" slot contains an environment containing
# the "expr" varaibale holding the expression matrix.
# The "expr" variable in the "assayData" environment can be accessed by
# eset@assayData$expr
# Note the "@" used for object slot access, and "$" used for accessing a
# variable from environment.
# Ref: https://www.r-bloggers.com/2011/06/environments-in-r/



eset = MetaGxBreast::loadBreastEsets(loadString = c("METABRIC", "TCGA"))
# Note that according to MetaGxBreast maintainer, there is a bug in loadBreastEsets(),
# that is if the load string doesn't contain a vector of two elements it will
# throw error, except for the single element vector "majority".
# In short, minimu two datasets should be specified in loadString

# loadBreastEsets() log file:
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# snapshotDate(): 2019-10-22
# see ?MetaGxBreast and browseVignettes('MetaGxBreast') for documentation
# downloading 1 resources
# retrieving 1 resource
# |====================================================================| 100%
# loading from cache
# see ?MetaGxBreast and browseVignettes('MetaGxBreast') for documentation
# loading from cache
# Clean up the esets.
# including experiment hub dataset METABRIC
# including experiment hub dataset TCGA
# Ids with missing data: METABRIC
# Warning messages:
#   1: `select_()` is deprecated as of dplyr 0.7.0.
# Please use `select()` instead.
# 2: `filter_()` is deprecated as of dplyr 0.7.0.
# Please use `filter()` instead.
# See vignette('programming') for more help


# metabric <- eset$esets$METABRIC
# tcga <- eset$esets$TCGA


metabric <- list()
tcga <- list()

metabric[["annot"]] <- eset$esets$METABRIC@featureData@data %>% tibble::as_tibble()
metabric[["clin"]] <- eset$esets$METABRIC@phenoData@data %>% tibble::as_tibble()
metabric[["expr"]] <- eset$esets$METABRIC@assayData$exprs %>% tibble::as_tibble(rownames = "probeset") # assayData is an environment


tcga[["annot"]] <- eset$esets$TCGA@featureData@data %>% tibble::as_tibble()
tcga[["clin"]] <- eset$esets$TCGA@phenoData@data %>% tibble::as_tibble(rownames = "sample_name_expr")
tcga[["expr"]] <- eset$esets$TCGA@assayData$exprs %>% tibble::as_tibble(rownames = "probeset") # assayData is an environment


# reverting to basic classes
metabric$annot <- hablar::retype(metabric$annot) %>%
  dplyr::mutate(EntrezGene.ID = EntrezGene.ID %>% as.character())
tcga$annot <- hablar::retype(tcga$annot) %>%
  dplyr::mutate(EntrezGene.ID = EntrezGene.ID %>% as.character())


# max-var collapsing and Ncbi_gene_id mapping
# (Assuming that the annot$best_probe is the max-var/best probe.)
if(identical(metabric$annot$probeset, metabric$expr$probeset)){

  metabric$expr <- metabric$expr %>%
    dplyr::filter(metabric$annot$best_probe == 1)

  metabric$annot <- metabric$annot %>%
    dplyr::filter(best_probe == 1)

  metabric$expr <- metabric$expr %>%
    dplyr::mutate(probeset = str_c("ncbi_",metabric$annot$EntrezGene.ID)) %>%
    dplyr::rename( Ncbi_gene_id = "probeset") # Ncbi_gene_id: local convension
}

if(identical(tcga$annot$probeset, tcga$expr$probeset)){

  tcga$expr <- tcga$expr %>%
    dplyr::filter(tcga$annot$best_probe == 1)

  tcga$annot <- tcga$annot %>%
    dplyr::filter(best_probe == 1)

  tcga$expr <- tcga$expr %>%
    dplyr::mutate(probeset = str_c("ncbi_",tcga$annot$EntrezGene.ID)) %>%
    dplyr::rename( Ncbi_gene_id = "probeset") # Ncbi_gene_id: local convension

  tcga$annot <- tcga$annot %>%
    dplyr::rename(
      Ncbi_gene_id = "EntrezGene.ID",# Ncbi_gene_id: local convension
      Hugo_gene_symbol = "gene",# Hugo_gene_symbol: local convension
    )
}



# # Checking metabric sample congruence
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# identical(metabric$clin$sample_name, names(metabric$expr)[-1]) # TRUE
# # Metabric is congruant
#
#
# # Checking tcga sample congruence
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# identical(tcga$clin$sample_name_expr, names(tcga$expr)[-1]) # TRUE
#
# x <- tcga$clin$sample_name %>% str_replace_all("-",".")
# (x %in% names(tcga$expr)[-1]) %>% table()
# # TRUE
# # 1073
# idx <- purrr::map_int(x,function(xx, id){which(id == xx)}, id = names(tcga$expr)[-1])
# identical(names(tcga$expr)[-1][idx], x) # TRUE
# # TCGA is not congruant

# NOTE !!!!!!!!!!!!:
# Clinical dta and expression data from TCGA are not congruent.


save(metabric, file = str_c(out_data, "metabric.RData"))
save(tcga, file = str_c(out_data, "tcga.RData"))

rm(eset)

#
# ==============================================================================






# 3. Load, clean and summarize gene-modules.
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

# Consolidated and cleaned gene modules
write_tsv(x = module_consolidated,
          path = str_c(out_data, "module_consolidated.tsv"))


module_list <- purrr::map(
  module_consolidated$Module_id %>% unique(),
  function(nme, module_consolidated){

    module_consolidated %>%
      dplyr::filter(Module_id == nme & (!is.na(Ncbi_gene_id))) %>%
      dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE) %>%
      dplyr::select(Ncbi_gene_id, Hugo_gene_symbol, Direction)

  },
  module_consolidated
)
names(module_list) <- module_consolidated$Module_id %>% unique()



module_list_subset <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
  },
  genes = expr_neoadj$Ncbi_gene_id %>% str_replace("ncbi_", "")
)



module_stat <- tibble(
  Module_id = names(module_list),
  Original_size = purrr::map_int(module_list, nrow),
  Subset_size = purrr::map_int(module_list_subset, nrow)
) %>%
  dplyr::mutate(
    Subset_size_percent = (Subset_size/Original_size) * 100
  )

#
# ==============================================================================




# 4. Gene-module validation in TCGA/METABRIC.
# ==============================================================================




#
# ==============================================================================



# 5. Module-score computation of geo data.
# ==============================================================================

#
# ==============================================================================



# 6. Cell type estimation of geo data.
# ==============================================================================


#
# ==============================================================================





# 7. Update clin_neoadj with module score and celltype estmiates.
# ==============================================================================



#
# ==============================================================================


