# MTA papaer supplementary table preparation



# Supplementary table: Gene modules associated with TIL
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# The below detailed independent tables were consolidated into on table:
# results/table/Supplementary_til_associated_gene_modules.xlsx


# TILsig gene list:
#     Extracted from results/table/tilsig_consolidated_manually_for_DAVID_input.txt

# Enriched biological processes:
#     Extracted from results/table/david_results_fdr.05.xlsx.

# Gene modules representing common biological processes
# (extracted from ALL version of DAVID analysis results):
#     Computationally generated below using results/data/tilsig_bp_merged.RData
load("results/data/tilsig_bp_merged.RData")
# padding "" terms to make equal no. of rows in each list to apply bind_cols()
max_terms <- purrr::map(tilsig_bp_merged, nrow) %>% unlist() %>% max()
siglist <- purrr::map(tilsig_bp_merged,
                      function(x, max_terms){
                        x <- x$Ncbi_gene_id %>% str_replace("ncbi_", "")

                        if(length(x) < max_terms){
                          x <- c(x, rep("", times = max_terms - length(x)))
                          return(x)
                        } else {
                          return(x)
                        }

                      },max_terms)
siglist <- bind_cols(siglist)
write_xlsx(siglist,path = "results/tables/til_bp_gene_modules.xlsx")
#     Extracted from results/tables/til_bp_gene_modules.xlsx






# Supplementary table: published gene modules
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pubmod <- list()

# Gene module extraction procedure from original publications
pubmod[["Extraction_details"]] <- read_tsv(file = "results/tables/Supplementary_published_gene_module_references.csv")


# Reliability testing of published gene modules
load(str_c(out_data, "validation_neoadj_stat.RData"))
load(str_c(out_data, "validation_finher_stat.RData"))
pubmod[["Reliability_of_finher_subset"]] <- validation_finher_stat
pubmod[["Reliability_of_neoadj_subset"]] <- validation_neoadj_stat


# Original gene modules extracted and its subsets used in this analysis
load(str_c(out_data, "module_list.RData"))
load(str_c(out_data, "module_list_neoadj.RData"))
load(str_c(out_data, "module_list_finher.RData"))

module_original <- purrr::map(names(module_list),
                 function(nme,mlist){
                   x <- mlist[[nme]]
                   x <- x %>%
                     dplyr::mutate(Module_name = nme)
                   x <- x[,c("Module_name",names(x)[-ncol(x)])]
                   x
                 }, mlist = module_list)
names(module_original) <- names(module_list)
module_original <- bind_rows(module_original) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id, "ncbi_",""))
pubmod[["Extracted_modules"]] <- module_original


module_finher <- purrr::map(names(module_list_finher),
                            function(nme,mlist){
                              x <- mlist[[nme]]
                              x <- x %>%
                                dplyr::mutate(Module_name = nme)
                              x <- x[,c("Module_name",names(x)[-ncol(x)])]
                              x
                            }, mlist = module_list_finher)
names(module_finher) <- names(module_list_finher)
module_finher <- bind_rows(module_finher) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id, "ncbi_",""))
pubmod[["Finher_module_subsets"]] <- module_finher


module_neoadj <- purrr::map(names(module_list_neoadj),
                            function(nme,mlist){
                              x <- mlist[[nme]]
                              x <- x %>%
                                dplyr::mutate(Module_name = nme)
                              x <- x[,c("Module_name",names(x)[-ncol(x)])]
                              x
                            }, mlist = module_list_neoadj)
names(module_neoadj) <- names(module_list_neoadj)
module_neoadj <- bind_rows(module_neoadj) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id, "ncbi_",""))
pubmod[["Neoadj_module_subsets"]] <- module_neoadj


# Module name mapping
pubmod[["Module_name_mapping"]] <- tibble(
  Mapped_name = c(
    "TILsig",
    "TILsig_APP_Fc",
    "TILsig_Immune",
    "TILsig_IFNg",
    "TILsig_ECM",
    "TILsig_Adhesion",
    "Immune1",
    "Immune2",
    "Immune3",
    "Interferon1",
    "Interferon2",
    "Interferon3",
    "Cholesterol1",
    "Cholesterol2",
    "Cholesterol3",
    "Fibrosis1",
    "Fibrosis2",
    "Fibrosis3",
    "Proliferation1",
    "Proliferation2",
    "Proliferation3",
    "Tcell",
    "CLymphocyte",
    "Bcell",
    "NKcell",
    "Monocyte",
    "MDendritic",
    "Fibroblast"),

  Original_name = c(
    "TILsig",
    "TILsig_APP_Fc",
    "TILsig_Immune",
    "TILsig_IFNg",
    "TILsig_ECM",
    "TILsig_Adhesion",
    "Gruosso_2019_Immune.cdsig1",
    "Hamy_2016_Immune",
    "Yang_2018_Immune",
    "Gruosso_2019_Interferon.edsig2",
    "Hamy_2016_Interferon",
    "Nirmal_2018_Interferon",
    "Gruosso_2019_Cholesterol.edsig5",
    "Sorrentino_2014_Cholesterol.mevalonate",
    "Simigdala_2016_Cholesterol",
    "Gruosso_2019_Fibrosis.cdsig3",
    "Hamy_2016_Ecm",
    "Triulzi_2013_Ecm",
    "Desmedt_2008_Proliferation",
    "Yang_2018_Proliferation",
    "Nirmal_2018_Proliferation",
    "MCPcounter_T.Cells",
    "MCPcounter_Cytotoxic.Lymphocytes",
    "MCPcounter_B.Lineage",
    "MCPcounter_Nk.Cells",
    "MCPcounter_Monocytic.Lineage",
    "MCPcounter_Myeloid.Dendritic.Cells",
    "MCPcounter_Fibroblasts"
  )
)

write_xlsx(pubmod,
           path = "results/tables/Supplementary_published_gene_modules.xlsx")


