# s12_manuscript_tables.R


# MTA paper table preparation


# Table1: FinHER clinical summary
# ==============================================================================

load("results/data/clin_finher.RData")

clin <- read_tsv(file = "data/finher_clin.tsv") # Private clinical data

clin_finher <- clin_finher %>%
  dplyr::left_join(clin %>% dplyr::select("CEL_filename", "Grade", "Herceptin") %>%
                     dplyr::rename(Grade_original = "Grade",
                                   Herceptin_original = "Herceptin"),
                   by = "CEL_filename" )


table1 <- rbind(
  # Age
  clin_finher %>%
    dplyr::mutate(Age_bin = if_else(Age_Years <= 50, "Age<=50", "Age>50")) %>%
    clin_summary_table(xvar = "Age_bin"),
  # Grade
  clin_finher %>%
    dplyr::mutate(Grade = Grade_original) %>%
    clin_summary_table(xvar = "Grade"),
  # Node
  clin_finher %>%
    dplyr::mutate(Node_bin = Nodal_Status) %>%
    clin_summary_table(xvar = "Node_bin"),
  # Size
  # T stage
  # TX Primary tumor cannot be assessed
  # T0 No evidence of primary tumor
  # T1 size < 2cm
  # T2 2cm < size < 5cm
  # T3 size > 5cm
  # T4 Tumor of any size with direct extensionto the chest wall and/or to the skin
  # (ulceration or skin nodules).
  clin_finher %>%
    dplyr::mutate(Size_cat =  purrr::map_chr(
      Tumor_Size_mm/10,
      ~(case_when(.x < 2 ~ "T0-1",
                  (.x >= 2 & .x < 5) ~ "T2",
                  .x >= 5 ~ "T3",
                  TRUE ~ NA_character_))
    )) %>%
    clin_summary_table(xvar = "Size_cat"),
  # ER
  clin_finher %>%
    clin_summary_table(xvar = "ER_IHC"),
  # PR
  clin_finher %>%
    clin_summary_table(xvar = "PR_IHC"),
  # HER2
  clin_finher %>%
    clin_summary_table(xvar = "HER2_IHC_CISH"),
  # Subtype
  clin_finher %>%
    clin_summary_table(xvar = "Subtype_IHC_2"),
  # Chemo
  clin_finher %>%
    clin_summary_table(xvar = "Chemo"),
  # Herceptin
  clin_finher %>%
    dplyr::mutate(Herceptin = Herceptin_original) %>%
    clin_summary_table(xvar = "Herceptin"),
  # Hormone
  clin_finher %>%
    clin_summary_table(xvar = "Hormone"),
  # Event
  tibble(
    Variable = "Event",
    Values = c("RFS", "DDFS","OS"),
    All = c(clin_finher$RFS_Event %>% sum(),
          clin_finher$DDFS_Event %>% sum(),
          clin_finher$OS_Event %>% sum())
  )
)

write_tsv(table1, file = str_c(out_data, "Table1_finher_clinical_summary.tsv"))
# write_xlsx() is not writing NAs out as is.



#
# ==============================================================================




# Table2: GEO pooled dataset clinical summary
# ==============================================================================

load("results/data/clin_neoadj_finneo.RData")

clin_neoadj_finneo <- clin_neoadj_finneo %>%
  dplyr::filter(str_detect(Arm_consolidated, "AAA")) %>%
  dplyr::filter(Subtype_ihc != "HER2")

clin_neoadj_finneo %>% dim()
# [1] 988 134


table2 <- rbind(
  # Age
  clin_neoadj_finneo %>%
    dplyr::mutate(Age_bin_50 = if_else(Age_bin_50 == "old", ">50","<=50") %>%
                    factor(levels = c("<=50", ">50", NA))) %>%
    clin_summary_table(xvar = "Age_bin_50", yvar = "Series_matrix_accession"),
  # Grade
  clin_neoadj_finneo %>%
    dplyr::mutate(Grade = if_else(Grade == "G1/2", NA_character_,Grade) %>%
                    factor(levels = c("G1", "G2", "G3", NA))) %>%
    clin_summary_table(xvar = "Grade", yvar = "Series_matrix_accession"),
  # Node
  clin_neoadj_finneo %>%
    dplyr::mutate(Node = Node_bin %>%
                    factor(levels = c("neg", "pos", NA))) %>%
    clin_summary_table(xvar = "Node", yvar = "Series_matrix_accession"),
  # Size
  clin_neoadj_finneo %>%
    dplyr::mutate(Size = if_else((Size_cat == "T0" | Size_cat == "T1"), "T0-1", Size_cat) %>%
                    factor(levels = c("T0-1", "T2", "T3", "T4", NA))) %>%
    clin_summary_table(xvar = "Size", yvar = "Series_matrix_accession"),
  # ER
  clin_neoadj_finneo %>%
    dplyr::mutate(Er = Er %>%
                    factor(levels = c("neg", "pos", NA))) %>%
    clin_summary_table(xvar = "Er", yvar = "Series_matrix_accession"),
  # PR
  clin_neoadj_finneo %>%
    dplyr::mutate(Pr = Pr %>%
                    factor(levels = c("neg", "pos", NA))) %>%
    clin_summary_table(xvar = "Pr", yvar = "Series_matrix_accession"),
  # HER2
  clin_neoadj_finneo %>%
    dplyr::mutate(Her2 = Her2 %>%
                    factor(levels = c("neg", "pos", NA))) %>%
    clin_summary_table(xvar = "Her2", yvar = "Series_matrix_accession"),
  # Subtype
  clin_neoadj_finneo %>%
    dplyr::mutate(Subtype_ihc = Subtype_ihc %>%
                    factor(levels = c("HR", "HER2", "TN", NA))) %>%
    clin_summary_table(xvar = "Subtype_ihc", yvar = "Series_matrix_accession"),
  # Arm
  clin_neoadj_finneo %>%
    dplyr::mutate(Arm_consolidated = clin_neoadj_finneo$Arm_consolidated %>%
                    str_replace("///No_her2_agent///No_hormone_therapy///No_other_therapy","") %>%
                    str_replace("\\+noTaxane","") %>%
                    factor(levels = c("AAA", "AAA+Taxane", NA))) %>%
    clin_summary_table(xvar = "Arm_consolidated", yvar = "Series_matrix_accession"),
  # pCR
  clin_neoadj_finneo %>%
    dplyr::mutate(Response = if_else(Response == 1,"Yes", "No") %>%
                    factor(levels = c("No", "Yes", NA))) %>%
    clin_summary_table(xvar = "Response", yvar = "Series_matrix_accession")
)


# Order

x <- table2 %>%
  dplyr::filter(Variable == "Age_bin_50") %>%
  dplyr::select(-c("Variable", "Values", "All")) %>%
  purrr::map(~sum(.x)) %>% unlist() %>% sort(decreasing = T)

table2 <- table2 %>%
  dplyr::select(c("Variable", "Values", "All"), all_of(names(x)))


write_tsv(table2, file = str_c(out_data, "Table2_geo_clinical_summary.tsv"))
# write_xlsx() is not writing NAs out as is.

#
# ==============================================================================




# Table 3: Regimen summary
# ==============================================================================

table(clin_finher[,c("Chemo", "Subtype_IHC_2")])
#       Subtype_IHC_2
# Chemo HER2 TN
# DTX   87 60
# NVB   93 60

table(clin_neoadj_finneo[,c("Arm_consolidated", "Subtype_ihc")])
#                                                                       Subtype_ihc
# Arm_consolidated                                                      HR  TN
# AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91  83
# AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   521 293

#
# ==============================================================================




# Supplementary Table 1: TIL-BP gene modules
# ==============================================================================

load("results/data/tilsig_clean.RData")
load("results/data/tilsig_bp.RData")
load("results/data/tilsig_bp_merged.RData")
david_tnher2 <- readxl::read_excel(path = "results/data/tilsig_david_results_fdr.05.xlsx",
                                   sheet = "TN+HER2")
david_tn <- readxl::read_excel(path = "results/data/tilsig_david_results_fdr.05.xlsx",
                                   sheet = "TN")
david_her2 <- readxl::read_excel(path = "results/data/tilsig_david_results_fdr.05.xlsx",
                                   sheet = "HER2")

names(tilsig_bp)
# [1] "APP"           "Immune"        "Immune_reg."   "Fc_epsilon"    "IFN_gamma"
# [6] "Immune_innate" "Complement1"   "Complement2"   "ECM"           "Adhesion"
nme1 = c("APP", "Immune", "Immune_reg.", "Fc_epsilon", "IFN_gamma",
         "Immune_innate", "Complement1", "Complement2")
nme2 = c("ECM", "Adhesion")

tilsig_bp_immune <- purrr::map2(
  nme1, tilsig_bp[nme1],
  function(x,y){
    y$BP = x
    y
  }
)
tilsig_bp_immune <- Reduce(full_join, tilsig_bp_immune)


tilsig_bp_ecm <- purrr::map2(
  nme2, tilsig_bp[nme2],
  function(x,y){
    y$BP = x
    y
  }
)
tilsig_bp_ecm <- Reduce(full_join, tilsig_bp_ecm)



out = list(
TILsig_TN.HER2 = tilsig_clean$ALL %>%
  dplyr::mutate(
    Estimate = round(Estimate, digits = 3),
    Std_error = round(Std_error, digits = 3),
    T_value = round(T_value, digits = 3),
    P = round(P, digits = 5),
    P_adj = round(P_adj, digits = 5)),
TILsig_TN = tilsig_clean$TN %>%
  dplyr::mutate(
    Estimate = round(Estimate, digits = 3),
    Std_error = round(Std_error, digits = 3),
    T_value = round(T_value, digits = 3),
    P = round(P, digits = 5),
    P_adj = round(P_adj, digits = 5)),
TILsig_HER2 = tilsig_clean$HER2 %>%
  dplyr::mutate(
    Estimate = round(Estimate, digits = 3),
    Std_error = round(Std_error, digits = 3),
    T_value = round(T_value, digits = 3),
    P = round(P, digits = 5),
    P_adj = round(P_adj, digits = 5)),
DAVID_TN.HER2 = david_tnher2,
DAVID_TN = david_tn,
DAVID_HER2 = david_her2,
TILsig_BP_Imm = tilsig_bp_immune,
TILsig_BP_ECM = tilsig_bp_ecm
)

write_xlsx(out,path = str_c(out_data, "STable1_tilsig.xlsx"))

#
# ==============================================================================



# Supplementary Table 2: Published gene modules: Original and subset
# ==============================================================================


module_consolidated <- read_tsv(
  file = "data/gene_modules/Gene_modules_consolidated_clean_ncbi_hugo.tsv"
) %>%
  # unique row for module extraction details
  dplyr::select(Author_year, Pubmed_id, Module_id, Module_extraction_details) %>%
  dplyr::distinct(Module_id,.keep_all = T)


load(str_c(out_data, "module_list.RData"))
load(str_c(out_data, "module_list_sizecor_selected_finneo.RData"))
load(str_c(out_data, "module_list_sizecor_selected_neo.RData"))
load(str_c(out_data, "validation_neo_stat.RData"))
load(str_c(out_data, "validation_finneo_stat.RData"))







out = list(

  # Consolidated orginal module
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>

  Original_module = purrr::map2(
      names(module_list), module_list,
      function(x,y){
        y$Module_id2 = x
        y
      }
    ) %>%
    Reduce(f = full_join),


  # Module references and extraction details
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  Module_extraction = tibble(

    # The following module name map is extracted from s3_pub_gene_module_processing.R
    Module_id = c(
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
      # "Bergamaschi_2008_Ecm1",
      # The original paper contain ECM signature representing four subtypes,
      # it is not sure which one is associated with TIL.
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
      "Becht_2016_Fibroblasts"
    ),


    Module_id2 =  c(
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
      "Fibroblasts"
    )

  ) %>%
    dplyr::left_join(module_consolidated, by = "Module_id"),


  # FinGEO module subset
  # >>>>>>>>>>>>>>>>>>>>

  FinGEO_subset = purrr::map2(
    names(module_list_sizecor_selected_finneo),
    module_list_sizecor_selected_finneo,
    function(x,y){
      y$Module_id2 = x
      y
    }
  ) %>%
    Reduce(f = full_join),


  # GEO module subset
  # >>>>>>>>>>>>>>>>>

  GEO_subset = purrr::map2(
    names(module_list_sizecor_selected_neo),
    module_list_sizecor_selected_neo,
    function(x,y){
      y$Module_id2 = x
      y
    }
  ) %>%
    Reduce(f = full_join),


  # FinGEO module validation statistics
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  FinGEO_stat = validation_finneo_stat,



  # GEO module validation statistics
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  GEO_stat = validation_neo_stat

)


write_xlsx(out, path = str_c(out_data, "STable2_published_gene_modules.xlsx"))

#
# ==============================================================================

