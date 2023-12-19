
library(tidyverse)
library(pROC)
source(file = "src/r/functions.R")


# Data
# ====


load("data/metabric.RData")
load("data/tcga.RData")

# expr_metabric <- metabric$expr
# clin_metabric <- metabric$clin
#
# expr_tcga <- tcga$expr
# clin_tcga <- tcga$clin


# Convert expression matrix to samples x genes
# ! This format is needed for module score algorithm
#
metabric$expr <- metabric$expr %>% t_tibble(names_x_desc = "sample_name")
tcga$expr <- tcga$expr %>% t_tibble(names_x_desc = "sample_name")


# Discard genes with atleast one NA expression values in any samples
# ! NAs will create problems in module score algorithm
#
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




# Signatures
# ==========



# PAM50 centroids
pam50 <- list()
pam50[["centroid"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_centroids.txt")
pam50[["annot"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")

pam50 <- pam50$centroid %>%
  dplyr::rename(pcrID = "...1") %>%
  dplyr::left_join(pam50$annot, by = "pcrID") %>%
  dplyr::mutate(EntrezGene = str_c("ncbi_", EntrezGene)) %>%
  dplyr::filter(EntrezGene %in% names(tcga$expr)) %>%
  dplyr::filter(EntrezGene %in% names(metabric$expr)) %>%
  dplyr::rename(Ncbi_gene_id = "EntrezGene")
nrow(pam50) # 36

pam50_siglist <- list()
pam50_siglist[["LumA"]] <- tibble(Ncbi_gene_id = pam50$Ncbi_gene_id, Direction = dplyr::if_else(pam50$LumA < 0, -1, 1 ))
pam50_siglist[["LumB"]] <- tibble(Ncbi_gene_id = pam50$Ncbi_gene_id, Direction = dplyr::if_else(pam50$LumB < 0, -1, 1 ))
pam50_siglist[["Her2"]] <- tibble(Ncbi_gene_id = pam50$Ncbi_gene_id, Direction = dplyr::if_else(pam50$Her2 < 0, -1, 1 ))




# er/her2
erher2_siglist <- list()
erher2_siglist [["er_desmedt2008"]] <- readxl::read_excel(path = "data/gene_modules_erher2/Desmedt2008_Immune-Proliferation_Table_S1_clean_er.xlsx")
erher2_siglist [["her2_desmedt2008"]] <- readxl::read_excel(path = "data/gene_modules_erher2/Desmedt2008_Immune-Proliferation_Table_S1_clean_her2.xlsx")
erher2_siglist [["er_vantveer2002"]] <- readxl::read_excel(path = "data/gene_modules_erher2/vantVeer2002_er_41586_2002_BF415530a_MOESM10_ESM_clean.xlsx")
erher2_siglist [["her2_kalari2013"]] <- readxl::read_excel(path = "data/gene_modules_erher2/Kalari2013_her2_pone.0079298.s005_clean.xlsx")

erher2_siglist$er_desmedt2008 <- erher2_siglist$er_desmedt2008 %>%
  dplyr::mutate(
    Ncbi_gene_id = str_c("ncbi_", EntrezGene.ID),
    Direction = if_else(coefficient < 0 ,-1, 1)
    ) %>%
  dplyr::select(Ncbi_gene_id, Direction) %>%
  na.omit()


erher2_siglist$her2_desmedt2008 <- erher2_siglist$her2_desmedt2008 %>%
  dplyr::mutate(
    Ncbi_gene_id = str_c("ncbi_", EntrezGene.ID),
    Direction = if_else(coefficient < 0 ,-1, 1)
  ) %>%
  dplyr::select(Ncbi_gene_id, Direction) %>%
  na.omit()


erher2_siglist$er_vantveer2002 <- erher2_siglist$er_vantveer2002 %>%
  dplyr::filter(bioDBnet_Gene_ID != "-") %>%
  dplyr::mutate(
    Ncbi_gene_id = str_c("ncbi_", bioDBnet_Gene_ID),
    Direction = if_else(correlation < 0 ,-1, 1)
  ) %>%
  dplyr::select(Ncbi_gene_id, Direction) %>%
  na.omit()


erher2_siglist$her2_kalari2013 <- erher2_siglist$her2_kalari2013 %>%
  dplyr::filter(bioDBnet_Gene_ID != "-") %>%
  dplyr::rename(median_er = "Median ER+",
                median_her2 = "Median HER2+") %>%
  dplyr::mutate(
    Ncbi_gene_id = str_c("ncbi_", bioDBnet_Gene_ID),
    Direction = if_else(median_er > median_her2 ,-1, 1)
  ) %>%
  dplyr::select(Ncbi_gene_id, Direction) %>%
  na.omit()


# gene filtering
pam50_siglist <- purrr::map(
  pam50_siglist,
  function(x,id1,id2){
    x %>%
      dplyr::filter(Ncbi_gene_id %in% id1) %>%
      dplyr::filter(Ncbi_gene_id %in% id2) %>%
      dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE)
  },
  id1 = names(metabric$expr),
  id2 = names(tcga$expr))


erher2_siglist <- purrr::map(
  erher2_siglist,
  function(x,id1,id2){
    x %>%
      dplyr::filter(Ncbi_gene_id %in% id1) %>%
      dplyr::filter(Ncbi_gene_id %in% id2) %>%
      dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE)
  },
  id1 = names(metabric$expr),
  id2 = names(tcga$expr))


# Signature score
# ===============

# Algorithms
# Mean - not considered, as for PAM50 simple mean will be identical for all signature
#   due to matrix nature of signature (all genes are in all sig, only direction changes)
# Weighted average
# Weighted average-2




# Metabric
score <- list()
score[["wavg"]] <- get_module_score(x = metabric$expr,
                                    module_list = c(pam50_siglist, erher2_siglist),
                                    by = "Ncbi_gene_id")
score[["wavg2"]] <- get_module_score_2(x = metabric$expr,
                                       module_list = c(pam50_siglist, erher2_siglist),
                                       by = "Ncbi_gene_id")
score <- purrr::map(names(score),
                    function(nme,score){
                      xscore <- score[[nme]]
                      # 1st column sample name, needed intact for merging later
                      names(xscore)[-1] <- str_c(nme, "_", names(xscore)[-1])
                      xscore
                    },
                    score)
score <- left_join(score[[1]], score[[2]], by = "sample_name")

metabric$clin <- metabric$clin %>% left_join(score, by = "sample_name")


# TCGA
score <- list()
score[["wavg"]] <- get_module_score(x = tcga$expr,
                                    module_list = c(pam50_siglist, erher2_siglist),
                                    by = "Ncbi_gene_id")
score[["wavg2"]] <- get_module_score_2(x = tcga$expr,
                                       module_list = c(pam50_siglist, erher2_siglist),
                                       by = "Ncbi_gene_id")
score <- purrr::map(names(score),
                    function(nme,score){
                      xscore <- score[[nme]]
                      # 1st column sample name, needed intact for merging later
                      names(xscore)[-1] <- str_c(nme, "_", names(xscore)[-1])
                      xscore
                    },
                    score)
score <- left_join(score[[1]], score[[2]], by = "sample_name")


tcga$clin <- tcga$clin %>%
  left_join(score %>% dplyr::rename(sample_name_expr = "sample_name"),
            by = "sample_name_expr")



# Analysis
# ========

tcga$clin <- tcga$clin %>%
  dplyr::mutate(
    er_response = if_else(er == "positive", 1, 0) %>% # NA will propagate
      factor(levels = c(0,1)),
    her2_response = if_else(her2 == "positive", 1, 0) %>% # NA will propagate
      factor(levels = c(0,1))
  )

metabric$clin <- metabric$clin %>%
  dplyr::mutate(
    er_response = if_else(er == "positive", 1, 0) %>% # NA will propagate
      factor(levels = c(0,1)),

    her2_response = if_else(her2 == "positive", 1, 0) %>% # NA will propagate
      factor(levels = c(0,1))
  )


# TCGA ROC AUC comparison (pROC)
for(i in c("wavg", "wavg2")){

  luma_formula = str_c("er_response ~ ",i,"_LumA")
  lumb_formula = str_c("er_response ~ ",i,"_LumB")
  her2_formula = str_c("her2_response ~ ",i,"_Her2")
  erdesmedt_formula = str_c("er_response ~ ",i,"_er_desmedt2008")
  ervantveer_formula = str_c("er_response ~ ",i,"_er_vantveer2002 ")
  her2desmedt_formula = str_c("her2_response ~ ",i,"_her2_desmedt2008")
  her2kalari_formula = str_c("her2_response ~ ",i,"_her2_kalari2013")

  print(
    str_c(i, ";",
          "LumA AUC:", roc(formula = as.formula(luma_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "LumB AUC:", roc(formula = as.formula(lumb_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "HER2 AUC:",roc(formula = as.formula(her2_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "ERdesmedt AUC:", roc(formula = as.formula(erdesmedt_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "ERvantveer AUC:", roc(formula = as.formula(ervantveer_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "HER2desmedt AUC:", roc(formula = as.formula(her2desmedt_formula), data = tcga$clin) %>% auc() %>% round(digits=2), ",",
          "HER2kalari AUC:", roc(formula = as.formula(her2kalari_formula), data = tcga$clin) %>% auc() %>% round(digits=2))
  )
}
# [1] "wavg;LumA AUC:0.51,LumB AUC:0.53,HER2 AUC:0.53/n
# ERdesmedt AUC:0.51,ERvantveer AUC:0.51/n
# HER2desmedt AUC:0.5,HER2kalari AUC:0.54"

# [1] "wavg2;LumA AUC:0.51,LumB AUC:0.53,HER2 AUC:0.53/n
# ERdesmedt AUC:0.51,ERvantveer AUC:0.49/n
# HER2desmedt AUC:0.5,HER2kalari AUC:0.54"


cor(tcga$clin %>% dplyr::select("wavg_LumA", "wavg2_LumA",
                                "wavg_LumB", "wavg2_LumB",
                                "wavg_Her2", "wavg2_Her2"))
#             wavg_LumA wavg2_LumA  wavg_LumB wavg2_LumB  wavg_Her2 wavg2_Her2
# wavg_LumA   1.0000000  0.9991192 -0.2255264 -0.1670173 -0.7315388 -0.7142787
# wavg2_LumA  0.9991192  1.0000000 -0.2300646 -0.1767352 -0.7342911 -0.7195432
# wavg_LumB  -0.2255264 -0.2300646  1.0000000  0.9907692  0.7104382  0.7201121
# wavg2_LumB -0.1670173 -0.1767352  0.9907692  1.0000000  0.6718982  0.6897821
# wavg_Her2  -0.7315388 -0.7342911  0.7104382  0.6718982  1.0000000  0.9980433
# wavg2_Her2 -0.7142787 -0.7195432  0.7201121  0.6897821  0.9980433  1.0000000


# Metabric ROC AUC comparison (pROC)
for(i in c("wavg", "wavg2")){

    luma_formula = str_c("er_response ~ ",i,"_LumA")
    lumb_formula = str_c("er_response ~ ",i,"_LumB")
    her2_formula = str_c("her2_response ~ ",i,"_Her2")
    erdesmedt_formula = str_c("er_response ~ ",i,"_er_desmedt2008")
    ervantveer_formula = str_c("er_response ~ ",i,"_er_vantveer2002 ")
    her2desmedt_formula = str_c("her2_response ~ ",i,"_her2_desmedt2008")
    her2kalari_formula = str_c("her2_response ~ ",i,"_her2_kalari2013")

  print(
    str_c(i, ";",
          "LumA AUC:", roc(formula = as.formula(luma_formula), data = metabric$clin) %>% auc() %>% round(digits=2), ",",
          "LumB AUC:", roc(formula = as.formula(lumb_formula), data = metabric$clin) %>% auc() %>% round(digits=2), ",",
          "HER2 AUC:",roc(formula = as.formula(her2_formula), data = metabric$clin) %>% auc() %>% round(digits=2),",",
          "ERdesmedt AUC:", roc(formula = as.formula(erdesmedt_formula), data = metabric$clin) %>% auc() %>% round(digits=2), ",",
          "ERvantveer AUC:", roc(formula = as.formula(ervantveer_formula), data = metabric$clin) %>% auc() %>% round(digits=2), ",",
          "HER2desmedt AUC:", roc(formula = as.formula(her2desmedt_formula), data = metabric$clin) %>% auc() %>% round(digits=2), ",",
          "HER2kalari AUC:", roc(formula = as.formula(her2kalari_formula), data = metabric$clin) %>% auc() %>% round(digits=2))
  )
}
# [1] "wavg;LumA AUC:0.89,LumB AUC:0.74,HER2 AUC:0.82,
# ERdesmedt AUC:0.95,ERvantveer AUC:0.93,
# HER2desmedt AUC:0.87,HER2kalari AUC:0.59"
#
# [1] "wavg2;LumA AUC:0.89,LumB AUC:0.75,HER2 AUC:0.82,
# ERdesmedt AUC:0.95,ERvantveer AUC:0.94,
# HER2desmedt AUC:0.85,HER2kalari AUC:0.6"


cor(metabric$clin %>% dplyr::select("wavg_er_desmedt2008", "wavg2_er_desmedt2008",
                                    "wavg_er_vantveer2002", "wavg2_er_vantveer2002",
                                    "wavg_her2_desmedt2008", "wavg2_her2_desmedt2008",
                                    "wavg_her2_kalari2013", "wavg2_her2_kalari2013"))
#                        wavg_er_desmedt2008 wavg2_er_desmedt2008 wavg_er_vantveer2002 wavg2_er_vantveer2002
# wavg_er_desmedt2008              1.0000000            0.9991587            0.9453572             0.9562427
# wavg2_er_desmedt2008             0.9991587            1.0000000            0.9417112             0.9514247
# wavg_er_vantveer2002             0.9453572            0.9417112            1.0000000             0.9976216
# wavg2_er_vantveer2002            0.9562427            0.9514247            0.9976216             1.0000000
# wavg_her2_desmedt2008            0.2427807            0.2434679            0.2079980             0.2087766
# wavg2_her2_desmedt2008           0.2793802            0.2837855            0.2600448             0.2548882
# wavg_her2_kalari2013            -0.1933838           -0.1905944           -0.2711714            -0.2451825
# wavg2_her2_kalari2013           -0.2434136           -0.2405500           -0.3096919            -0.2841898
#                        wavg_her2_desmedt2008 wavg2_her2_desmedt2008 wavg_her2_kalari2013 wavg2_her2_kalari2013
# wavg_er_desmedt2008                0.2427807              0.2793802           -0.1933838            -0.2434136
# wavg2_er_desmedt2008               0.2434679              0.2837855           -0.1905944            -0.2405500
# wavg_er_vantveer2002               0.2079980              0.2600448           -0.2711714            -0.3096919
# wavg2_er_vantveer2002              0.2087766              0.2548882           -0.2451825            -0.2841898
# wavg_her2_desmedt2008              1.0000000              0.9439484           -0.1752559            -0.1940744
# wavg2_her2_desmedt2008             0.9439484              1.0000000           -0.2212225            -0.2359283
# wavg_her2_kalari2013              -0.1752559             -0.2212225            1.0000000             0.9883810
# wavg2_her2_kalari2013             -0.1940744             -0.2359283            0.9883810             1.0000000

cor(metabric$clin %>% dplyr::select("wavg_LumA", "wavg2_LumA",
                                    "wavg_LumB", "wavg2_LumB",
                                    "wavg_Her2", "wavg2_Her2"))



#
# ==============================================================================



















