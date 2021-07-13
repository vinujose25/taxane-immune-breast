

# 1. Load and prepare data (geo, finher, tcga, metabric)
# ==============================================================================

# clinical
load("results/data/clin_neoadj.RData")
load("results/data/clin_finher.RData")

# expression (to subset original gene modules and validate in tcga/metabric)
load("results/data/expr_neoadj.RData")
load("results/data/expr_finher.RData")
# TCGA/METABRIC from metagxbreast R-package
# Ref: inhouse project "tcga-metabric-metagxbreast" (https://osf.io/jb8za/)
# load("data/tcga.RData") # clincal and expression data are not congruent
load("data/metabric.RData")

dim(expr_neoadj) # 9184 1500
dim(expr_finher) # 3350 301

# dim(tcga$expr) # [1] 19405  1074
dim(metabric$expr) # [1] 24924  2115

expr_metabric <- metabric$expr
clin_metabric <- metabric$clin


# Convert expression matrix to samples x genes (for module score algorithm)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
expr_neoadj <- expr_neoadj %>% t_tibble(names_x_desc = "Sample_geo_accession")
expr_finher <- expr_finher %>% t_tibble(names_x_desc = "CEL_filename")
# tcga <- tcga$expr %>% t_tibble(names_x_desc = "Sample_id")
expr_metabric <- expr_metabric %>% t_tibble(names_x_desc = "sample_name")



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

# tcga %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# # FALSE
# # 19406
# # tcga: No genes with NAs, No genes to discard

expr_metabric %>% purrr::map_lgl(~any(is.na(.x))) %>% table() # includes Non-NA Sample_id
# FALSE  TRUE
# 24918     7
# metabric: 7 genes with NAs, 7 genes to discard
idx <- expr_metabric %>% purrr::map_lgl(~!any(is.na(.x))) %>% which()
expr_metabric <- expr_metabric[,idx]

#
# ==============================================================================



# 2. Load and, clean gene-modules.
# ==============================================================================

# gene modules (vectorised)
module_consolidated <- readxl::read_xlsx(
  path = "data/gene_modules/original/Desmedt2008_Immune-Proliferation_Table_S1.xlsx",
  sheet = 1) %>%
  dplyr::filter(! (is.na(module) | is.na(EntrezGene.ID) | is.na(coefficient))) %>%
  dplyr::mutate(Direction = if_else(coefficient < 0, -1, 1)) %>%
  dplyr::rename(Module_id = "module",
                Ncbi_gene_id = "EntrezGene.ID",
                Hugo_gene_symbol = "HUGO.gene.symbol") %>%
  dplyr::select(Module_id, Ncbi_gene_id, Hugo_gene_symbol, Direction)




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



module_list_metabric <- purrr::map(
  module_list,
  function(sig, genes){
    sig %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes))
    # Gene-ids are already formatted to ncbi_*
  },
  genes = names(expr_metabric)[-1]
)


purrr::map(module_list_neoadj,~(.x$Direction %>% table()))
# $ESR1
# -1   1
# 150 192
# $ERBB2
# -1  1
# 3 17
purrr::map(module_list_finher,~(.x$Direction %>% table()))
# $ESR1
# -1  1
# 65 82
# $ERBB2
# -1  1
# 1 14
purrr::map(module_list_metabric,~(.x$Direction %>% table()))
# $ESR1
# -1   1
# 193 272
# $ERBB2
# -1  1
# 5 23

#
# ==============================================================================



# 3. Module score algorithm omparison in neadj and finher
# ==============================================================================

# Subset neoadj and finher with ER and HER2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clin_neoadj %>%
  dplyr::group_by(Er, Her2) %>%
  dplyr::summarise(N = n())
#   Er    Her2      N
# 1 neg   neg     569
# 2 neg   pos     148
# 3 pos   neg     708
# 4 pos   pos      74

# dim(clin_neoadj)
# [1] 1499  161
# All samples have ER/HER2 status available


clin_finher %>%
  dplyr::group_by(ER_IHC, HER2_IHC_CISH) %>%
  dplyr::summarise(N = n())
#   ER_IHC   HER2_IHC_CISH     N
# 1 Negative Negative        120
# 2 Negative Positive         93
# 3 Positive Positive         87

# dim(clin_finher)
# [1] 300 118
# All samples have ER/HER2 status available


clin_metabric %>%
  dplyr::group_by(er, her2) %>%
  dplyr::summarise(N = n())
#   er       her2         N
# 1 negative negative   134
# 2 negative positive    72
# 3 negative NA         230
# 4 positive negative   530
# 5 positive positive    74
# 6 positive NA         897
# 7 NA       negative    11
# 8 NA       positive     2
# 9 NA       NA         164

# dim(clin_metabric)
# [1] 2114 26
# Not all samples have ER/HER2 status available
# Metabric NA filtering
clin_metabric <- clin_metabric %>%
  dplyr::filter( !(is.na(er)| is.na(her2)) )
clin_metabric %>%
  dplyr::group_by(er, her2) %>%
  dplyr::summarise(N = n())
#   er       her2         N
# 1 negative negative   134
# 2 negative positive    72
# 3 positive negative   530
# 4 positive positive    74



# Compute module score (avg, wavg, wavg2)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

score <- list()
score[["avg"]] <- get_module_score(x = expr_finher,
                                   module_list = module_list_finher %>% purrr::map(~{.x$Direction = 1; .x}),
                                   by = "Ncbi_gene_id")
score[["wavg"]] <- get_module_score(x = expr_finher,
                                    module_list = module_list_finher,
                                    by = "Ncbi_gene_id")
score[["wavg2"]] <- get_module_score_2(x = expr_finher,
                                       module_list = module_list_finher,
                                       by = "Ncbi_gene_id")
score <- purrr::map(names(score),
                    function(nme,score){
                      xscore <- score[[nme]] %>%
                        dplyr::select(CEL_filename, ESR1, ERBB2)
                      names(xscore) <- str_c(names(xscore), "_", nme)
                      xscore
                      },
                    score)

score_finher <- bind_cols(score) %>%
  dplyr::select(-c(CEL_filename_wavg,CEL_filename_wavg2)) %>%
  dplyr::rename(CEL_filename = "CEL_filename_avg")

clin_finher <- clin_finher %>%
  dplyr::left_join(score_finher, by = "CEL_filename")



score <- list()
score[["avg"]] <- get_module_score(x = expr_neoadj,
                                   module_list = module_list_neoadj %>% purrr::map(~{.x$Direction = 1; .x}),
                                   by = "Ncbi_gene_id")
score[["wavg"]] <- get_module_score(x = expr_neoadj,
                                    module_list = module_list_neoadj,
                                    by = "Ncbi_gene_id")
score[["wavg2"]] <- get_module_score_2(x = expr_neoadj,
                                       module_list = module_list_neoadj,
                                       by = "Ncbi_gene_id")
score <- purrr::map(names(score),
                    function(nme,score){
                      xscore <- score[[nme]] %>%
                        dplyr::select(Sample_geo_accession, ESR1, ERBB2)
                      names(xscore) <- str_c(names(xscore), "_", nme)
                      xscore
                    },
                    score)

score_neoadj <- bind_cols(score) %>%
  dplyr::select(-c(Sample_geo_accession_wavg,Sample_geo_accession_wavg2)) %>%
  dplyr::rename(Sample_geo_accession = "Sample_geo_accession_avg")

clin_neoadj <- clin_neoadj %>%
  dplyr::left_join(score_neoadj, by = "Sample_geo_accession")




score <- list()
score[["avg"]] <- get_module_score(x = expr_metabric,
                                   module_list = module_list_metabric %>% purrr::map(~{.x$Direction = 1; .x}),
                                   by = "Ncbi_gene_id")
score[["wavg"]] <- get_module_score(x = expr_metabric,
                                    module_list = module_list_metabric,
                                    by = "Ncbi_gene_id")
score[["wavg2"]] <- get_module_score_2(x = expr_metabric,
                                       module_list = module_list_metabric,
                                       by = "Ncbi_gene_id")
score <- purrr::map(names(score),
                    function(nme,score){
                      xscore <- score[[nme]] %>%
                        dplyr::select(sample_name, ESR1, ERBB2)
                      names(xscore) <- str_c(names(xscore), "_", nme)
                      xscore
                    },
                    score)

score_metabric <- bind_cols(score) %>%
  dplyr::select(-c(sample_name_wavg,sample_name_wavg2)) %>%
  dplyr::rename(sample_name = "sample_name_avg")

clin_metabric <- clin_metabric %>%
  dplyr::left_join(score_metabric, by = "sample_name")




# ROC AUC comparison (pROC)
# >>>>>>>>>>>>>>>>>>

# FinHER
for(i in c("avg", "wavg", "wavg2")){

  er_formula = str_c("ER_IHC ~ ESR1_",i)
  her2_formula = str_c("HER2_IHC_CISH ~ ERBB2_",i)
  print(
    str_c(i,"; ER AUC:", roc(formula = as.formula(er_formula), data = clin_finher) %>% auc() %>% round(digits=2)
          ,", HER2 AUC:",roc(formula = as.formula(her2_formula), data = clin_finher) %>% auc() %>% round(digits=2))
  )
}
# [1] "avg; ER AUC:0.56, HER2 AUC:0.85"
# [1] "wavg; ER AUC:0.84, HER2 AUC:0.88"
# [1] "wavg2; ER AUC:0.86, HER2 AUC:0.81"
# wavg2: improved ER discrimination,
#         poor HER2 discrimination due to single gene in one direction

cor(clin_finher %>% dplyr::select(ESR1_avg,ESR1_wavg,ESR1_wavg2))
#             ESR1_avg ESR1_wavg ESR1_wavg2
# ESR1_avg   1.0000000 0.3775780  0.2649933
# ESR1_wavg  0.3775780 1.0000000  0.9929299
# ESR1_wavg2 0.2649933 0.9929299  1.0000000
cor(clin_finher %>% dplyr::select(ERBB2_avg,ERBB2_wavg,ERBB2_wavg2))
#             ERBB2_avg ERBB2_wavg ERBB2_wavg2
# ERBB2_avg   1.0000000  0.9818995   0.3713654
# ERBB2_wavg  0.9818995  1.0000000   0.5405013
# ERBB2_wavg2 0.3713654  0.5405013   1.0000000

module_list_finher$ESR1$Direction %>% table()
# -1  1
# 65 82 # minimum 65 genes in single direction

# Simialr number of up and down genes implies
# (here 44% down and 56% up genes)
#   1) significant difference between avg and wavg/wavg2
#   2) highly correlated wavg and wavg2
#   (Though wavg2 normalize seperatly for up and down genes,
#   due to similar number of up and down genes both wavg and wavg2 gives
#   highly correlated scores)

module_list_finher$ERBB2$Direction %>% table()
# -1  1
# 1 14 # minimum 1 gene in single direction

# Majority of gene in single direction imples
# (here > 90% up genes)
#   1) highly correlated avg and wavg
#   2) significant difference between wavg and wavg2
#   (As wavg2 normalize seperatly for up and down genes)




# Neoadj
for(i in c("avg", "wavg", "wavg2")){

  er_formula = str_c("Er ~ ESR1_",i)
  her2_formula = str_c("Her2 ~ ERBB2_",i)
  print(
    str_c(i,"; ER AUC:", roc(formula = as.formula(er_formula), data = clin_neoadj) %>% auc() %>% round(digits=2)
          ,", HER2 AUC:",roc(formula = as.formula(her2_formula), data = clin_neoadj) %>% auc() %>% round(digits=2))
  )
}
# [1] "avg; ER AUC:0.68, HER2 AUC:0.61"
# [1] "wavg; ER AUC:0.91, HER2 AUC:0.68"
# [1] "wavg2; ER AUC:0.91, HER2 AUC:0.72"
# wavg2: Simialr ER discrimination (due to similar number of up and down genes),
#         improved HER2 discrimination (minimum 3 genes in single direction)


cor(clin_neoadj %>% dplyr::select(ESR1_avg,ESR1_wavg,ESR1_wavg2))
#             ESR1_avg ESR1_wavg ESR1_wavg2
# ESR1_avg   1.0000000 0.3116655  0.2539653
# ESR1_wavg  0.3116655 1.0000000  0.9981905
# ESR1_wavg2 0.2539653 0.9981905  1.0000000
cor(clin_neoadj %>% dplyr::select(ERBB2_avg,ERBB2_wavg,ERBB2_wavg2))
#             ERBB2_avg ERBB2_wavg ERBB2_wavg2
# ERBB2_avg   1.0000000  0.9046988   0.4442825
# ERBB2_wavg  0.9046988  1.0000000   0.7836359
# ERBB2_wavg2 0.4442825  0.7836359   1.0000000


module_list_neoadj$ESR1$Direction %>% table()
# -1   1
# 150 192 # minimum 150 genes in single direction

# Simialr number of up and down genes implies
# (here 44% down and 56% up genes)
#   1) significant difference between avg and wavg/wavg2
#   2) highly correlated wavg and wavg2
#   (Though wavg2 normalize seperatly for up and down genes,
#   due to similar number of up and down genes both wavg and wavg2 gives
#   highly correlated scores)

module_list_neoadj$ERBB2$Direction %>% table()
# -1  1
# 3 17 # minimum 3 genes in single direction

# Majority of gene in single direction imples
# (here > 90% up genes)
#   1) highly correlated avg and wavg
#   2) significant difference between wavg and wavg2
#   (As wavg2 normalize seperatly for up and down genes)




# Metabric
for(i in c("avg", "wavg", "wavg2")){

  er_formula = str_c("er ~ ESR1_",i)
  her2_formula = str_c("her2 ~ ERBB2_",i)
  print(
    str_c(i,"; ER AUC:", roc(formula = as.formula(er_formula), data = clin_metabric) %>% auc() %>% round(digits=2)
          ,", HER2 AUC:",roc(formula = as.formula(her2_formula), data = clin_metabric) %>% auc() %>% round(digits=2))
  )
}
# [1] "avg; ER AUC:0.86, HER2 AUC:0.88"
# [1] "wavg; ER AUC:0.95, HER2 AUC:0.87"
# [1] "wavg2; ER AUC:0.95, HER2 AUC:0.85"
# wavg2: Simialr ER discrimination (due to similar number of up and down genes),
#         worse HER2 discrimination (minimum 3 genes in single direction)


cor(clin_metabric %>% dplyr::select(ESR1_avg,ESR1_wavg,ESR1_wavg2))
#             ESR1_avg ESR1_wavg ESR1_wavg2
# ESR1_avg   1.0000000 0.6830630  0.6546959
# ESR1_wavg  0.6830630 1.0000000  0.9992717
# ESR1_wavg2 0.6546959 0.9992717  1.0000000
cor(clin_metabric %>% dplyr::select(ERBB2_avg,ERBB2_wavg,ERBB2_wavg2))
#             ERBB2_avg ERBB2_wavg ERBB2_wavg2
# ERBB2_avg   1.0000000  0.9518865   0.8035897
# ERBB2_wavg  0.9518865  1.0000000   0.9473208
# ERBB2_wavg2 0.8035897  0.9473208   1.0000000


module_list_metabric$ESR1$Direction %>% table()
# -1   1
# 193 272 # minimum 193 genes in single direction

# Simialr number of up and down genes implies
# (here 42% down and 58% up genes)
#   1) significant difference between avg and wavg/wavg2
#   2) highly correlated wavg and wavg2
#   (Though wavg2 normalize seperatly for up and down genes,
#   due to similar number of up and down genes both wavg and wavg2 gives
#   highly correlated scores)

module_list_metabric$ERBB2$Direction %>% table()
# -1  1
# 5 23 # minimum 5 genes in single direction

# Majority of gene in single direction imples
# (here > 80% up genes)
#   1) highly correlated avg and wavg
#   2) significant difference between wavg and wavg2
#   !!!!!!!!!!!!!! here wavg and wavg2 are correlated !!!!!!!!!!!!!!!!!
#   (As wavg2 normalize seperatly for up and down genes)



# Summary
# >>>>>>>
# wavg and wavg2 gave similar results. However, independent normalization
# of up and down regulated genes in wavg2 is biologically more interpretable than
# wavg. i.e wavg2 score can be interpreted as the difference between
# average up and average down regulated genes.
# Hence wavg2 is used in this analysis.


#
# ==============================================================================

