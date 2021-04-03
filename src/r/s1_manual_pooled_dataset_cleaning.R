# s1_manual_pooled_dataset_cleaning.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Filter pooled dataset based on the criteria specified in the text.
# The pooled dataset is filtered using the following criteria;
# a) Discard samples with missing IHC status for HR and HER2
#
# For each IHC subtype (HR,HER2,TN):
# b) Discard samples with missing treatment data
# c) Discard samples with missing clinical endpoint (pCR/DFS) data
# d) Discard matched samples from longitudnal studies (series matrices)
#   (i.e. consider only pre-treatment expression profiles),
# e) Discard treatment regimens with samples only from a single study,
# f) Discard treatment regimens with sample size less than 50, and
# g) Discard studies from a treatment regimen with less than 10 samples.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and prepare data
# 2. Clinical data filtering (filter expression data later)
# 3. Summarize neoadj/adj dataset
# 4. Expression data filtering
# 5. Additional formating of clinincal data to aid in analysis
# 6. Save filtered clinical and expression data



# 1. Load and prepare data
# ==============================================================================

load("results/data/geo_clin_curated_locally.RData")
clin <- geo_clin_curated_locally
rm(geo_clin_curated_locally)

expr <- read_tsv(file = "data/geo_expr_Quantile.tsv")


dim(clin)
# [1] 3736   96
dim(expr)
# [1] 9184 3737


# Making expression and clinical data congruent
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cat(
  "Congruent expression and clinical data? -",
  if_else(
    identical(names(expr)[-1], clin$Sample_geo_accession),
    "Yes.\n",
    "No.\n"
  )
)

if(! identical(names(expr)[-1], clin$Sample_geo_accession)){

  cat("Making expression and clinical data congruent.\n")
  expr <- expr %>%
    dplyr::select(1, clin$Sample_geo_accession)
}


# Review clinical endpoint data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin %>%
  dplyr::group_by(Response = !is.na(Response),
                  Response_path = !is.na(Response_pathological),
                  Response_clin = !is.na(Response_clinical)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Response Response_path Response_clin     N
# 1 FALSE    FALSE         FALSE           848 No response data
# 2 TRUE     FALSE         TRUE            242
# 3 TRUE     TRUE          FALSE          2646

# Response is clean
# Both clinical and pathological response is integrated into Response.



clin %>%
  dplyr::group_by(Event_dfs = !is.na(Event_dfs),
                  Time_dfs = !is.na(Time_dfs)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Event_dfs Time_dfs     N
# 1 FALSE     FALSE     2344  No followup data
# 2 FALSE     TRUE         7 ! Missing event
# 3 TRUE      FALSE      108 ! Missing time
# 4 TRUE      TRUE      1277

# !!!!!!!!!!!!!!!!!
# For this analysis
# Set event as NA for samples with missing time data and vice-versa.

clin <- clin %>%
  dplyr::mutate(
    Event_dfs = if_else(is.na(Time_dfs), NA_real_, Event_dfs),
    Time_dfs = if_else(is.na(Event_dfs), NA_real_, Time_dfs)
  )

clin %>%
  dplyr::group_by(Event_dfs = !is.na(Event_dfs),
                  Time_dfs = !is.na(Time_dfs)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Event_dfs Time_dfs     N
# 1 FALSE     FALSE     2459
# 2 TRUE      TRUE      1277



# Consolidate Arm for easy filtering
# (Arm_chemo + Arm_her2 + Arm_hr + Arm_other)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin %>%
  dplyr::group_by(Arm_chemo = !is.na(Arm_chemo),
                  Arm_her2 = !is.na(Arm_her2),
                  Arm_hormone = !is.na(Arm_hormone),
                  Arm_other = !is.na(Arm_other)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Arm_chemo Arm_her2 Arm_hormone Arm_other     N
# 1 FALSE     FALSE    FALSE       FALSE       448 ! No treatment data
# 2 TRUE      TRUE     TRUE        TRUE       3288

# Arm is clean

clin <- clin %>%
  dplyr::mutate(
    Arm_consolidated = str_c(clin$Arm_chemo, "///",
                             clin$Arm_her2, "///",
                             clin$Arm_hormone, "///",
                             clin$Arm_other)
  )

clin$Arm_consolidated %>% table() %>% sum()  # 3288  !!! NAs present

#
# ==============================================================================




# 2. Clinical data filtering (filter expression data later)
# ==============================================================================


# To implement per subtype sample/study filtering
# Later in the script will merge this list to a single tibble
clin_neoadj <- vector(mode = "list", length = 3)
clin_adj <- vector(mode = "list", length = 3)

names(clin_neoadj) <- c("HR","HER2","TN")
names(clin_adj) <- c("HR","HER2","TN")



clin %>%
  dplyr::group_by(Regimen_updated, Subtype_ihc) %>%
  dplyr::summarise(N= n(), .groups = "keep")
#   Regimen_updated Subtype_ihc     N
# 1 adj             HER2          109
# 2 adj             HR            198
# 3 adj             TN             63
# 4 adj             NA            384
# 5 neoadj          HER2          647
# 6 neoadj          HR           1034
# 7 neoadj          TN            787
# 8 neoadj          NA            514


# neoadj filtering
#
# a) Discard samples with missing IHC status for HR, HER2, and TN
# For each IHC subtype (HR,HER2,TN):
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_neoadj$HR <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "neoadj" & Subtype_ihc == "HR"),
  type = "neoadj"
)
# neoadj : Started with 1034 no.of samples.
# neoadj : Discarding 102 samples with missing treatment data.
# neoadj : Discarding 18 samples with missing clinical response data.
# neoadj : Discarding 0 matched samples from longitudnal studies.
# neoadj : Discarding 71 samples from 13 series matrices with <20 samples per treatment regimen.
# neoadj : Discarding 32 samples from 1 treatment regimens with <50 samples.
# neoadj : Discarding 60 samples from 1 treatment regimens unique to a single study.
# neoadj : End with 751 no.of samples.



clin_neoadj$HER2 <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "neoadj" & Subtype_ihc == "HER2"),
  type = "neoadj"
)
# neoadj : Started with 647 no.of samples.
# neoadj : Discarding 68 samples with missing treatment data.
# neoadj : Discarding 16 samples with missing clinical response data.
# neoadj : Discarding 0 matched samples from longitudnal studies.
# neoadj : Discarding 103 samples from 16 series matrices with <20 samples per treatment regimen.
# neoadj : Discarding 188 samples from 6 treatment regimens with <50 samples.
# neoadj : Discarding 50 samples from 1 treatment regimens unique to a single study.
# neoadj : End with 222 no.of samples.


clin_neoadj$TN <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "neoadj" & Subtype_ihc == "TN"),
  type = "neoadj"
)
# neoadj : Started with 787 no.of samples.
# neoadj : Discarding 52 samples with missing treatment data.
# neoadj : Discarding 37 samples with missing clinical response data.
# neoadj : Discarding 0 matched samples from longitudnal studies.
# neoadj : Discarding 55 samples from 11 series matrices with <20 samples per treatment regimen.
# neoadj : Discarding 58 samples from 2 treatment regimens with <50 samples.
# neoadj : Discarding 59 samples from 1 treatment regimens unique to a single study.
# neoadj : End with 526 no.of samples.





# adj filtering
#
# a) Discard samples with missing IHC status for HR, HER2 and TN
# For each IHC subtype (HR,HER2,TN):
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


clin_adj$HR <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "adj" & Subtype_ihc == "HR"),
  type = "adj"
)
# adj : Started with 198 no.of samples.
# adj : Discarding 4 samples with missing treatment data.
# adj : Discarding 52 samples with missing clinical follow-up data.
# adj : Discarding 0 matched samples from longitudnal studies.
# adj : Discarding 97 samples from 20 series matrices with <20 samples per treatment regimen.
# adj : Discarding 45 samples from 1 treatment regimens with <50 samples.
# adj : Discarding 0 samples from 0 treatment regimens unique to a single study.
# adj : End with 0 no.of samples.


clin_adj$HER2 <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "adj" & Subtype_ihc == "HER2"),
  type = "adj"
)
# adj : Started with 109 no.of samples.
# adj : Discarding 7 samples with missing treatment data.
# adj : Discarding 30 samples with missing clinical follow-up data.
# adj : Discarding 0 matched samples from longitudnal studies.
# adj : Discarding 52 samples from 18 series matrices with <20 samples per treatment regimen.
# adj : Discarding 20 samples from 1 treatment regimens with <50 samples.
# adj : Discarding 0 samples from 0 treatment regimens unique to a single study.
# adj : End with 0 no.of samples.


clin_adj$TN <- filter_data(
  x = clin %>%
    dplyr::filter(Regimen_updated == "adj" & Subtype_ihc == "TN"),
  type = "adj"
)
# adj : Started with 63 no.of samples.
# adj : Discarding 4 samples with missing treatment data.
# adj : Discarding 20 samples with missing clinical follow-up data.
# adj : Discarding 0 matched samples from longitudnal studies.
# adj : Discarding 39 samples from 15 series matrices with <20 samples per treatment regimen.
# adj : Discarding 0 samples from 0 treatment regimens with <50 samples.
# adj : Discarding 0 samples from 0 treatment regimens unique to a single study.
# adj : End with 0 no.of samples.




# Merging the clin_neoadj/clin_adj (list) as a single tibble.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_neoadj <- bind_rows(clin_neoadj)
# clin_adj <- bind_rows(clin_adj)
# !!!!!!! No dataset from adjuvant setting satisfying filterng criteria


#
#===============================================================================




# 3. Summarize neoadj/adj dataset
# ==============================================================================

# Note !!!
# No adj dataset cleared filtering criteria
# Summary of neoadj dataset is given below



# Clinical-vars, subtype, treatment and response data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

xsum <- clin_neoadj %>%
  dplyr::group_by(Age = !is.na(Age_bin_50),
                  Grade = !is.na(Grade_bin),
                  Node = !is.na(Node_bin),
                  Size = !is.na(Size_bin),
                  Subtype = !is.na(Subtype_ihc),
                  Arm = !is.na(Arm_consolidated),
                  pCR = !is.na(Response)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n(), .groups = "keep")
#   Age   Grade Node  Size  Subtype Arm   pCR       N
# 1 FALSE TRUE  FALSE TRUE  TRUE    TRUE  TRUE      1
# 2 TRUE  FALSE FALSE FALSE TRUE    TRUE  TRUE      2
# 3 TRUE  FALSE FALSE TRUE  TRUE    TRUE  TRUE    229
# 4 TRUE  FALSE TRUE  FALSE TRUE    TRUE  TRUE      1
# 5 TRUE  FALSE TRUE  TRUE  TRUE    TRUE  TRUE     82
# 6 TRUE  TRUE  FALSE FALSE TRUE    TRUE  TRUE    138
# 7 TRUE  TRUE  FALSE TRUE  TRUE    TRUE  TRUE     50
# 8 TRUE  TRUE  TRUE  TRUE  TRUE    TRUE  TRUE    996


xsum$N %>% sum() # 1499; 1499 - 996 = 503
# !!!!!!!!!!!!!!!!!!!
# 996 samples got all relevant clinical-pathological variables

xsum <- clin_neoadj %>%
  dplyr::filter(!is.na(Age_bin_50) &
                  !is.na(Grade_bin) &
                  !is.na(Node_bin) &
                  !is.na(Size_bin) &
                  !is.na(Subtype_ihc) &
                  !is.na(Arm_consolidated) &
                  !is.na(Response)) %>%
  dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n(), .groups = "keep") %>%
  dplyr::filter(N >= 50)
#   Arm_consolidated                                                     Subtype_ihc     N
# 1 AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy HR             58
# 2 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HER2           82
# 3 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HR            486
# 4 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   TN            269

# Potential analyses !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# If age, grade, node, size, subtype, arm, and pCR needs to be used,
# the analysis that can be performed are
# 1) Interaction between immune, subtype and pCR in AAA+Taxane regimen.
# 2) Prognosis of immune in AAA+Taxane (pan subtype)

# 3) Interaction between immune, AAA+/-taxane and pCR in HR subtype.
# 4) Prognosis of immune in HR (pan tretment regimen)

# 5) Prognosis of immune in BC (pooled dataset)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# Subtype, treatment and response data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

xsum <-clin_neoadj %>%
  dplyr::group_by(Subtype = !is.na(Subtype_ihc),
                  Arm = !is.na(Arm_consolidated),
                  pCR = !is.na(Response)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n(), .groups = "keep")
#   Subtype Arm   pCR       N
# 1 TRUE    TRUE  TRUE   1499

# !!!!!!!!!!!!!!!!!!!
# 1499 samples got all subtype, arm and pCR data

xsum <- clin_neoadj %>%
  dplyr::filter(!is.na(Subtype_ihc) &
                  !is.na(Arm_consolidated) &
                  !is.na(Response)) %>%
  dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
  # dplyr::group_by(Arm_consolidated, Subtype_ihc, Hr) %>% # new
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n(), .groups = "keep") %>%
  dplyr::filter(N >= 50)
#   Arm_consolidated                                                     Subtype_ihc     N
# 1 A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HR            139
# 2 A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   TN            150

# 3 AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy HR             91
# 4 AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy TN             83

# 5 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HER2          135
# 6 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HR            521
# 7 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   TN            293
# 8 AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy     HER2           87


# Potential analyses !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# If subtype, arm and pCR needs to be used,
# the analysis that can be performed are

# TN
# 1) Interaction between immune, AAA(+/-taxane) and pCR in TN subtype.
# 2) Interaction between immune, Taxane+(AAA/A0A) and pCR in TN subtype.
# 3) Prognosis of immune in TN (pan treatment)

# HER2 (Note include HR stratification in the anlysis)
# 4) Interaction between immune, AAA+Taxane(+/-TRA). and pCR in HER2 subtype.
# 5) Prognosis of immune in HER2 (pan treatment)

# HR
# 6) Interaction between immune, AAA+-Taxane and pCR in HR
# 7) Interaction between immune, Taxane+(AAA/A0A) and pCR in HR
# 8) Prognosis of immune in HR (pan treatment)

# All
# 9) Interaction between immune, subtype and pCR in AAA+Taxane, AAA+noTaxane, and AOA+Taxane regimen.
# 10) Prognosis of immune in breast cancer (pan subtype, pan treatment)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




# Neoadjuvant analysis summary
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

# TN
# immune * Taxane interaction on pCR in AAA regimen
# immune * AAA/A0A interaction on pCR in Taxane containing regimen
# prognostic value of immune irrespective of treatment regimens (pan treatment)

# HER2 (Note include HR stratification in the anlysis)
# immune * Trastuzumab interaction on pCR in AAA+Taxane regimen
# prognostic value of immune irrespective of treatment regimens (pan treatment)

# HR
# immune * Taxane interaction on pCR in AAA regimen
# immune * AAA/A0A interaction on pCR in Taxane containing regimen
# prognostic value of immune irrespective of treatment regimens (pan treatment)

# All
# immune * subtype interacton on pCR in AAA+Taxane, AAA+noTaxane, A0A+Taxane
# prognosis of immune in breast cancer (pan subtype, pan treatment)


# Summary of HR stratified HER2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xsum <- clin_neoadj %>%
  dplyr::filter(!is.na(Subtype_ihc) &
                  !is.na(Arm_consolidated) &
                  !is.na(Response)) %>%
  dplyr::filter(Subtype_ihc == "HER2") %>%
  dplyr::group_by(Arm_consolidated, Subtype_ihc, Hr) %>%
  dplyr::summarise(N = n(), .groups = "keep")
#   Arm_consolidated                                              Subtype_ihc Hr        N
# 1 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_t… HER2        neg      80
# 2 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_t… HER2        pos      55
# 3 AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_the… HER2        neg      58
# 4 AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_the… HER2        pos      29

# Note !!!
# Include transtuzumab interaction analysis in each HR+ and HR- HER2 amplified
# subtype as supplementary.

#
#===============================================================================



# # Old version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # In new version no adjuvant dataset cleared filtering criteria
# # 4. Summarize adj dataset
# # ==============================================================================
#
#
# # Clinical-vars, subtype, treatment and response data
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# xsum <- clin_adj %>%
#   dplyr::group_by(Age = !is.na(Age_bin_50),
#                   Grade = !is.na(Grade_bin),
#                   Node = !is.na(Node_bin),
#                   Size = !is.na(Size_bin),
#                   Subtype = !is.na(Subtype_ihc),
#                   Arm = !is.na(Arm_consolidated),
#                   Event = !is.na(Event_dfs)) %>%
#   # TRUE = data, FALSE = no data, NA = no data
#   dplyr::summarise(N = n(), .groups = "keep")
# #   Age   Grade Node  Size  Subtype Arm   Event     N
# # 1 TRUE  TRUE  TRUE  TRUE  TRUE    TRUE  TRUE     55
#
# # 55 samples got all relevant clinical-pathological variables
#
# xsum <- clin_adj %>%
#   dplyr::filter(!is.na(Age_bin_50) &
#                   !is.na(Grade_bin) &
#                   !is.na(Node_bin) &
#                   !is.na(Size_bin) &
#                   !is.na(Subtype_ihc) &
#                   !is.na(Arm_consolidated) &
#                   !is.na(Event_dfs)) %>%
#   dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
#   # TRUE = data, FALSE = no data, NA = no data
#   dplyr::summarise(N = n(), .groups = "keep") %>%
#   dplyr::filter(N >= 50)
# #   Arm_consolidated                                            Subtype_ihc     N
# # 1 000+noTaxane///No_her2_agent///Tamoxifen///No_other_therapy HR             55
#
# # Potential analyses !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # If age, grade, node, size, subtype, arm, and followup needs to be used,
# # the analysis that can be performed are
# # 1) Prognosis in Tamoxifen arm in HR
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#
# # Subtype, treatment and response data
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# xsum <- clin_adj %>%
#   dplyr::group_by(Subtype = !is.na(Subtype_ihc),
#                   Arm = !is.na(Arm_consolidated),
#                   Event = !is.na(Event_dfs)) %>%
#   # TRUE = data, FALSE = no data, NA = no data
#   dplyr::summarise(N = n(), .groups = "keep")
# #   Subtype Arm   Event     N
# # 1 TRUE    TRUE  TRUE     55
#
# # 55 samples got all subtype, arm and followup data
#
# xsum <- clin_adj %>%
#   dplyr::filter(!is.na(Subtype_ihc) &
#                   !is.na(Arm_consolidated) &
#                   !is.na(Event_dfs)) %>%
#   dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
#   # TRUE = data, FALSE = no data, NA = no data
#   dplyr::summarise(N = n(), .groups = "keep") %>%
#   dplyr::filter(N >= 50)
# #   Arm_consolidated                                            Subtype_ihc     N
# # 1 000+noTaxane///No_her2_agent///Tamoxifen///No_other_therapy HR             55
#
#
# # Potential analyses !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # If subtype, arm and followup needs to be used,
# # the analysis that can be performed are
# # 1) Prognosis in Tamoxifen arm in HR
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#
# # Adjuvant analysis summary
# # >>>>>>>>>>>>>>>>>>>>>>>>>
# # Discard Adjuvant dataset, as we cannot check the interaction between
# # treditional immunologic/immunomodulatory treatment and tumor-immune-respone on
# # clinical response.
#
# #
# #===============================================================================




# 4. Expression data filtering
# ==============================================================================

expr_neoadj <- expr %>%
      dplyr::select(1, all_of(clin_neoadj$Sample_geo_accession))


#
# ==============================================================================



# 5. Additional formating of clinincal data to aid in analysis
# ==============================================================================

clin_neoadj <- clin_neoadj %>%
  dplyr::mutate(
    # Recoding "Response"
    Response = if_else(Response == "pCR", 1, 0),

    # Defining starta as the combination of series matrix accession and arm.
    # Clean strata by removing redundant arm information
    # noTaxane, No_her2_agent, No_hormone_therapy, No_other_therapy to ""
    Strata = str_c(Series_matrix_accession,"///", Arm_consolidated) %>%
      str_replace("\\+noTaxane", "") %>%
      str_replace("///No_her2_agent", "") %>%
      str_replace("///No_hormone_therapy", "") %>%
      str_replace("///No_other_therapy", "") %>%
      str_replace("///Trastuzumab", "+TRA") %>%
      str_replace("\\+Taxane", "+T") %>%
      str_replace("///", ":"), # GSE..///AAA > GSE...:AAA

    # Added on 3/4/21
    # Integrate HR status in the the treatment strata for HER2
    # Objective:
    # 1) To asses HR status based difference in response
    # 2) Usually in clinics HR+ve samples were given hormone therapy
    #     irrespective of HER2 status.
    Strata = if_else(Subtype_ihc == "HER2", str_c(Strata,":",Hr), Strata) %>%
      str_replace("pos","HR+") %>%
      str_replace("neg", "HR-")
  )

clin_neoadj %>%
  dplyr::group_by(Subtype_ihc, Hr, Strata) %>%
  dplyr::summarise(N = n(), Pcr = Response %>% sum()) %>%
  as.data.frame()
#    Subtype_ihc  Hr                 Strata   N Pcr
# 1         HER2 neg     GSE20194:AAA+T:HR-  25  13
# 2         HER2 neg     GSE32646:AAA+T:HR-  18   9
# 3         HER2 neg GSE42822:AAA+T+TRA:HR-  15  10
# 4         HER2 neg     GSE50948:AAA+T:HR-  37   9
# 5         HER2 neg GSE50948:AAA+T+TRA:HR-  43  25
# 6         HER2 pos     GSE20194:AAA+T:HR+  25   4
# 7         HER2 pos     GSE32646:AAA+T:HR+  16   3
# 8         HER2 pos GSE42822:AAA+T+TRA:HR+   9   2
# 9         HER2 pos     GSE50948:AAA+T:HR+  14   4
# 10        HER2 pos GSE50948:AAA+T+TRA:HR+  20   6
# 11          HR pos         GSE20194:AAA+T 141   9
# 12          HR pos           GSE20271:AAA  49   3
# 13          HR pos         GSE20271:AAA+T  44   3
# 14          HR pos           GSE22093:AAA  42  10
# 15          HR pos         GSE23988:AAA+T  32   7
# 16          HR pos         GSE25066:A0A+T  35   3
# 17          HR pos         GSE25066:AAA+T 194  22
# 18          HR pos         GSE32646:AAA+T  55   5
# 19          HR pos         GSE41998:A0A+T 104  12
# 20          HR pos         GSE42822:AAA+T  29   8
# 21          HR pos         GSE50948:AAA+T  26   4
# 22          TN neg         GSE20194:AAA+T  64  25
# 23          TN neg           GSE20271:AAA  28   4
# 24          TN neg         GSE20271:AAA+T  31   9
# 25          TN neg           GSE22093:AAA  55  18
# 26          TN neg         GSE23988:AAA+T  29  13
# 27          TN neg         GSE25066:A0A+T  25   6
# 28          TN neg         GSE25066:AAA+T 119  43
# 29          TN neg         GSE32646:AAA+T  26  10
# 30          TN neg         GSE41998:A0A+T 125  51
# 31          TN neg         GSE42822:AAA+T  24  12
#
# ==============================================================================




# 6. Save filtered clinical and expression data
# ========================================================================= =====


save(clin_neoadj, file = str_c(out_data,"clin_neoadj.RData"))
save(expr_neoadj, file = str_c(out_data,"expr_neoadj.RData"))

rm(clin, clin_adj, expr)
#
# ==============================================================================

