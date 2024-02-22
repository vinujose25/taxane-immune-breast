# s0_treatment_curated_breast_dataset.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Dealing with some caveats in the manual curation of sample characteristics
# (clinical data) in the treatment-curated-breast-dataset project.
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# These caveats should be addressed in the original project before publishing.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load clinical data
# 2. Format clinical data
#   (age, grade, node, size,
#   response, follow-up,
#   arm,
#   subtype_ihc (er, pr, hr, her2),
#   subtype_pam50)
# 3. Save the curated clinical data under results/data !!!!




# 1. Load data
# ==============================================================================

clin <- read_tsv(file = "data/geo_clin_Quantile.tsv",
                 guess_max = 3700)
# guess_max =  3736 ?
# Due to large missingness in some variables, large no.of rows needs to read in
# to correctly guess the variable data type.

#
# ==============================================================================





# 2. Format clinical data
#
#   (age, grade, node, size,
#   response, follow-up,
#   arm,
#   subtype_ihc (er, pr, hr, her2),
#   subtype_pam50)
# ==============================================================================


# Age
# >>>>

clin %>%
  dplyr::group_by(Age = !is.na(Age), Age_bin, Age_detailed) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Age   Age_bin Age_detailed                 N
# 1 FALSE old     young:<=50,old>50           50 !!! integrate Age_bin and Age
# 2 FALSE young   young:<=50,old>50           70 !!! integrate Age_bin and Age
# 3 FALSE NA      NA                         791 No age data
# 4 TRUE  old     Young:age<=40;Old:age>40    27 ignore Age has data
# 5 TRUE  young   Young:age<=40;Old:age>40    17 ignore Age has data
# 6 TRUE  NA      NA                        2781 ignore Age has data


# Use a 50 cutoff with Age and by introduce a new variable
# Age_bin_50 to maximaize the samples with Age data
# (Note: The 120 Age_bin samples to be integrated with Age uses 50 cutoff)

clin <- clin %>%
  dplyr::mutate(
    Age_bin_50 = if_else(Age<=50, "young", "old"),
    Age_bin_50 = if_else(is.na(Age) & str_detect(Age_detailed,"50"),
                         Age_bin, Age_bin_50)
  )

clin$Age_bin_50 %>% table()
# old young
# 1374  1571


# Verification !!!!!!!!!!!

# The original Age_bin contains binarized age data using different cutoffs.
#
# Hence to avoid ambiguity, add a new variable Age_bin_50 that categorize age
# using 50 years as cutoff.
# Age_bin_50 to maximize the samples with Age data

# End of verification !!!!!!!!!!!





# Grade
# >>>>>

clin$Grade %>% is.na() %>% table()
# FALSE  TRUE
# 2284  1452
clin$Grade %>%  table()
# G1 G1/2   G2   G3
# 161   23  869 1231


# Integrate  Grade 1 and 2 to maximize usable data with Grade

clin <- clin %>%
  dplyr::mutate(
    Grade_bin = purrr::map_chr(
      Grade,
      ~switch(
        .x,
        "G1" = "G12",
        "G1/2" = "G12",
        "G2" = "G12",
        "G3" = "G3",
        NA
      )
    )
  )

clin$Grade_bin %>%  table()
# G12   G3
# 1053 1231


# Verification !!!!!!!!!!!!

# No need to introduce Grade_bin, leave it to the user

# End of verification !!!!!!!!!!!!



# Node
# >>>>

clin %>%
  dplyr::group_by(Node = !is.na(Node), Node_bin, Node_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Node  Node_bin Node_cat     N
# 1 FALSE neg      N0         537 ignore Node_bin has data
# 2 FALSE neg      NA         260 ignore Node_bin has data
# 3 FALSE pos      N1         664 ignore Node_bin has data
# 4 FALSE pos      N2         215 ignore Node_bin has data
# 5 FALSE pos      N3         142 ignore Node_bin has data
# 6 FALSE pos      NA         278 ignore Node_bin has data
# 7 FALSE NA       NA        1358 No node data
# 8 TRUE  neg      N0          49 ignore Node_bin has data
# 9 TRUE  neg      N1          55 ignore Node_bin has data
# 10 TRUE  neg      NA          22 ignore Node_bin has data
# 11 TRUE  pos      N0          37 ignore Node_bin has data
# 12 TRUE  pos      N1          66 ignore Node_bin has data
# 13 TRUE  pos      N2          17 ignore Node_bin has data
# 14 TRUE  pos      NA          33 ignore Node_bin has data
# 15 TRUE  NA       NA           3 !!! Explore Node counts must be categorized

clin$Node[ (!is.na(clin$Node)) & is.na(clin$Node_bin)]
# [1] "N-" "N-" "N-"

# Node_bin contains maximum number of nodal information.
# Explore why the three samples with Node counts are not categorized.


table(clin$Node_bin)
# neg  pos
# 923 1452

# Verification !!!!!!!!!!!!!

# This anomaly is from GSE4779.
# In the original file N- means NA

# End of verification !!!!!!!!!!!!!



# Size
# >>>>

clin %>%
  dplyr::group_by(Size = !is.na(Size), Size_bin, Size_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Size  Size_bin Size_cat     N
# 1 FALSE large    T2        1095 ignore Size_bin has data
# 2 FALSE large    T3         453 ignore Size_bin has data
# 3 FALSE large    T4         238 ignore Size_bin has data
# 4 FALSE samall   T0           3 ignore Size_bin has data
# 5 FALSE samall   T1         131 ignore Size_bin has data
# 6 FALSE small    T0           6 ignore Size_bin has data
# 7 FALSE small    T1          62 ignore Size_bin has data
# 8 FALSE NA       NA         745 No size data
# 9 TRUE  large    T2         461 ignore Size_bin has data
# 10 TRUE  large    T3         251 ignore Size_bin has data
# 11 TRUE  large    NA          35 ignore Size_bin has data
# 12 TRUE  small    T1         233 ignore Size_bin has data
# 13 TRUE  small    NA          20 ignore Size_bin has data
# 14 TRUE  NA       NA           3 !!! Explore Size must be categorized


clin$Size[ (!is.na(clin$Size)) & is.na(clin$Size_bin)]
# [1] "T-" "T-" "T-"

# Size_bin cotains maximum number of nodal information.
# Explore why the three samples with Size are not categorized.
# Correct the clerical error in Size_bin

clin <- clin %>%
  dplyr::mutate(
    Size_bin = str_replace(Size_bin, "samall", "small")
  )

table(clin$Size_bin)
# large small
# 2533   455
# small = T0-1, large = T2-3-4 (>T1)

# Verification !!!!!!!!!!!!!

# This anomaly is from GSE4779.
# In the original file T- means NA

# End of verification !!!!!!!!!!!!!



# Response
# >>>>>>>>

clin %>%
  dplyr::group_by(Response = !is.na(Response),
                  Response_path = !is.na(Response_pathological),
                  Response_clin = !is.na(Response_clinical)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Response Response_path Response_clin     N
# 1 FALSE    FALSE         FALSE           848 No response data
# 2 FALSE    TRUE          FALSE           488 !!! Integrate Response path to Response
# 3 TRUE     FALSE         TRUE            242 ignore Response has data
# 4 TRUE     TRUE          FALSE          2158 ignore Response has data


clin$Response %>% table()
# npCR  pCR
# 1755  645

clin$Response_pathological[is.na(clin$Response) &!is.na(clin$Response_pathological)] %>%
  table()
# npCR  pCR
# 389   99
# Sum = 488,

# Integrate Response_pathological to missing Response data without recoding
clin <- clin %>%
  dplyr::mutate(
    Response = if_else(is.na(Response), Response_pathological, Response)
  )

clin$Response %>% table()
# npCR  pCR
# 2144  744

# Verification !!!!!!!!!!!!!

# The 488 samples with Response = NA and Response_path = TRUE is due to
# clerical error in the manual curation step of GSE25066.
# The manual curation step failed to update Response = Response_path in GSE25066.
# This is updated in the manual curation code related to GSE25066 in the following file;
# ~/projects/on/treatment-curated-breast-dataset/src/r/s3_manual_sample_characterisitic_curation.R

# End of verification !!!!!!!!!!!!!




# Follow up
# >>>>>>>>>

clin %>%
  dplyr::group_by(Event_dfs = !is.na(Event_dfs),
                  Time_dfs = !is.na(Time_dfs)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Event_dfs Time_dfs     N
# 1 FALSE     FALSE     2344 ! No followup data
# 2 FALSE     TRUE         7 ! Missing event #  "GSE16391"; Verified; Missing in original file.
# 3 TRUE      FALSE      108 ! Missing time # "GSE16446"; Verified; 7 samples have no time data
#                                             "GSE19615"; Verified; 101 samples have no time data
# 4 TRUE      TRUE      1277 Consider only samples with both Event and Time data


# Consider only samples with both Event and Time data
# !!! Explore why Time data has no matching Event data and vice-versa


# Verification !!!!!!!!!!!!!!!

# 2 FALSE     TRUE         7 ! Missing event #  "GSE16391"; Verified; Missing in original file.
# 3 TRUE      FALSE      108 ! Missing time # "GSE16446"; Verified; 7 samples have no time data
#                                             "GSE19615"; Verified; 101 samples have no time data

# End of verification !!!!!!!!!




# Arm
# >>>

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



# Subtype_ihc
# >>>>>>>>>>>

clin %>%
  dplyr::group_by(Er = !is.na(Er),
                  Pr = !is.na(Pr),
                  Hr = !is.na(Hr)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Er    Pr    Hr        N
# 1 FALSE FALSE FALSE   539 No hormone receptor status
# 2 TRUE  FALSE TRUE    325 !!! Hormon reeptor status without Pr status
# 3 TRUE  TRUE  TRUE   2872 Clean hormone receptor status


# Exploring: Hormon reeptor status without Pr status
clin[is.na(clin$Pr) & !is.na(clin$Hr), ] %>%
  dplyr::group_by(Er, Pr, Hr) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Er    Pr    Hr        N
# 1 neg   NA    neg     244
# 2 pos   NA    pos      81

# Identical Er and Hr status suggests that for some samples Hr status is
# defined purely based on Er status.
# Of note, some studies only have Hr status, not Er and Pr status.
# !!!!!! Verify this in treatment-curated-breast-dataset !!!!!!!!!!!!


clin %>%
  dplyr::group_by(Hr = !is.na(Hr),
                  Her2 = !is.na(Her2),
                  Subtype_ihc = !is.na(Subtype_ihc)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Hr    Her2  Subtype_ihc     N
# 1 FALSE FALSE FALSE         444
# 2 FALSE TRUE  FALSE          95
# 3 TRUE  FALSE FALSE         359
# 4 TRUE  TRUE  TRUE         2838 Clean Subtype IHC

# Subtype_ihc is defined for samples with both Hr and Her2 status


# Verification !!!!!!!!!!!!!!!!

x <- clin[!is.na(clin$Er) & is.na(clin$Pr),]
dim(x)  # 325 93
x$Series_matrix_accession %>% table()
# GSE16446 GSE22093 GSE23988 GSE25066 GSE41998 GSE45255  GSE4779  GSE6861 GSE69031
# 120       98       61        1        1        5        2       36        1

x <- clin[clin$Series_matrix_accession %in% c("GSE16446", "GSE22093", "GSE23988", "GSE25066", "GSE41998", "GSE45255", "GSE4779",  "GSE6861", "GSE69031"), ]
x$Series_matrix_accession %>% table()
# GSE16446 GSE22093 GSE23988 GSE25066 GSE41998 GSE45255  GSE4779  GSE6861 GSE69031
# 120      103       61      508      279      139      102      161      118


# Datasets in which 325 samples have only Er status and no Pr status available.

# GSE16446: PR status not available
# GSE22093: PR status not available; ER available for 98
# GSE23988: PR status not available
# GSE25066: PR set as NA for one more sample than taht for ER status (Intermediate score = NA)
# GSE41998: PR unknown for one sample
# GSE45255: PR unknown for five sample
# GSE4779: PR unknown for two sample
# GSE6861: PR unknown for 36 samples; Er=neg for all samples from literature
# GSE69031: PR unknown for two samples; 12 samples were removed during
#           expression data integration due to NA expression (130-12=118),
#           which may include one sample with PR = NA.

# End verification !!!!!!!!!!!!!!!!




# Subtype_pam50
# >>>>>>>>>>>>>

clin$Subtype_pam50 %>% is.na() %>% table()
# FALSE
# 3736

clin$Subtype_pam50 %>% table()
# Basal         Her2         LumA         LumB       Normal
# 1047          636         1066          566          401
# Unclassified
# 20

# Clean Subtype_pam50.
# Subtype_pam50 is defined for all samples.



# Subtype_ihc vs Subtype_pam50
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin %>%
  dplyr::group_by(Subtype_ihc = !is.na(Subtype_ihc),
                  Subtype_pam50 = !is.na(Subtype_pam50)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Subtype_ihc Subtype_pam50     N
# 1 FALSE       TRUE            898 Missing Subtype_ihc
# 2 TRUE        TRUE           2838


#
# ==============================================================================





# 3. Save the curated clinical data under results/data !!!!
# ==============================================================================

geo_clin_curated_locally <- clin
save(geo_clin_curated_locally,
     file = str_c(out_data, "geo_clin_curated_locally.RData"))

#
# ==============================================================================


