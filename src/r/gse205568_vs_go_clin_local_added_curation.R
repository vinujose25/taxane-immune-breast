

# This script check the agreement between the public version of pooled clinical data
# and the in house version used in the taxane-immune-breast project.

# in house version
# >>>>>>>>>>>>>>>>
load("results/data/geo_clin_curated_locally.RData")
clin_mta <- geo_clin_curated_locally

# public version
# >>>>>>>>>>>>>>>
clin <- read_tsv(file = "data/tmp/GSE205568_quantile_corrected/quantile_corrected/geo_clin_Quantile.tsv",
                 guess_max = 3700)
clin_gse205568 <- clin


# Comparison
# >>>>>>>>>>>


glimpse(clin_mta)
glimpse(clin_gse205568)


dim(clin_mta) # 3736   95
dim(clin_gse205568) # 3736   93; missing Age_bin_50 and Grade_bin, both irrelevant, users can recode it


# Size
# ====

clin_mta %>%
  dplyr::group_by(Size = is.na(Size), Size_bin, Size_cat) %>%
  dplyr::summarise(N = n())
#   Size  Size_bin Size_cat     N
# 1 FALSE large    T2         461
# 2 FALSE large    T3         251
# 3 FALSE large    NA          35
# 4 FALSE small    T1         233
# 5 FALSE small    NA          20
# 6 FALSE NA       NA           3
# 7 TRUE  large    T2        1095
# 8 TRUE  large    T3         453
# 9 TRUE  large    T4         238
# 10 TRUE  small    T0           9
# 11 TRUE  small    T1         193
# 12 TRUE  NA       NA         745


clin_gse205568 %>%
  dplyr::group_by(Size = is.na(Size), Size_bin, Size_cat) %>%
  dplyr::summarise(N = n())
#   Size  Size_bin Size_cat     N
# 1 FALSE large    T2         461
# 2 FALSE large    T3         251
# 3 FALSE large    NA          35
# 4 FALSE small    T1         233
# 5 FALSE small    NA          20
# 6 FALSE NA       NA           3
# 7 TRUE  large    T2        1095
# 8 TRUE  large    T3         453
# 9 TRUE  large    T4         238
# 10 TRUE  samall   T0           3
# 11 TRUE  samall   T1         131
# 12 TRUE  small    T0           6
# 13 TRUE  small    T1          62
# 14 TRUE  NA       NA         745



# Verdict: NOT clean !!!!  Problem with size_bin = small vs samall (clerical error)


# Node
# ====

clin_mta %>%
  dplyr::group_by(Node = is.na(Node), Node_bin, Node_cat) %>%
  dplyr::summarise(N = n())
#   Node  Node_bin Node_cat     N
# 1 FALSE neg      N0          49
# 2 FALSE neg      N1          55 !!!!! wrong, correct node_bin=pos
# 3 FALSE neg      NA          22
# 4 FALSE pos      N0          37 !!!!! wrong, correct node_bin=neg
# 5 FALSE pos      N1          66
# 6 FALSE pos      N2          17
# 7 FALSE pos      NA          33
# 8 FALSE NA       NA           3
# 9 TRUE  neg      N0         537
# 10 TRUE  neg      NA         260
# 11 TRUE  pos      N1         664
# 12 TRUE  pos      N2         215
# 13 TRUE  pos      N3         142
# 14 TRUE  pos      NA         278
# 15 TRUE  NA       NA        1358

clin_gse205568 %>%
  dplyr::group_by(Node = is.na(Node), Node_bin, Node_cat) %>%
  dplyr::summarise(N = n())
#   Node  Node_bin Node_cat     N
# 1 FALSE neg      N0          86
# 2 FALSE neg      NA          22
# 3 FALSE pos      N1         121
# 4 FALSE pos      N2          17
# 5 FALSE pos      NA          33
# 6 FALSE NA       NA           3
# 7 TRUE  neg      N0         537
# 8 TRUE  neg      NA         260
# 9 TRUE  pos      N1         664
# 10 TRUE  pos      N2         215
# 11 TRUE  pos      N3         142
# 12 TRUE  pos      NA         278
# 13 TRUE  NA       NA        1358


# Verdict: clean gse205568 !!!!  In house version has wrong categorization !!!


# Response
# ========
clin_mta %>%
  dplyr::group_by(Response = !is.na(Response),
                  Response_path = !is.na(Response_pathological),
                  Response_clin = !is.na(Response_clinical)) %>%
  dplyr::summarise(N = n())
# Response Response_path Response_clin     N
# 1 FALSE    FALSE         FALSE           848
# 2 TRUE     FALSE         TRUE            242
# 3 TRUE     TRUE          FALSE          2646

clin_gse205568 %>%
  dplyr::group_by(Response = !is.na(Response),
                  Response_path = !is.na(Response_pathological),
                  Response_clin = !is.na(Response_clinical)) %>%
  dplyr::summarise(N = n())
#   Response Response_path Response_clin     N
# 1 FALSE    FALSE         FALSE           848
# 2 TRUE     FALSE         TRUE            242
# 3 TRUE     TRUE          FALSE          2646

clin_mta$Response %>% table()
# npCR  pCR
# 2144  744
clin_gse205568$Response %>% table()
# npCR  pCR
# 2144  744


# Verdict: clean !!!!


# Follow up
# ========

clin_mta %>%
  dplyr::group_by(Event_dfs = !is.na(Event_dfs),
                  Time_dfs = !is.na(Time_dfs)) %>%
  dplyr::summarise(N = n())
#   Event_dfs Time_dfs     N
# 1 FALSE     FALSE     2344
# 2 FALSE     TRUE         7
# 3 TRUE      FALSE      108
# 4 TRUE      TRUE      1277

clin_gse205568 %>%
  dplyr::group_by(Event_dfs = !is.na(Event_dfs),
                  Time_dfs = !is.na(Time_dfs)) %>%
  dplyr::summarise(N = n())
#   Event_dfs Time_dfs     N
# 1 FALSE     FALSE     2344
# 2 FALSE     TRUE         7
# 3 TRUE      FALSE      108
# 4 TRUE      TRUE      1277


# Verdict: clean !!!!


# Arm
# ========

clin_mta %>%
  dplyr::group_by(Arm_chemo = !is.na(Arm_chemo),
                  Arm_her2 = !is.na(Arm_her2),
                  Arm_hormone = !is.na(Arm_hormone),
                  Arm_other = !is.na(Arm_other)) %>%
  dplyr::summarise(N = n())
#   Arm_chemo Arm_her2 Arm_hormone Arm_other     N
# 1 FALSE     FALSE    FALSE       FALSE       448
# 2 TRUE      TRUE     TRUE        TRUE       3288

clin_gse205568 %>%
  dplyr::group_by(Arm_chemo = !is.na(Arm_chemo),
                  Arm_her2 = !is.na(Arm_her2),
                  Arm_hormone = !is.na(Arm_hormone),
                  Arm_other = !is.na(Arm_other)) %>%
  dplyr::summarise(N = n())
#   Arm_chemo Arm_her2 Arm_hormone Arm_other     N
# 1 FALSE     FALSE    FALSE       FALSE       448
# 2 TRUE      TRUE     TRUE        TRUE       3288


# Verdict: clean !!!!


# Subtype ihc
# ===========

clin_mta %>%
  dplyr::group_by(Er = !is.na(Er),
                  Pr = !is.na(Pr),
                  Hr = !is.na(Hr),
                  Her2 = !is.na(Her2),
                  Subtype_ihc = !is.na(Subtype_ihc)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Er    Pr    Hr    Her2  Subtype_ihc     N
# 1 FALSE FALSE FALSE FALSE FALSE         444
# 2 FALSE FALSE FALSE TRUE  FALSE          95
# 3 TRUE  FALSE TRUE  FALSE FALSE          68
# 4 TRUE  FALSE TRUE  TRUE  TRUE          257
# 5 TRUE  TRUE  TRUE  FALSE FALSE         291
# 6 TRUE  TRUE  TRUE  TRUE  TRUE         2581

clin_gse205568 %>%
  dplyr::group_by(Er = !is.na(Er),
                  Pr = !is.na(Pr),
                  Hr = !is.na(Hr),
                  Her2 = !is.na(Her2),
                  Subtype_ihc = !is.na(Subtype_ihc)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Er    Pr    Hr    Her2  Subtype_ihc     N
# 1 FALSE FALSE FALSE FALSE FALSE         444
# 2 FALSE FALSE FALSE TRUE  FALSE          95
# 3 TRUE  FALSE TRUE  FALSE FALSE          68
# 4 TRUE  FALSE TRUE  TRUE  TRUE          257
# 5 TRUE  TRUE  TRUE  FALSE FALSE         291
# 6 TRUE  TRUE  TRUE  TRUE  TRUE         2581

# Verdict: clean !!!!


# Subtype pam50
# =============

clin_mta$Subtype_pam50 %>% table()
# Basal         Her2         LumA         LumB       Normal Unclassified
# 1047          636         1066          566          401           20
clin_gse205568$Subtype_pam50 %>% table()
# Basal         Her2         LumA         LumB       Normal Unclassified
# 1047          636         1066          566          401           20

# Verdict: clean !!!!

