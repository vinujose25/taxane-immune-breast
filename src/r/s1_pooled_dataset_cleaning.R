# s1_pooled_dataset_cleaning.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Filter pooled dataset based on the criteria specified in the manuscript.
# The pooled dataset is filtered using the following criteria;
# a) Discard samples with missing IHC status for HR and HER2
#
# For each IHC subtype (HR,HER2,TN):
# b) Discard samples with missing treatment data
# c) Discard samples with missing clinical endpoint (pCR/DFS) data
# d) Discard matched samples from longitudinal studies (series matrices)
#   (i.e. consider only pre-treatment expression profiles),
# e) Discard treatment regimens with samples only from a single study,
# f) Discard treatment regimens with sample size less than 50, and
# g) Discard studies from a treatment regimen with less than 10 samples.



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and prepare data
# 2. Clinical data filtering (filter expression data later)
# 3. Expression data filtering
# 4. Additional formatting of clinical data to aid in analysis
# 5. Summarize neoadj/adj dataset
# 6. Save filtered clinical and expression data



# 1. Load and prepare data
# ==============================================================================

# load("results/data/geo_clin_curated_locally.RData")
# clin <- geo_clin_curated_locally
# rm(geo_clin_curated_locally)

clin <- read_tsv(file = "data/gse205568/GSE205568_quantile_corrected/quantile_corrected/geo_clin_Quantile.tsv", guess_max = 3700)
# guess_max =  3736 ?
# Due to large missingness in some variables, large no.of rows needs to read in
# to correctly guess the variable data type.

expr <- read_tsv(file = "data/gse205568/GSE205568_quantile_corrected/quantile_corrected/geo_expr_Quantile.tsv")

dim(clin)
# [1] 3736   93
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



# Review age, grade, node, size, subtype, clinical endpoint data, treatment arm
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



clin %>%
  dplyr::select("Age", "Age_bin", "Age_detailed", "Grade",
                "Node", "Node_bin", "Node_cat",
                "Size", "Size_bin", "Size_cat",
                "Hr", "Er", "Pr", "Her2", "Subtype_ihc",
                "Subtype_pam50", "Subtype_pam50_2") %>%
  glimpse()
# Rows: 3,736
# Columns: 17
# $ Age             <chr> "37.8", "45.8", "40.7", "40.8", "35.5", … !!!!! character
# $ Age_bin         <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
# $ Age_detailed    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
# $ Grade           <chr> "G2", "G2", "G3", NA, "G2", "G3", NA, NA…
# $ Node            <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, … !!!!! character
# $ Node_bin        <chr> "pos", "pos", "neg", "pos", "pos", "neg"…
# $ Node_cat        <chr> "N1", "N1", "N0", "N1", "N1", "N0", "N0"…
# $ Size            <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, … !!!!! character
# $ Size_bin        <chr> "large", "large", "large", "large", "lar…
# $ Size_cat        <chr> "T2", "T3", "T3", "T2", "T3", "T4", "T2"…
# $ Hr              <chr> "pos", "pos", "neg", "neg", "pos", "neg"…
# $ Er              <chr> "pos", "pos", "neg", "neg", "pos", "neg"…
# $ Pr              <chr> "pos", "pos", "neg", "neg", "pos", "neg"…
# $ Her2            <chr> "neg", "neg", "neg", "neg", "neg", "neg"…
# $ Subtype_ihc     <chr> "HR", "HR", "TN", "TN", "HR", "TN", NA, …
# $ Subtype_pam50   <chr> "LumA", "Normal", "Basal", "Basal", "Lum…
# $ Subtype_pam50_2 <chr> "HR", NA, "TN", "TN", "HR", "TN", "TN",



# Age >>>>>>>>>>>
#


# Age is character !!!! Explore further
x <- clin %>%
  dplyr::mutate(Age2 = as.integer(Age)) %>%
  dplyr::filter((!is.na(Age))&is.na(Age2))

x$Age %>% table()
# -- >40 ≤40
# 1  27  17
# Age contains these characters: --, >40, ≤40. Remove them !!!


# Explore the characters in Age further
x %>%
  dplyr::group_by(Age = !is.na(Age), Age2 = !is.na(Age2), Age_bin, Age_detailed) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Age   Age2  Age_bin Age_detailed                 N
#   <lgl> <lgl> <chr>   <chr>                    <int>
# 1 TRUE  FALSE old     Young:age<=40;Old:age>40    27
# 2 TRUE  FALSE young   Young:age<=40;Old:age>40    17
# 3 TRUE  FALSE NA      NA                           1


# Summarize Age
clin %>%
  dplyr::group_by(Age = !is.na(Age), Age_bin, Age_detailed) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Age   Age_bin Age_detailed                 N
# 1 FALSE old     young:<=50,old>50           50
# 2 FALSE young   young:<=50,old>50           70
# 3 FALSE NA      NA                         791
# 4 TRUE  old     Young:age<=40;Old:age>40    27
# 5 TRUE  young   Young:age<=40;Old:age>40    17
# 6 TRUE  NA      NA                        2781


# !!!! Set Age as as.integer() will remove categories from continuous variable
# Will introduce warnings !!!!
# Use 50 as age cutoff and introduce a new variable Age_bin_50
# to maximize the sample size with Age
# (Note: The 120 samples with Age_bin has different age cutoffs)

# Update clin
clin <- clin %>%
  dplyr::mutate(
    # Remove characters (--, >40, ≤40) from Age
    # Age = as.integer(Age),
    Age = purrr::map_chr(
      Age,
      ~switch(
        .x,
        "--" = NA_character_,
        ">40" = NA_character_,
        "≤40" = NA_character_,
        .x
      )
    ), # to suppress  warning NAs introduced by coercion
    Age = as.integer(Age),
    Age_bin_50 = if_else(Age <= 50, "young", "old"),
    Age_bin_50 = if_else(is.na(Age) & str_detect(Age_detailed,"50"),
                         Age_bin, Age_bin_50)
  )


# Verification
clin %>%
  dplyr::group_by(Age = !is.na(Age), Age_bin_50, Age_bin, Age_detailed) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Age   Age_bin_50 Age_bin Age_detailed                 N
# 1 FALSE old        old     young:<=50,old>50           50
# 2 FALSE young      young   young:<=50,old>50           70
# 3 FALSE NA         old     Young:age<=40;Old:age>40    27
# 4 FALSE NA         young   Young:age<=40;Old:age>40    17
# 5 FALSE NA         NA      NA                         792 # one more NA due to Age = "-"
# 6 TRUE  old        NA      NA                        1309
# 7 TRUE  young      NA      NA                        1471


clin$Age_bin_50 %>% table()
# old young
# 1359 1541



# Grade >>>>>>>>>>>>>
#


clin %>%
  dplyr::group_by(Grade) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Grade     N
#   <chr> <int>
# 1 G1      161
# 2 G1/2     23
# 3 G2      869
# 4 G3     1231
# 5 NA     1452



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


# Verification
clin %>%
  dplyr::group_by(Grade, Grade_bin) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Grade Grade_bin     N
# 1 G1    G12         161
# 2 G1/2  G12          23
# 3 G2    G12         869
# 4 G3    G3         1231
# 5 NA    NA         1452

clin$Grade_bin %>%  table()
# G12   G3
# 1053 1231




# Node >>>>>>>>>>>>>>>>>
#


# Node is character !!!! Explore further
x <- clin %>%
  dplyr::mutate(Node2 = as.integer(Node)) %>%
  dplyr::filter((!is.na(Node))&is.na(Node2))

x$Node %>% table()
# N-  N0  N1  N2
# 3  86 121  17
# Node contains these characters: N-  N0  N1  N2. Remove them or merge to Node_cat !!!


# Explore the characters in Node further
x %>%
  dplyr::group_by(Node = !is.na(Node), Node2 = !is.na(Node2), Node_bin, Node_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Node  Node2 Node_bin Node_cat     N
# 1 TRUE  FALSE neg      N0          86
# 2 TRUE  FALSE pos      N1         121
# 3 TRUE  FALSE pos      N2          17
# 4 TRUE  FALSE NA       NA           3 # N-


# Summarize Node
clin %>%
  dplyr::group_by(Node = !is.na(Node), Node_bin, Node_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Node  Node_bin Node_cat     N
# 1 FALSE neg      N0         537
# 2 FALSE neg      NA         260
# 3 FALSE pos      N1         664
# 4 FALSE pos      N2         215
# 5 FALSE pos      N3         142
# 6 FALSE pos      NA         278
# 7 FALSE NA       NA        1358
# 8 TRUE  neg      N0          86
# 9 TRUE  neg      NA          22
# 10 TRUE  pos      N1         121
# 11 TRUE  pos      N2          17
# 12 TRUE  pos      NA          33
# 13 TRUE  NA       NA           3


# Explore Node
clin %>%
  dplyr::filter(!is.na(Node)) %>% dplyr::select(Node, Node_bin, Node_cat) %>%
  as.data.frame()

# Node is intended to be continuous with no. of positive node,
# but contains ONLY redundant Node_bin or Node_cat data.
# Set Node = NA
clin <- clin %>%
  dplyr::mutate(Node = NA)


# Verification !!!!!!!!!!!!!
clin %>%
  dplyr::group_by(Node = !is.na(Node), Node_bin, Node_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Node  Node_bin Node_cat     N
# 1 FALSE neg      N0         623
# 2 FALSE neg      NA         282
# 3 FALSE pos      N1         785
# 4 FALSE pos      N2         232
# 5 FALSE pos      N3         142
# 6 FALSE pos      NA         311
# 7 FALSE NA       NA        1361


table(clin$Node_bin)
# neg  pos
# 905 1470



# Size >>>>>>>>>>>>>
#


# Size is character !!!! Explore further
x <- clin %>%
  dplyr::mutate(Size2 = as.integer(Size)) %>%
  dplyr::filter((!is.na(Size))&is.na(Size2))

x$Size %>% table()
# T-  T1  T2  T3
# 3   5 139  80
# Size contains these characters: T-  T1  T2  T3 . Remove them or merge to Size_cat !!!


# Explore the characters in Size further
x %>%
  dplyr::group_by(Size = !is.na(Size), Size2 = !is.na(Size2), Size_bin, Size_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Size  Size2 Size_bin Size_cat     N
# 1 TRUE  FALSE large    T2         139
# 2 TRUE  FALSE large    T3          80
# 3 TRUE  FALSE small    T1           5
# 4 TRUE  FALSE NA       NA           3 # T-


# Summarize Size
clin %>%
  dplyr::group_by(Size = !is.na(Size), Size_bin, Size_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())
#   Size  Size_bin Size_cat     N
# 1 FALSE large    T2        1095
# 2 FALSE large    T3         453
# 3 FALSE large    T4         238
# 4 FALSE samall   T0           3
# 5 FALSE samall   T1         131
# 6 FALSE small    T0           6
# 7 FALSE small    T1          62
# 8 FALSE NA       NA         745
# 9 TRUE  large    T2         461
# 10 TRUE  large    T3         251
# 11 TRUE  large    NA          35
# 12 TRUE  small    T1         233
# 13 TRUE  small    NA          20
# 14 TRUE  NA       NA           3 # T- ~ NA


# Explore Node
x <- clin %>%
  dplyr::filter(!is.na(Size)) %>% dplyr::select(Size, Size_bin, Size_cat) %>%
  as.data.frame()
dim(x) # 1003 3

# Size is intended to be continuous with no. of positive node,
# but contains redundant Size_bin or Size_cat data.
# Remove redundant values with "T
clin <- clin %>%
  dplyr::mutate(
    Size = if_else(str_detect(Size, "T"), NA_character_, Size),
    Size_bin = str_replace(Size_bin, "samall", "small")
  )

clin %>%
  dplyr::filter(!is.na(Size)) %>% dplyr::select(Size, Size_bin, Size_cat) %>%
  dim() # 776 3



# Summarize Size
clin %>%
  dplyr::group_by(Size = !is.na(Size), Size_bin, Size_cat) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n())

#   Size  Size_bin Size_cat     N
# 1 FALSE large    T2        1234
# 2 FALSE large    T3         533
# 3 FALSE large    T4         238
# 4 FALSE small    T0           9
# 5 FALSE small    T1         198
# 6 FALSE NA       NA         748
# 7 TRUE  large    T2         322
# 8 TRUE  large    T3         171
# 9 TRUE  large    NA          35
# 10 TRUE  small    T1         228
# 11 TRUE  small    NA          20


table(clin$Size_bin)
# large small
# 2533   455


# Subtype >>>>>>>>>>>
#


clin %>%
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

# Subtype_ihc is defined for samples with both Hr and Her2 status



# Clinical endpoint data >>>>>>>>>>>>>>>>>>>>>
#



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




# Treatment arm >>>>>>>>>>>>>>>>>>>>>>>>>
#


# Consolidate Arm for easy filtering;
# Arm_chemo + Arm_her2 + Arm_hr + Arm_other


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




# 3. Expression data filtering
# ==============================================================================

expr_neoadj <- expr %>%
      dplyr::select(1, all_of(clin_neoadj$Sample_geo_accession))


#
# ==============================================================================



# 4. Additional formatting of clinical data to aid in analysis
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

    # # Added on 3/4/21
    # # Integrate HR status in the the treatment strata for HER2
    # # Objective:
    # # 1) To asses HR status based difference in response
    # # 2) Usually in clinics HR+ve samples were given hormone therapy
    # #     irrespective of HER2 status.
    # Strata = if_else(Subtype_ihc == "HER2", str_c(Strata,":",Hr), Strata) %>%
    #   str_replace("pos","HR+") %>%
    #   str_replace("neg", "HR-")

    # Added on 5/4/21
    # Reclassify HER2 subtype as HR-HER2+ and HR+HER2+
    # Objective:
    # 1) To asses HR status based difference in response
    # 2) Usually in clinics HR+ve samples were given hormone therapy
    #     irrespective of HER2 status.
    Subtype_ihc2 = if_else(Subtype_ihc == "HER2",
                           str_c(Hr, Subtype_ihc, "+"),
                           Subtype_ihc) %>%
      str_replace("pos","HR+") %>%
      str_replace("neg", "HR-")
  )

clin_neoadj %>%
  dplyr::group_by(Subtype_ihc, Subtype_ihc2, Strata) %>%
  dplyr::summarise(N = n(), Pcr = Response %>% sum()) %>%
  as.data.frame()
#    Subtype_ihc Subtype_ihc2             Strata   N Pcr
# 1         HER2     HR-HER2+     GSE20194:AAA+T  25  13
# 2         HER2     HR-HER2+     GSE32646:AAA+T  18   9
# 3         HER2     HR-HER2+ GSE42822:AAA+T+TRA  15  10
# 4         HER2     HR-HER2+     GSE50948:AAA+T  37   9
# 5         HER2     HR-HER2+ GSE50948:AAA+T+TRA  43  25
# 6         HER2     HR+HER2+     GSE20194:AAA+T  25   4
# 7         HER2     HR+HER2+     GSE32646:AAA+T  16   3
# 8         HER2     HR+HER2+ GSE42822:AAA+T+TRA   9   2
# 9         HER2     HR+HER2+     GSE50948:AAA+T  14   4
# 10        HER2     HR+HER2+ GSE50948:AAA+T+TRA  20   6
# 11          HR           HR     GSE20194:AAA+T 141   9
# 12          HR           HR       GSE20271:AAA  49   3
# 13          HR           HR     GSE20271:AAA+T  44   3
# 14          HR           HR       GSE22093:AAA  42  10
# 15          HR           HR     GSE23988:AAA+T  32   7
# 16          HR           HR     GSE25066:A0A+T  35   3
# 17          HR           HR     GSE25066:AAA+T 194  22
# 18          HR           HR     GSE32646:AAA+T  55   5
# 19          HR           HR     GSE41998:A0A+T 104  12
# 20          HR           HR     GSE42822:AAA+T  29   8
# 21          HR           HR     GSE50948:AAA+T  26   4
# 22          TN           TN     GSE20194:AAA+T  64  25
# 23          TN           TN       GSE20271:AAA  28   4
# 24          TN           TN     GSE20271:AAA+T  31   9
# 25          TN           TN       GSE22093:AAA  55  18
# 26          TN           TN     GSE23988:AAA+T  29  13
# 27          TN           TN     GSE25066:A0A+T  25   6
# 28          TN           TN     GSE25066:AAA+T 119  43
# 29          TN           TN     GSE32646:AAA+T  26  10
# 30          TN           TN     GSE41998:A0A+T 125  51
# 31          TN           TN     GSE42822:AAA+T  24  12
#
# ==============================================================================



# 5. Summarize neoadj/adj dataset
# ==============================================================================

# Note !!!
# No adj dataset cleared filtering criteria
# Summary of neoadj dataset is given below



# Clinical-vars, subtype, treatment and response data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_neoadj %>%
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


# 996 samples got all relevant clinical-pathological variables
# Total: 1499;
# 1499 - 996 = 503 have at least one relevant clinical-pathological variable missing


clin_neoadj %>%
  dplyr::filter(!is.na(Age_bin_50) &
                  !is.na(Grade_bin) &
                  !is.na(Node_bin) &
                  !is.na(Size_bin) &
                  !is.na(Subtype_ihc) &
                  !is.na(Arm_consolidated) &
                  !is.na(Response)) %>%
  dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
  dplyr::summarise(N = n(), .groups = "keep") %>%
  dplyr::filter(N >= 50)
#   Arm_consolidated                                                     Subtype_ihc     N
# 1 AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy HR             58
# 2 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HER2           82
# 3 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   HR            486
# 4 AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy   TN            269

# If age, grade, node, size, subtype, arm, and pCR needs to be used,
# AAA+/-Taxane * immune interaction is possible in HR subtype



# Subtype, treatment and response data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_neoadj %>%
  dplyr::group_by(Subtype = !is.na(Subtype_ihc),
                  Arm = !is.na(Arm_consolidated),
                  pCR = !is.na(Response)) %>%
  # TRUE = data, FALSE = no data, NA = no data
  dplyr::summarise(N = n(), .groups = "keep")
#   Subtype Arm   pCR       N
# 1 TRUE    TRUE  TRUE   1499

# !!!!!!!!!!!!!!!!!!!
# 1499 samples got all subtype, arm and pCR data

clin_neoadj %>%
  dplyr::filter(!is.na(Subtype_ihc) &
                  !is.na(Arm_consolidated) &
                  !is.na(Response)) %>%
  dplyr::group_by(Arm_consolidated, Subtype_ihc) %>%
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


# If subtype, arm and pCR needs to be used,
# AAA+/-Taxane * immune interaction is possible in HR and TN subtype


# !!! Limit analysis to AAA+/-Taxane containing regime to make the treatment regimen
# congruent to FinHER FEC+MTA regimen

#
#===============================================================================




# 6. Save filtered clinical and expression data
# ==============================================================================

save(clin_neoadj, file = str_c(out_data,"clin_neoadj.RData"))
save(expr_neoadj, file = str_c(out_data,"expr_neoadj.RData"))

#
# ==============================================================================


# Clean memory
# ==============================================================================

rm(clin, clin_adj, clin_neoadj,
   expr,  expr_neoadj)

# ==============================================================================
