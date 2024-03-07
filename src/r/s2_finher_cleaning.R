# s2_finher_cleaning.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Clean finher clinical and expression data



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and prepare data
# 2. Clinical data cleaning
# 3. Summarize clinical dataset
# 4. Save filtered clinical and expression data


# 1. Load and prepare data
# ==============================================================================

clin <- read_tsv(file = "data/finher_clin.tsv") # Private clinical data
metadata <- read_tsv(file = "data/gse47994/gse47994_metadata.tsv") # inhouse - geo-id mapping
expr <- read_tsv(file = "data/gse47994/GSE47994_reprocessed_sample_filtered/mas/expr_maxvar_combat_after_sample_and_pset_filtering.tsv")
annot <- read_tsv(file = "data/gse47994/GSE47994_reprocessed_sample_filtered/mas/annot_maxvar_after_sample_and_pset_filtering.tsv")



clin <- dplyr::left_join(
  clin,
  metadata %>%
    dplyr::select(Sample_title, Sample_geo_accession, Filename) %>%
    dplyr::mutate(CEL_filename = str_c(Filename, ".CEL.gz")),
  by = "CEL_filename")


# expr ordering sample names
expr <- expr %>%
  dplyr::rename(Probeset_Id = "ID_REF") %>%
  dplyr::select(Probeset_Id, all_of(clin$Sample_geo_accession)) %>%
  dplyr::rename_with(
    ~purrr::map_chr(.x,
                function(.x,map){
                  if(.x %in% map$Sample_geo_accession){
                    return(map$CEL_filename[which(.x == map$Sample_geo_accession)])
                  }else(
                    return(.x)
                  )
                }, map = clin)
  )


identical(names(expr)[-1], clin$CEL_filename) # TRUE

# expr ordering gene names
expr <- left_join(expr %>% dplyr::select(Probeset_Id), annot, by = "Probeset_Id") %>%
  left_join(expr, by = "Probeset_Id") %>%
  dplyr::select(-c("Probeset_Id", "Gene_Symbol")) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))

dim(clin)
# [1] 300 80
dim(expr)
# [1] 3350  301
dim(annot)
# [1] 3350  3

#
# ==============================================================================



# 2. Clinical data cleaning
# ==============================================================================

# TIL
# age, grade, node, size, er, pr, hr, her2, subtype_ihc,
# rfs, ddfs, os
# arm(chemo, trastuzumab, hormone)

glimpse(clin)

# InTuLy_Mean, StrLy_Mean, TotalLy_Mean
clin %>%
  dplyr::select(InTuLy_Mean, StrLy_Mean, TotalLy_Mean) %>%
  summary()
# InTuLy_Mean      StrLy_Mean     TotalLy_Mean
# Min.   : 0.00   Min.   : 0.00   Min.   :  0.00
# 1st Qu.: 0.50   1st Qu.:10.00   1st Qu.: 10.50
# Median : 2.50   Median :17.50   Median : 21.00
# Mean   : 4.45   Mean   :23.95   Mean   : 28.40
# 3rd Qu.: 5.50   3rd Qu.:37.50   3rd Qu.: 42.62
# Max.   :70.00   Max.   :85.00   Max.   :155.00
# NA's   :12      NA's   :12      NA's   :12


# Age_Years, Grade, Nodal_Status, Axill_Nodes_Investigated, Axill_Nodes_With_Cancer
# Tumor_Size_mm
clin %>%
  dplyr::select(Age_Years, Axill_Nodes_Investigated,
                Axill_Nodes_With_Cancer,Tumor_Size_mm) %>%
  summary()
#   Age_Years     Axill_Nodes_Investigated Axill_Nodes_With_Cancer Tumor_Size_mm
# Min.   :25.45   Min.   : 1.00            Min.   : 0.000          Min.   :  6.00
# 1st Qu.:43.60   1st Qu.: 9.00            1st Qu.: 1.000          1st Qu.: 19.00
# Median :51.28   Median :12.00            Median : 2.000          Median : 24.00
# Mean   :49.94   Mean   :13.04            Mean   : 3.257          Mean   : 27.95
# 3rd Qu.:56.42   3rd Qu.:16.00            3rd Qu.: 4.000          3rd Qu.: 30.00
# Max.   :65.81   Max.   :52.00            Max.   :29.000          Max.   :160.00
#                                                                  NA's   :1
clin %>%
  dplyr::select(Grade, Nodal_Status) %>%
  purrr::map(~(table(.x)))
# $Grade
# Grade1 Grade2 Grade3
# 5     79    208
# $Nodal_Status
# Negative Positive
# 64      236
# Suggested change !!!!!!!!!!!!!!!!!!!!!
# Map Grade1 to Grade2 !!!!!!!!!!!!!!!!!


# ER_IHC, PR_IHC, HR_IHC, HER2_IHC_CISH, Subtype_IHC, Subtype_IHC_2
clin %>%
  dplyr::group_by(ER_IHC, PR_IHC, HR_IHC, HER2_IHC_CISH, Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarise(N = n())
#   ER_IHC   PR_IHC   HR_IHC   HER2_IHC_CISH Subtype_IHC Subtype_IHC_2     N
# 1 Negative Negative Negative Negative      TN          TN              120
# 2 Negative Negative Negative Positive      HR-HER2+    HER2             89
# 3 Negative Positive Positive Positive      HR+HER2+    HER2              4
# 4 Positive Negative Positive Positive      HR+HER2+    HER2             27
# 5 Positive Positive Positive Positive      HR+HER2+    HER2             60
# Total 300, No NAs !!!!

# RFS_Event, DDFS_Event, OS_Event
clin %>%
  dplyr::group_by(RFS_Event, DDFS_Event, OS_Event) %>%
  dplyr::summarise(N = n())
#   RFS_Event DDFS_Event OS_Event     N
# 1         0          0        0   220
# 2         1          0        0     8
# 3         1          1        0    24
# 4         1          1        1    48
# Total 300, No NAs !!!!


# RFS_Time_Years, DDFS_Time_Years, OS_Time_Years
clin %>%
  dplyr::select(RFS_Time_Years, DDFS_Time_Years, OS_Time_Years) %>%
  summary()
# RFS_Time_Years    DDFS_Time_Years  OS_Time_Years
# Min.   :0.00548   Min.   :0.1287   Min.   :0.1287
# 1st Qu.:4.01711   1st Qu.:4.0931   1st Qu.:4.3799
# Median :4.89802   Median :4.9090   Median :5.1006
# Mean   :4.48520   Mean   :4.5925   Mean   :4.9150
# 3rd Qu.:5.87885   3rd Qu.:5.9055   3rd Qu.:5.9377
# Max.   :6.76523   Max.   :6.7652   Max.   :6.7652


# Chemo, Herceptin, HR_IHC
clin %>%
  dplyr::group_by(Chemo, Herceptin, HR_IHC) %>%
  dplyr::summarise(N = n())
#   Chemo Herceptin HR_IHC       N
# 1 DTX   No        Negative    27
# 2 DTX   No        Positive    18
# 3 DTX   Yes       Negative    26
# 4 DTX   Yes       Positive    16
# 5 DTX   NA        Negative    60 # TN
# 6 NVB   No        Negative    19
# 7 NVB   No        Positive    27
# 8 NVB   Yes       Negative    17
# 9 NVB   Yes       Positive    30
# 10 NVB   NA        Negative    60 # TN
# Suggested change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Chemo: keep intact
# Herceptin: No/NA as noTRA, NA are set to noTRA for later grouping of samples
# introduce Hormone: as HR for Positive and noHR for Negative
# Obsolete: Integrate Chemo Herceptin HR_IHC into a single variable Arm !!!!!!!!!!!!!!!!!


# Implement suggested changes to clin
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clin <- clin %>%
  dplyr::mutate(
    Grade = if_else(Grade == "Grade1", "Grade2", Grade),
    Herceptin = if_else( (is.na(Herceptin) | Herceptin == "No"), "noTRA", "TRA"),
    Hormone = if_else( HR_IHC == "Negative", "noHR", "HR"),
  )

clin$Grade %>% table()
# Grade2 Grade3
# 84    208

clin %>%
  dplyr::group_by(Chemo, Herceptin, Hormone) %>%
  dplyr::summarise(N = n())
#   Chemo Herceptin Hormone     N
# 1 DTX   noTRA     HR         18
# 2 DTX   noTRA     noHR       87 # TN + HER2
# 3 DTX   TRA       HR         16
# 4 DTX   TRA       noHR       26
# 5 NVB   noTRA     HR         27
# 6 NVB   noTRA     noHR       79 # TN + HER2
# 7 NVB   TRA       HR         30
# 8 NVB   TRA       noHR       17

#
# ==============================================================================



# 3. Summarize clinical dataset
# ==============================================================================

clin %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2, Arm = str_c(Chemo,"+",Herceptin,"+",Hormone)) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC Subtype_IHC_2 Arm                N
# 1 HR-HER2+    HER2          DTX+noTRA+noHR    27
# 2 HR-HER2+    HER2          DTX+TRA+noHR      26
# 3 HR-HER2+    HER2          NVB+noTRA+noHR    19
# 4 HR-HER2+    HER2          NVB+TRA+noHR      17

# 5 HR+HER2+    HER2          DTX+noTRA+HR      18
# 6 HR+HER2+    HER2          DTX+TRA+HR        16
# 7 HR+HER2+    HER2          NVB+noTRA+HR      27
# 8 HR+HER2+    HER2          NVB+TRA+HR        30

# 9 TN          TN            DTX+noTRA+noHR    60
# 10 TN          TN            NVB+noTRA+noHR    60


clin %>%
  dplyr::group_by(Subtype_IHC_2, Arm = str_c(Chemo)) %>%
  dplyr::summarise(N = n())

#   Subtype_IHC_2 Arm       N
# 1 HER2          DTX      87
# 2 HER2          NVB      93
# 3 TN            DTX      60
# 4 TN            NVB      60


#
# ==============================================================================



# 4. Save filtered clinical and expression data
# ==============================================================================

clin_finher <- clin
expr_finher <- expr

save(clin_finher, file = str_c(out_data,"clin_finher.RData"))
save(expr_finher, file = str_c(out_data,"expr_finher.RData"))



#
# ==============================================================================



# Clear memory
# ==============================================================================

rm(clin, clin_finher,
   expr, expr_finher,
   metadata, annot)

#
# ==============================================================================
