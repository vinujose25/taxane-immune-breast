# s4_exploring_interaction.R



# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the interaction between
#   1) Arm, Signature, and Survival within each Subtype, and
#   2) Subtype, Signature and Survival within each Arm.
# Prepare data necessary to plot interaction plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data
# 2. Structure analysis and define strata and interacting variables
# 3. MTA interaction
# 4. Trastuzumab interaction
# 5. Subtype interaction between TN and HR-HER2+ (no-hormone + noTRA) interaction
# 6. Save Robjects of interaction summaries




# 1. Load data
# ==============================================================================

load("results/data/clin_finher.RData")
dim(clin_finher)
# 300 118


# clean Arm and factor levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
clin_finher$Arm %>% factor() %>% levels()
# [1] "FEC_DTX"             "FEC_DTX_Hormone"     "FEC_DTX_TRA"         "FEC_DTX_TRA_Hormone"
# [5] "FEC_NVB"             "FEC_NVB_Hormone"     "FEC_NVB_TRA"         "FEC_NVB_TRA_Hormone"

clin_finher <- clin_finher %>%
  dplyr::mutate(
    Arm =  Arm %>%
      str_replace_all("_", "+")%>%
      factor(levels = c(
        "FEC+DTX", "FEC+DTX+Hormone", "FEC+DTX+TRA", "FEC+DTX+TRA+Hormone",
        "FEC+NVB", "FEC+NVB+Hormone", "FEC+NVB+TRA", "FEC+NVB+TRA+Hormone"
      ))
  )

clin_finher$Arm %>% levels()
# [1] "FEC+DTX"             "FEC+DTX+Hormone"     "FEC+DTX+TRA"         "FEC+DTX+TRA+Hormone"
# [5] "FEC+NVB"             "FEC+NVB+Hormone"     "FEC+NVB+TRA"         "FEC+NVB+TRA+Hormone"



# set subtype ihc levels
# >>>>>>>>>>>>>>>>>>>>>>
clin_finher$Subtype_IHC %>% factor() %>% levels()
# [1] "HR-HER2+" "HR+HER2+" "TN"
clin_finher$Subtype_IHC_2 %>% factor() %>% levels()
# [1] "HER2" "TN"

clin_finher <- clin_finher %>%
  dplyr::mutate(
    Subtype_IHC =  Subtype_IHC %>%
      factor(levels = c("TN", "HR-HER2+", "HR+HER2+")),
    Subtype_IHC_2 =  Subtype_IHC_2 %>%
      factor(levels = c("TN", "HER2"))
  )

clin_finher$Subtype_IHC %>% levels()
# [1] "TN" "HR-HER2+" "HR+HER2+"
clin_finher$Subtype_IHC_2 %>% levels()
# [1] "TN" "HER2"


#
# ==============================================================================




# 2. Structure analysis and define strata and interacting variables
# ==============================================================================


# Summarize clin_finher
# >>>>>>>>>>>>>>>>>>>>>

# All adj regimen

clin_finher %>%
  dplyr::group_by(
    Subtype = Subtype_IHC, # Subtype_IHC_2
    Arm = Arm
  ) %>%
  dplyr::summarise(N = n(),
                   RFS = which(RFS_Event == 1) %>% length(),
                   DDFS = which(DDFS_Event == 1) %>% length(),
                   OS = which(OS_Event == 1) %>% length()) %>%
  dplyr::mutate(
    RFS_perc = (RFS/N)*100,
    DDFS_perc = (DDFS/N)*100,
    OS_perc = (OS/N)*100
  ) %>%
  as.data.frame()

#     Subtype                 Arm  N RFS DDFS OS RFS_perc DDFS_perc   OS_perc
# 1        TN             FEC+DTX 60  14   13  9 23.33333  21.66667 15.000000
# 2        TN             FEC+NVB 60  21   19 14 35.00000  31.66667 23.333333

# 3  HR-HER2+             FEC+DTX 27   7    7  5 25.92593  25.92593 18.518519
# 4  HR-HER2+         FEC+DTX+TRA 26   3    3  2 11.53846  11.53846  7.692308
# 5  HR-HER2+             FEC+NVB 19   6    6  6 31.57895  31.57895 31.578947
# 6  HR-HER2+         FEC+NVB+TRA 17   6    5  3 35.29412  29.41176 17.647059

# 7  HR+HER2+     FEC+DTX+Hormone 18   4    4  2 22.22222  22.22222 11.111111
# 8  HR+HER2+ FEC+DTX+TRA+Hormone 16   2    0  0 12.50000   0.00000  0.000000
# 9  HR+HER2+     FEC+NVB+Hormone 27   7    6  3 25.92593  22.22222 11.111111
# 10 HR+HER2+ FEC+NVB+TRA+Hormone 30  10    9  4 33.33333  30.00000 13.333333




# TN analysis
# >>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival


# HR-HER2 analysis
# >>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival + strata(TRA|noTRA)

# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival + strata(DTX|NVB)


# HR+HER2 analysis
# >>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen (! RFS only)
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival + strata(TRA|noTRA)

# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen (! RFS only)
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival + strata(DTX|NVB)


# HER2 analysis (! stratified HR- and  HR+)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival strata(TRA|noTRA) regimen

# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival strata(DTX|NVB) regimen


# Pan-Subtype (TN + HR-HER2+) analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. Subtype interaction in DTX + noTRA regimen: TN/HR-HER2+ * Sig * Survival
# 2. Subtype interaction in NVB + noTRA regimen: TN/HR-HER2+ * Sig * Survival
# 3. Subtype interaction in noTRA regimen: TN/HR-HER2+ * Sig * Survival + strata(DTX|NVB)


# Analysis note: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# 1. Script structure:
#     Create inter_sum (interaction summary) list object.
#     Each elemet represent anlysis done using each endpoint on each subtype.




# Introduce strata variable for
# 1) DTX|NVB interaction (HR_IHC * Herceptin) and
# 2) TRA|noTRA interaction analysis (HR_IHC * Chemo)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_finher <- clin_finher %>%
  dplyr::mutate(
    # HR_IHC * Herceptin stratum for DTX|NVB interaction analysis
    Stratum_hr_tra = str_c(HR_IHC, "_", Herceptin), # NA propogation, Herceptin = NA for TN
    Stratum_hr_tra = if_else(is.na(Stratum_hr_tra), "TN", Stratum_hr_tra),
    # Note that the dummy stratum "TN" is included to simplify coding

    # HR_IHC * Chemo stratum for TRA|noTRA interaction analysis
    Stratum_hr_chemo = str_c(HR_IHC, "_", Chemo)
  )

clin_finher$Stratum_hr_chemo %>% table()
# Negative_DTX Negative_NVB Positive_DTX Positive_NVB
#          113           96           34           57
clin_finher$Stratum_hr_tra %>% table()
# Negative_No Negative_Yes  Positive_No Positive_Yes           TN
#          46           43           45           46          120




# Interacting variables
# >>>>>>>>>>>>>>>>>>>>>

# Chemo (DTX|NVB)
# Herceptin (TRA|noTRA)
# Subtype_IHC_2 (HER2|TN)



#
# ==============================================================================




# 3. MTA interaction
# ==============================================================================



# TN analysis
# >>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival


# HR-HER2 analysis
# >>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival + strata(TRA|noTRA)


# HR+HER2 analysis
# >>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen (! RFS only)
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival + strata(TRA|noTRA)


# HER2 analysis (! stratified HR- and  HR+)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. DTX|NVB interaction: DTX|NVB * Sig * Survival in no TRA regimen
# 2. DTX|NVB interaction: DTX|NVB * Sig * Survival in TRA regimen
# !!! main analysis
# 3. DTX|NVB interaction: DTX|NVB * Sig * Survival strata(TRA|noTRA) regimen


# Define variable
# >>>>>>>>>>>>>>>

inter_sum_finher_mta <- list()


# MTA interaction stratified on HR*TRA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_mta[["all"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

    function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("TN", "HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                     subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Chemo",
              strata_variable = "Stratum_hr_tra",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("TN", "HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher
)

names(inter_sum_finher_mta[["all"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_mta later in the script




# MTA interaction in HER2 + noTRA stratified on HR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_mta[["notra"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Chemo",
              strata_variable = "Stratum_hr_tra",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher %>% dplyr::filter(Herceptin == "No") # for TN Herceptin == NA
)

names(inter_sum_finher_mta[["notra"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_notra_mta later in the script




# MTA interaction in HER2 + TRA stratified on HR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_mta[["tra"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Chemo",
              strata_variable = "Stratum_hr_tra",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher %>% dplyr::filter(Herceptin == "Yes") # for TN Herceptin == NA
)
# Convergence issue warnings !!!!!!!!
# There were 50 or more warnings (use warnings() to see the first 50)
#   "1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#     Loglik converged before variable  2 ; coefficient may be infinite. "
# Warnings are due to no-event strata DTX+TRA+Hormone in HR+HER2+
# Hence coeffcients from the strata DTX+TRA+Hormone are invalid.
#
# Note
# The above script did not look for individual term p values, but likilihood pvalues
# Hence pvalues are unaffected by convergence issue.
# See in ?coxph, the section "Convergence".

names(inter_sum_finher_mta[["tra"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_tra_mta later in the script


#
# ==============================================================================




# 4. Trastuzumab interaction
# ==============================================================================


# HR-HER2 analysis
# >>>>>>>>>>>>>>>>
# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival + strata(DTX|NVB)


# HR+HER2 analysis
# >>>>>>>>>>>>>>>>
# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen (! RFS only)
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival + strata(DTX|NVB)


# HER2 analysis (! stratified HR- and  HR+)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in DTX regimen
# 5. TRA|noTRA interaction: TRA|noTRA * Sig * Survival in NVB regimen
# !!! main analysis
# 6. TRA|noTRA interaction: TRA|noTRA * Sig * Survival strata(DTX|NVB) regimen



# Define variable
# >>>>>>>>>>>>>>>

inter_sum_finher_tra <- list()



# Trastuzumab interaction in HER2 stratified on HR*Chemotherapy
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_tra[["all"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))

        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Herceptin",
              strata_variable = "Stratum_hr_chemo",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher
)

names(inter_sum_finher_tra[["all"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_tra later in the script




# Trastuzumab interaction in HER2 + DTX stratified on HR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_tra[["dtx"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))

        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Herceptin",
              strata_variable = "Stratum_hr_chemo",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher %>% dplyr::filter(Chemo == "DTX")
)
# Convergence issue warnings !!!!!!!!
# There were 50 or more warnings (use warnings() to see the first 50)
#   "1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#     Loglik converged before variable  2 ; coefficient may be infinite. "
# Warnings are due to no-event strata DTX+TRA+Hormone in HR+HER2+
# Hence coeffcients from the strata DTX+TRA+Hormone are invalid.
#
# Note
# The above script did not look for individual term p values, but likilihood pvalues
# Hence pvalues are unaffected by convergence issue.
# See in ?coxph, the section "Convergence".

names(inter_sum_finher_tra[["dtx"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_dtx_tra later in the script




# Trastuzumab interaction in HER2 + NVB stratified on HR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


inter_sum_finher_tra[["nvb"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        # subtype= "TN"
        # event_type = "DDFS"
        # clin = clin_finher
        # sig = "Immune1"


        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))

        # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Herceptin",
              strata_variable = "Stratum_hr_chemo",
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH"),
            Subtype = subtype
          )

      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },
  clin = clin_finher %>% dplyr::filter(Chemo == "NVB")
)

names(inter_sum_finher_tra[["nvb"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_nvb_tra later in the script

#
# ==============================================================================




# 5. Subtype interaction between TN and HR-HER2+ (no-hormone + noTRA) interaction
# ==============================================================================



# Pan-Subtype (TN + HR-HER2+) analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. Subtype interaction in DTX + noTRA regimen: TN/HR-HER2+ * Sig * Survival
# 2. Subtype interaction in NVB + noTRA regimen: TN/HR-HER2+ * Sig * Survival
# 3. Subtype interaction in noTRA regimen: TN/HR-HER2+ * Sig * Survival + strata(DTX|NVB)

clin_finher %>%
  dplyr::group_by(Subtype_IHC, Herceptin) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC Herceptin     N
# 1 TN          NA          120
# 2 HR-HER2+    No           46
# 3 HR-HER2+    Yes          43
# 4 HR+HER2+    No           45
# 5 HR+HER2+    Yes          46

clin_finher %>%
  dplyr::filter(Subtype_IHC == "TN" | Subtype_IHC == "HR-HER2+") %>%
  dplyr::filter(is.na(Herceptin) | Herceptin == "No") %>%
  dplyr::group_by(Subtype_IHC, Herceptin) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC Herceptin     N
# 1 TN          NA          120
# 2 HR-HER2+    No           46

clin_finher %>%
  dplyr::filter(Subtype_IHC == "TN" | Subtype_IHC == "HR-HER2+") %>%
  dplyr::filter(is.na(Herceptin) | Herceptin == "No") %>%
  dplyr::group_by(Subtype_IHC_2, Subtype_IHC, Herceptin, Stratum_hr_chemo) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC_2 Subtype_IHC Herceptin Stratum_hr_chemo     N
# 1 TN            TN          NA        Negative_DTX        60
# 2 TN            TN          NA        Negative_NVB        60
# 3 HER2          HR-HER2+    No        Negative_DTX        27
# 4 HER2          HR-HER2+    No        Negative_NVB        19



# Define variable
# >>>>>>>>>>>>>>>

inter_sum_finher_subtype <- list()



# Subtype interaction in noTRA + noHormone stratified on chemo(DTX|NVB)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_subtype[["all"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

      # Per sig mta interaction
        yy <- purrr::map(

          c("Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_inter(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              interaction_variable = "Subtype_IHC_2",
              strata_variable = "Stratum_hr_chemo",
              xdata = xdata
            )

          },

          event_type,
          xdata = clin
        )


        # Consolidate per sig stat and format
        bind_rows(yy) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH")
          )

  },
  clin = clin_finher %>%
    dplyr::filter(Subtype_IHC == "TN" | Subtype_IHC == "HR-HER2+") %>%
    dplyr::filter(is.na(Herceptin) | Herceptin == "No")
)

names(inter_sum_finher_subtype[["all"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_subtype later in the script




# Subtype interaction in noTRA + noHormone + DTX
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_subtype[["dtx"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)


    # Per sig mta interaction
    yy <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1",
        "Cholesterol2",
        # "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", #"CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        "Fibroblast"),

      function(sig, event_type, xdata){ # subtype,

        get_adj_inter(
          event_variable = str_c(event_type, "_Event"),
          time_variable = str_c(event_type, "_Time_Years"),
          biomarker = sig,
          interaction_variable = "Subtype_IHC_2",
          strata_variable = "Stratum_hr_chemo",
          xdata = xdata
        )

      },

      event_type,
      xdata = clin
    )


    # Consolidate per sig stat and format
    bind_rows(yy) %>%
      dplyr::mutate(
        P_adj = p.adjust(p = P, method = "BH")
      )

  },
  clin = clin_finher %>%
    dplyr::filter(Subtype_IHC == "TN" | Subtype_IHC == "HR-HER2+") %>%
    dplyr::filter(is.na(Herceptin) | Herceptin == "No") %>%
    dplyr::filter(Chemo == "DTX")
)

names(inter_sum_finher_subtype[["dtx"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_dtx_subtype later in the script




# Subtype interaction in noTRA + noHormone + NVB
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inter_sum_finher_subtype[["nvb"]] <- purrr::map(
  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)


    # Per sig mta interaction
    yy <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1",
        "Cholesterol2",
        # "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", #"CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        "Fibroblast"),

      function(sig, event_type, xdata){ # subtype,

        get_adj_inter(
          event_variable = str_c(event_type, "_Event"),
          time_variable = str_c(event_type, "_Time_Years"),
          biomarker = sig,
          interaction_variable = "Subtype_IHC_2",
          strata_variable = "Stratum_hr_chemo",
          xdata = xdata
        )

      },

      event_type,
      xdata = clin
    )


    # Consolidate per sig stat and format
    bind_rows(yy) %>%
      dplyr::mutate(
        P_adj = p.adjust(p = P, method = "BH")
      )

  },
  clin = clin_finher %>%
    dplyr::filter(Subtype_IHC == "TN" | Subtype_IHC == "HR-HER2+") %>%
    dplyr::filter(is.na(Herceptin) | Herceptin == "No") %>%
    dplyr::filter(Chemo == "NVB")
)

names(inter_sum_finher_subtype[["nvb"]]) <- c("RFS", "DDFS", "OS")
# !!!! Save inter_sum_finher_nvb_subtype later in the script


#
# ==============================================================================




# 6. Save Robjects of interaction summaries
# ==============================================================================


# FinHER interaction summary objects
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# "inter_sum_finher_mta"
# "inter_sum_finher_tra"
# "inter_sum_finher_subtype"

save(inter_sum_finher_mta, file = str_c(out_data,"inter_sum_finher_mta.RData"))
save(inter_sum_finher_tra, file = str_c(out_data,"inter_sum_finher_tra.RData"))
save(inter_sum_finher_subtype, file = str_c(out_data,"inter_sum_finher_subtype.RData"))


#
# ==============================================================================





