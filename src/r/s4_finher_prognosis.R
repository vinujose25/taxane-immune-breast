# s3_exploring_prognosis.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the prognostic value of TIL, TIL associated signatures and Immune celltypes.
# Account for arm heterogenity.
# Prepare data necessary to plot prognostic plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis.
# 3. Save Robjects



# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")


# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

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
# 1  HR-HER2+             FEC_DTX 27   7    7  5 25.92593  25.92593 18.518519
# 2  HR-HER2+         FEC_DTX_TRA 26   3    3  2 11.53846  11.53846  7.692308
# 3  HR-HER2+             FEC_NVB 19   6    6  6 31.57895  31.57895 31.578947
# 4  HR-HER2+         FEC_NVB_TRA 17   6    5  3 35.29412  29.41176 17.647059
# 5  HR+HER2+     FEC_DTX_Hormone 18   4    4  2 22.22222  22.22222 11.111111
# 6  HR+HER2+ FEC_DTX_TRA_Hormone 16   2    0  0 12.50000   0.00000  0.000000
# 7  HR+HER2+     FEC_NVB_Hormone 27   7    6  3 25.92593  22.22222 11.111111
# 8  HR+HER2+ FEC_NVB_TRA_Hormone 30  10    9  4 33.33333  30.00000 13.333333
# 9        TN             FEC_DTX 60  14   13  9 23.33333  21.66667 15.000000
# 10       TN             FEC_NVB 60  21   19 14 35.00000  31.66667 23.333333


# Analysis
# 1. Per subtype prognosis (arm heterogenity)
# 2. Pan-subtype prognosis (subtype-arm heterogenity)

# Note: By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to clim a significant treatment
# interaction.


# Note !!!!!
# The event rate in HR+ subtype is low

#
# ==============================================================================





# 2. Explore prognosis.
# ==============================================================================

# Per subtype prognosis (Arm (~ Strata) heterogenity)


# Reference
# >>>>>>>>>

# altmeta::metahet(y, s2, n.resam = 1000)
# y	: a numeric vector indicating the observed effect sizes in the collected studies; they are assumed to be normally distributed.
# s2	: a numeric vector indicating the within-study variances.#
# n.resam	: a positive integer indicating the number of resampling iterations for calculating p-values of test statistics and 95% confidence interval of heterogeneity measures.

# Reference: Negative value of I2 (I^2) set as zero
# https://www.researchgate.net/post/How-should-I-consider-a-negative-value-of-the-I-squared-during-the-heterogeneity-assessment-in-meta-analyses




# Analysis/plot structure
# >>>>>>>>>>>>>>>>>>>>>>>

# Basic prognosis plot strucutre ( per subtype - per sig - arm heterogenity)
# (!!! Note : This plot will be generated in the NEXT SECTION of the script)
#
# TN - sig1 - arm1 - prognosis
# TN - sig1 - arm2 - prognosis
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# TN - sig1 - armX - prognosis - arm heterogenity (global)
# Note:
# In prognosis forest plot include the global prognosis (pooled arm) for sig1 per subtype.
# The global analysis account for arm heterogenity,
# and it also asseses the heterogenity level.


# Higher level prognosis plot strucutre
# TN/HER2 (event summary for each subtype, event/size: e/n)
# (!!! Note : This plot will be generated in the THIS SECTION of the script)
#
# immune sig1 - armX - prognosis - arm heterogenity (global)
# immune sig2 - armX - prognosis - arm heterogenity (global)
# immune sig3 - armX - prognosis - arm heterogenity (global)
# interferon sig1 - armX - prognosis - arm heterogenity (global)
# interferon sig2 - armX - prognosis - arm heterogenity (global)
# interferon sig3 - armX - prognosis - arm heterogenity (global)
# fibrosis sig1 - armX - prognosis - arm heterogenity (global)
# fibrosis sig2 - armX - prognosis - arm heterogenity (global)
# fibrosis sig3 - armX - prognosis - arm heterogenity (global)
# cholesterol sig1 - armX - prognosis - arm heterogenity (global)
# cholesterol sig2 - armX - prognosis - arm heterogenity (global)
# cholesterol sig3 - armX - prognosis - arm heterogenity (global)
# proliferation sig1 - armX - prognosis - arm heterogenity (global)
# proliferation sig2 - armX - prognosis - arm heterogenity (global)
# proliferation sig3 - armX - prognosis - arm heterogenity (global)



# Get prognostic effect and heterogenity (interaction) measurement
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Prognosis summary
prog_sum_finher <- purrr::map(
  c("RFS", "DDFS", "OS"),
  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("TN", "HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)
        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")

        # xdata <- filter_no_event_strata(
        #   x = subset(x = clin,
        #              subset = (clin[, subtype_varible] == subtype)),
        #   strata_variable = "Arm",
        #   event_variable = str_c(event_type, "_Event"))

        # The idea behind filtering no-event strata is to ignore convergense issue
        # warning due to no-event strata during interaction test(heterogenity test).
        # Convergense issue will be problematic if we use coeffcient and
        # wald test results of no-event interacting strata as it will be unreliable.
        # However p-value from likilihood ratio test is unaffected.
        # See ?coxph under Convergense section
        #
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

        xdata <- subset(x = clin,
                     subset = (clin[, subtype_varible] == subtype))


        # Per sig prognosis + heterogenity
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

            get_adj_prog(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
              biomarker = sig,
              strata_variable = "Arm",
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
# There were 50 or more warnings (use warnings() to see the first 50)
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#      Loglik converged before variable  2 ; coefficient may be infinite.

names(prog_sum_finher) <- c("RFS", "DDFS", "OS")


# !!!! Save prog_sum_finher later in the script

#
# ==============================================================================



# 3. Save Robjects
# ==============================================================================

# prog_sum : Prognosis summary
save(prog_sum_finher, file = str_c(out_data,"prog_sum_finher.RData"))

#
# ==============================================================================

