# finher_her2_ph_fix.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Acounting for non-PH in Denovo_TILsig * Chemo interaction models on RFS in HER2

# Note on PH fixation
# >>>>>>>>>>>>>>>>>>>
# Step-function ref:
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# http://www.sthda.com/english/wiki/cox-model-assumptions



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore non-PH in TRA interaction models


# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher_finneo.RData")



# Format clinical data
# >>>>>>>>>>>>>>>>>>>>

clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10)

clin_finher_finneo  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 HR-HER2+    HER2             89
# 2 HR+HER2+    HER2             91

#
# ==============================================================================




# 2. Explore non-PH in scaled_Denovo_TILsig * Chemo interaction on RFS in HER2
# ==============================================================================


xdata <- clin_finher_finneo  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")



# RFS HER2 chemo interaction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosis_rfs_her2_chemo_interaction_regular_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Denovo_TILsig")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(RFS_Time_Years, RFS_Event) ~", i,
                               "* Chemo + strata(Hormone, Herceptin)")),
    data = xdata
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_Denovo_TILsig: summary"
#                                    coef exp(coef)  se(coef)          z   Pr(>|z|)
# scaled_Denovo_TILsig           1.419365 4.1344924 0.9240416  1.5360396 0.12452866
# ChemoNVB                       1.135172 3.1117087 0.6262598  1.8126215 0.06989022
# scaled_Denovo_TILsig:ChemoNVB -1.172689 0.3095334 1.1993055 -0.9778069 0.32816985
# [1] "scaled_Denovo_TILsig: ph test"
#                            chisq df    p
# scaled_Denovo_TILsig        2.44  1 0.12
# Chemo                       0.12  1 0.73
# scaled_Denovo_TILsig:Chemo  3.86  1 0.05
# GLOBAL                      6.19  3 0.10




# RFS HER2 chemo interaction PH corrected
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Split time-to-event data to satisfy PH assumption
# Cut time points are derived from cox.zph() plot of scaled
# Schonfield residueals from the non-PH variable (scaled_General_Immune) and
# transformed time of the regular model.

xformula <- str_c("Surv(RFS_Time_Years, RFS_Event) ~ ",
                  "scaled_Denovo_TILsig + Herceptin + Hormone + Chemo")
xdata_ddfs_split <- survSplit(formula = as.formula(xformula),
                         data = clin_finher_finneo  %>%
                           dplyr::filter(Subtype_IHC_2 == "HER2"),
                         cut = c(1.43, 3), # approx median
                         # cut = c(2,4),
                         start = "RFS_Time_Start", # New variable introduced
                         end = "RFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                         zero = 0,
                         episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))



pdf("results/figures/PH_diagnosis_rfs_her2_chemo_interaction_phfixed_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Denovo_TILsig")){

  # TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo)
  # DDFS_Time_Start, DDFS_Time_End, DDFS_Event

  m1 <- coxph(
    formula = as.formula(str_c("Surv(RFS_Time_Start, RFS_Time_End, RFS_Event) ~", i,
                               "* Chemo * Interval_Id - Interval_Id + strata(Hormone, Herceptin)")),
    data = xdata_ddfs_split
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_Denovo_TILsig: summary"
#                                                   coef    exp(coef)  se(coef)           z   Pr(>|z|)
# scaled_Denovo_TILsig                        0.98361104   2.67409510 1.5448640  0.63669750 0.52432190
# ChemoNVB                                    1.94879425   7.02021787 0.9586494  2.03285393 0.04206728
# scaled_Denovo_TILsig:ChemoNVB              -2.98768222   0.05040413 2.0375245 -1.46632948 0.14255855
# scaled_Denovo_TILsig:Interval_Id2           1.86457911   6.45321919 2.1819737  0.85453785 0.39280704
# scaled_Denovo_TILsig:Interval_Id3          -0.59719144   0.55035518 2.4021319 -0.24860893 0.80366330
# ChemoNVB:Interval_Id2                      -0.49734617   0.60814242 1.5315240 -0.32473939 0.74537832
# ChemoNVB:Interval_Id3                      -2.34450327   0.09589482 1.5419410 -1.52048830 0.12838830
# scaled_Denovo_TILsig:ChemoNVB:Interval_Id2  0.07003944   1.07255048 3.0241902  0.02315973 0.98152286
# scaled_Denovo_TILsig:ChemoNVB:Interval_Id3  5.53679475 253.86299858 3.0554313  1.81211561 0.06996834
# [1] "scaled_Denovo_TILsig: ph test"
#                                         chisq df    p
# scaled_Denovo_TILsig                    0.738  1 0.39
# Chemo                                   1.845  1 0.17
# scaled_Denovo_TILsig:Chemo              0.639  1 0.42
# scaled_Denovo_TILsig:Interval_Id        0.288  2 0.87
# Chemo:Interval_Id                       1.316  2 0.52
# scaled_Denovo_TILsig:Chemo:Interval_Id  4.581  2 0.10
# GLOBAL                                 10.924  9 0.28



# Summary
# In HER2 chemo interaction on RFS, the Denovo_TILsig signature generated from
# FinHER TIL-H&E counts failed PH assumption.

# The PH fixed model using the splitted time to event data at cutoffs 1.43 and 3 years
# did not showed any significant interaction. Hence the PH failed model's estimates are
# used in the paper for consistency.

#
# ==============================================================================




