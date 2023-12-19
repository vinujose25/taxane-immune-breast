# finher_her2_ph_fix.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Acounting for non-PH in Chemo interaction models in HER2
# The following models failed PH asumptions in both prognostic and interaction modelling;
# General_Immune(DDFS and RFS), Control_Proliferation (OS).

# Note on PH fixation
# >>>>>>>>>>>>>>>>>>>
# Step-function ref:
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# http://www.sthda.com/english/wiki/cox-model-assumptions



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore non-PH in TRA interaction models (TIL, Fibrosis1-3, and Fibroblast)


# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")



# Format clinical data
# >>>>>>>>>>>>>>>>>>>>

clin_finher <- clin_finher %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10)

clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 HR-HER2+    HER2             89
# 2 HR+HER2+    HER2             91

#
# ==============================================================================





# 2. Explore non-PH on DDFS in HER2 prognosis
# ==============================================================================


xdata <- clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")



# DDFS regular HER2 prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosis_ddfs_her2_prognosis_regular_model2.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i,
                               "+ Hormone + strata(Herceptin, Chemo)")),
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

# [1] "scaled_General_Immune: summary"
# coef exp(coef)  se(coef)          z  Pr(>|z|)
# scaled_General_Immune -0.4797089 0.6189635 0.5751409 -0.8340721 0.4042404
# HormonenoHR            0.3923379 1.4804378 0.3301817  1.1882485 0.2347355
# [1] "scaled_General_Immune: ph test"
# chisq df      p
# scaled_General_Immune  8.75  1 0.0031
# Hormone                3.36  1 0.0667
# GLOBAL                10.22  2 0.0060



# DDFS PH correction using step function in HER2 prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Step-function ref:
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# http://www.sthda.com/english/wiki/cox-model-assumptions


xformula <- str_c("Surv(DDFS_Time_Years, DDFS_Event) ~ ",
                  "scaled_General_Immune + scaled_Control_Proliferation",
                  "+ Herceptin + Hormone + Chemo")
xdata_ddfs_split <- survSplit(formula = as.formula(xformula),
                              data = clin_finher  %>%
                                dplyr::filter(Subtype_IHC_2 == "HER2"),
                              cut = c(1, 2.8),
                              # Cutpoint derived from "PH_diagnosis_ddfs_her2_prognosis_regular_model.pdf"
                              # Only one point fluctuation.
                              start = "DDFS_Time_Start", # New variable introduced
                              end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                              zero = 0,
                              episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))


pdf("results/figures/PH_diagnosis_ddfs_her2_prognosis_phfixed_model2.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  # TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo)
  # DDFS_Time_Start, DDFS_Time_End, DDFS_Event
  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~", i,
                               "* Interval_Id - Interval_Id + Hormone + strata(Herceptin, Chemo)")),
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

# [1] "scaled_General_Immune: summary"
#                                          coef  exp(coef)  se(coef)         z   Pr(>|z|)
# scaled_General_Immune               0.9138589 2.49392771 0.8765107  1.042610 0.29712896
# scaled_General_Immune:Interval_Id2 -2.7497280 0.06394525 1.2107882 -2.271023 0.02314558

# [1] "scaled_General_Immune: ph test"
#                                   chisq df     p
# scaled_General_Immune              3.11  1 0.078
# scaled_General_Immune:Interval_Id  4.98  1 0.026
# GLOBAL                             5.00  2 0.082

# !!!!!!!!!!!!!!!!!!!!!!
# Note:
# Step functions (by time split data using survsplit()) did not fix PH violation.

# Seperate analysis for HR+ve and HR-ve as these have profound influence on
# breast cancer patient survival.



# DDFS regular HER2+HR-ve prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosis_ddfs_her2+hr-_prognosis_regular_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i,
                               "+ strata(Hormone, Herceptin, Chemo)")),
    data = xdata %>%
      dplyr::filter(Subtype_IHC == "HR-HER2+")
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_General_Immune: summary"
#                            coef exp(coef)  se(coef)         z  Pr(>|z|)
# scaled_General_Immune 0.3287524  1.389234 0.7984622 0.4117319 0.6805359
# [1] "scaled_General_Immune: ph test"
#                       chisq df    p
# scaled_General_Immune 0.941  1 0.33
# GLOBAL                0.941  1 0.33

# In HER2+HR-ve, General_Immune model follow PH assumption and
# has no significant prognostic value within this subgroup.



# DDFS regular HER2+HR+ve prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosis_ddfs_her2+hr+_prognosis_regular_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i,
                               "+ strata(Hormone, Herceptin, Chemo)")),
    data = xdata %>%
      dplyr::filter(Subtype_IHC == "HR+HER2+")
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_General_Immune: summary"
#                            coef exp(coef)  se(coef)         z   Pr(>|z|)
# scaled_General_Immune -1.515807 0.2196308 0.8831161 -1.716431 0.08608327
# [1] "scaled_General_Immune: ph test"
#                       chisq df     p
# scaled_General_Immune   6.1  1 0.014
# GLOBAL                  6.1  1 0.014

# Failed PH
# Survsplit in HR+HER2+ using cut points from "PH_diagnosis_ddfs_her2+hr+_prognosis_regular_model.pdf


# Survsplit HR+HER2+ for fixing prognostic PH
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

xformula <- str_c("Surv(DDFS_Time_Years, DDFS_Event) ~ ",
                  "scaled_General_Immune + scaled_Control_Proliferation",
                  "+ Herceptin + Hormone + Chemo")
xdata_ddfs_split <- survSplit(formula = as.formula(xformula),
                              data = clin_finher  %>%
                                dplyr::filter(Subtype_IHC == "HR+HER2+"),
                              cut = c(1.7, 3.8),
                              # Cutpoints derived from "PH_diagnosis_ddfs_her2+hr+_prognosis_regular_model.pdf"
                              # Only one point fluctuation.
                              start = "DDFS_Time_Start", # New variable introduced
                              end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                              zero = 0,
                              episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))


# DDFS PH fixed HER2+HR+ve prognosis
pdf("results/figures/PH_diagnosis_ddfs_her2+hr+_prognosis_phfixed_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~", i,
                               "* Interval_Id - Interval_Id + strata(Hormone, Herceptin, Chemo)")),
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

# [1] "scaled_General_Immune: summary"
#                                         coef   exp(coef) se(coef)          z   Pr(>|z|)
# scaled_General_Immune               1.089427 2.972570078 1.845288  0.5903831 0.55493384
# scaled_General_Immune:Interval_Id2 -2.224379 0.108134515 2.272017 -0.9790327 0.32756382
# scaled_General_Immune:Interval_Id3 -5.902404 0.002732866 2.754754 -2.1426245 0.03214326 * Significant

# [1] "scaled_General_Immune: ph test"
# chisq df    p
# scaled_General_Immune              1.24  1 0.27
# scaled_General_Immune:Interval_Id  3.44  2 0.18
# GLOBAL                             5.02  3 0.17


# Summary
# In HER2 prognosis on DDFS/RFS (both with similar events) among different immune
# signatures derived from multiple ways, the General_immune signature generated by
# pooling published immune signatures failed PH assumption.

# On interrogation, it is revealed that, the signature followed PH assumption in
# HER2+HR-ve subgroup (but insignificant General_immune signature),
# and PH failure is limited to HER2+HR+ve subgroup.
# The fitted line in the scatter plot between scaled Schoenfeld residuals and
# transformed time suggests two cut points (curve fluctuations) to split the
# time to event data into three intervals into which to fit a step function for interaction.

# This PH fixed model followed PH assumption and a significant difference in
# survival data in the last (3rd interval) interval in HER2+HR+ subgroup.

# However, given that among the multiple immune signatures derived from
# different strategies, only the pooled published immune signature has got the
# PH assumption issue. We ignored this PH failure, as it could be
# due to the confounding effect (by non-immune signal) due to pooling from multiple studies.

#
# ==============================================================================






# 3. Explore non-PH on DDFS in HER2 chemo interaction
# ==============================================================================


xdata <- clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")



# DDFS HER2 chemo interaction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosis_ddfs_her2_chemo_interaction_regular_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i,
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

# [1] "scaled_General_Immune: summary"
#                                      coef exp(coef)  se(coef)          z   Pr(>|z|)
# scaled_General_Immune           0.5741963 1.7757028 1.0620171  0.5406657 0.58873801
# ChemoNVB                        1.5124865 4.5380003 0.8023419  1.8850897 0.05941772
# scaled_General_Immune:ChemoNVB -1.4464162 0.2354124 1.2586019 -1.1492246 0.25046338
# [1] "scaled_General_Immune: ph test"
#                             chisq df     p
# scaled_General_Immune       5.987  1 0.014
# Chemo                       0.100  1 0.752
# scaled_General_Immune:Chemo 0.671  1 0.413
# GLOBAL                      6.664  3 0.083




# DDFS HER2 chemo interaction PH corrected
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Split time-to-event data to satisfy PH assumption
# Cut time points are derived from cox.zph() plot of scaled
# Schonfield residueals from the non-PH variable (scaled_General_Immune) and
# transformed time of the regular model.

xformula <- str_c("Surv(DDFS_Time_Years, DDFS_Event) ~ ",
                  "scaled_General_Immune + scaled_Control_Proliferation",
                  "+ Herceptin + Hormone + Chemo")
xdata_ddfs_split <- survSplit(formula = as.formula(xformula),
                         data = clin_finher  %>%
                           dplyr::filter(Subtype_IHC_2 == "HER2"),
                         cut = c(1.3), # approx median
                         # cut = c(2,4),
                         start = "DDFS_Time_Start", # New variable introduced
                         end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                         zero = 0,
                         episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))



pdf("results/figures/PH_diagnosis_ddfs_her2_chemo_interaction_phfixed_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_General_Immune")){

  # TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo)
  # DDFS_Time_Start, DDFS_Time_End, DDFS_Event

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~", i,
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

# 1] "scaled_General_Immune: summary"
#                                                   coef exp(coef) se(coef)          z  Pr(>|z|)
# scaled_General_Immune                        2.2805289 9.7818525 1.865489  1.2224829 0.2215251
# ChemoNVB                                     1.8827153 6.5713237 1.622115  1.1606544 0.2457825
# scaled_General_Immune:ChemoNVB              -1.8525150 0.1568422 2.197332 -0.8430745 0.3991868
# scaled_General_Immune:Interval_Id2          -2.6229129 0.0725911 2.275653 -1.1525977 0.2490755
# ChemoNVB:Interval_Id2                       -0.5437291 0.5805792 1.859881 -0.2923461 0.7700220
# scaled_General_Immune:ChemoNVB:Interval_Id2  0.6346289 1.8863221 2.693598  0.2356064 0.8137381
# [1] "scaled_General_Immune: ph test"
#                                           chisq df     p
# scaled_General_Immune                   2.94974  1 0.086
# Chemo                                   0.51649  1 0.472
# scaled_General_Immune:Chemo             0.00117  1 0.973
# scaled_General_Immune:Interval_Id       3.11239  1 0.078
# Chemo:Interval_Id                       0.62029  1 0.431
# scaled_General_Immune:Chemo:Interval_Id 0.07990  1 0.777
# GLOBAL                                  8.49005  6 0.204



# Summary
# In HER2 chemo interaction on DDFS/RFS (both with similar events) among different immune
# signatures derived from multiple ways, the General_immune signature generated by
# pooling published immune signatures failed PH assumption.

# The PH fixed model using the splitted time to event adta at cut of 1.3 years
# did not showed any significant interaction. Hence the PH failed model estimates are
# used in the paper for consistency.

#
# ==============================================================================





# 4. Explore non-PH on OS in HER2 prognosis
# ==============================================================================


xdata <- clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")



# OS regular and PH corrected models
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# OS regular HER2 prognosis
pdf("results/figures/PH_diagnosis_os_her2_prognosis_regular_model2.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Control_Proliferation")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(OS_Time_Years, OS_Event) ~", i,
                               "+ Hormone + strata(Herceptin, Chemo)")),
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

# [1] "scaled_Control_Proliferation: summary"
#                                  coef exp(coef)  se(coef)        z  Pr(>|z|)
# scaled_Control_Proliferation 1.150943  3.161172 0.8760245 1.313825 0.1889052
# [1] "scaled_Control_Proliferation: ph test"
#                              chisq df      p
# scaled_Control_Proliferation  7.48  1 0.0062
# GLOBAL                        7.48  1 0.0062




# OS PH correction using step function in HER2 prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Step-function ref:
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# http://www.sthda.com/english/wiki/cox-model-assumptions


xformula <- str_c("Surv(OS_Time_Years, OS_Event) ~ ",
                  "scaled_General_Immune + scaled_Control_Proliferation",
                  "+ Herceptin + Hormone + Chemo")
xdata_os_split <- survSplit(formula = as.formula(xformula),
                              data = clin_finher  %>%
                                dplyr::filter(Subtype_IHC_2 == "HER2"),
                              cut = c(2.3, 3.5),
                              # Cutpoint derived from "PH_diagnosis_os_her2_prognosis_regular_model.pdf"
                              start = "OS_Time_Start", # New variable introduced
                              end = "OS_Time_End", # Renaming of OS_Time_Years to avoid confusion
                              zero = 0,
                              episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))



# OS PH corrected HER2 prognosis

pdf("results/figures/PH_diagnosis_os_her2_prognosis_phfixed_model2.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Control_Proliferation")){

  # TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo)
  # DDFS_Time_Start, DDFS_Time_End, DDFS_Event
  m1 <- coxph(
    formula = as.formula(str_c("Surv(OS_Time_Start, OS_Time_End, OS_Event) ~", i,
                               "* Interval_Id - Interval_Id + Hormone + strata(Herceptin, Chemo)")),
    data = xdata_os_split
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_Control_Proliferation: summary"
#                                                coef    exp(coef) se(coef)         z    Pr(>|z|)
# scaled_Control_Proliferation              -2.680072   0.06855825 1.734507 -1.545149 0.122310333
# scaled_Control_Proliferation:Interval_Id2  4.904171 134.85105190 2.493038  1.967147 0.049166294 * significant
# scaled_Control_Proliferation:Interval_Id3  6.070132 432.73768932 2.252611  2.694709 0.007045018 ** significant

# [1] "scaled_Control_Proliferation: ph test"
#                                          chisq df    p
# scaled_Control_Proliferation             0.457  1 0.50
# scaled_Control_Proliferation:Interval_Id 4.093  2 0.13
# GLOBAL                                   4.242  3 0.24


# Summary
# Proliferation is significantly associated with survival on interval 2 and 3
# in the PH fixed model

#
# ==============================================================================




# 5. Explore non-PH on OS in HER2 chemo interaction
# ==============================================================================



xdata <- clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")



# OS regular HER2 chemo interaction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/PH_diagnosiss_os_her2_chemo_interaction_regular_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Control_Proliferation")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(OS_Time_Years, OS_Event) ~", i,
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

# [1] "scaled_Control_Proliferation: summary"
#                                              coef exp(coef) se(coef)           z  Pr(>|z|)
# scaled_Control_Proliferation           0.09936608 1.1044705 1.449247  0.06856391 0.9453367
# ChemoNVB                              -0.08290406 0.9204395 1.020652 -0.08122654 0.9352618
# scaled_Control_Proliferation:ChemoNVB  1.75554363 5.7865927 1.802455  0.97397356 0.3300697
# [1] "scaled_Control_Proliferation: ph test"
#                                    chisq df      p
# scaled_Control_Proliferation       8.445  1 0.0037
# Chemo                              0.205  1 0.6511
# scaled_Control_Proliferation:Chemo 0.891  1 0.3452
# GLOBAL                             9.035  3 0.0288



# OS PH correction using step function in HER2 prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Step-function ref:
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# http://www.sthda.com/english/wiki/cox-model-assumptions


xformula <- str_c("Surv(OS_Time_Years, OS_Event) ~ ",
                  "scaled_General_Immune + scaled_Control_Proliferation",
                  "+ Herceptin + Hormone + Chemo")
xdata_os_split <- survSplit(formula = as.formula(xformula),
                            data = clin_finher  %>%
                              dplyr::filter(Subtype_IHC_2 == "HER2"),
                            cut = c(2.2, 3.3),
                            # Cutpoint derived from "PH_diagnosiss_os_her2_chemo_interaction_regular_model.pdf"
                            start = "OS_Time_Start", # New variable introduced
                            end = "OS_Time_End", # Renaming of OS_Time_Years to avoid confusion
                            zero = 0,
                            episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))




# OS HER2 chemo interaction PH corrected

pdf("results/figures/PH_diagnosiss_os_her2_chemo_interaction_phfixed_model.pdf",
    width = 5, height = 5, onefile = T)

for(i in c("scaled_Control_Proliferation")){

  # TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo)
  # DDFS_Time_Start, DDFS_Time_End, DDFS_Event

  m1 <- coxph(
    formula = as.formula(str_c("Surv(OS_Time_Start, OS_Time_End, OS_Event) ~", i,
                               "* Chemo * Interval_Id - Interval_Id + strata(Hormone, Herceptin)")),
    data = xdata_os_split
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "scaled_Control_Proliferation: summary"
#                                                          coef    exp(coef) se(coef)          z   Pr(>|z|)
# scaled_Control_Proliferation                       -4.7998645 8.230862e-03 2.658728 -1.8053238 0.07102403
# ChemoNVB                                           -0.8777973 4.156975e-01 1.348706 -0.6508442 0.51514705
# scaled_Control_Proliferation:ChemoNVB               4.6205856 1.015535e+02 3.340292  1.3832881 0.16657654
# scaled_Control_Proliferation:Interval_Id2           6.6316410 7.587262e+02 4.206919  1.5763651 0.11494168
# scaled_Control_Proliferation:Interval_Id3           8.1160307 3.347706e+03 3.577053  2.2689152 0.02327349
# ChemoNVB:Interval_Id2                               2.4337242 1.140126e+01 2.691405  0.9042580 0.36585862
# ChemoNVB:Interval_Id3                               0.9900094 2.691260e+00 2.482496  0.3987960 0.69004354
# scaled_Control_Proliferation:ChemoNVB:Interval_Id2 -5.3439324 4.777049e-03 5.051030 -1.0579887 0.29006058
# scaled_Control_Proliferation:ChemoNVB:Interval_Id3 -3.6803039 2.521531e-02 4.525107 -0.8133075 0.41604175
# [1] "scaled_Control_Proliferation: ph test"
#                                                chisq df    p
# scaled_Control_Proliferation                   1.444  1 0.23
# Chemo                                          0.274  1 0.60
# scaled_Control_Proliferation:Chemo             0.579  1 0.45
# scaled_Control_Proliferation:Interval_Id       2.990  2 0.22
# Chemo:Interval_Id                              0.148  2 0.93
# scaled_Control_Proliferation:Chemo:Interval_Id 1.272  2 0.53
# GLOBAL                                         7.696  9 0.57



# Summary
# In HER2 chemo interaction on OS the General_proliferation signature generated by
# pooling published immune signatures failed PH assumption.

# The PH fixed model using the splitted time to event data at cut of 2.2 & 3.3 years
# did not showed any significant interaction. Hence the PH failed model estimates are
# used in the paper for consistency.



#
# ==============================================================================

