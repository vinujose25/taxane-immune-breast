# s8_finher_plots_mta_paper.R


load("results/data/finher_tn.RData")
load("results/data/finher_her2.RData")
load("results/data/clin_finher_finneo.RData")


# Finher MTA interaction sample size
# ==============================================================================
clin_finher_finneo  %>%
  dplyr::group_by(Chemo, HR_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())
#   Chemo HR_IHC   Subtype_IHC_2     N
# 1 DTX   Negative TN               60
# 2 DTX   Negative HER2             53
# 3 DTX   Positive HER2             34
# 4 NVB   Negative TN               60
# 5 NVB   Negative HER2             36
# 6 NVB   Positive HER2             57

#
# ==============================================================================




# FinHER KM curves
# ==============================================================================


# Main plot
# >>>>>>>>>

# TN subset
dat_km <- clin_finher_finneo %>%
  dplyr::filter(Subtype_IHC == "TN") %>%
  dplyr::rename(TIL = "StrLy_Mean")


# Relevant medians
dat_km_medians <- dat_km %>%
  dplyr::select(TIL,
                scaled_Denovo_TILsig,
                scaled_Denovo_Immune,
                scaled_Denovo_ECM) %>%
  purrr::map(~median(.x, na.rm = T))

dat_km_medians
# $TIL
# [1] 25
# $scaled_Denovo_TILsig
# [1] 0.5016433
# $scaled_Denovo_Immune
# [1] 0.5901034
# $scaled_Denovo_ECM
# [1] 0.5066598

# Grouping data
dat_km <- bind_rows(

  dat_km %>%
    dplyr::mutate(
      Subgroup = "TIL",
      Score = if_else(TIL < dat_km_medians$TIL, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Denovo_TILsig",
      Score = if_else(scaled_Denovo_TILsig < dat_km_medians$scaled_Denovo_TILsig, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Denovo_Immune",
      Score = if_else(scaled_Denovo_Immune < dat_km_medians$scaled_Denovo_Immune, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Denovo_ECM",
      Score = if_else(scaled_Denovo_ECM < dat_km_medians$scaled_Denovo_ECM, "Lo", "Hi")
    )
)


dat_km <- dat_km %>%  dplyr::mutate(
  # Subgroup = str_c(Chemo, "*", Subgroup)
  # Cat = str_c(Chemo,"*",Subgroup2, "-", Score)
  Cat = str_c(Chemo,"*Score-", Score," ")
)

dat_km %>%
  group_by(Cat, Score) %>%
  summarise(n = n())
#   Cat             Score     n
# 1 "DTX*Score-Hi " Hi      118
# 2 "DTX*Score-Lo " Lo      121
# 3 "NVB*Score-Hi " Hi      123
# 4 "NVB*Score-Lo " Lo      116
# 5  NA             NA        2



plot_finher_km_main(
  dat_km = dat_km,
  event_var = "DDFS_Event",
  time_var = "DDFS_Time_Years",
  prefix = "DDFS")

plot_finher_km_main(
  dat_km = dat_km,
  event_var = "RFS_Event",
  time_var = "RFS_Time_Years",
  prefix = "RFS")

plot_finher_km_main(
  dat_km = dat_km,
  event_var = "OS_Event",
  time_var = "OS_Time_Years",
  prefix = "OS")





# Supplementary plot
# >>>>>>>>>>>>>>>>>>

# TN subset
dat_km <- clin_finher_finneo %>%
  dplyr::filter(Subtype_IHC == "TN") %>%
  dplyr::rename(TIL = "StrLy_Mean")


# Relevant medians
dat_km_medians <- dat_km %>%
  dplyr::select(
    "scaled_Hamy2016_Immune", "scaled_Yang2018_Immune",
    "scaled_Gruosso2019_Interferon", "scaled_Farmer2009_MX1",
    "scaled_Hamy2016_Interferon", "scaled_Nirmal2018_Interferon",
    "scaled_Hamy2016_Ecm", "scaled_Naba2014_Ecmcore",
    "scaled_Triulzi2013_Ecm", "scaled_Sorrentino2014_Chol",
    "scaled_Pooled_Immune", "scaled_Pooled_Interferon",
    "scaled_Pooled_Fibrosis", "scaled_Pooled_Cholesterol"
  ) %>%
  purrr::map(~median(.x, na.rm = T))



# Grouping data
dat_km <- bind_rows(

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Hamy2016_Immune",
      Score = if_else(scaled_Hamy2016_Immune < dat_km_medians$scaled_Hamy2016_Immune, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Yang2018_Immune",
      Score = if_else(scaled_Yang2018_Immune < dat_km_medians$scaled_Yang2018_Immune, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Gruosso2019_Interferon",
      Score = if_else(scaled_Gruosso2019_Interferon < dat_km_medians$scaled_Gruosso2019_Interferon, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Farmer2009_MX1",
      Score = if_else(scaled_Farmer2009_MX1 < dat_km_medians$scaled_Farmer2009_MX1, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Hamy2016_Interferon",
      Score = if_else(scaled_Hamy2016_Interferon < dat_km_medians$scaled_Hamy2016_Interferon, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Nirmal2018_Interferon",
      Score = if_else(scaled_Nirmal2018_Interferon < dat_km_medians$scaled_Nirmal2018_Interferon, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Hamy2016_Ecm",
      Score = if_else(scaled_Hamy2016_Ecm < dat_km_medians$scaled_Hamy2016_Ecm, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Naba2014_Ecmcore",
      Score = if_else(scaled_Naba2014_Ecmcore < dat_km_medians$scaled_Naba2014_Ecmcore, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Triulzi2013_Ecm",
      Score = if_else(scaled_Triulzi2013_Ecm < dat_km_medians$scaled_Triulzi2013_Ecm, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Sorrentino2014_Chol",
      Score = if_else(scaled_Sorrentino2014_Chol < dat_km_medians$scaled_Sorrentino2014_Chol, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Pooled_Immune",
      Score = if_else(scaled_Pooled_Immune < dat_km_medians$scaled_Pooled_Immune, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Pooled_Interferon",
      Score = if_else(scaled_Pooled_Interferon < dat_km_medians$scaled_Pooled_Interferon, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Pooled_Fibrosis",
      Score = if_else(scaled_Pooled_Fibrosis < dat_km_medians$scaled_Pooled_Fibrosis, "Lo", "Hi")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Pooled_Cholesterol",
      Score = if_else(scaled_Pooled_Cholesterol < dat_km_medians$scaled_Pooled_Cholesterol, "Lo", "Hi")
    )
)


dat_km <- dat_km %>%  dplyr::mutate(
  # Subgroup = str_c(Chemo, "*", Subgroup)
  Cat = str_c(Chemo,"*Score-", Score," ")
)

dat_km %>%
  group_by(Subgroup, Score) %>%
  summarise(n = n())
#   Subgroup               Score     n
# 1 Farmer2009_MX1         Hi       60
# 2 Farmer2009_MX1         Lo       60
# 3 Gruosso2019_Interferon Hi       60
# 4 Gruosso2019_Interferon Lo       60
# 5 Hamy2016_Ecm           Hi       60
# 6 Hamy2016_Ecm           Lo       60
# 7 Hamy2016_Immune        Hi       60
# 8 Hamy2016_Immune        Lo       60
# 9 Hamy2016_Interferon    Hi       60
# 10 Hamy2016_Interferon    Lo       60
# ...


plot_finher_km_supp(
  dat_km = dat_km,
  event_var = "DDFS_Event",
  time_var = "DDFS_Time_Years",
  prefix = "DDFS")

plot_finher_km_supp(
  dat_km = dat_km,
  event_var = "RFS_Event",
  time_var = "RFS_Time_Years",
  prefix = "RFS")

plot_finher_km_supp(
  dat_km = dat_km,
  event_var = "OS_Event",
  time_var = "OS_Time_Years",
  prefix = "OS")


#
# ==============================================================================




# Finher MTA interaction
# ==============================================================================


# Individual sigs

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_individual.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_individual.sig),
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_individual.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_individual.sig),
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_individual.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_individual.sig),
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)




# Pooled sigs

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_pooled.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_pooled.sig),
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_pooled.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_pooled.sig),
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)

plot_finher_mta_interaction(
  tn_cox_sum = list(inter.chemo = finher_tn$inter.chemo_pooled.sig),
  her2_cox_sum =  list(inter.chemo = finher_her2$inter.chemo_pooled.sig),
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)


#
# ==============================================================================



# Finher prognosis
# ==============================================================================



# Individual sig

plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_individual.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_individual.sig),
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)


plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_individual.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_individual.sig),
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)


plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_individual.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_individual.sig),
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_individual-sig"),
  width = 7.5,
  height = 6.5)



# Pooled.sig

plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_pooled.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_pooled.sig),
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)

plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_pooled.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_pooled.sig),
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)

plot_finher_prognosis(
  tn_cox_sum = list(prog = finher_tn$prog_pooled.sig),
  her2_cox_sum =  list(prog = finher_her2$prog_pooled.sig),
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_pooled-sig"),
  width = 7.5,
  height = 4.25)

#
# ==============================================================================



# Clear memory
# ==============================================================================

rm(clin_finher_finneo, dat_km, dat_km_medians)

#
# ==============================================================================


