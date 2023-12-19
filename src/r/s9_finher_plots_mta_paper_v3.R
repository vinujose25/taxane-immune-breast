# finher_plots_mta_paper.R


load("results/data/finher_tn_v2.RData")
load("results/data/finher_her2_v2.RData")

load("results/data/clin_finher_finneo.RData")

# For MTA interaction PH is valid for all her2 models. Hence time interaction
# analysis in combination with survsplit() is not needed
# load("results/data/finher_her2_split.RData")
# load("results/data/clin_finher_finneo.RData") # for sig correlations and km plots


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
#   Subgroup           Score     n
#   <chr>              <chr> <int>
# 1 Denovo_ECM-High    High     60
# 2 Denovo_ECM-Low     Low      60
# 3 Denovo_Immune-High High     60
# 4 Denovo_Immune-Low  Low      60
# 5 Denovo_TILsig-High High     60
# 6 Denovo_TILsig-Low  Low      60
# 7 TIL-High           High     61
# 8 TIL-Low            Low      57
# 9 NA                 NA        2



plot_finher_km_main <- function(
  dat_km,
  event_var,
  time_var,
  prefix){

  # dat_km = dat_km
  # event_var = "DDFS_Event"
  # time_var = "DDFS_Time_Years"
  # prefix = "DDFS"

  # Model fitting
  km_fit = survminer::surv_fit(
    # as.formula(str_c("Surv(", time_var, ",", event_var, ") ~ Score")),
    # as.formula(str_c("Surv(", time_var, ",", event_var, ") ~ Chemo")),
    as.formula(str_c("Surv(", time_var, ",", event_var, ") ~ Cat")),
    group.by = "Subgroup", # var in data_km
    data = as.data.frame(dat_km)
  )

  # names(km_fit) will be used as title by ggsurvplot()
  names(km_fit) = str_split_fixed(
    str_split_fixed(names(km_fit), "::", 2)[ ,1],
    "\\.",
    2)[ ,2]



  # col_pal <- list()
  # col_pal[["Dark_red"]] <- "#e41a1c" #"#67001f"
  # col_pal[["Light_red"]] <- "#fb9a99" #"#e41a1c"
  # col_pal[["Dark_blue"]] <- "#377eb8" #"#053061"
  # col_pal[["Light_blue"]] <- "#a6cee3" #"#377eb8"




  km_fig <- ggsurvplot(
    fit = km_fit,
    pval = T,
    # pval.method = T,
    risk.table = T,
    newpage = F,
    # palette  = c(col_pal$Dark_red, col_pal$Light_red,
    #              col_pal$Dark_blue, col_pal$Light_blue),
    palette  = c(
      # "#e41a1c", #"#000000", #"black" DTX - high
      # "#67001f", #"red" DTX - low
      # "#377eb8", # blue NVB - High
      # "#053061" #"#bdbdbd" #"gray" NVB - low
      "#e41a1c", #"#000000", #"black" DTX - high
      "#f4a582", #"red" DTX - low
      "#377eb8", # blue NVB - High
      "#92c5de" #"#bdbdbd" #"gray" NVB - low
    ),
    xlab = "Years",
    ylab = str_c(prefix, " probability"),
    break.time.by = 1,
    xlim = c(0,6),
    legend.title = "",

    risk.table.title = element_blank(),
    risk.table.fontsize = 3,
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = T,     # show bars instead of names in text annotations
                               # in legend of risk table
    risk.tables.col = "Score",

    ggtheme = theme_bw() +
      theme(panel.grid = element_blank(),
            plot.margin = margin(l=1, r=1, unit = "lines"),
            panel.spacing.y = unit(0,"lines"),
            axis.title.x.bottom = element_text(margin = margin()),
            # axis.title.y.left = element_text(vjust = -4),
            axis.title.y.left = element_text(vjust = 0.5),
            plot.title = element_text(vjust = 0, hjust = .5, face = "bold"),
            legend.direction = "vertical"),

    tables.theme = theme_cleantable() +
      theme(
        axis.title.x.bottom = element_blank()
        # axis.title.x.bottom = element_text(color = "White")
      )
  )


  km_fig_legend <- get_legend(km_fig[[1]])
  # grid.draw(km_fig_legend)

  # Removing legends
  km_fig <- purrr::map(km_fig,
                       function(x){
                         x + guides(color = "none")
                       })

  # arrange_ggsurvplots() only support arranging by column.
  # Rather than setting factor levels it easy to rearrange the list
  # km_fig <- km_fig[ c("Immune1", "Tcell","Interferon1", "TIL") ]

  # km_fig <- km_fig[c("DTX*Immune1", "NVB*Immune1",
  #                    "DTX*Interferon1","NVB*Interferon1",
  #                    "DTX*Tcell", "NVB*Tcell",
  #                    "DTX*TIL", "NVB*TIL" ) ]

  # km_fig <- km_fig[c("DTX*TIL", "DTX*Denovo_TILsig",
  #                    "DTX*Denovo_Immune", "DTX*Denovo_ECM",
  #
  #                    "NVB*TIL", "NVB*Denovo_TILsig",
  #                    "NVB*Denovo_Immune", "NVB*Denovo_ECM") ]

  # km_fig <- km_fig[c("TIL-High", "Denovo_TILsig-High",
  #                    "Denovo_Immune-High", "Denovo_ECM-High",
  #
  #                    "TIL-Low", "Denovo_TILsig-Low",
  #                    "Denovo_Immune-Low", "Denovo_ECM-Low")]

  km_fig <- km_fig[c("TIL", "Denovo_Immune",
                     "Denovo_TILsig","Denovo_ECM")]

  zz <- arrange_ggsurvplots(
    x = km_fig, print = F,
    ncol = 2,
    # nrow = 4,
    nrow = 2,
    surv.plot.height = .75,
    risk.table.height = .25,
    newpage = F,
    byrow =T
  )
  # pdf(file = paste0(out_figures, "finher_km_", prefix, "_v3.pdf"),
  #     width = 7.5, height = 8)
  png(filename = paste0(out_figures, "finher_km_", prefix, "_v3.png"),
      width = 7.5, height = 7, units = "in", res = 300)
  print(zz)
  dev.off()


  # Risk table is text annotated; no need of legend

  # pdf(file = paste0(out_figures, "finher_km_", prefix, "legend.pdf"),
  #     width = 1.85, height = 1.25)
  # grid.draw(km_fig_legend)
  # dev.off()

}


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

plot_finher_km_supp <- function(
  dat_km,
  event_var,
  time_var,
  prefix){

  # dat_km = dat_km
  # event_var = "DDFS_Event"
  # time_var = "DDFS_Time_Years"
  # prefix = "DDFS"


  # Model fitting
  km_fit = survminer::surv_fit(
    as.formula(str_c("Surv(", time_var, ",", event_var, ") ~ Cat")),
    group.by = "Subgroup", # var in data_km
    data = as.data.frame(dat_km)
  )

  # names(km_fit) will be used as title by ggsurvplot()
  names(km_fit) = str_split_fixed(
    str_split_fixed(names(km_fit), "::", 2)[ ,1],
    "\\.",
    2)[ ,2]



  # col_pal <- list()
  # col_pal[["Dark_red"]] <- "#e41a1c" #"#67001f"
  # col_pal[["Light_red"]] <- "#fb9a99" #"#e41a1c"
  # col_pal[["Dark_blue"]] <- "#377eb8" #"#053061"
  # col_pal[["Light_blue"]] <- "#a6cee3" #"#377eb8"




  km_fig <- ggsurvplot(
    fit = km_fit,
    pval = T,
    # pval.method = T,
    risk.table = T,
    newpage = F,
    # palette  = c(col_pal$Dark_red, col_pal$Light_red,
    #              col_pal$Dark_blue, col_pal$Light_blue),
    palette  = c(
      # "#000000", #"black",
      # "#404040", #"gray"
      # "#ca0020" # red
      # "#bdbdbd" #"gray"
      "#e41a1c", #"#000000", #"black" DTX - high
      "#f4a582", #"red" DTX - low
      "#377eb8", # blue NVB - High
      "#92c5de" #"#bdbdbd" #"gray" NVB - low
    ),
    xlab = "Years",
    ylab = str_c(prefix, " probability"),
    break.time.by = 1,
    xlim = c(0,6),
    legend.title = "",

    risk.table.title = element_blank(),
    risk.table.fontsize = 3,
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = T,     # show bars instead of names in text annotations
    # in legend of risk table
    risk.tables.col = "Score",

    ggtheme = theme_bw() +
      theme(panel.grid = element_blank(),
            plot.margin = margin(l=1, r=1, unit = "lines"),
            panel.spacing.y = unit(0,"lines"),
            axis.title.x.bottom = element_text(margin = margin()),
            # axis.title.y.left = element_text(vjust = -4),
            axis.title.y.left = element_text(vjust = 0.5),
            plot.title = element_text(vjust = 0, hjust = .5, face = "bold"),
            legend.direction = "vertical"),

    tables.theme = theme_cleantable() +
      theme(
        axis.title.x.bottom = element_blank()
        # axis.title.x.bottom = element_text(color = "White")
      )
  )


  km_fig_legend <- get_legend(km_fig[[1]])
  # grid.draw(km_fig_legend)

  # Removing legends
  km_fig <- purrr::map(km_fig,
                       function(x){
                         x + guides(color = "none")
                       })

  # arrange_ggsurvplots() only support arranging by column.
  # Rather than setting factor levels it easy to rearrange the list

  km_fig <- km_fig[
    c( # page 1
      "Pooled_Immune", "Pooled_Fibrosis",
      "Pooled_Interferon", "Pooled_Cholesterol",

      # page 2
      "Hamy2016_Immune", "Gruosso2019_Interferon",
      "Yang2018_Immune","Farmer2009_MX1",

      #page 3
      "Hamy2016_Interferon", "Hamy2016_Ecm",
      "Nirmal2018_Interferon", "Naba2014_Ecmcore",

      # page 4
      "Triulzi2013_Ecm", "Pooled_Immune",
      "Sorrentino2014_Chol", "Pooled_Immune"
    )
  ]

  zz <- arrange_ggsurvplots(
    x = km_fig, print = F,
    ncol = 2,
    nrow = 2,
    surv.plot.height = .75,
    risk.table.height = .25,
    newpage = T,
    byrow =T
  )
  pdf(file = paste0(out_figures, "finher_km_", prefix, "_supp_v3.pdf"),
      width = 7.5, height = 8)
  print(zz)
  dev.off()
}


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

