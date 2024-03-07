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
  pdf(file = paste0(out_figures, "finher_km_", prefix, "-MTA-immune_interaction.pdf"),
      width = 7.5, height = 8)
  # png(filename = paste0(out_figures, "finher_km_", prefix, ".png"),
  #     width = 7.5, height = 7, units = "in", res = 300)
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
  pdf(file = paste0(out_figures, "finher_km_", prefix, "-MTA-immune_interaction_supp.pdf"),
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




# Finher MTA interaction
# ==============================================================================

plot_finher_mta_interaction <- function(
    tn_cox_sum,
    her2_cox_sum,
    event_type,
    inner_text_size,
    outer_text_size,
    color = NULL,
    prefix,
    width,
    height){


  # Function for plotting MTA interaction plot (A:TN, B:HER2)


  #   # Test data
  #   tn_cox_sum <- finher_tn
  #   her2_cox_sum <- finher_her2
  #   event_type <- "DDFS"
  #   color = NULL
  #   inner_text_size = 2.75
  #   outer_text_size = 10
  #   prefix = str_c(out_figures,"finher_mta_inter")
  #   width = 7.5
  #   height = 8.5

  if(is.null(color)){
    color <- c(NVB = "#377eb8", # blue
               DTX = "#e41a1c") # red

    # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
  }



  # Common theme
  th_finher <- theme(
    text = element_text(size = outer_text_size),
    axis.text.x.bottom = element_text(color = "black"),
    axis.title.y = element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    strip.text.x = element_blank(),
    panel.grid = element_blank(),
    # to adjust space between facet panels
    panel.spacing.y = unit(x = 2, units = "pt")
  )



  cox_sum <- purrr::map(
    list(
      tn = tn_cox_sum[["inter.chemo"]][[event_type]]$main_effect,
      her2 = her2_cox_sum[["inter.chemo"]][[event_type]]$HER2$main_effect
    ),
    function(x){
      x %>%
        dplyr::mutate(
          HR = HR %>% round(digits=2),
          Low95 = exp(Low95) %>% round(digits=2),
          Up95 = exp(Up95) %>% round(digits=2),
          HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
          P_text = str_c(
            if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
            "(",
            if_else(is.na(P_adj),
                    "NA",
                    if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
            ),
            ")"),
          P_text = if_else(is.na(P_text), "", P_text),


          Variable = str_split_fixed(string = Variable, pattern = ":", n = 2)[, 1] %>%
            str_replace("scaled_",""),

          # grouping modules
          Variable_class = case_when(
            Variable == "TIL" ~ "H&E",

            str_detect(Variable, "Denovo") ~ "De-novo",
            str_detect(Variable, "Pooled") ~ "Published-pooled",
            str_detect(Variable, "MCPcounter") ~ "Cell",

            Variable == "Hamy2016_Immune" ~ "Immune",
            Variable == "Yang2018_Immune" ~ "Immune",

            Variable == "Gruosso2019_Interferon" ~ "Interferon",
            Variable == "Farmer2009_MX1" ~ "Interferon",
            Variable == "Hamy2016_Interferon" ~ "Interferon",
            Variable == "Nirmal2018_Interferon" ~ "Interferon",


            Variable == "Hamy2016_Ecm" ~ "Fibrosis",
            Variable == "Naba2014_Ecmcore" ~ "Fibrosis",
            Variable == "Triulzi2013_Ecm" ~ "Fibrosis",

            Variable == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

            TRUE ~ "Unknown"
          ) %>%
            factor(
              levels = c(
                "H&E", "De-novo",
                "Immune", "Interferon", "Chol.", "Fibrosis",
                "Published-pooled", "Cell")
            ),

          Variable =  Variable %>%
            str_replace("Denovo","") %>%
            str_replace("Pooled","") %>%
            str_replace("MCPcounter","") %>%

            str_replace("Hamy2016_Immune","Hamy2016") %>%
            str_replace("Yang2018_Immune","Yang2018") %>%

            str_replace("Gruosso2019_Interferon","Gruosso2019") %>%
            str_replace("Farmer2009_MX1","Farmer2009") %>%
            str_replace("Hamy2016_Interferon","Hamy2016") %>%
            str_replace("Nirmal2018_Interferon","Nirmal2018") %>%

            str_replace("Hamy2016_Ecm","Hamy2016") %>%
            str_replace("Naba2014_Ecmcore","Naba2014") %>%
            str_replace("Triulzi2013_Ecm","Triulzi2013") %>%

            str_replace("Sorrentino2014_Chol","Sorrentino2014") %>%

            str_replace("Tcell","T.Cell") %>%
            str_replace("B.Lineage","B.Cell") %>%
            str_replace("Monocytic.Lineage","Monocyte") %>%
            str_replace("Myeloid.Dendritic","M.Dendritic") %>%
            # Endothelial
            str_replace("Fibroblasts","Fibroblast") %>%

            str_replace("_","") %>%

            factor(levels = c(

              # H&E
              "TIL",

              # denovo: "TILsig", "Immune","ECM",
              # pooled: "Immune", "Interferon", "Cholesterol", "Fibrosis",

              # denovo + pooled
              "TILsig", "Immune","ECM", "Interferon", "Cholesterol", "Fibrosis",

              # immune: "Hamy2016","Yang2018",
              # interferon: "Gruosso2019", "Farmer2009", "Hamy2016", "Nirmal2018",
              # cholesterol:"Sorrentino2014",
              # fibrosis: "Hamy2016", "Naba2014", "Triulzi2013",

              # immune + interferon + cholesterol + fibrosis
              "Gruosso2019", "Farmer2009", "Hamy2016", "Nirmal2018",
              "Yang2018",
              "Naba2014", "Triulzi2013",
              "Sorrentino2014",

              # cell
              "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
              "Monocyte", "M.Dendritic", "Neutrophil",
              "Endothelial", "Fibroblast"

            ) %>% rev()
            ),


          Therapy_class = str_split_fixed(Therapy, ":", n = 2)[ , 2],
          Therapy_class = str_replace(Therapy_class, "Chemo", "") %>%
            factor(levels = c("NVB","DTX"))
        ) %>%
        dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class, Therapy_class, Subtype)
    }
  )

  p_forest <- purrr::map(
    cox_sum,
    function(x){
      x %>%
        ggplot() +
        geom_point(aes(x = log(HR), y = Variable, color = Therapy_class, group = Therapy_class),
                   position = ggstance::position_dodgev(height = .75)) +
        geom_errorbarh(aes(xmin = log(Low95), xmax = log(Up95), y = Variable, color = Therapy_class, group = Therapy_class),
                       position = ggstance::position_dodgev(height = .75),
                       height = .1) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
        scale_color_manual(values = color) +
        guides(color = "none") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_grid(facets = Variable_class ~ 1,
                   switch = "y",
                   scales = "free_y",
                   space = "free_y") +
        theme_bw() +
        theme(
          # to manage space in ggarrage()
          plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
          # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
        )
    })



  p_annot <- purrr::map(
    cox_sum,
    function(x){
      x %>%
        tidyr::gather("key", "value", HR_text, P_text) %>%
        # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
        dplyr::mutate(key = purrr::map_chr(key,
                                           ~switch(.x,
                                                   HR_text = "HR^a(95%CI)",
                                                   # HR_text = expression("HR"^a*"(95%CI) "),
                                                   P_text = "P(Padj^b)",
                                                   Het_text= "Q(P)/I^2")) %>%
                        factor(levels = c("HR^a(95%CI)","P(Padj^b)","Q(P)/I^2"),
                               # labels = c(str_c("HR",supsc('a'),"(95%CI)"),
                               #            str_c("P(Padj",supsc('b'),")"),
                               #            str_c("Q(P)/I",supsc('2')))
                               labels = c(expression("HR"^a*"(95%CI) "),
                                          expression("P(Padj"^b*") "),
                                          expression("Q(P)/I"^2))
                               )
                      ) %>%
        ggplot() +
        geom_text(aes(x = key, y = Variable, label = value, color = Therapy_class, group = Therapy_class),
                  position = ggstance::position_dodgev(height = .75),
                  size = inner_text_size) +
        scale_color_manual(values = color) +
        scale_x_discrete(labels = ~parse(text = .x)) +
        # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
        guides(color = "none") +
        facet_grid(facets = Variable_class ~ 1,
                   switch = "y",
                   scales = "free_y",
                   space = "free_y") +
        theme_bw() +
        theme(
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          # to manage space in ggarrage()
          plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
          # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
        )
    })



  p <- ggpubr::ggarrange(
    p_forest$tn + th_finher + labs(x = expression("Log(HR"^a*") ")) +  # *close superscript, space is for margin
      theme (axis.title.x.bottom = element_text(size = outer_text_size-1)),
    p_annot$tn + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    p_forest$her2 + th_finher + labs(x = expression("Log(HR"^a*") ")) +
      theme(axis.text.y = element_blank(),
            axis.title.x.bottom = element_text(size = outer_text_size-1)),
    p_annot$her2 + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    ncol = 4,
    nrow = 1,
    # widths = c(.16,.375,.09,.375), #old
    # widths = c(.21,.37,.14,.37),
    widths = c(.23,.36,.13,.36),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix, "_", event_type, "-MTA-immune_interaction.pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()
}



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


plot_finher_prognosis <- function(
    tn_cox_sum,
    her2_cox_sum,
    event_type,
    inner_text_size,
    outer_text_size,
    color = NULL,
    prefix,
    width,
    height){


  # Function for prognosis plot (A:TN, B:HER2)


  # # Test data
  # tn_cox_sum <- finher_tn
  # her2_cox_sum <- finher_her2
  # event_type <- "DDFS"
  # color = NULL
  # inner_text_size = 2.75
  # outer_text_size = 10
  # prefix = str_c(out_figures,"finher_prognosis")
  # width = 7.5
  # height = 8.5


  # Common theme
  th_finher <- theme(
    text = element_text(size = outer_text_size),
    axis.text.x.bottom = element_text(color = "black"),
    axis.title.y = element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    strip.text.x = element_blank(),
    panel.grid = element_blank(),
    # to adjust space between facet panels
    panel.spacing.y = unit(x = 2, units = "pt")
  )



  cox_sum <- purrr::map(
    list(
      tn = tn_cox_sum[["prog"]][[event_type]]$main_effect,
      her2 = her2_cox_sum[["prog"]][[event_type]]$HER2$main_effect
    ),
    function(x){
      x %>%
        dplyr::mutate(
          HR = HR %>% round(digits=2),
          Low95 = Low95 %>% round(digits=2),
          Up95 = Up95 %>% round(digits=2),
          HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
          P_text = str_c(
            if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
            "(",
            if_else(is.na(P_adj),
                    "NA",
                    if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
            ),
            ")"),


          Variable = str_split_fixed(string = Variable, pattern = ":", n = 2)[, 1] %>%
            str_replace("scaled_",""),

          # grouping modules
          Variable_class = case_when(
            Variable == "TIL" ~ "H&E",

            str_detect(Variable, "Denovo") ~ "De-novo",
            str_detect(Variable, "Pooled") ~ "Published-pooled",
            str_detect(Variable, "MCPcounter") ~ "Cell",

            Variable == "Hamy2016_Immune" ~ "Immune",
            Variable == "Yang2018_Immune" ~ "Immune",

            Variable == "Gruosso2019_Interferon" ~ "Interferon",
            Variable == "Farmer2009_MX1" ~ "Interferon",
            Variable == "Hamy2016_Interferon" ~ "Interferon",
            Variable == "Nirmal2018_Interferon" ~ "Interferon",


            Variable == "Hamy2016_Ecm" ~ "Fibrosis",
            Variable == "Naba2014_Ecmcore" ~ "Fibrosis",
            Variable == "Triulzi2013_Ecm" ~ "Fibrosis",

            Variable == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

            TRUE ~ "Unknown"
          ) %>%
            factor(
              levels = c(
                "H&E", "De-novo",
                "Immune", "Interferon", "Chol.", "Fibrosis",
                "Published-pooled", "Cell")
            ),

          Variable =  Variable %>%
            str_replace("Denovo","") %>%
            str_replace("Pooled","") %>%
            str_replace("MCPcounter","") %>%

            str_replace("Hamy2016_Immune","Hamy2016") %>%
            str_replace("Yang2018_Immune","Yang2018") %>%

            str_replace("Gruosso2019_Interferon","Gruosso2019") %>%
            str_replace("Farmer2009_MX1","Farmer2009") %>%
            str_replace("Hamy2016_Interferon","Hamy2016") %>%
            str_replace("Nirmal2018_Interferon","Nirmal2018") %>%

            str_replace("Hamy2016_Ecm","Hamy2016") %>%
            str_replace("Naba2014_Ecmcore","Naba2014") %>%
            str_replace("Triulzi2013_Ecm","Triulzi2013") %>%

            str_replace("Sorrentino2014_Chol","Sorrentino2014") %>%

            str_replace("Tcell","T.Cell") %>%
            str_replace("B.Lineage","B.Cell") %>%
            str_replace("Monocytic.Lineage","Monocyte") %>%
            str_replace("Myeloid.Dendritic","M.Dendritic") %>%
            # Endothelial
            str_replace("Fibroblasts","Fibroblast") %>%

            str_replace("_","") %>%

            factor(levels = c(

              # H&E
              "TIL",

              # denovo: "TILsig", "Immune","ECM",
              # pooled: "Immune", "Interferon", "Cholesterol", "Fibrosis",

              # denovo + pooled
              "TILsig", "Immune","ECM", "Interferon", "Cholesterol", "Fibrosis",

              # immune: "Hamy2016","Yang2018",
              # interferon: "Gruosso2019", "Farmer2009", "Hamy2016", "Nirmal2018",
              # cholesterol:"Sorrentino2014",
              # fibrosis: "Hamy2016", "Naba2014", "Triulzi2013",

              # immune + interferon + cholesterol + fibrosis
              "Gruosso2019", "Farmer2009", "Hamy2016", "Nirmal2018",
              "Yang2018",
              "Naba2014", "Triulzi2013",
              "Sorrentino2014",

              # cell
              "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
              "Monocyte", "M.Dendritic", "Neutrophil",
              "Endothelial", "Fibroblast"

            ) %>% rev())



        ) %>%
        dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class, Subtype)

    }
  )


  p_forest <- purrr::map(
    cox_sum,
    function(x){
      x %>%
        ggplot() +
        geom_point(aes(x = log(HR), y = Variable)) +
        geom_errorbarh(aes(xmin = log(Low95), xmax = log(Up95), y = Variable), height = .1) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
        facet_grid(facets = Variable_class ~ 1,
                   switch = "y",
                   scales = "free_y",
                   space = "free_y") +
        theme_bw() +
        theme(
          # to manage space in ggarrage()
          plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
          # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
        )
    })


  p_annot <- purrr::map(
    cox_sum,
    function(x){
      x %>%
        tidyr::gather("key", "value", HR_text, P_text) %>%
        dplyr::mutate(key = purrr::map_chr(
          key,
          ~switch(.x,
                  # HR_text = "HR(95%CI)",
                  # P_text = "P(Padj*)",
                  # Het_text = "Q(P)/I^2"
                  HR_text = "HR^a(95%CI)",
                  P_text = "P(Padj^b)",
                  Het_text= "Q(P)/I^2")
        ) %>%
          factor(
            # levels = c("HR(95%CI)","P(Padj*)", "Q(P)/I^2")
            levels = c("HR^a(95%CI)","P(Padj^b)","Q(P)/I^2"),
            labels = c(expression("HR"^a*"(95%CI) "),
                       expression("P(Padj"^b*") "),
                       expression("Q(P)/I"^2))
          )
        ) %>%
        ggplot() +
        geom_text(aes(x = key, y = Variable, label = value), size = inner_text_size) +
        # geom_vline(xintercept = 1.6, linetype = "solid", color = "gray") +
        # geom_vline(xintercept = 2.4, linetype = "solid", color = "gray") +
        facet_grid(facets = Variable_class ~ 1,
                   switch = "y",
                   scales = "free_y",
                   space = "free_y") +
        scale_x_discrete(labels = ~parse(text = .x)) +
        # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
        theme_bw() +
        theme(
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          # to manage space in ggarrage()
          plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
          # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
        )
    })



  p <- ggpubr::ggarrange(
    p_forest$tn + th_finher + # labs(x = "Log-HR") +
      labs(x = expression("Log(HR"^a*") ")) +  # *close superscript, space is for margin
      theme(axis.title.x.bottom = element_text(size = outer_text_size-1)),
    p_annot$tn + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    p_forest$her2 + th_finher + # labs(x = "Log-HR") +
      labs(x = expression("Log(HR"^a*") ")) +  # *close superscript, space is for margin
      theme(axis.text.y = element_blank(),
            axis.title.x.bottom = element_text(size = outer_text_size-1)),
    p_annot$her2 + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    ncol = 4,
    nrow = 1,
    # widths = c(.16,.375,.09,.375), #old
    # widths = c(.21,.37,.14,.37),
    widths = c(.23,.36,.13,.36),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix, "_", event_type, "_prognosis.pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()
}




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


