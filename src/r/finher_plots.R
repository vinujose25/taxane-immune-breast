# finher plots


load("results/data/finher_tn.RData")
load("results/data/finher_her2.RData")
load("results/data/finher_her2_split.RData")

load("results/data/clin_finher.RData") # for sig correlations
# load("results/data/expr_finher.RData") # for TILsig computation
# load("results/data/tilsig_clean.RData") # for TILsig computation
#
# # Update clin_finher with TILsig score
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sig <- tilsig_clean$ALL %>%
#   dplyr::select(Direction, Ncbi_gene_id1)
#
# # Compute module score and update clin_finher
# score <- get_module_score_2(
#   x = expr_finher %>% t_tibble(names_x_desc = "CEL_filename"),
#   module_list = list(TILsig_imm = sig %>%
#                        dplyr::filter(Direction == 1) %>%
#                        dplyr::mutate(Direction = 1), # average imm
#                      TILsig_fib = sig %>%
#                        dplyr::filter(Direction == -1) %>%
#                        dplyr::mutate(Direction = 1), # average fib
#                      TILsig = sig), # weighted average
#   by = "Ncbi_gene_id1"
# ) %>%
#   dplyr::mutate(TILsig_imm = TILsig_imm %>% genefu::rescale(q=0.05),
#                 TILsig_fib = TILsig_fib %>% genefu::rescale(q=0.05),
#                 TILsig = TILsig %>% genefu::rescale(q=0.05))
#
# clin_finher <- clin_finher %>% left_join(score, by = "CEL_filename")




# Finher-TN prognosis and chemo teraction
# ==============================================================================

plot_finher_tn <- function(x, event_type,
                          inner_text_size, outer_text_size, color = NULL,
                          prefix, width, height){

  # function specific for plotting finher_tn results


  # # # # Test data
  # x <- finher_tn
  # event_type <- "OS"
  # color = NULL
  # inner_text_size = 2.75
  # outer_text_size = 10

  if(is.null(color)){
    color <- c(DTX = "#e41a1c", # red
               NVB = "#377eb8") # blue
    # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
  }


  # TN prognosis plot
  x_prog <- x[["prog"]][[event_type]]$main_effect %>%
    dplyr::mutate(
      HR = HR %>% round(digits=2),
      Low95 = Low95 %>% round(digits=2),
      Up95 = Up95 %>% round(digits=2),
      HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
      P_text = str_c(
        if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
        "/",
        if_else(is.na(P_adj),
                "NA",
                if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
        )),
      # Het_text = "0.0(.960)-552.5", # to test text size

      Variable = str_replace(Variable, "_scaled","") %>%
        str_replace("_Fc","") %>%
        str_replace("Innate","Inn") %>%
        str_replace("Adhesion","Adh.") %>%
        str_replace("Immune","Imm") %>%
        str_replace("Interferon","Ifn") %>%
        str_replace("Cholesterol","Chl") %>%
        str_replace("Fibrosis","Fib") %>%
        str_replace("Proliferation","Prolif") %>%
        str_replace("Fibroblast","F.blast") %>%
        str_replace("_","."),

      Variable_class = case_when(
        str_detect(Variable, "TIL") ~ "TIL",
        TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
      ) %>%
        factor(levels = c(
          "TIL",
          "Celltype / Proliferation / TIL-localization signatures"
        )),

      Variable = factor(Variable, levels = c(
        # "TIL", "TILsig", "TILsig_imm","TILsig_fib",
        # "Immune1", "Immune2", "Immune3",
        # "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1", "Cholesterol2", "Cholesterol3",
        # "Fibrosis1", "Fibrosis2", "Fibrosis3",
        # "Proliferation1", "Proliferation2", "Proliferation3",
        # "Tcell", "CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        # "Fibroblast"
        # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
        "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
        "TILsig.ECM", "TILsig.Adh.",
        "Imm1", "Imm2", "Imm3",
        "Ifn1", "Ifn2", "Ifn3",
        "Chl1", "Chl2", "Chl3",
        "Fib1", "Fib2", "Fib3",
        "Prolif1", "Prolif2", "Prolif3",
        "Tcell", "CLymph", "Bcell",
        "NKcell", "Monocyte", "MDendri.",
        "F.blast"
        ) %>%
          rev())
    ) %>%
    dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class)
  # dplyr::select(HR, Low95, Up95, HR_text, P_text, Het_text, Variable, Variable_class)

  p_prog <- x_prog %>%
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

  p_prog_annot <- x_prog %>%
    # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
    tidyr::gather("key", "value", HR_text, P_text) %>%
    dplyr::mutate(key = purrr::map_chr(key,
                                       ~switch(.x,
                                               HR_text = "HR(95%CI)",
                                               P_text = "P/Padj*",
                                               Het_text = "Q(P)/I^2")) %>%
                    factor(levels = c("HR(95%CI)","P/Padj*", "Q(P)/I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Variable, label = value), size = inner_text_size) +
    # geom_vline(xintercept = 1.6, linetype = "solid", color = "gray") +
    # geom_vline(xintercept = 2.4, linetype = "solid", color = "gray") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  # Common theme
  th_finher <- theme(
    text = element_text(size = outer_text_size),
    axis.text.x.bottom = element_text(color = "black"),
    axis.title.y = element_blank(),
    # strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    # to adjust space between facet panels
    panel.spacing.y = unit(x = 2, units = "pt")
  )



  # TN chemo interaction plot
  x_chemo <- x[["inter.chemo"]][[event_type]]$main_effect %>%
    dplyr::mutate(
      HR = HR %>% round(digits=2),
      Low95 = exp(Low95) %>% round(digits=2),
      Up95 = exp(Up95) %>% round(digits=2),
      HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
      P_text = str_c(
        if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
        "/",
        if_else(is.na(P_adj),
                "NA",
                if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
        )),
      P_text = if_else(is.na(P_text), "", P_text),
      # Het_text = "0.0(.960)-552.5", # to test text size

      Variable = str_replace(Variable, "_scaled","") %>%
        str_replace("_Fc","") %>%
        str_replace("Innate","Inn") %>%
        str_replace("Adhesion","Adh.") %>%
        str_replace("Immune","Imm") %>%
        str_replace("Interferon","Ifn") %>%
        str_replace("Cholesterol","Chl") %>%
        str_replace("Fibrosis","Fib") %>%
        str_replace("Proliferation","Prolif") %>%
        str_replace("Fibroblast","F.blast") %>%
        str_replace("_","."),

      Variable_class = case_when(
        str_detect(Variable, "TIL") ~ "TIL",
        TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
      ) %>%
        factor(levels = c(
          "TIL",
          "Celltype / Proliferation / TIL-localization signatures"
        )),

      Variable = factor(Variable, levels = c(
        # "TIL", "TILsig", "TILsig_imm","TILsig_fib",
        # "Immune1", "Immune2", "Immune3",
        # "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1", "Cholesterol2", "Cholesterol3",
        # "Fibrosis1", "Fibrosis2", "Fibrosis3",
        # "Proliferation1", "Proliferation2", "Proliferation3",
        # "Tcell", "CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        # "Fibroblast"
        # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
        "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
        "TILsig.ECM", "TILsig.Adh.",
        "Imm1", "Imm2", "Imm3",
        "Ifn1", "Ifn2", "Ifn3",
        "Chl1", "Chl2", "Chl3",
        "Fib1", "Fib2", "Fib3",
        "Prolif1", "Prolif2", "Prolif3",
        "Tcell", "CLymph", "Bcell",
        "NKcell", "Monocyte", "MDendri.",
        "F.blast"
      ) %>%
        rev()),
      Therapy_class = str_split_fixed(Therapy, ":", n = 2)[ , 2],
      Therapy_class = str_replace(Therapy_class, "Chemo", "") %>%
        factor(levels = c("DTX","NVB"))
    ) %>%
    dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class, Therapy_class)
  # dplyr::select(HR, Low95, Up95, HR_text, P_text, Het_text, Variable, Variable_class, Therapy_class)



  p_chemo <- x_chemo %>%
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

  p_chemo_annot <- x_chemo %>%
    tidyr::gather("key", "value", HR_text, P_text) %>%
    # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
    dplyr::mutate(key = purrr::map_chr(key,
                                       ~switch(.x,
                                               HR_text = "HR(95%CI)",
                                               P_text = "P/Padj*",
                                               Het_text= "Q(P)/I^2")) %>%
                    factor(levels = c("HR(95%CI)","P/Padj*","Q(P)/I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Variable, label = value, color = Therapy_class, group = Therapy_class),
              position = ggstance::position_dodgev(height = .75),
              size = inner_text_size) +
    scale_color_manual(values = color) +
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
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  p <- ggpubr::ggarrange(
    p_prog + th_finher + labs(x = "Log-HR"),
    p_prog_annot + th_finher + labs(x = ""),

    p_chemo + th_finher + labs(x = "Log-HR") + theme(axis.text.y = element_blank()),
    p_chemo_annot + th_finher + labs(x = ""),

    ncol = 4,
    nrow = 1,
    widths = c(.16,.375,.09,.375),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix, "_", event_type, ".pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()
}


plot_finher_tn(x = finher_tn,
              event_type = "DDFS",
              inner_text_size = 2.75,
              outer_text_size = 10,
              color =  c(DTX = "#e41a1c", # red
                         NVB = "#377eb8"), # blue,
              prefix = str_c(out_figures,"finher_tn"),
              width = 7.5,
              height = 7)
plot_finher_tn(x = finher_tn,
              event_type = "OS",
              inner_text_size = 2.75,
              outer_text_size = 10,
              color =  c(DTX = "#e41a1c", # red
                         NVB = "#377eb8"), # blue,
              prefix = str_c(out_figures,"finher_tn"),
              width = 7.5,
              height = 7)


#
# ==============================================================================




# Finher-HER2 prognosis and chemo/tra interaction
# ==============================================================================

plot_finher_her2 <- function(x, xsplit, event_type,
                            inner_text_size, outer_text_size,
                            color_chemo = NULL, color_tra = NULL,
                            prefix, width, height){

  # function specific for plotting finher_her2 results


  # # Test data
  # x = finher_her2
  # xsplit = finher_her2_split
  # event_type = "DDFS"
  # inner_text_size = 2.75
  # outer_text_size = 7.5
  # color_chemo = NULL
  # color_tra = NULL

  if(is.null(color_chemo)){
    color_chemo <- c(DTX = "#e41a1c", # red
                     NVB = "#377eb8") # blue
    # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
  }

  if(is.null(color_tra)){
    color_tra <- c(TRA = "#a65628", # brown
                   noTRA = "#984ea3") # purple
    # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
  }





  # HER2 prognosis plot
  # >>>>>>>>>>>>>>>>>>>

  x_prog <- x[["prog"]][[event_type]][["HER2"]]$main_effect

  # Updating only P from survsplit models, if ph of full model failed
  # The HR and CI from full model is used, as it can be interpreted as an average
  # HR for entire followup period.
  x_prog2 <- xsplit[["prog"]][[event_type]][["HER2"]]$main_effect1
  x_prog$P <- purrr::pmap_dbl(list(x_prog$P, x_prog2$P, x_prog$Ph_ok, x_prog2$Ph_ok),
                              ~(if_else(..3,..1,if_else(..4,..2,100))))
  x_prog$P_adj[!is.na(x_prog$P_adj)] <- p.adjust(p = x_prog$P[!is.na(x_prog$P_adj)],
                                                 method = "BH")

  # if( any(x_prog$Ph_ok == F) ){
  #   print(str_c("Using 3yr splitted survival analysis for prog-", event_type,"-HER2-main_effect1"))
  #   x_prog <- xsplit[["prog"]][[event_type]][["HER2"]]$main_effect1
  # }
  #
  #
  # if( any(x_prog$Ph_ok == F) ){
  #   stop("Error PH failed in prog after survSplit.")
  # }

  x_prog <- x_prog %>%
    dplyr::mutate(
      HR = HR %>% round(digits=2),
      Low95 = Low95 %>% round(digits=2),
      Up95 = Up95 %>% round(digits=2),

      HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
      # for PH failed models HR can be interpreted as average HR over followup period
      HR_text = if_else(Ph_ok, HR_text, str_c(HR_text,"~")),

      P_text = str_c(
        if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
        "/",
        if_else(is.na(P_adj),
                "NA",
                if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
        )),
      # for PH failed models P is extracted from the survSplit model
      P_text = if_else(Ph_ok, P_text, str_c(P_text,"#")),
      P_text = if_else(is.na(P_text), "", P_text),

      # Het_text = "0.0(.960)-552.5", # to test text size

      Variable = str_replace(Variable, "_scaled","") %>%
        str_replace("_Fc","") %>%
        str_replace("Innate","Inn") %>%
        str_replace("Adhesion","Adh.") %>%
        str_replace("Immune","Imm") %>%
        str_replace("Interferon","Ifn") %>%
        str_replace("Cholesterol","Chl") %>%
        str_replace("Fibrosis","Fib") %>%
        str_replace("Proliferation","Prolif") %>%
        str_replace("Fibroblast","F.blast") %>%
        str_replace("_","."),

      Variable_class = case_when(
        str_detect(Variable, "TIL") ~ "TIL",
        TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
      ) %>%
        factor(levels = c(
          "TIL",
          "Celltype / Proliferation / TIL-localization signatures"
        )),

      Variable = factor(Variable, levels = c(
        # "TIL", "TILsig", "TILsig_imm","TILsig_fib",
        # "Immune1", "Immune2", "Immune3",
        # "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1", "Cholesterol2", "Cholesterol3",
        # "Fibrosis1", "Fibrosis2", "Fibrosis3",
        # "Proliferation1", "Proliferation2", "Proliferation3",
        # "Tcell", "CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        # "Fibroblast"
        # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
        "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
        "TILsig.ECM", "TILsig.Adh.",
        "Imm1", "Imm2", "Imm3",
        "Ifn1", "Ifn2", "Ifn3",
        "Chl1", "Chl2", "Chl3",
        "Fib1", "Fib2", "Fib3",
        "Prolif1", "Prolif2", "Prolif3",
        "Tcell", "CLymph", "Bcell",
        "NKcell", "Monocyte", "MDendri.",
        "F.blast"
      ) %>%
        rev())
    ) %>%
    dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class)
  # dplyr::select(HR, Low95, Up95, HR_text, P_text, Het_text, Variable, Variable_class)


  p_prog <- x_prog %>%
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
      # axis.text.y.left = element_text(angle = 60),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )

  p_prog_annot <- x_prog %>%
    # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
    tidyr::gather("key", "value", HR_text, P_text) %>%
    dplyr::mutate(key = purrr::map_chr(key,
                                       ~switch(.x,
                                               HR_text = "HR(95%CI)",
                                               P_text = "P/Padj*",
                                               Het_text = "Q(P)/I^2")) %>%
                    factor(levels = c("HR(95%CI)","P/Padj*", "Q(P)/I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Variable, label = value), size = inner_text_size) +
    # geom_vline(xintercept = 1.6, linetype = "solid", color = "gray") +
    # geom_vline(xintercept = 2.4, linetype = "solid", color = "gray") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )





  # Common theme
  # >>>>>>>>>>>>

  th_finher <- theme(
    text = element_text(size = outer_text_size),
    axis.text.x.bottom = element_text(color = "black"),
    axis.title.y = element_blank(),
    # strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    # to adjust space between facet panels
    panel.spacing.y = unit(x = 2, units = "pt")
  )




  # HER2 chemo interaction plot
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>

  x_chemo <- x[["inter.chemo"]][[event_type]][["HER2"]]$main_effect

  # Updating only P from survsplit models, if ph of full model failed
  # The HR and CI from full model is used, as it can be interpreted as an average
  # HR for entire followup period.
  x_chemo2 <- xsplit[["inter.chemo"]][[event_type]][["HER2"]]$main_effect1

  # Setting interaction pvale on DTX row (to make NAs congruent b/w x_chemo and x_chemo2)
  tmp_p <- rep(NA,nrow(x_chemo2))
  tmp_p[which(is.na(x_chemo2$P))] <- x_chemo2$P[which(is.na(x_chemo2$P))+1]
  x_chemo2$P <- tmp_p

  x_chemo$P <- purrr::pmap_dbl(list(x_chemo$P, x_chemo2$P, x_chemo$Ph_ok, x_chemo2$Ph_ok),
                              ~(if_else(..3,..1,if_else(..4,..2,100))))
  x_chemo$P_adj[!is.na(x_chemo$P_adj)] <- p.adjust(p = x_chemo$P[!is.na(x_chemo$P_adj)],
                                                 method = "BH")

  # if( any(x_chemo$Ph_ok == F) ){
  #   print(str_c("Using 3yr splitted survival analysis for intr.chemo-", event_type,"-HER2-main_effect1"))
  #   x_chemo <- xsplit[["inter.chemo"]][[event_type]][["HER2"]]$main_effect1
  # }
  #
  # if( any(x_chemo$Ph_ok == F) ){
  #  stop("Error PH failed in inter.chemo after survSplit.")
  # }

  x_chemo <- x_chemo %>%
    dplyr::mutate(
      HR = HR %>% round(digits=2),
      Low95 = exp(Low95) %>% round(digits=2),
      Up95 = exp(Up95) %>% round(digits=2),

      HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
      # for PH failed models HR can be interpreted as average HR over followup period
      HR_text = if_else(Ph_ok, HR_text, str_c(HR_text,"~")),

      P_text = str_c(
        if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
        "/",
        if_else(is.na(P_adj),
                "NA",
                if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
        )),
      # for PH failed models P is extracted from the survSplit model
      P_text = if_else(Ph_ok, P_text, str_c(P_text,"#")),
      P_text = if_else(is.na(P_text), "", P_text),

      # Het_text = "0.0(.960)-552.5", # to test text size

      Variable = str_replace(Variable, "_scaled","") %>%
        str_replace("_Fc","") %>%
        str_replace("Innate","Inn") %>%
        str_replace("Adhesion","Adh.") %>%
        str_replace("Immune","Imm") %>%
        str_replace("Interferon","Ifn") %>%
        str_replace("Cholesterol","Chl") %>%
        str_replace("Fibrosis","Fib") %>%
        str_replace("Proliferation","Prolif") %>%
        str_replace("Fibroblast","F.blast") %>%
        str_replace("_","."),

      Variable_class = case_when(
        str_detect(Variable, "TIL") ~ "TIL",
        TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
      ) %>%
        factor(levels = c(
          "TIL",
          "Celltype / Proliferation / TIL-localization signatures"
        )),

      Variable = factor(Variable, levels = c(
        # "TIL", "TILsig", "TILsig_imm","TILsig_fib",
        # "Immune1", "Immune2", "Immune3",
        # "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1", "Cholesterol2", "Cholesterol3",
        # "Fibrosis1", "Fibrosis2", "Fibrosis3",
        # "Proliferation1", "Proliferation2", "Proliferation3",
        # "Tcell", "CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        # "Fibroblast"
        # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
        "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
        "TILsig.ECM", "TILsig.Adh.",
        "Imm1", "Imm2", "Imm3",
        "Ifn1", "Ifn2", "Ifn3",
        "Chl1", "Chl2", "Chl3",
        "Fib1", "Fib2", "Fib3",
        "Prolif1", "Prolif2", "Prolif3",
        "Tcell", "CLymph", "Bcell",
        "NKcell", "Monocyte", "MDendri.",
        "F.blast"
      ) %>%
        rev()),
      Therapy = str_replace(Therapy, ":Interval_Id1",""), # considering only main_effect1
      Therapy_class = str_split_fixed(Therapy, ":", n = 2)[ , 2],
      Therapy_class = str_replace(Therapy_class, "Chemo", "") %>%
        factor(levels = c("DTX","NVB"))
    ) %>%
    dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class, Therapy_class)
  # dplyr::select(HR, Low95, Up95, HR_text, P_text, Het_text, Variable, Variable_class, Therapy_class)



  p_chemo <- x_chemo %>%
    ggplot() +
    geom_point(aes(x = log(HR), y = Variable, color = Therapy_class, group = Therapy_class),
               position = ggstance::position_dodgev(height = .75)) +
    geom_errorbarh(aes(xmin = log(Low95), xmax = log(Up95), y = Variable, color = Therapy_class, group = Therapy_class),
                   position = ggstance::position_dodgev(height = .75),
                   height = .1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_color_manual(values = color_chemo) +
    guides(color = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      # axis.text.y.left = element_text(angle = 60),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  p_chemo_annot <- x_chemo %>%
    tidyr::gather("key", "value", HR_text, P_text) %>%
    # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
    dplyr::mutate(key = purrr::map_chr(key,
                                       ~switch(.x,
                                               HR_text = "HR(95%CI)",
                                               P_text = "P/Padj*",
                                               Het_text= "Q(P)/I^2")) %>%
                    factor(levels = c("HR(95%CI)","P/Padj*","Q(P)/I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Variable, label = value, color = Therapy_class, group = Therapy_class),
              position = ggstance::position_dodgev(height = .75),
              size = inner_text_size) +
    scale_color_manual(values = color_chemo) +
    guides(color = "none") +
    # geom_vline(xintercept = 1.6, linetype = "solid", color = "gray") +
    # geom_vline(xintercept = 2.4, linetype = "solid", color = "gray") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )




  # HER2 tra interaction plot
  # >>>>>>>>>>>>>>>>>>>>>>>>>

  x_tra <- x[["inter.tra"]][[event_type]][["HER2"]]$main_effect


  # Updating only P from survsplit models, if ph of full model failed
  # The HR and CI from full model is used, as it can be interpreted as an average
  # HR for entire followup period.
  x_tra2 <- xsplit[["inter.tra"]][[event_type]][["HER2"]]$main_effect1

  # Setting interaction pvale on noTRA row (to make NAs congruent b/w x_tra and x_tra2)
  tmp_p <- rep(NA,nrow(x_tra2))
  tmp_p[which(is.na(x_tra2$P))] <- x_tra2$P[which(is.na(x_tra2$P))+1]
  x_tra2$P <- tmp_p

  x_tra$P <- purrr::pmap_dbl(list(x_tra$P, x_tra2$P, x_tra$Ph_ok, x_tra2$Ph_ok),
                               ~(if_else(..3,..1,if_else(..4,..2,100))))
  x_tra$P_adj[!is.na(x_tra$P_adj)] <- p.adjust(p = x_tra$P[!is.na(x_tra$P_adj)],
                                                   method = "BH")

  # # test
  # tibble(x=x_tra$P,y=x_tra2$P,p=x_tra$Ph_ok,q=x_tra2$Ph_ok) %>% as.data.frame()

  # if( any(x_tra$Ph_ok == F) ){
  #   print(str_c("Using 3yr splitted survival analysis for intr.tra-", event_type,"-HER2-main_effect1"))
  #   x_tra <- xsplit[["inter.tra"]][[event_type]][["HER2"]]$main_effect1
  # }
  #
  # if( any(x_tra$Ph_ok == F) ){
  #   stop("Error PH failed in inter.tra after survSplit.")
  # }

  x_tra <- x_tra %>%
    dplyr::mutate(
      HR = HR %>% round(digits=2),
      Low95 = exp(Low95) %>% round(digits=2),
      Up95 = exp(Up95) %>% round(digits=2),

      HR_text = str_c(HR,"(",Low95," to ",Up95,")"),
      # for PH failed models HR can be interpreted as average HR over followup period
      HR_text = if_else(Ph_ok, HR_text, str_c(HR_text,"~")),

      P_text = str_c(
        if_else(P < 0.001,"<.001", round(P, digits=3) %>% as.character()),
        "/",
        if_else(is.na(P_adj),
                "NA",
                if_else(P_adj < 0.001,"<.001", round(P_adj,digits=3) %>% as.character())
        )),
      # for PH failed models P is extracted from the survSplit model
      P_text = if_else(Ph_ok, P_text, str_c(P_text,"#")),
      P_text = if_else(is.na(P_text), "", P_text),


      # Het_text = "0.0(.960)-552.5", # to test text size

      Variable = str_replace(Variable, "_scaled","") %>%
        str_replace("_Fc","") %>%
        str_replace("Innate","Inn") %>%
        str_replace("Adhesion","Adh.") %>%
        str_replace("Immune","Imm") %>%
        str_replace("Interferon","Ifn") %>%
        str_replace("Cholesterol","Chl") %>%
        str_replace("Fibrosis","Fib") %>%
        str_replace("Proliferation","Prolif") %>%
        str_replace("Fibroblast","F.blast") %>%
        str_replace("_","."),

      Variable_class = case_when(
        str_detect(Variable, "TIL") ~ "TIL",
        TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
      ) %>%
        factor(levels = c(
          "TIL",
          "Celltype / Proliferation / TIL-localization signatures"
        )),

      Variable = factor(Variable, levels = c(
        # "TIL", "TILsig", "TILsig_imm","TILsig_fib",
        # "Immune1", "Immune2", "Immune3",
        # "Interferon1", "Interferon2", "Interferon3",
        # "Cholesterol1", "Cholesterol2", "Cholesterol3",
        # "Fibrosis1", "Fibrosis2", "Fibrosis3",
        # "Proliferation1", "Proliferation2", "Proliferation3",
        # "Tcell", "CLymphocyte", "Bcell",
        # "NKcell", "Monocyte", "MDendritic",
        # "Fibroblast"
        # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
        "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
        "TILsig.ECM", "TILsig.Adh.",
        "Imm1", "Imm2", "Imm3",
        "Ifn1", "Ifn2", "Ifn3",
        "Chl1", "Chl2", "Chl3",
        "Fib1", "Fib2", "Fib3",
        "Prolif1", "Prolif2", "Prolif3",
        "Tcell", "CLymph", "Bcell",
        "NKcell", "Monocyte", "MDendri.",
        "F.blast"
      ) %>%
        rev()),
      Therapy = str_replace(Therapy, ":Interval_Id1",""), # considering only main_effect1
      Therapy_class = str_split_fixed(Therapy, ":", n = 2)[ , 2],
      Therapy_class = str_replace(Therapy_class, "Herceptin", "") %>%
        factor(levels = c("TRA","noTRA"))
    ) %>%
    dplyr::select(HR, Low95, Up95, HR_text, P_text, Variable, Variable_class, Therapy_class)
  # dplyr::select(HR, Low95, Up95, HR_text, P_text, Het_text, Variable, Variable_class, Therapy_class)


  p_tra <- x_tra %>%
    ggplot() +
    geom_point(aes(x = log(HR), y = Variable, color = Therapy_class, group = Therapy_class),
               position = ggstance::position_dodgev(height = .75)) +
    geom_errorbarh(aes(xmin = log(Low95), xmax = log(Up95), y = Variable, color = Therapy_class, group = Therapy_class),
                   position = ggstance::position_dodgev(height = .75),
                   height = .1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_color_manual(values = color_tra) +
    guides(color = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      # axis.text.y.left = element_text(angle = 60),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  p_tra_annot <- x_tra %>%
    tidyr::gather("key", "value", HR_text, P_text) %>%
    # tidyr::gather("key", "value", HR_text, P_text, Het_text) %>%
    dplyr::mutate(key = purrr::map_chr(key,
                                       ~switch(.x,
                                               HR_text = "HR(95%CI)",
                                               P_text = "P/Padj*",
                                               Het_text= "Q(P)/I^2")) %>%
                    factor(levels = c("HR(95%CI)","P/Padj*","Q(P)/I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Variable, label = value, color = Therapy_class, group = Therapy_class),
              position = ggstance::position_dodgev(height = .75),
              size = inner_text_size) +
    scale_color_manual(values = color_tra) +
    guides(color = "none") +
    # geom_vline(xintercept = 1.6, linetype = "solid", color = "gray") +
    # geom_vline(xintercept = 2.4, linetype = "solid", color = "gray") +
    facet_grid(facets = Variable_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  # interaction plots
  # >>>>>>>>>>>>>>>>>

  p <- ggpubr::ggarrange(

    p_tra + th_finher + labs(x = "Log-HR"),
    p_tra_annot + th_finher + labs(x = ""),

    p_chemo + th_finher + labs(x = "Log-HR") + theme(axis.text.y = element_blank()),
    p_chemo_annot + th_finher + labs(x = ""),

    ncol = 4,
    nrow = 1,
    widths = c(.16,.375,.09,.375),
    # widths = c(.18,.36,.1,.36),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix,"_inter_",event_type,".pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()



  # Prognosis plot
  # >>>>>>>>>>>>>>

  p <- ggpubr::ggarrange(

    p_prog + th_finher + labs(x = "Log-HR"),
    p_prog_annot + th_finher + labs(x = ""),

    ncol = 2,
    nrow = 1,
    widths = c(.35,.65)
  )

  pdf(file = str_c(prefix,"_prog_",event_type,".pdf") %>% str_to_lower(),
      width = width/2, height = height)
  print(p)
  dev.off()

}


plot_finher_her2(x = finher_her2,
                xsplit = finher_her2_split,
                event_type = "DDFS",
                inner_text_size = 2.75,
                outer_text_size = 10,
                color_chemo =  c(DTX = "#e41a1c", # red
                                 NVB = "#377eb8"), # blue,
                color_tra =  c(TRA = "#a65628", # brown
                               noTRA = "#984ea3"), # purple
                prefix = str_c(out_figures,"finher_her2"),
                width = 7.5,
                height = 7)

plot_finher_her2(x = finher_her2,
                xsplit = finher_her2_split,
                event_type = "OS",
                inner_text_size = 2.75,
                outer_text_size = 10,
                color_chemo =  c(DTX = "#e41a1c", # red
                                 NVB = "#377eb8"), # blue,
                color_tra =  c(TRA = "#a65628", # brown
                               noTRA = "#984ea3"), # purple
                prefix = str_c(out_figures,"finher_her2"),
                width = 7.5,
                height = 7)


#
# ==============================================================================




# FinHER TILsig-TIL localization sig correlation plots
# ==============================================================================

nme <- c("Immune1", "Immune2", "Immune3",
         "Interferon1", "Interferon2", "Interferon3",
         "Cholesterol2",
         "Fibrosis1", "Fibrosis2", "Fibrosis3",
         "Proliferation1", "Proliferation2", "Proliferation3",
         "Tcell","Fibroblast",
         # "TIL", "TILsig", "TILsig_imm", "TILsig_fib")
         "TIL", "TILsig_scaled", "TILsig_APP_Fc", "TILsig_Immune",
         "TILsig_IFNg", "TILsig_Innate",
         "TILsig_ECM", "TILsig_Adhesion")

plot_corr <- purrr::map(

  c("ALL","HER2","TN"),

  function(subtype, clin, nme){

    # subtype = "TN"
    # clin <- clin_finher

    xclin <- switch(subtype,
                    "ALL" = clin,
                    "TN" = clin %>%
                      dplyr::filter(Subtype_IHC_2 == "TN"),
                    "HER2" = clin %>%
                      dplyr::filter(Subtype_IHC_2 == "HER2"))

    xcor <- cor(
      xclin %>%
        dplyr::filter(!is.na(StrLy_Mean)) %>%
        dplyr::rename(TIL = "StrLy_Mean") %>%
        dplyr::select(all_of(nme)),
      method = "spearman")

    as_tibble(xcor, rownames = "Rows") %>%
      tidyr::gather(key = "Columns", value = "Correlation",all_of(nme)) %>%
      dplyr::mutate(Subtype = subtype)
  },
  clin = clin_finher,
  nme
)

plot_corr <- bind_rows(plot_corr) %>%
  dplyr::mutate(
    Rows = factor(Rows, levels =  nme %>% rev()),
    Columns = factor(Columns, levels =  nme %>% rev()),
    Subtype = factor(Subtype, levels =  c("ALL","HER2","TN")),
    Cor_text = Correlation %>% round(digits = 1) %>% str_replace("0\\.", ".")
  )



p <- plot_corr %>%
  ggplot(aes(x = Columns, y= Rows)) +
  geom_raster(aes(fill = Correlation)) +
  geom_text(aes(label = Cor_text), size = 2) +
  scale_fill_gradient2(high = "#b2182b",# red
                       low = "#2166ac",# blue
                       limits = c(-1,1)) +
  guides(fill = guide_colorbar(title = "Spearman\ncorrelation", title.vjust = 1)) +
  facet_wrap(facets = ~Subtype , nrow = 2) +
  theme(legend.position = "bottom",
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust=.5),
        axis.title = element_blank())

pdf(file = "results/figures/finher_correlations.pdf",
    width = 7.5, height = 8) # 4.5)
print(p)
dev.off()



#
# ==============================================================================
