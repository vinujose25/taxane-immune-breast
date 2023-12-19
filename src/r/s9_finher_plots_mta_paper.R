# finher_plots_mta_paper.R


load("results/data/finher_tn.RData")
load("results/data/finher_her2.RData")
# For MTA interaction PH is valid for all her2 models. Hence time interaction
# analysis in combination with survsplit() is not needed
# load("results/data/finher_her2_split.RData")
# load("results/data/clin_finher.RData") # for sig correlations and km plots


# Finher MTA interaction sample size
# ==============================================================================
clin_finher  %>%
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
    color <- c(DTX = "#e41a1c", # red
               NVB = "#377eb8") # blue
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


          # Variable = str_replace(Variable, "_scaled","") %>%
          #   str_replace("_Fc","") %>%
          #   str_replace("Innate","Inn") %>%
          #   str_replace("Adhesion","Adh.") %>%
          #   str_replace("Immune","Imm") %>%
          #   str_replace("Interferon","Ifn") %>%
          #   str_replace("Cholesterol","Chl") %>%
          #   str_replace("Fibrosis","Fib") %>%
          #   str_replace("Proliferation","Prolif") %>%
          #   str_replace("Fibroblast","F.blast") %>%
          #   str_replace("_","."),
          #
          # Variable_class = case_when(
          #   str_detect(Variable, "TIL") ~ "TIL",
          #   TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
          # ) %>%
          #   factor(levels = c(
          #     "TIL",
          #     "Celltype / Proliferation / TIL-localization signatures"
          #   )),
          #
          # Variable = factor(Variable, levels = c(
          #   "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
          #   "TILsig.ECM", "TILsig.Adh.",
          #   "Imm1", "Imm2", "Imm3",
          #   "Ifn1", "Ifn2", "Ifn3",
          #   "Chl1", "Chl2", "Chl3",
          #   "Fib1", "Fib2", "Fib3",
          #   "Prolif1", "Prolif2", "Prolif3",
          #   "Tcell", "CLymph", "Bcell",
          #   "NKcell", "Monocyte", "MDendri.",
          #   "F.blast"
          # ) %>%
          #   rev()),

          Variable = str_split_fixed(string = Variable, pattern = ":", n = 2)[, 1] %>%
            str_replace("scaled_",""),

          # grouping modules
          Variable_class = case_when(
            Variable == "TIL" ~ "H&E",
            str_detect(Variable, "Denovo") ~ "De-novo",
            str_detect(Variable, "Gruosso2019") ~ "TIL-loc.",
            str_detect(Variable, "General") ~ "General",
            str_detect(Variable, "MCPcounter") ~ "MCPcounter",
            str_detect(Variable, "Control") ~ "Control",
            TRUE ~ "Unknown"
          ) %>%
            factor(levels = c(
              "H&E", "De-novo", "TIL-loc.", "General", "MCPcounter", "Control"
            )),

          Variable =  Variable %>%
            str_replace("Denovo","") %>%
            str_replace("Gruosso2019","") %>%
            str_replace("General","") %>%
            str_replace("MCPcounter","") %>%
            str_replace("Control","") %>%
            str_replace("_","") %>%
            str_replace("Tcell","T.Cell") %>%
            str_replace("Cyto.Lymphocyte","Cyto.Lympho") %>%
            str_replace("B.Lineage","B.Cell") %>%
            str_replace("NK.Cells","NK.Cell") %>%
            str_replace("Monocytic.Lineage","Monocyte") %>%
            str_replace("Myeloid.Dendritic","M.Dendritic") %>%
            str_replace("Neutrophils","Neutrophil") %>%
            str_replace("Fibroblasts","Fibroblast") %>%
            str_replace("Non.breast.tissue","Non.Breast") %>%
            str_replace("Human_Behaviour","Behavioural") %>%
            factor(levels = c(
              "TIL",
              "TILsig", #"Immune","ECM", # Denovo
              "Immune", "Interferon", "Cholesterol", "Fibrosis", #  grusso
              "ECM", # "Immune", "Interferon", "Cholesterol",  # general
              "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
              "Monocyte", "M.Dendritic", "Neutrophil",
              "Endothelial", "Fibroblast",
              "Proliferation","Non.Breast", "Behavioural" # control
            ) %>% rev()
            ),


          Therapy_class = str_split_fixed(Therapy, ":", n = 2)[ , 2],
          Therapy_class = str_replace(Therapy_class, "Chemo", "") %>%
            factor(levels = c("DTX","NVB"))
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
                                                   HR_text = "HR(95%CI)",
                                                   P_text = "P(Padj*)",
                                                   Het_text= "Q(P)/I^2")) %>%
                        factor(levels = c("HR(95%CI)","P(Padj*)","Q(P)/I^2"))) %>%
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
    })



  p <- ggpubr::ggarrange(
    p_forest$tn + th_finher + labs(x = "Log-HR"),
    p_annot$tn + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    p_forest$her2 + th_finher + labs(x = "Log-HR") + theme(axis.text.y = element_blank()),
    p_annot$her2 + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    ncol = 4,
    nrow = 1,
    # widths = c(.16,.375,.09,.375), #old
    # widths = c(.21,.37,.14,.37),
    widths = c(.23,.36,.13,.36),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix, "_", event_type, ".pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()
}


plot_finher_mta_interaction(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_mta_interaction"),
  width = 7.5,
  height = 8.5)

plot_finher_mta_interaction(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_mta_interaction"),
  width = 7.5,
  height = 8.5)

plot_finher_mta_interaction(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  c(DTX = "#e41a1c", # red
             NVB = "#377eb8"), # blue,
  prefix = str_c(out_figures,"finher_mta_interaction"),
  width = 7.5,
  height = 8.5)

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
          # Het_text = "0.0(.960)-552.5", # to test text size

          # Variable = str_replace(Variable, "_scaled","") %>%
          #   str_replace("_Fc","") %>%
          #   str_replace("Innate","Inn") %>%
          #   str_replace("Adhesion","Adh.") %>%
          #   str_replace("Immune","Imm") %>%
          #   str_replace("Interferon","Ifn") %>%
          #   str_replace("Cholesterol","Chl") %>%
          #   str_replace("Fibrosis","Fib") %>%
          #   str_replace("Proliferation","Prolif") %>%
          #   str_replace("Fibroblast","F.blast") %>%
          #   str_replace("_","."),
          #
          # Variable_class = case_when(
          #   str_detect(Variable, "TIL") ~ "TIL",
          #   TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
          # ) %>%
          #   factor(levels = c(
          #     "TIL",
          #     "Celltype / Proliferation / TIL-localization signatures"
          #   )),
          #
          # Variable = factor(Variable, levels = c(
          #   "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg","TILsig.Inn",
          #   "TILsig.ECM", "TILsig.Adh.",
          #   "Imm1", "Imm2", "Imm3",
          #   "Ifn1", "Ifn2", "Ifn3",
          #   "Chl1", "Chl2", "Chl3",
          #   "Fib1", "Fib2", "Fib3",
          #   "Prolif1", "Prolif2", "Prolif3",
          #   "Tcell", "CLymph", "Bcell",
          #   "NKcell", "Monocyte", "MDendri.",
          #   "F.blast"
          # ) %>%
          #   rev())

          Variable = str_split_fixed(string = Variable, pattern = ":", n = 2)[, 1] %>%
            str_replace("scaled_",""),

          # grouping modules
          Variable_class = case_when(
            Variable == "TIL" ~ "H&E",
            str_detect(Variable, "Denovo") ~ "De-novo",
            str_detect(Variable, "Gruosso2019") ~ "TIL-loc.",
            str_detect(Variable, "General") ~ "General",
            str_detect(Variable, "MCPcounter") ~ "MCPcounter",
            str_detect(Variable, "Control") ~ "Control",
            TRUE ~ "Unknown"
          ) %>%
            factor(levels = c(
              "H&E", "De-novo", "TIL-loc.", "General", "MCPcounter", "Control"
            )),

          Variable =  Variable %>%
            str_replace("Denovo","") %>%
            str_replace("Gruosso2019","") %>%
            str_replace("General","") %>%
            str_replace("MCPcounter","") %>%
            str_replace("Control","") %>%
            str_replace("_","") %>%
            str_replace("Tcell","T.Cell") %>%
            str_replace("Cyto.Lymphocyte","Cyto.Lympho") %>%
            str_replace("B.Lineage","B.Cell") %>%
            str_replace("NK.Cells","NK.Cell") %>%
            str_replace("Monocytic.Lineage","Monocyte") %>%
            str_replace("Myeloid.Dendritic","M.Dendritic") %>%
            str_replace("Neutrophils","Neutrophil") %>%
            str_replace("Fibroblasts","Fibroblast") %>%
            str_replace("Non.breast.tissue","Non.Breast") %>%
            str_replace("Human_Behaviour","Behavioural") %>%
            factor(levels = c(
              "TIL",
              "TILsig", #"Immune","ECM", # Denovo
              "Immune", "Interferon", "Cholesterol", "Fibrosis", #  grusso
              "ECM", # "Immune", "Interferon", "Cholesterol",  # general
              "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
              "Monocyte", "M.Dendritic", "Neutrophil",
              "Endothelial", "Fibroblast",
              "Proliferation","Non.Breast", "Behavioural" # control
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
        dplyr::mutate(key = purrr::map_chr(key,
                                           ~switch(.x,
                                                   HR_text = "HR(95%CI)",
                                                   P_text = "P(Padj*)",
                                                   Het_text = "Q(P)/I^2")) %>%
                        factor(levels = c("HR(95%CI)","P(Padj*)", "Q(P)/I^2"))) %>%
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
    })



  p <- ggpubr::ggarrange(
    p_forest$tn + th_finher + labs(x = "Log-HR"),
    p_annot$tn + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    p_forest$her2 + th_finher + labs(x = "Log-HR") + theme(axis.text.y = element_blank()),
    p_annot$her2 + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

    ncol = 4,
    nrow = 1,
    # widths = c(.16,.375,.09,.375), #old
    # widths = c(.21,.37,.14,.37),
    widths = c(.23,.36,.13,.36),
    labels = c("A","","B",""),
    hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
  )

  pdf(file = str_c(prefix, "_", event_type, ".pdf") %>% str_to_lower(),
      width = width, height = height)
  print(p)
  dev.off()
}


plot_finher_prognosis(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "DDFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_prognosis"),
  width = 7.5,
  height = 8.5)


plot_finher_prognosis(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "RFS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_prognosis"),
  width = 7.5,
  height = 8.5)


plot_finher_prognosis(
  tn_cox_sum = finher_tn,
  her2_cox_sum = finher_her2,
  event_type = "OS",
  inner_text_size = 2.75,
  outer_text_size = 10,
  color =  NULL,
  prefix = str_c(out_figures,"finher_prognosis"),
  width = 7.5,
  height = 8.5)


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
         "TIL", "TILsig", "TILsig_APP_Fc", "TILsig_Immune",
         "TILsig_IFNg", "TILsig_Innate",
         "TILsig_ECM", "TILsig_Adhesion")

plot_corr <- purrr::map(

  c("TN","HER2","ALL"),

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
        dplyr::rename(TIL = "StrLy_Mean",
                      TILsig_unscaled = "TILsig",
                      TILsig = "TILsig_scaled") %>%
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
    Subtype = factor(Subtype, levels =  c("TN","HER2","ALL")),
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





# FinHER KM curves
# ==============================================================================


# Prepare data
# ============

# TN subset
dat_km <- clin_finher %>%
  dplyr::filter(Subtype_IHC == "TN") %>%
  dplyr::rename(TIL = "StrLy_Mean")


# Relevant medians
dat_km_medians <- dat_km %>%
  dplyr::select(TILsig_Immune, TILsig_IFNg, Immune1, Interferon1, TIL, Tcell) %>%
  purrr::map(~median(.x, na.rm = T))


# Grouping data
dat_km <- bind_rows(

  dat_km %>%
    dplyr::mutate(
      Subgroup = "TIL",
      Score = if_else(TIL < dat_km_medians$TIL, "Low", "High")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "Tcell",
      Score = if_else(Tcell < dat_km_medians$Tcell, "Low", "High")
    ),

  # dat_km %>%
  #   dplyr::mutate(
  #     Subgroup = "Immune1",
  #     Score = if_else(Immune1 < dat_km_medians$Immune1, "Low", "High")
  #   ),
  #
  # dat_km %>%
  #   dplyr::mutate(
  #     Subgroup = "Interferon1",
  #     Score = if_else(Interferon1 < dat_km_medians$Interferon1, "Low", "High")
  #   ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "TILsig_Immune",
      Score = if_else(TILsig_Immune < dat_km_medians$TILsig_Immune, "Low", "High")
    ),

  dat_km %>%
    dplyr::mutate(
      Subgroup = "TILsig_IFNg",
      Score = if_else(TILsig_IFNg < dat_km_medians$TILsig_IFNg, "Low", "High")
    )
)


dat_km <- dat_km %>%  dplyr::mutate(
  Subgroup = str_c(Chemo, "*", Subgroup)
)

dat_km %>%
  group_by(Subgroup, Score) %>%
  summarise(n = n())
#   Subgroup          Score     n
# 1 DTX * Immune1     High     29
# 2 DTX * Immune1     Low      31
# 3 DTX * Interferon1 High     31
# 4 DTX * Interferon1 Low      29
# 5 DTX * Tcell       High     29
# 6 DTX * Tcell       Low      31
# 7 DTX * TIL         High     27
# 8 DTX * TIL         Low      32
# 9 DTX * TIL         NA        1
# 10 NVB * Immune1     High     31
# 11 NVB * Immune1     Low      29
# 12 NVB * Interferon1 High     29
# 13 NVB * Interferon1 Low      31
# 14 NVB * Tcell       High     31
# 15 NVB * Tcell       Low      29
# 16 NVB * TIL         High     34
# 17 NVB * TIL         Low      25
# 18 NVB * TIL         NA        1

#   Subgroup          Score     n
# 1 DTX*Tcell         High     29
# 2 DTX*Tcell         Low      31
# 3 DTX*TIL           High     27
# 4 DTX*TIL           Low      32
# 5 DTX*TIL           NA        1
# 6 DTX*TILsig_IFNg   High     31
# 7 DTX*TILsig_IFNg   Low      29
# 8 DTX*TILsig_Immune High     32
# 9 DTX*TILsig_Immune Low      28
# 10 NVB*Tcell         High     31
# 11 NVB*Tcell         Low      29
# 12 NVB*TIL           High     34
# 13 NVB*TIL           Low      25
# 14 NVB*TIL           NA        1
# 15 NVB*TILsig_IFNg   High     29
# 16 NVB*TILsig_IFNg   Low      31
# 17 NVB*TILsig_Immune High     28
# 18 NVB*TILsig_Immune Low      32



plot_finher_km <- function(
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
    as.formula(str_c("Surv(", time_var, ",", event_var, ") ~ Score")),
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
      "#000000", #"black",
      "#bdbdbd" #"gray"
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
    # risk.tables.col = "strata",

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

  km_fig <- km_fig[c("DTX*TILsig_Immune", "DTX*TILsig_IFNg",
                     "DTX*Tcell", "DTX*TIL",

                     "NVB*TILsig_Immune", "NVB*TILsig_IFNg",
                     "NVB*Tcell", "NVB*TIL") ]





  zz <- arrange_ggsurvplots(
    x = km_fig, print = F,
    ncol = 2,
    nrow = 4,
    surv.plot.height = .75,
    risk.table.height = .25,
    newpage = F,
    byrow =T
  )
  pdf(file = paste0(out_figures, "finher_km_", prefix, ".pdf"),
      width = 7.5, height = 8)
  print(zz)
  dev.off()


  # Risk table is text annotated; no need of legend

  # pdf(file = paste0(out_figures, "finher_km_", prefix, "legend.pdf"),
  #     width = 1.85, height = 1.25)
  # grid.draw(km_fig_legend)
  # dev.off()

}


plot_finher_km(
  dat_km = dat_km,
  event_var = "DDFS_Event",
  time_var = "DDFS_Time_Years",
  prefix = "DDFS")

plot_finher_km(
  dat_km = dat_km,
  event_var = "RFS_Event",
  time_var = "RFS_Time_Years",
  prefix = "RFS")

plot_finher_km(
  dat_km = dat_km,
  event_var = "OS_Event",
  time_var = "OS_Time_Years",
  prefix = "OS")



#
# ==============================================================================

