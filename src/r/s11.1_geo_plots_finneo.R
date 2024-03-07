# s11.1_geo plots_finneo.R


load("results/data/geo_prog_finneo.RData")
load("results/data/geo_inter_finneo.RData")

load("results/data/clin_neoadj_finneo.RData")


geo_prog <- geo_prog_finneo
geo_inter <- geo_inter_finneo

# Color scheme
# >>>>>>>>>>>>



# color_chemo <- c(DTX = "#e41a1c", # red
#                  NVB = "#377eb8") # blue
# # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
#
# color_tra <- c(TRA = "#a65628", # brown
#                noTRA = "#984ea3") # purple
# # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9

inner_text_size = 2.75
outer_text_size = 10

custom_colours <- c(
  # https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9

  "#e41a1c", # red for AAA+T (DTX)
  # "#377eb8", # blue for AAA (NVB)
  "#008B45FF", # green for AAA (NVB)
  "#a65628", # brown for AAA+T+TRA (TRA)

  "#008B45FF", # green for HR
  "#631879FF", # purple for HER2
  "#A20056FF", # majenta for TN

  "#ff7f00" # orange for A0A+T
)


nme <- c(
  # # nejm/material-design colurs from ggsci package
  "AAA+T", # "#BC3C29FF"
  "AAA", # "#0072B5FF"
  "AAA+T+TRA", # "#E18727FF"

  "HR", # "#20854EFF"
  "HER2", # "#7876B1FF"
  "TN", # "#EE4C97FF"

  # # aaas/nejm colurs from ggsci package
  "A0A+T" # "#1B1919FF"
)

names(custom_colours) <- nme



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



width = 7.5
height = 8.5



# GEO interaction subtype specific regimen summary
# ==============================================================================
clin_neoadj_finneo  %>%
  dplyr::group_by(Subtype_ihc2, Arm_consolidated) %>%
  dplyr::summarize(N = n())
#   Subtype_ihc2 Arm_consolidated                                                         N
# 1 HR           A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     139
# 2 HR           AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy    91
# 3 HR           AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     521
# 4 HR-HER2+     AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy      80
# 5 HR-HER2+     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy        58
# 6 HR+HER2+     AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy      55
# 7 HR+HER2+     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy        29
# 8 TN           A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     150
# 9 TN           AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy    83
# 10 TN           AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     293




#
# ==============================================================================




# GEO prognosis plots
# ==============================================================================

# 1. Consolidate subtype specific prognostic summary
# 2. Format consolidated data for forest plot and plot annotation



# Individual.sig >>>>>>>>>>>>>>>>>>>>>>>

geo_prog_merged <-

  # preparing geo_prog for consolidation
  purrr::map(

    names(geo_prog$individual.sig), # per subtype

    function(nme, x){

      x[[nme]] %>%
        dplyr::mutate(Subtype = nme) # updating subtype info in each per subtype summary
    },

    x = geo_prog$individual.sig
  ) %>%

  # Consolidation
  bind_rows() %>%

  dplyr::mutate(

    # Formatting prognosis annotation
    Text_or = str_c(
      # Estimate in Log scale; exp(Estimate) in original scale(ratio)
      exp(Estimate) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      exp(Lower_ci) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      exp(Upper_ci) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      ")"
    ),
    Text_p = str_c(
      if_else(
        P < .001,
        "<.001",
        formatC(x = P, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      )
      ,
      "(",
      if_else(
        P_adj < .001,
        "<.001",
        formatC(x = P_adj, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      ),")"
    ),
    Text_het = if_else(
      is.na(Q) & is.na(Q_p) & is.na(I2),
      "",
      str_c(
        Q %>% formatC(digits = 1, width = 1, format = "f"),
        "(",
        if_else(
          Q_p < .001,
          "<.001",
          formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        ),
        "),",
        I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
      )
    ),


    Module_name = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 1] %>%
      str_replace("scaled_",""),

    # grouping modules
    Module_class = case_when(
      Module_name == "TIL" ~ "H&E",

      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Pooled") ~ "Published-pooled",
      str_detect(Module_name, "MCPcounter") ~ "Cell",

      Module_name == "Hamy2016_Immune" ~ "Immune",
      Module_name == "Yang2018_Immune" ~ "Immune",

      Module_name == "Gruosso2019_Interferon" ~ "Interferon",
      Module_name == "Farmer2009_MX1" ~ "Interferon",
      Module_name == "Hamy2016_Interferon" ~ "Interferon",
      Module_name == "Nirmal2018_Interferon" ~ "Interferon",


      Module_name == "Hamy2016_Ecm" ~ "Fibrosis",
      Module_name == "Naba2014_Ecmcore" ~ "Fibrosis",
      Module_name == "Triulzi2013_Ecm" ~ "Fibrosis",

      Module_name == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

      TRUE ~ "Unknown"
    ) %>%
      factor(
        levels = c(
          "H&E", "De-novo",
          "Immune", "Interferon", "Chol.", "Fibrosis",
          "Published-pooled", "Cell")
      ),

    Module_name =  Module_name %>%
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
      )
  )


# forest plotting
p_prog <- p_prog_annot <- list() # template

for(subtype in c("TN", "HR")){

  # subtype = "TN"
  print(subtype)

  # forest plot
  p_prog[[subtype]] <- geo_prog_merged %>%
    dplyr::filter(Subtype == subtype) %>%
    ggplot() +
    geom_point(aes(x = Estimate, y = Module_name)) +
    geom_errorbarh(aes(xmin = Lower_ci, xmax = Upper_ci, y = Module_name), height = .1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    facet_grid(facets = Module_class ~ 1, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  # forest plot annotation
  p_prog_annot[[subtype]] <- geo_prog_merged %>%
    tidyr::gather("key", "value", "Text_or","Text_p","Text_het") %>%
    dplyr::filter(Subtype == subtype) %>%
    dplyr::mutate(
      key = purrr::map_chr(key, ~switch(
        .x,
                                        Text_or = "OR(95%CI)", # "Log-OR (95% CI)"
                                        Text_p = "P(Padj)", # "Pval (Pval-adj)"
                                        Text_het = "Q(P),I^2" # "Cochran's Q (Q-Pval), I^2"
      )) %>% factor(levels = c("OR(95%CI)", "P(Padj)", "Q(P),I^2"),
                    labels =  c(expression("OR"^a*"(95%CI) "),
                                "P(Padj)",
                                expression("Q(P),I"^2))
                    )
      ) %>%
    ggplot() +
    geom_text(aes(x = key, y = Module_name, label = value),
              size = inner_text_size) +
    facet_grid(facets = Module_class ~ 1, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    scale_x_discrete(labels = ~parse(text = .x)) +
    # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}

p <- ggpubr::ggarrange(
  p_prog$TN + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*")  ")) +
    theme(axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_prog_annot$TN + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_prog$HR + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*")  ")) +
    theme(axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_prog_annot$HR + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.22,.37,.13,.37),
  widths = c(.24,.36,.12,.36),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

# pdf(file = str_c(out_figures,"geo_prognosis_individual.sig_v2.1.pdf") %>% str_to_lower(),
#     width = 7.5, height = 6.5)
pdf(file = str_c(out_figures,"geo_individual-sig_finneo_prognosis.pdf") %>% str_to_lower(),
    width = 7.5, height = 6.5)
print(p)
dev.off()




# Pooled.sig >>>>>>>>>>>>>>>>>>>>>>>

geo_prog_merged <-

  # preparing geo_prog for consolidation
  purrr::map(

    names(geo_prog$pooled.sig), # per subtype

    function(nme, x){

      x[[nme]] %>%
        dplyr::mutate(Subtype = nme) # updating subtype info in each per subtype summary
    },

    x = geo_prog$pooled.sig
  ) %>%

  # Consolidation
  bind_rows() %>%

  dplyr::mutate(

    # Formatting prognosis annotation
    Text_or = str_c(
      # Estimate in Log scale; exp(Estimate) in original scale(ratio)
      exp(Estimate) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      exp(Lower_ci) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      exp(Upper_ci) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      ")"
    ),
    Text_p = str_c(
      if_else(
        P < .001,
        "<.001",
        formatC(x = P, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      )
      ,
      "(",
      if_else(
        P_adj < .001,
        "<.001",
        formatC(x = P_adj, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      ),")"
    ),
    Text_het = if_else(
      is.na(Q) & is.na(Q_p) & is.na(I2),
      "",
      str_c(
        Q %>% formatC(digits = 1, width = 1, format = "f"),
        "(",
        if_else(
          Q_p < .001,
          "<.001",
          formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        ),
        "),",
        I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
      )
    ),


    Module_name = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 1] %>%
      str_replace("scaled_",""),

    # grouping modules
    Module_class = case_when(
      Module_name == "TIL" ~ "H&E",

      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Pooled") ~ "Published-pooled",
      str_detect(Module_name, "MCPcounter") ~ "Cell",

      Module_name == "Hamy2016_Immune" ~ "Immune",
      Module_name == "Yang2018_Immune" ~ "Immune",

      Module_name == "Gruosso2019_Interferon" ~ "Interferon",
      Module_name == "Farmer2009_MX1" ~ "Interferon",
      Module_name == "Hamy2016_Interferon" ~ "Interferon",
      Module_name == "Nirmal2018_Interferon" ~ "Interferon",


      Module_name == "Hamy2016_Ecm" ~ "Fibrosis",
      Module_name == "Naba2014_Ecmcore" ~ "Fibrosis",
      Module_name == "Triulzi2013_Ecm" ~ "Fibrosis",

      Module_name == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

      TRUE ~ "Unknown"
    ) %>%
      factor(
        levels = c(
          "H&E", "De-novo",
          "Immune", "Interferon", "Chol.", "Fibrosis",
          "Published-pooled", "Cell")
      ),

    Module_name =  Module_name %>%
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
      )
  )


# forest plotting
p_prog <- p_prog_annot <- list() # template

for(subtype in c("TN", "HR")){

  # subtype = "TN"
  print(subtype)

  # forest plot
  p_prog[[subtype]] <- geo_prog_merged %>%
    dplyr::filter(Subtype == subtype) %>%
    ggplot() +
    geom_point(aes(x = Estimate, y = Module_name)) +
    geom_errorbarh(aes(xmin = Lower_ci, xmax = Upper_ci, y = Module_name), height = .1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    facet_grid(facets = Module_class ~ 1, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )


  # forest plot annotation
  p_prog_annot[[subtype]] <- geo_prog_merged %>%
    tidyr::gather("key", "value", "Text_or","Text_p","Text_het") %>%
    dplyr::filter(Subtype == subtype) %>%
    dplyr::mutate(
      key = purrr::map_chr(key, ~switch(
        .x,
        Text_or = "OR(95%CI)", # "Log-OR (95% CI)"
        Text_p = "P(Padj)", # "Pval (Pval-adj)"
        Text_het = "Q(P),I^2" # "Cochran's Q (Q-Pval), I^2"
      )) %>% factor(levels = c("OR(95%CI)", "P(Padj)", "Q(P),I^2"),
                    labels =  c(expression("OR"^a*"(95%CI) "),
                                "P(Padj)",
                                expression("Q(P),I"^2))
                    )
      ) %>%
    ggplot() +
    geom_text(aes(x = key, y = Module_name, label = value),
              size = inner_text_size) +
    facet_grid(facets = Module_class ~ 1, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    scale_x_discrete(labels = ~parse(text = .x)) +
    # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}

p <- ggpubr::ggarrange(
  p_prog$TN + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_prog_annot$TN + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_prog$HR + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_prog_annot$HR + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.22,.37,.13,.37),
  widths = c(.225,.36,.135,.36),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

# pdf(file = str_c(out_figures,"geo_prognosis_pooled.sig_v2.1.pdf") %>% str_to_lower(),
#     width = 7.5, height = 4.25)
pdf(file = str_c(out_figures,"geo_pooled-sig_finneo_prognosis.pdf") %>% str_to_lower(),
    width = 7.5, height = 4.25)
print(p)
dev.off()


#
# ==============================================================================



# GEO interaction plots
# ==============================================================================


# 1. Consolidate subtype-analysis specific interaction summary
# 2. Format consolidated data for forest plot and plot annotation


# individual.sig >>>>>>>>>>>>>>>>>>

geo_inter_merged <-

  # preparing geo_inter for consolidation
  purrr::map(

    names(geo_inter$individual.sig), # per subtype

    function(subtype, x){


      purrr::map(

        names(x[[subtype]]), # per subtype per analysis

        function(analysis, subtype, xx){
          xx[[analysis]] %>%
            dplyr::mutate(Subtype = subtype, Analysis = analysis)
          # updating subtype and analysis(per-subtype-analysis) info in each per subtype summary
        },

        subtype = subtype,
        xx = x[[subtype]]
      ) %>%

        # Per subtype, analysis level consolidation
        bind_rows()

    },

    x = geo_inter$individual.sig
  ) %>%

  # Subtype level consolidation
  bind_rows()




geo_inter_merged <-

  geo_inter_merged %>%
  # subsettting interaction terms
  dplyr::filter(str_detect(Module_name, ":")) %>%

  dplyr::mutate(

    # # format sample size and pcr per arm
    # Text_arm = str_c(
    #   Interaction_var_levels,
    #   ": ",
    #   "N-", N_patients,
    #   ", pCR-", N_response,"(", Percent_repsonse %>% round(digits = 2), "%)"
    # ),

    # Formatting prognosis annotation
    Text_or = str_c(
      # Estimate in Log scale; exp(Estimate) in original scale(ratio)
      exp(Estimate) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      exp(l95) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      exp(u95) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      ")"
    ),
    Text_p = str_c(
      if_else(
        P < .001,
        "<.001",
        formatC(x = P, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      )
      ,
      "(",
      if_else(
        P_inter_adj < .001,
        "<.001",
        formatC(x = P_inter_adj, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      ),
      ")"
    ),
    Text_p = if_else(is.na(Text_p), "", Text_p),

    # manually aligning non-MTA OR text
    Text_or_size = purrr::map_int(Text_or,str_length),
    Text_or = format(x = Text_or, width = max(Text_or_size), justify = "left"),
    Text_or = str_pad(Text_or, max(Text_or_size) + 12.5, "left", " "),

    Text_het = if_else(
      is.na(Q) & is.na(Q_p) & is.na(I2),
      "",
      str_c(
        Q %>% formatC(digits = 1, width = 1, format = "f"),
        "(",
        if_else(
          Q_p < .001,
          "<.001",
          formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        ),
        "),",
        I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
      )
    ),

    Arm_or_subtype =  str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
      str_replace("Subtype_ihc", "") %>% # for subtype * sig interaction per arm
      str_replace("Arm_consolidated", "") %>%
      #   str_replace("\\+noTaxane", "") %>%
      #   str_replace("\\+Taxane", "+T") %>%
      #   str_replace("///Trastuzumab", "+TRA") %>%
      #   str_replace("///No_her2_agent", "") %>%
      #   str_replace("///No_hormone_therapy", "") %>%
      #   str_replace("///No_other_therapy", "") %>%
      factor(levels = c("AAA+T+TRA", "AAA", "A0A+T", "AAA+T", "HR", "HER2", "TN")),


    Module_name = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 1] %>%
      str_replace("scaled_",""),

    # grouping modules
    Module_class = case_when(
      Module_name == "TIL" ~ "H&E",

      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Pooled") ~ "Published-pooled",
      str_detect(Module_name, "MCPcounter") ~ "Cell",

      Module_name == "Hamy2016_Immune" ~ "Immune",
      Module_name == "Yang2018_Immune" ~ "Immune",

      Module_name == "Gruosso2019_Interferon" ~ "Interferon",
      Module_name == "Farmer2009_MX1" ~ "Interferon",
      Module_name == "Hamy2016_Interferon" ~ "Interferon",
      Module_name == "Nirmal2018_Interferon" ~ "Interferon",


      Module_name == "Hamy2016_Ecm" ~ "Fibrosis",
      Module_name == "Naba2014_Ecmcore" ~ "Fibrosis",
      Module_name == "Triulzi2013_Ecm" ~ "Fibrosis",

      Module_name == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

      TRUE ~ "Unknown"
    ) %>%
      factor(
        levels = c(
          "H&E", "De-novo",
          "Immune", "Interferon", "Chol.", "Fibrosis",
          "Published-pooled", "Cell")
      ),

    Module_name =  Module_name %>%
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

    Analysis_group = str_c(Subtype, ":", Analysis)

  )



# forest plotting
p_inter <- p_inter_annot <- list() # template

for(anagroup in geo_inter_merged$Analysis_group %>% unique()){

  # anagroup = "TN:AAA.T_or_NoT"
  print(anagroup)

  # forest plot
  p_inter[[anagroup]] <- geo_inter_merged %>%
    dplyr::filter(Analysis_group == anagroup) %>%
    ggplot() +
    geom_point(aes(x = Estimate, y = Module_name, color = Arm_or_subtype, group = Arm_or_subtype),
               position = ggstance::position_dodgev(height = .75)) +
    geom_errorbarh(aes(xmin = l95, xmax = u95, y = Module_name, color = Arm_or_subtype, group = Arm_or_subtype),
                   position = ggstance::position_dodgev(height = .75),
                   height = .1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_color_manual(values = custom_colours) +
    guides(color = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(facets = Module_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )



  p_inter_annot[[anagroup]] <- geo_inter_merged %>%
    dplyr::filter(Analysis_group == anagroup) %>%
    dplyr::mutate(Subtype = Arm_or_subtype) %>%
    # !!! geo_inter_merged already contains a subtype variable
    dplyr::rename("OR(95%CI)" = "Text_or",
                  "P(Padj)" = "Text_p",
                  "Q(P),I^2" = "Text_het") %>%
    tidyr::gather("key", "value",
                  "OR(95%CI)", "P(Padj)", "Q(P),I^2") %>%
    dplyr::mutate(key = key %>%
                    factor(levels = c("OR(95%CI)", "P(Padj)", "Q(P),I^2"),
                           labels =  c(expression("OR"^a*"(95%CI) "),
                                       "P(Padj)",
                                       expression("Q(P),I"^2)))
                  ) %>%
    ggplot() +
    geom_text(aes(x = key, y = Module_name, label = value, color = Arm_or_subtype, group = Arm_or_subtype),
              position = ggstance::position_dodgev(height = .75),
              size = inner_text_size) +
    scale_color_manual(values = custom_colours) +
    guides(color = "none") +
    facet_grid(facets = Module_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    scale_x_discrete(labels = ~parse(text = .x)) +
    # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}



p <- ggpubr::ggarrange(
  p_inter$`TN:AAA.T_or_NoT` + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_inter_annot$`TN:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_inter$`HR:AAA.T_or_NoT` + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_inter_annot$`HR:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.21,.37,.14,.37),
  widths = c(.24,.36,.12,.36),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

# pdf(file = str_c(out_figures,"geo_mta_interaction_individual.sig_v2.1.pdf") %>% str_to_lower(),
#     width = 7.5, height = 6.5)
pdf(file = str_c(out_figures,"geo_individual-sig_finneo_pcr-mta-immune_interactionv.pdf") %>% str_to_lower(),
    width = 7.5, height = 6.5)
print(p)
dev.off()


# pooled.sig >>>>>>>>>>>>>>>>>>

geo_inter_merged <-

  # preparing geo_inter for consolidation
  purrr::map(

    names(geo_inter$pooled.sig), # per subtype

    function(subtype, x){


      purrr::map(

        names(x[[subtype]]), # per subtype per analysis

        function(analysis, subtype, xx){
          xx[[analysis]] %>%
            dplyr::mutate(Subtype = subtype, Analysis = analysis)
          # updating subtype and analysis(per-subtype-analysis) info in each per subtype summary
        },

        subtype = subtype,
        xx = x[[subtype]]
      ) %>%

        # Per subtype, analysis level consolidation
        bind_rows()

    },

    x = geo_inter$pooled.sig
  ) %>%

  # Subtype level consolidation
  bind_rows()




geo_inter_merged <-

  geo_inter_merged %>%
  # subsettting interaction terms
  dplyr::filter(str_detect(Module_name, ":")) %>%

  dplyr::mutate(

    # # format sample size and pcr per arm
    # Text_arm = str_c(
    #   Interaction_var_levels,
    #   ": ",
    #   "N-", N_patients,
    #   ", pCR-", N_response,"(", Percent_repsonse %>% round(digits = 2), "%)"
    # ),

    # Formatting prognosis annotation
    Text_or = str_c(
      # Estimate in Log scale; exp(Estimate) in original scale(ratio)
      exp(Estimate) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      exp(l95) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      exp(u95) %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      ")"
    ),
    Text_p = str_c(
      if_else(
        P < .001,
        "<.001",
        formatC(x = P, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      )
      ,
      "(",
      if_else(
        P_inter_adj < .001,
        "<.001",
        formatC(x = P_inter_adj, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      ),
      ")"
    ),
    Text_p = if_else(is.na(Text_p), "", Text_p),

    # manually aligning non-MTA OR text
    Text_or_size = purrr::map_int(Text_or,str_length),
    Text_or = format(x = Text_or, width = max(Text_or_size), justify = "left"),
    Text_or = str_pad(Text_or, max(Text_or_size) + 9, "left", " "),

    Text_het = if_else(
      is.na(Q) & is.na(Q_p) & is.na(I2),
      "",
      str_c(
        Q %>% formatC(digits = 1, width = 1, format = "f"),
        "(",
        if_else(
          Q_p < .001,
          "<.001",
          formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        ),
        "),",
        I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
      )
    ),

    Arm_or_subtype =  str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
      str_replace("Subtype_ihc", "") %>% # for subtype * sig interaction per arm
      str_replace("Arm_consolidated", "") %>%
      #   str_replace("\\+noTaxane", "") %>%
      #   str_replace("\\+Taxane", "+T") %>%
      #   str_replace("///Trastuzumab", "+TRA") %>%
      #   str_replace("///No_her2_agent", "") %>%
      #   str_replace("///No_hormone_therapy", "") %>%
      #   str_replace("///No_other_therapy", "") %>%
      factor(levels = c("AAA+T+TRA", "AAA", "A0A+T", "AAA+T", "HR", "HER2", "TN")),


    Module_name = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 1] %>%
      str_replace("scaled_",""),

    # grouping modules
    Module_class = case_when(
      Module_name == "TIL" ~ "H&E",

      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Pooled") ~ "Published-pooled",
      str_detect(Module_name, "MCPcounter") ~ "Cell",

      Module_name == "Hamy2016_Immune" ~ "Immune",
      Module_name == "Yang2018_Immune" ~ "Immune",

      Module_name == "Gruosso2019_Interferon" ~ "Interferon",
      Module_name == "Farmer2009_MX1" ~ "Interferon",
      Module_name == "Hamy2016_Interferon" ~ "Interferon",
      Module_name == "Nirmal2018_Interferon" ~ "Interferon",


      Module_name == "Hamy2016_Ecm" ~ "Fibrosis",
      Module_name == "Naba2014_Ecmcore" ~ "Fibrosis",
      Module_name == "Triulzi2013_Ecm" ~ "Fibrosis",

      Module_name == "Sorrentino2014_Chol" ~ "Chol.", #"Cholesterol"

      TRUE ~ "Unknown"
    ) %>%
      factor(
        levels = c(
          "H&E", "De-novo",
          "Immune", "Interferon", "Chol.", "Fibrosis",
          "Published-pooled", "Cell")
      ),

    Module_name =  Module_name %>%
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

    Analysis_group = str_c(Subtype, ":", Analysis)

  )



# forest plotting
p_inter <- p_inter_annot <- list() # template

for(anagroup in geo_inter_merged$Analysis_group %>% unique()){

  # anagroup = "TN:AAA.T_or_NoT"
  print(anagroup)

  # forest plot
  p_inter[[anagroup]] <- geo_inter_merged %>%
    dplyr::filter(Analysis_group == anagroup) %>%
    ggplot() +
    geom_point(aes(x = Estimate, y = Module_name, color = Arm_or_subtype, group = Arm_or_subtype),
               position = ggstance::position_dodgev(height = .75)) +
    geom_errorbarh(aes(xmin = l95, xmax = u95, y = Module_name, color = Arm_or_subtype, group = Arm_or_subtype),
                   position = ggstance::position_dodgev(height = .75),
                   height = .1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_color_manual(values = custom_colours) +
    guides(color = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(facets = Module_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    theme(
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 0.1, b = 5.5, l = 5.5, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )



  p_inter_annot[[anagroup]] <- geo_inter_merged %>%
    dplyr::filter(Analysis_group == anagroup) %>%
    dplyr::mutate(Subtype = Arm_or_subtype) %>%
    # !!! geo_inter_merged already contains a subtype variable
    dplyr::rename("OR(95%CI)" = "Text_or",
                  "P(Padj)" = "Text_p",
                  "Q(P),I^2" = "Text_het") %>%
    tidyr::gather("key", "value",
                  "OR(95%CI)", "P(Padj)", "Q(P),I^2") %>%
    dplyr::mutate(key = key %>%
                    factor(levels = c("OR(95%CI)", "P(Padj)", "Q(P),I^2"),
                           labels =  c(expression("OR"^a*"(95%CI) "),
                                       "P(Padj)",
                                       expression("Q(P),I"^2)))
                  ) %>%
    ggplot() +
    geom_text(aes(x = key, y = Module_name, label = value, color = Arm_or_subtype, group = Arm_or_subtype),
              position = ggstance::position_dodgev(height = .75),
              size = inner_text_size) +
    scale_color_manual(values = custom_colours) +
    guides(color = "none") +
    facet_grid(facets = Module_class ~ 1,
               switch = "y",
               scales = "free_y",
               space = "free_y") +
    theme_bw() +
    scale_x_discrete(labels = ~parse(text = .x)) +
    # Ref: https://stackoverflow.com/questions/73014834/how-to-create-subscripts-in-the-names-of-variables-in-r
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 4, l = 0.1, unit = "pt") # b=4 to accommodate subscript
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}



p <- ggpubr::ggarrange(
  p_inter$`TN:AAA.T_or_NoT` + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_inter_annot$`TN:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_inter$`HR:AAA.T_or_NoT` + th_finher + # labs(x = "Log-OR") +
    labs(x = expression("Log(OR"^a*") ")) +
    theme(axis.text.y = element_blank(),
          axis.title.x.bottom = element_text(size = outer_text_size-1)),
  p_inter_annot$`HR:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.21,.37,.14,.37),
  widths = c(.225,.36,.135,.36),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

# pdf(file = str_c(out_figures,"geo_mta_interaction_pooled.sig_v2.1.pdf") %>% str_to_lower(),
#     width = 7.5, height = 4.25)
pdf(file = str_c(out_figures,"geo_pooled-sig_finneo_pcr-mta-immune_interaction.pdf") %>% str_to_lower(),
    width = 7.5, height = 4.25)
print(p)
dev.off()



#
# ==============================================================================


# Clear memory
# ==============================================================================


rm(geo_prog, geo_prog_finneo,
   geo_inter, geo_inter_finneo)

#
# ==============================================================================
