# geo plots


load("results/data/geo_prog.RData")
load("results/data/geo_inter.RData")


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
  "#377eb8", # blue for AAA (NVB)
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



# GEO interaction subtype specific regimen summary
# ==============================================================================
clin_neoadj  %>%
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

geo_prog_merged <-

  # preparing geo_prog for consolidation
  purrr::map(

    names(geo_prog), # per subtype

    function(nme, x){

      x[[nme]] %>%
        dplyr::mutate(Subtype = nme) # updating subtype info in each per subtype summary
    },

    x = geo_prog
  ) %>%

  # Consolidation
  bind_rows() %>%

  dplyr::mutate(

    # Formatting prognosis annotation
    Text_or = str_c(
      Estimate %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      Lower_ci %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      Upper_ci %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
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
      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Gruosso2019") ~ "TIL-loc.",
      str_detect(Module_name, "General") ~ "General",
      str_detect(Module_name, "MCPcounter") ~ "MCPcounter",
      str_detect(Module_name, "Control") ~ "Control",
      TRUE ~ "Unknown"
    ) %>%
      factor(levels = c(
        "De-novo", "TIL-loc.", "General", "MCPcounter", "Control"
      )),

    Module_name =  Module_name %>%
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
        "TILsig", #"Immune","ECM", # Denovo
        "Immune", "Interferon", "Cholesterol", "Fibrosis", #  grusso
        "ECM", # "Immune", "Interferon", "Cholesterol",  # general
        "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
        "Monocyte", "M.Dendritic", "Neutrophil",
        "Endothelial", "Fibroblast",
        "Proliferation","Non.Breast", "Behavioural" # control
      ) %>% rev()
      )
  )


# forest plotting
p_prog <- p_prog_annot <- list() # template

for(subtype in c("TN", "HER2", "HR", "ALL")){

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
                                        Text_or = "Log-OR(95%CI)", # "Log-OR (95% CI)"
                                        Text_p = "P(Padj)", # "Pval (Pval-adj)"
                                        Text_het = "Q(P),I^2" # "Cochran's Q (Q-Pval), I^2"
      )) %>% factor(levels = c("Log-OR(95%CI)", "P(Padj)", "Q(P),I^2"))) %>%
    ggplot() +
    geom_text(aes(x = key, y = Module_name, label = value),
              size = inner_text_size) +
    facet_grid(facets = Module_class ~ 1, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}



#
# ==============================================================================



# GEO interaction plots
# ==============================================================================


# 1. Consolidate subtype-analysis specific interaction summary
# 2. Format consolidated data for forest plot and plot annotation


geo_inter_merged <-

  # preparing geo_inter for consolidation
  purrr::map(

    names(geo_inter), # per subtype

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

    x = geo_inter
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
      Estimate %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      "(",
      l95 %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " to ",
      u95 %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
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
      str_detect(Module_name, "Denovo") ~ "De-novo",
      str_detect(Module_name, "Gruosso2019") ~ "TIL-loc.",
      str_detect(Module_name, "General") ~ "General",
      str_detect(Module_name, "MCPcounter") ~ "MCPcounter",
      str_detect(Module_name, "Control") ~ "Control",
      TRUE ~ "Unknown"
    ) %>%
      factor(levels = c(
        "De-novo", "TIL-loc.", "General", "MCPcounter", "Control"
      )),

    Module_name =  Module_name %>%
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
        "TILsig", #"Immune","ECM", # Denovo
        "Immune", "Interferon", "Cholesterol", "Fibrosis", #  grusso
        "ECM", # "Immune", "Interferon", "Cholesterol",  # general
        "T.Cell", "Cyto.Lympho", "B.Cell", "NK.Cell",
        "Monocyte", "M.Dendritic", "Neutrophil",
        "Endothelial", "Fibroblast",
        "Proliferation","Non.Breast", "Behavioural" # control
      ) %>% rev()
      ),



    # Module_name = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 1] %>%
    #
    #   str_replace("_scaled","") %>%
    #   str_replace("_Fc","") %>%
    #   str_replace("Adhesion","Adh.") %>%
    #   str_replace("Immune","Imm") %>%
    #   str_replace("Interferon","Ifn") %>%
    #   str_replace("Cholesterol","Chl") %>%
    #   str_replace("Fibrosis","Fib") %>%
    #   str_replace("Proliferation","Prolif") %>%
    #   str_replace("CLymphocyte","C.Lymph") %>%
    #   str_replace("MDendritic","M.Dendri") %>%
    #   # str_replace("Fibroblast","F.blast") %>%
    #   str_replace_all("_",".") %>%
    #
    #   factor(levels = c(
    #     # "TIL", "TILsig", "TILsig.imm","TILsig.fib",
    #     "TIL", "TILsig", "TILsig.APP", "TILsig.Imm","TILsig.IFNg",
    #     "TILsig.ECM", "TILsig.Adh.",
    #     "Imm1", "Imm2", "Imm3",
    #     "Ifn1", "Ifn2", "Ifn3",
    #     "Chl1", "Chl2", "Chl3",
    #     "Fib1", "Fib2", "Fib3",
    #     "Prolif1", "Prolif2", "Prolif3",
    #     "Tcell", "C.Lymph", "Bcell",
    #     "NKcell", "Monocyte", "M.Dendri",
    #     "Fibroblast"
    #   ) %>% rev()
    #   ),
    #
    # # grouping modules
    # Module_class = case_when(
    #   str_detect(Module_name, "TIL") ~ "TIL",
    #   TRUE ~ "Celltype / Proliferation / TIL-localization signatures"
    # ) %>%
    #   factor(levels = c(
    #     "TIL",
    #     "Celltype / Proliferation / TIL-localization signatures"
    #   )),

    Analysis_group = str_c(Subtype, ":", Analysis)

  )


# response_sum_inter <- geo_inter_merged %>%
#   dplyr::group_by(Analysis_group) %>%
#   dplyr::summarise(Title = str_c(Text_arm %>% unique(), collapse = " | ")) %>%
#   dplyr::mutate(
#     Analysis_group2 = if_else(str_detect(Analysis_group, "ALL"),
#                               str_split_fixed(str_split_fixed(Analysis_group, ":", 2)[,2],
#                                               "\\.", 2)[, 1],
#                               str_split_fixed(Analysis_group,":",2)[, 1]),
#     Analysis_group2 = Analysis_group2 %>% str_replace_all("_plus_","+"),
#     Title = str_c(
#       Analysis_group2,
#       "  [", Title, "]"
#     )
#   )


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
    dplyr::rename("Log-OR(95%CI)" = "Text_or",
                  "P(Padj)" = "Text_p",
                  "Q(P),I^2" = "Text_het") %>%
    tidyr::gather("key", "value",
                  "Log-OR(95%CI)", "P(Padj)", "Q(P),I^2") %>%
    dplyr::mutate(key = key %>%
                    factor(levels = c("Log-OR(95%CI)", "P(Padj)", "Q(P),I^2"))) %>%
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
    theme(
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      # to manage space in ggarrage()
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.1, unit = "pt")
      # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt") # default margin
    )
}


#
# ==============================================================================




# # Print plots combined
# # ==============================================================================
#
# # Common theme
# th_finher <- theme(
#   text = element_text(size = outer_text_size),
#   axis.text.x.bottom = element_text(color = "black"),
#   axis.title.y = element_blank(),
#   # strip.background = element_blank(),
#   strip.text = element_blank(),
#   panel.grid = element_blank(),
#   # to adjust space between facet panels
#   panel.spacing.y = unit(x = 2, units = "pt")
# )
#
# width = 7.5
# height = 7
#
#
# # TN plot
# p <- ggpubr::ggarrange(
#   p_prog$TN + th_finher + labs(x = "Log-OR"),
#   p_prog_annot$TN + th_finher + labs(x = ""),
#
#   p_inter$`TN:AAA.T_or_NoT` + th_finher + labs(x = "Log-OR") + theme(axis.text.y = element_blank()),
#   p_inter_annot$`TN:AAA.T_or_NoT` + th_finher + labs(x = ""),
#
#   ncol = 4,
#   nrow = 1,
#   widths = c(.16,.375,.09,.375),
#   labels = c("A","","B",""),
#   hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
# )
#
# pdf(file = str_c(out_figures,"geo_tn_prog_plus_inter_taxane.pdf") %>% str_to_lower(),
#     width = width, height = height)
# print(p)
# dev.off()
#
#
# # HER2 plot
# p <- ggpubr::ggarrange(
#   p_prog$HER2 + th_finher + labs(x = "Log-OR"),
#   p_prog_annot$HER2 + th_finher + labs(x = ""),
#
#   p_inter$`HER2:AAA_plus_T.TRA_or_NoTRA` + th_finher + labs(x = "Log-OR") + theme(axis.text.y = element_blank()),
#   p_inter_annot$`HER2:AAA_plus_T.TRA_or_NoTRA` + th_finher + labs(x = ""),
#
#   ncol = 4,
#   nrow = 1,
#   widths = c(.16,.375,.09,.375),
#   labels = c("A","","B",""),
#   hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
# )
#
# pdf(file = str_c(out_figures,"geo_her2_prog_plus_inter_tra.pdf") %>% str_to_lower(),
#     width = width, height = height)
# print(p)
# dev.off()
#
#
# # HR plot
# p <- ggpubr::ggarrange(
#   p_prog$HR + th_finher + labs(x = "Log-OR"),
#   p_prog_annot$HR + th_finher + labs(x = ""),
#
#   p_inter$`HR:AAA.T_or_NoT` + th_finher + labs(x = "Log-OR") + theme(axis.text.y = element_blank()),
#   p_inter_annot$`HR:AAA.T_or_NoT` + th_finher + labs(x = ""),
#
#   ncol = 4,
#   nrow = 1,
#   widths = c(.16,.375,.09,.375),
#   labels = c("A","","B",""),
#   hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
# )
#
# pdf(file = str_c(out_figures,"geo_hr_prog_plus_inter_taxane.pdf") %>% str_to_lower(),
#     width = width, height = height)
# print(p)
# dev.off()
#
#
#
# # All plot
# p <- ggpubr::ggarrange(
#   p_inter$`ALL:AAA_plus_T.TN_or_HER2_or_HR` + th_finher + labs(x = "Log-OR"),
#   p_inter_annot$`ALL:AAA_plus_T.TN_or_HER2_or_HR` + th_finher + labs(x = ""),
#   ncol = 2,
#   nrow = 1,
#   widths = c(.35,.65)
#   # labels = c("A","","B",""),
#   # hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
# )
#
# pdf(file = str_c(out_figures,"geo_all_aaa_t_inter_subtype.pdf") %>% str_to_lower(),
#     width = width/1.5, height = height+3)
# print(p)
# dev.off()
#
# #
# # ==============================================================================



# Print MTA paper plots
# ==============================================================================

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


# MTA interaction plot
# >>>>>>>>>>>>>>>>>>>>

p <- ggpubr::ggarrange(
  p_inter$`TN:AAA.T_or_NoT` + th_finher + labs(x = "Log-OR"),
  p_inter_annot$`TN:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_inter$`HR:AAA.T_or_NoT` + th_finher + labs(x = "Log-OR") + theme(axis.text.y = element_blank()),
  p_inter_annot$`HR:AAA.T_or_NoT` + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.21,.37,.14,.37),
  widths = c(.23,.36,.13,.36),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

pdf(file = str_c(out_figures,"geo_mta_interaction.pdf") %>% str_to_lower(),
    width = width, height = height)
print(p)
dev.off()


# Prognosis plot
# >>>>>>>>>>>>>>

p <- ggpubr::ggarrange(
  p_prog$TN + th_finher + labs(x = "Log-OR"),
  p_prog_annot$TN + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  p_prog$HR + th_finher + labs(x = "Log-OR") + theme(axis.text.y = element_blank()),
  p_prog_annot$HR + th_finher + labs(x = "") + theme(strip.text.y = element_blank()),

  ncol = 4,
  nrow = 1,
  # widths = c(.21,.37,.14,.37),
  widths = c(.22,.37,.13,.37),
  labels = c("A","","B",""),
  hjust = c(-0.5,-0.5,0.5,-0.5) #default = -0.5
)

pdf(file = str_c(out_figures,"geo_prognosis.pdf") %>% str_to_lower(),
    width = width, height = height)
print(p)
dev.off()


#
# ==============================================================================

#
#
# # GEO TILsig-TIL localization sig correlation plots
# # ==============================================================================
#
# nme <- c("Immune1", "Immune2", "Immune3",
#          "Interferon1", "Interferon2", "Interferon3",
#          "Cholesterol1","Cholesterol2","Cholesterol3",
#          "Fibrosis1", "Fibrosis2", "Fibrosis3",
#          "Proliferation1", "Proliferation2", "Proliferation3",
#          "Tcell", "CLymphocyte", "Bcell", "NKcell", "Monocyte",
#          "MDendritic", "Fibroblast",
#          "TILsig", "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg",
#          "TILsig_ECM", "TILsig_Adhesion")
#
#
# plot_corr <- purrr::map(
#
#   c("ALL","HR","HER2", "TN"),
#
#   function(subtype, clin, nme){
#
#     # subtype = "TN"
#     # clin <- clin_finher
#
#     xclin <- switch(subtype,
#                     "ALL" = clin,
#                     "TN" = clin %>%
#                       dplyr::filter(Subtype_ihc == "TN"),
#                     "HER2" = clin %>%
#                       dplyr::filter(Subtype_ihc == "HER2"),
#                     "HR" = clin %>%
#                       dplyr::filter(Subtype_ihc == "HR"))
#
#     xcor <- cor(
#       xclin %>%
#         dplyr::rename(TILsig_unscaled = "TILsig",
#                       TILsig = "TILsig_scaled") %>%
#         dplyr::select(all_of(nme)),
#       method = "spearman")
#
#     as_tibble(xcor, rownames = "Rows") %>%
#       tidyr::gather(key = "Columns", value = "Correlation",all_of(nme)) %>%
#       dplyr::mutate(Subtype = subtype)
#   },
#   clin = clin_neoadj,
#   nme
# )
#
# plot_corr <- bind_rows(plot_corr) %>%
#   dplyr::mutate(
#     Rows = factor(Rows, levels =  nme %>% rev()),
#     Columns = factor(Columns, levels =  nme %>% rev()),
#     Subtype = factor(Subtype, levels =  c("TN","HER2","HR","ALL")),
#     Cor_text = Correlation %>% round(digits = 1) %>% str_replace("0\\.", ".")
#   )
#
#
#
# p <- plot_corr %>%
#   ggplot(aes(x = Columns, y= Rows)) +
#   geom_raster(aes(fill = Correlation)) +
#   geom_text(aes(label = Cor_text), size = 2) +
#   scale_fill_gradient2(high = "#b2182b",# red
#                        low = "#2166ac",# blue
#                        limits = c(-1,1)) +
#   guides(fill = guide_colorbar(title = "Spearman\ncorrelation", title.vjust = 1)) +
#   facet_wrap(facets = ~Subtype , nrow = 2) +
#   theme(legend.position = "bottom",
#         axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust=.5),
#         axis.title = element_blank())
#
# pdf(file = "results/figures/geo_correlations.pdf",
#     width = 7.5, height = 8.5)
# print(p)
# dev.off()
#
#
#
# #
# # ==============================================================================
#
#





