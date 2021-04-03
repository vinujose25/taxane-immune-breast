# s5_figures_tables_data.R



# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Generate journal specific plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data
# 2. Gene module list as table output
# 3. Prepare common dataframes to use.
# 4. Prognosis plot
# 5. Prognosis heterogenity plots
# 6. Interaction plot




# 1. Load data.
# ==============================================================================

load("results/data/clin_neoadj.RData") # clinical data celaned

load("results/data/module_consolidated.RData") # original module
load("results/data/module_list_neoadj.RData") # module subset from neoadj genes
load("results/data/module_list_finher.RData") # module subset from finher genes

load("results/data/prog_sum.RData") # prognosis summary
load("results/data/prog_het_sum.RData") # prognosis heterogenity summary
load("results/data/inter_sum.RData") # interaction summary

#
# ==============================================================================



# 2. Gene module list as table output
# ==============================================================================

# Consolidated and cleaned gene modules
write_tsv(x = module_consolidated,
          path = str_c(out_tables, "module_consolidated.tsv"))

# module_list_neoadj
x <- purrr::map(
  names(module_list_neoadj),
  function(nme, lst){

    lst[[nme]] %>%
      dplyr::mutate(Module_id = nme)

  },
  lst = module_list_neoadj
) %>% bind_rows()

write_tsv(x = x,
          path = str_c(out_tables, "module_consolidated_neoadj.tsv"))


# module_list_finher
x <- purrr::map(
  names(module_list_finher),
  function(nme, lst){

    lst[[nme]] %>%
      dplyr::mutate(Module_id = nme)

  },
  lst = module_list_finher
) %>% bind_rows()

write_tsv(x = x,
          path = str_c(out_tables, "module_consolidated_finher.tsv"))


#
# ==============================================================================




# 3. Prepare common dataframes to use.
# ==============================================================================

# Get event summary per subtype
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

response_sum  <- clin_neoadj %>%
  dplyr::group_by(Subtype_ihc) %>%
  dplyr::summarise(N_response = which(Response == 1) %>% length(),
                   N_patients = n(),
                   Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2))

# adding row for All samples
response_sum  <- bind_rows(
  response_sum,
  tibble(
    Subtype_ihc = "ALL",
    N_response = response_sum$N_response %>% sum(),
    N_patients = response_sum$N_patients %>% sum(),
    Percent_repsonse = (N_response/N_patients) * 100
  )
)


# formatting
response_sum <- response_sum %>%
  dplyr::mutate(
    Subtype_ihc = factor(Subtype_ihc, levels = c("TN", "HER2", "HR", "ALL")),
    # Subtype_label = str_c(Subtype_ihc,
    #                       " (n.Sample: ", N_patients,
    #                       ", n.pCR: ", N_response," ~ ", Percent_repsonse %>% round(digits = 2), "%)")
    Subtype_label = str_c(Subtype_ihc,
                          "  [N: ", N_patients,
                          ", pCR: ", N_response," (", Percent_repsonse %>% round(digits = 2), "%)]")
  )





# Colors
# >>>>>>>>>>>>>>
# library("scales")

pdf(file = str_c(out_figures,"palatte_ggsci.pdf"), onefile = T)

# par(mfrow=c(2,2)) # by row

show_col(pal_npg("nrc")(10))
title(main = "NPG: palatte type - nrc")

show_col(pal_aaas("default")(10))
title(main = "AAAS: palatte type - default")

show_col(pal_nejm("default")(8))
title(main = "NEJM: palatte type - default")

show_col(pal_lancet("lanonc")(9))
title(main = "Lancet: palatte type - lanonc")

show_col(pal_jama("default")(7))
title(main = "JAMA: palatte type - default")

show_col(pal_jco("default")(10))
title(main = "JCO: palatte type - default")

show_col(pal_ucscgb("default")(26))
title(main = "UCSCGB: palatte type - default")

show_col(pal_d3("category10")(10))
title(main = "D3: palatte type - category10")
show_col(pal_d3("category20")(20))
title(main = "D3: palatte type - category20")
show_col(pal_d3("category20b")(20))
title(main = "D3: palatte type - category20b")
show_col(pal_d3("category20c")(20))
title(main = "D3: palatte type - category20c")

show_col(pal_locuszoom("default")(7))
title(main = "LocusZoom: palatte type - default")

show_col(pal_igv("default")(51))
title(main = "IGV: palatte type - default")
show_col(pal_igv("alternating")(2))
title(main = "IGV: palatte type - alternating")

show_col(pal_uchicago("default")(9))
title(main = "UChicago: palatte type - default")
show_col(pal_uchicago("light")(9))
title(main = "UChicago: palatte type - light")
show_col(pal_uchicago("dark")(9))
title(main = "UChicago: palatte type - dark")

show_col(pal_startrek("uniform")(7))
title(main = "Star Trek: palatte type - uniform")

show_col(pal_tron("legacy")(7))
title(main = "Tron Legacy: palatte type - legacy")

show_col(pal_futurama("planetexpress")(12))
title(main = "Futurama: palatte type - planetexpress")

show_col(pal_rickandmorty("schwifty")(12))
title(main = "Rick and Morty: palatte type - schwifty")

show_col(pal_simpsons("springfield")(16))
title(main = "The Simpsons: palatte type - springfield")

show_col(pal_gsea("default")(12))
title(main = "GSEA: palatte type - default")

show_col(pal_material("red")(10))
title(main = "Material Design: palatte type - red")
show_col(pal_material("pink")(10))
title(main = "Material Design: palatte type - pink")
show_col(pal_material("purple")(10))
title(main = "Material Design: palatte type - purple")
show_col(pal_material("deep-purple")(10))
title(main = "Material Design: palatte type - deep-purple")
show_col(pal_material("indigo")(10))
title(main = "Material Design: palatte type - indigo")
show_col(pal_material("blue")(10))
title(main = "Material Design: palatte type - blue")
show_col(pal_material("light-blue")(10))
title(main = "Material Design: palatte type - light-blue")
show_col(pal_material("cyan")(10))
title(main = "Material Design: palatte type - cyan")
show_col(pal_material("teal")(10))
title(main = "Material Design: palatte type - teal")
show_col(pal_material("green")(10))
title(main = "Material Design: palatte type - green")
show_col(pal_material("light-green")(10))
title(main = "Material Design: palatte type - light-green")
show_col(pal_material("lime")(10))
title(main = "Material Design: palatte type - lime")
show_col(pal_material("yellow")(10))
title(main = "Material Design: palatte type - yellow")
show_col(pal_material("amber")(10))
title(main = "Material Design: palatte type - amber")
show_col(pal_material("orange")(10))
title(main = "Material Design: palatte type - orange")
show_col(pal_material("deep-orange")(10))
title(main = "Material Design: palatte type - deep-orange")
show_col(pal_material("brown")(10))
title(main = "Material Design: palatte type - brown")
show_col(pal_material("grey")(10))
title(main = "Material Design: palatte type - gray")
show_col(pal_material("blue-grey")(10))
title(main = "Material Design: palatte type - blue-gray", cex = .75)

dev.off()


custom_colours_set1 <- c(

  # # nejm colurs from ggsci package
  "#BC3C29FF", # red for AAA+T
  "#0072B5FF", # blue for AAA
  "#E18727FF", # orange for AAA+T+TRA

  "#20854EFF", # green for HR
  "#7876B1FF", # purple for HER2
  "#EE4C97FF", # majenta for TN

  # # aaas colurs from ggsci package
  "#808180FF" # gray
  # "#1B1919FF"  # black for A0A+T
)

custom_colours_set2 <- c(

  # # Material Design: palatte type − blue−gray, colurs from ggsci package
  "#90A4ADFF", # light-gray for AAA+T
  "#536D79FF", # gray for AAA
  "#263238FF", # dark-gray for AAA+T+TRA

  "#90A4ADFF", # light-gray for HR
  "#536D79FF", # gray for HER2
  "#263238FF", # dark-gray for TN

  # # nejm colurs from ggsci package
  "#EE4C97FF" # majenta for A0A+T
)


custom_colours_set3 <- c(

  # # AAAS: palatte type − default, colurs from ggsci package
  "#1B1919FF", # black for AAA+T
  "#BB0021FF", # red for AAA
  "#008280FF", # cyan for AAA+T+TRA

  "#008B45FF", # green for HR
  "#631879FF", # purple for HER2
  "#A20056FF", # majenta for TN

  "#5F559BFF" # blue for A0A+T
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

names(custom_colours_set1) <- nme
names(custom_colours_set2) <- nme
names(custom_colours_set3) <- nme

#
# ==============================================================================



# 4. Prognosis plot
# ==============================================================================


# 1. Consolidate subtype specific prognostic summary
# 2. Format consolidated data for forest plot and plot annotation

  prog_sum_merged <-

    # preparing prog_sum for consolidation
    purrr::map(

      names(prog_sum), # per subtype

      function(nme, x){

        x[[nme]] %>%
          dplyr::mutate(Subtype = nme) # updating subtype info in each per subtype summary
      },

      x = prog_sum
    ) %>%

    # Consolidation
    bind_rows() %>%

    # integrate with response summary
    left_join(response_sum %>%
                dplyr::rename(Subtype = "Subtype_ihc"),
              by = "Subtype") %>%

    dplyr::mutate(

      # Formatting prognosis annotation
      Text_or = str_c(
        Estimate %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
        " (",
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
        "/",
        if_else(
          P_adj < .001,
          "<.001",
          formatC(x = P_adj, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        )
      ),
      Text_het = if_else(
        is.na(Q) & is.na(Q_p) & is.na(I2),
        "",
        str_c(
          Q %>% formatC(digits = 1, width = 1, format = "f"),
          "/",
          if_else(
            Q_p < .001,
            "<.001",
            formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
              str_replace("0\\.", ".") # str_replace("0.", " .")
          ),
          "/",
          I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
        )
      ),

      Module_name = factor(Module_name, levels = c(
        "Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        "Cholesterol1", "Cholesterol2", "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", "CLymphocyte", "Bcell",
        "NKcell", "Monocyte", "MDendritic",
        "Fibroblast") %>% rev()
      ),

      # grouping modules
      Module_class = case_when(
        str_detect(Module_name, "Immune") ~ "Immune",
        str_detect(Module_name, "Interferon") ~ "Interferon",
        str_detect(Module_name, "Cholesterol") ~ "Cholesterol",
        str_detect(Module_name, "Fibrosis") ~ "Fibrosis",
        str_detect(Module_name, "Proliferation") ~ "Proliferation",
        TRUE ~ "Celltype",
      ) %>%
        factor(levels = c(
          "Immune", "Interferon", "Cholesterol", "Fibrosis",
          "Proliferation", "Celltype")
        )

    )


  # forest plotting
  p1 <- p2 <- list() # template

  for(subtype in c("TN", "HER2", "HR", "ALL")){

    print(subtype)

    # forest plot
    p1[[subtype]] <- prog_sum_merged %>%
      dplyr::filter(Subtype == subtype) %>%
      ggplot() +
      geom_point(aes(x = Estimate, y = Module_name)) +
      geom_errorbarh(aes(xmin = Lower_ci, xmax = Upper_ci, y = Module_name), height = .1) +
      facet_grid(facets = Module_class ~ Subtype_label, switch = "y",
                 scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            panel.grid = element_blank()
      )


    # forest plot annotation
    p2[[subtype]] <- prog_sum_merged %>%
      tidyr::gather("key", "value", "Text_or","Text_p","Text_het") %>%
      dplyr::filter(Subtype == subtype) %>%
      dplyr::mutate(
        key = purrr::map_chr(key, ~switch(
          .x,
          Text_or = "Log-OR (95%CI)", # "Log-OR (95% CI)"
          Text_p = "P/P-adj", # "Pval (Pval-adj)"
          Text_het = "Q/Q-P/I^2" # "Cochran's Q (Q-Pval), I^2"
        )),
        key = factor(key, levels = c("Log-OR (95%CI)",
                                     "P/P-adj",
                                     "Q/Q-P/I^2"))
      ) %>%
      ggplot() +
      geom_text(aes(x = key, y = Module_name, label = value)) +
      facet_grid(facets = Module_class ~ Subtype_label, switch = "y",
                 scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.grid = element_blank()
      )
  }



# Printing plot
pdf(file = str_c(out_figures, "forest_prognosis.pdf"),
    width = 7.5, height = 8.5, onefile = T)

for(subtype in c("TN", "HER2", "HR", "ALL")){

  print(subtype)

  p <- ggpubr::ggarrange(
    p1[[subtype]] + labs(x = "Log-OR") + geom_vline(xintercept = 0, linetype = "dashed") ,
    p2[[subtype]] + labs(x = "") + theme(axis.text.x.bottom = element_text(color = "black")),
    ncol = 2,
    nrow = 1,
    widths = c(.375, .625)
  ) %>%
  # Seting a global title for arranged plot
  ggpubr::annotate_figure(
    top = text_grob(label = response_sum %>%
                      dplyr::filter(Subtype_ihc == subtype) %>%
                      dplyr::select(Subtype_label) %>%
                      tibble::deframe(),
                    face = "bold", size = 10),
    fig.lab.pos = "bottom"#,
    #fig.lab.size = 3
  )
  print(p)
}

dev.off()


# # Select individual pdf pages
# pdftools::pdf_subset(
#   input = "plot_prognosis.pdf",
#   pages = 2,
#   output = "plot_prognosis_her2.pdf"
# )

#
# ==============================================================================



# 5. Prognosis heterogenity plots
# ==============================================================================


# forest plotting
p1 <- p2 <- list()  # template


# Per-subtype, per-sig study+arm (Strata) heterogenity exploraion

for(subtype in c("TN", "HER2", "HR", "ALL")){

  print(subtype)

  # To plot similar sigs in a single page
  for(module_class in c("Immune", "Interferon","Cholesterol","Fibrosis",
                        "Proliferation","Celltype1","Celltype2")) {

    print(module_class)
    # subtype = "TN"
    # module_class = "Immune"

    # grouping of similar sigs
    module_names <- switch (
      module_class,
      Immune = c("Immune1", "Immune2", "Immune3"),
      Interferon = c("Interferon1", "Interferon2", "Interferon3"),
      Cholesterol = c("Cholesterol1", "Cholesterol2", "Cholesterol3"),
      Fibrosis = c("Fibrosis1", "Fibrosis2", "Fibrosis3"),
      Proliferation = c("Proliferation1", "Proliferation2", "Proliferation3"),
      Celltype1 = c("Tcell", "CLymphocyte", "Bcell", "NKcell"),
      Celltype2 = c("Monocyte", "MDendritic", "Fibroblast")
    )

    # Consolidate prognosis summary of similar sigs for plotting
    x <- prog_het_sum[[subtype]][module_names] %>%
      bind_rows()
    # x <- bind_rows(x)


    # Format consolidated prognosis summary for plotting
    x <- x %>%
      dplyr::mutate(

        # Formatting prognosis annotation
        Text_or = str_c(
          Estimate %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
          " (",
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
          "/",
          if_else(
            P_adj < .001,
            "<.001",
            formatC(x = P_adj, digits = 3, width = 1, format = "f") %>%
              str_replace("0\\.", ".") # str_replace("0.", " .")
          )
        ),
        Text_het = if_else(
          is.na(Q) & is.na(Q_p) & is.na(I2),
          "",
          str_c(
            Q %>% formatC(digits = 1, width = 1, format = "f"),
            "/",
            if_else(
              Q_p < .001,
              "<.001",
              formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
                str_replace("0\\.", ".") # str_replace("0.", " .")
            ),
            "/",
            I2 %>% formatC(digits = 1, width = 1, format = "f") # , flag = " "
          )
        ),
        Text_size = str_c(
          N_patients,
          "/",
          N_response
        ),

        Module_name = factor(Module_name, levels = module_names),

        Subtype = subtype # for facetting
      )


    # forest plot
    p1[[subtype]][[module_class]] <- x %>%
      ggplot() +
      geom_point(aes(x = Estimate, y = Strata)) +
      geom_errorbarh(aes(xmin = Lower_ci, xmax = Upper_ci, y = Strata), height = .1) +
      facet_grid(facets = Module_name ~ Subtype, switch = "y",
                 scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            panel.grid = element_blank()
      )

    # forest plot annotation
    p2[[subtype]][[module_class]]  <- x %>%
      tidyr::gather("key", "value", "Text_size", "Text_or","Text_p", "Text_het") %>%
      dplyr::mutate(
        key = purrr::map_chr(key, ~switch(
          .x,
          Text_size = "N/pCR", # |pCR%
          Text_or = "Log-OR (95%CI)", # "Log-OR (95% CI)"
          Text_p = "P/P-adj", # "Pval (Pval-adj)"
          Text_het = "Q/Q-P/I^2" # "Cochran's Q (Q-Pval), I^2"
        )),
        key = factor(key, levels = c("N/pCR",
                                     "Log-OR (95%CI)",
                                     "P/P-adj",
                                     "Q/Q-P/I^2"))
      ) %>%
      ggplot() +
      geom_text(aes(x = key, y = Strata, label = value)) +
      facet_grid(facets = Module_name ~ Subtype, switch = "y",
                 scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.grid = element_blank()
      )

  }

}



# printing plot
pdf(file = str_c(out_figures, "forest_prognosis_hetrogenity.pdf"),
    width = 7.5, height = 8.5, onefile = T)

for(subtype in c("TN", "HER2", "HR", "ALL")){

  print(subtype)

  for(module_class in c("Immune", "Interferon","Cholesterol","Fibrosis",
                        "Proliferation","Celltype1","Celltype2")){

    # subtype = "ALL"
    # module_class = "Celltype1"

    print(module_class)

    p <- ggpubr::ggarrange(
      p1[[subtype]][[module_class]] + labs(x = "Log-OR") +
        geom_vline(xintercept = 0, linetype = "dashed"),
      p2[[subtype]][[module_class]] + labs(x = "") +
        theme(axis.text.x.bottom = element_text(color = "black")),
      ncol = 2,
      nrow = 1,
      widths = c(.375, .625)
    ) %>%
      ggpubr::annotate_figure(
        top = text_grob(label = response_sum %>%
                          dplyr::filter(Subtype_ihc == subtype) %>%
                          dplyr::select(Subtype_label) %>%
                          tibble::deframe(),
                        face = "bold", size = 10),
        fig.lab.pos = "bottom"#,
        #fig.lab.size = 3
      )
      # ggpubr::annotate_figure(
      #   top = response_sum %>%
      #     dplyr::filter(Subtype_ihc == subtype) %>%
      #     dplyr::select(Subtype_label) %>%
      #     tibble::deframe(),
      #   fig.lab.pos = "bottom",
      #   fig.lab.size = 4
      # )
    print(p)
  }
}

dev.off()


# # Select individual pdf pages
# pdftools::pdf_subset(
#   input = "plot_prognosis.pdf",
#   pages = 2,
#   output = "plot_prognosis_her2.pdf"
# )


#
# ==============================================================================



# 6. Interaction plot
# ==============================================================================

# custom_colors <- c(
#   "#e41a1c", # red
#   "#377eb8", # blue
#   "#4daf4a", # green
#   "#984ea3", # purple
#   "#ff7f00"  # orange
# )
# Ref : https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#:~:text=In%20this%20case%2C%20we'll,that%20are%20color%2Dblind%20friendly.


# 1. Consolidate subtype-analysis specific interaction summary
# 2. Format consolidated data for forest plot and plot annotation


inter_sum_merged <-

  # preparing inter_sum for consolidation
  purrr::map(

    names(inter_sum), # per subtype

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
      )%>%

        # Per subtype, analysis level consolidation
        bind_rows()

    },

    x = inter_sum
  ) %>%

  # Subtype level consolidation
  bind_rows()




inter_sum_merged <-

  inter_sum_merged %>%
  # subsettting interaction terms
  dplyr::filter(str_detect(Module_name, ":")) %>%

  # # integrate with response summary
  # left_join(response_sum %>%
  #             dplyr::rename(Subtype = "Subtype_ihc"),
  #           by = "Subtype") %>%

  dplyr::mutate(

    # format sample size and pcr per arm
    Text_arm = str_c(
      Interaction_var_levels,
      ": ",
      "N-", N_patients,
      ", pCR-", N_response,"(", Percent_repsonse %>% round(digits = 2), "%)"
    ),

    # Formatting prognosis annotation
    Text_or = str_c(
      Estimate %>% formatC(digits = 1, width = 1, format = "f"), # , flag = "  "
      " (",
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
      "/",
      if_else(
        P_inter_adj < .001,
        "<.001",
        formatC(x = P_inter_adj, digits = 3, width = 1, format = "f") %>%
          str_replace("0\\.", ".") # str_replace("0.", " .")
      )
    ),
    Text_het = if_else(
      is.na(Q) & is.na(Q_p) & is.na(I2),
      "",
      str_c(
        Q %>% formatC(digits = 1, width = 1, format = "f"),
        "/",
        if_else(
          Q_p < .001,
          "<.001",
          formatC(x = Q_p, digits = 3, width = 1, format = "f") %>%
            str_replace("0\\.", ".") # str_replace("0.", " .")
        ),
        "/",
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
      factor(levels = c(
      "Immune1", "Immune2", "Immune3",
      "Interferon1", "Interferon2", "Interferon3",
      "Cholesterol1", "Cholesterol2", "Cholesterol3",
      "Fibrosis1", "Fibrosis2", "Fibrosis3",
      "Proliferation1", "Proliferation2", "Proliferation3",
      "Tcell", "CLymphocyte", "Bcell",
      "NKcell", "Monocyte", "MDendritic",
      "Fibroblast") %>% rev()
    ),

    # grouping modules
    Module_class = case_when(
      str_detect(Module_name, "Immune") ~ "Immune",
      str_detect(Module_name, "Interferon") ~ "Interferon",
      str_detect(Module_name, "Cholesterol") ~ "Cholesterol",
      str_detect(Module_name, "Fibrosis") ~ "Fibrosis",
      str_detect(Module_name, "Proliferation") ~ "Proliferation",
      TRUE ~ "Celltype",
    ) %>%
      factor(levels = c(
        "Immune", "Interferon", "Cholesterol", "Fibrosis",
        "Proliferation", "Celltype")
      ),

    Analysis_group = str_c(Subtype, ":", Analysis)

  )


# response_sum_inter <- inter_sum_merged %>%
#   dplyr::group_by(Analysis_group, Analysis, Interaction_var_levels, Text_arm) %>%
#   dplyr::summarise() %>%
#   dplyr::mutate(Subtype = str_split_fixed(Analysis_group,":",2)[,1])

response_sum_inter <- inter_sum_merged %>%
  dplyr::group_by(Analysis_group) %>%
  dplyr::summarise(Title = str_c(Text_arm %>% unique(), collapse = " | ")) %>%
  dplyr::mutate(
    Analysis_group2 = if_else(str_detect(Analysis_group, "ALL"),
                              str_split_fixed(str_split_fixed(Analysis_group, ":", 2)[,2],
                                              "\\.", 2)[, 1],
                              str_split_fixed(Analysis_group,":",2)[, 1]),
    Analysis_group2 = Analysis_group2 %>% str_replace_all("_plus_","+"),
    Title = str_c(
      Analysis_group2,
      "  [", Title, "]"
    )
  )


# forest plotting
p1 <- p2 <- list() # template

for(anagroup in inter_sum_merged$Analysis_group %>% unique()){

  # anagroup = "ALL:AAA_plus_T.TN_or_HER2_or_HR"
  print(anagroup)

  # forest plot
  p1[[anagroup]] <- inter_sum_merged %>%
    dplyr::filter(Analysis_group == anagroup) %>%
    ggplot() +
    geom_point(aes(y = Estimate, x = Module_name, color = Arm_or_subtype, group = Arm_or_subtype),
               position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = l95, ymax = u95, x = Module_name, color = Arm_or_subtype),
                  position = position_dodge(width = 1), width = .5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # scale_color_lancet() +
    # scale_color_manual(values = custom_colours_set1) +
    # scale_color_manual(values = custom_colours_set2) +
    scale_color_manual(values = custom_colours_set3) +
    # facet_grid(facets = Module_class ~ Subtype_label, switch = "y",
    #            scales = "free_y", space = "free_y") +
    facet_grid(facets = Module_class ~ Analysis_group, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          panel.grid = element_blank()) +
    labs(y = "Log-OR") +
    coord_flip()



  # forest plot annotation
  if (str_detect(anagroup, "ALL")) {
    xx <- inter_sum_merged %>%
      dplyr::filter(Analysis_group == anagroup) %>%
      dplyr::mutate(Subtype = Arm_or_subtype) %>%
      # !!! Inter_sum_merged already contains a subtype variable
      dplyr::rename("Log-OR (95%CI)" = "Text_or",
                    "P/P-adj" = "Text_p",
                    "Q/Q-P/I^2" = "Text_het") %>%
      tidyr::gather("key", "value",
                    "Log-OR (95%CI)", "P/P-adj", "Q/Q-P/I^2", "Subtype") %>%
      dplyr::mutate(key = key %>%
                      factor(levels = c("Subtype", "Log-OR (95%CI)", "P/P-adj", "Q/Q-P/I^2")))

  } else {
    xx <- inter_sum_merged %>%
      dplyr::filter(Analysis_group == anagroup) %>%
      dplyr::mutate(Arm = Arm_or_subtype) %>%
      dplyr::rename("Log-OR (95%CI)" = "Text_or",
                    "P/P-adj" = "Text_p",
                    "Q/Q-P/I^2" = "Text_het") %>%
      tidyr::gather("key", "value",
                    "Log-OR (95%CI)", "P/P-adj", "Q/Q-P/I^2", "Arm") %>%
      dplyr::mutate(key = key %>%
                      factor(levels = c("Arm", "Log-OR (95%CI)", "P/P-adj", "Q/Q-P/I^2")))

  }


  p2[[anagroup]] <- xx %>%
    ggplot() +
    geom_text(aes(y = key, x = Module_name, label = value,
                  color = Arm_or_subtype, group = Arm_or_subtype),
              position = position_dodge(width = 1),
              size = 3) +
    # scale_color_lancet() +
    # scale_color_manual(values = custom_colours_set1) +
    # scale_color_manual(values = custom_colours_set2) +
    scale_color_manual(values = custom_colours_set3) +
    # facet_grid(facets = Module_class ~ Subtype_label, switch = "y",
    #            scales = "free_y", space = "free_y") +
    facet_grid(facets = Module_class ~ Analysis_group, switch = "y",
               scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank()) +
    labs(y = "") +
    coord_flip()
}



# Printing plot
pdf(file = str_c(out_figures, "forest_interaction.pdf"),
    width = 7.5, height = 8.5, onefile = T)

for(anagroup in inter_sum_merged$Analysis_group %>% unique()){

  print(anagroup)

  p <- ggpubr::ggarrange(
    p1[[anagroup]] + guides(color = "none"), #  geom_vline(xintercept = 0, linetype = "dashed") +
    p2[[anagroup]] + theme(axis.text.x.bottom = element_text(color = "black")) + guides(color = "none"),
    # p1 + geom_vline(xintercept = 0, linetype = "dashed") + guides(color = "none") ,
    # p2 + theme(axis.text.x.bottom = element_text(color = "black")) + guides(color = "none"),
    ncol = 2,
    nrow = 1,
    widths = c(.375, .625)
  ) %>%
    # Seting a global title for arranged plot
    ggpubr::annotate_figure(
      top = text_grob(label = response_sum_inter %>%
                        dplyr::filter(Analysis_group == anagroup) %>%
                        dplyr::select(Title) %>%
                        tibble::deframe(),
                      face = "bold", size = 10),
      fig.lab.pos = "bottom"#,
      #fig.lab.size = 3
    )
  print(p)
}

dev.off()


# # Select individual pdf pages
# pdftools::pdf_subset(
#   input = "plot_prognosis.pdf",
#   pages = 2,
#   output = "plot_prognosis_her2.pdf"
# )

#
# ==============================================================================



