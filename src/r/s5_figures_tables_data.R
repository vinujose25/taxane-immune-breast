# s5_figures_tables_data.R



# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Generate journal specific plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data
# 2. Prepare common dataframes to use.
# 3. Prognosis plot
# 4. Prognosis heterogenity plots
# 5. Interaction plot



# 1. Load data.
# ==============================================================================

load("results/data/clin_neoadj.RData") # clinical data celaned
load("results/data/prog_sum.RData") # prognosis summary
load("results/data/prog_het_sum.RData") # prognosis heterogenity summary

#
# ==============================================================================



# 2. Prepare common dataframes to use.
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




#
# ==============================================================================



# 3. Prognosis plot
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
      top = response_sum %>%
        dplyr::filter(Subtype_ihc == subtype) %>%
        dplyr::select(Subtype_label) %>%
        tibble::deframe(),
      fig.lab.pos = "bottom",
      fig.lab.size = 4
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



# 4. Prognosis heterogenity plots
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
        top = response_sum %>%
          dplyr::filter(Subtype_ihc == subtype) %>%
          dplyr::select(Subtype_label) %>%
          tibble::deframe(),
        fig.lab.pos = "bottom",
        fig.lab.size = 4
      )
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



# 5. Interaction plot
# ==============================================================================

#
# ==============================================================================
