# finher_her2_3yr_3year.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Prognostic, chemo interaction and herceptin interaction with TIL, TIL associated signatures and Immune celltypes.
# Account for arm heterogenity.
# Prepare data for plotting.
# Generate sample plot for HER2



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis, chemo interaction, and herceptin interaction in HER2.
# 3. Save Robjects
# 4. Explore non-PH in ddfs


# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")
load("results/data/expr_finher.RData") # for TILsig computation
load("results/data/tilsig_clean.RData") # for TILsig computation

# Update clin_finher with TILsig score
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sig <- tilsig_clean$ALL %>%
  dplyr::select(Direction, Ncbi_gene_id1)

# Compute module score and update clin_finher
score <- get_module_score(
  x = expr_finher %>% t_tibble(names_x_desc = "CEL_filename"),
  module_list = list(TILsig_imm = sig %>%
                       dplyr::filter(Direction == 1) %>%
                       dplyr::mutate(Direction = 1), # average imm
                     TILsig_fib = sig %>%
                       dplyr::filter(Direction == -1) %>%
                       dplyr::mutate(Direction = 1), # average fib
                     TILsig = sig), # weighted average
  by = "Ncbi_gene_id1"
) %>%
  dplyr::mutate(TILsig_imm = TILsig_imm %>% genefu::rescale(q=0.05),
                TILsig_fib = TILsig_fib %>% genefu::rescale(q=0.05),
                TILsig = TILsig %>% genefu::rescale(q=0.05))

clin_finher <- clin_finher %>% left_join(score, by = "CEL_filename")



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



# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

clin_finher %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(
    Subtype = Subtype_IHC, # Subtype_IHC_2
    Arm
  ) %>%
  dplyr::summarise(N = n(),
                   RFS = which(RFS_Event == 1) %>% length(),
                   DDFS = which(DDFS_Event == 1) %>% length(),
                   OS = which(OS_Event == 1) %>% length()) %>%
  dplyr::mutate(
    RFS_perc = (RFS/N)*100,
    DDFS_perc = (DDFS/N)*100,
    OS_perc = (OS/N)*100
  ) %>%
  as.data.frame()

#   Subtype            Arm  N RFS DDFS OS RFS_perc DDFS_perc   OS_perc
# 1 HR-HER2+ DTX.noTRA.noHR 27   7    7  5 25.92593  25.92593 18.518519
# 2 HR-HER2+ NVB.noTRA.noHR 19   6    6  6 31.57895  31.57895 31.578947
# 3 HR-HER2+   DTX.TRA.noHR 26   3    3  2 11.53846  11.53846  7.692308
# 4 HR-HER2+   NVB.TRA.noHR 17   6    5  3 35.29412  29.41176 17.647059
# 5 HR+HER2+   DTX.noTRA.HR 18   4    4  2 22.22222  22.22222 11.111111
# 6 HR+HER2+   NVB.noTRA.HR 27   7    6  3 25.92593  22.22222 11.111111
# 7 HR+HER2+     DTX.TRA.HR 16   2    0  0 12.50000   0.00000  0.000000
# 8 HR+HER2+     NVB.TRA.HR 30  10    9  4 33.33333  30.00000 13.333333

# Average size of strata, c(27,19,26,17,18,27,16,30) %>% mean() = 22.5

# Analysis
# Note that there are *eight* stratas in HER2 BC !!!!!!!!!!!!!!!
# 1. Sig's prognostic value in HER2 (with arm heterogenity; interaction with *eight* stratas)
# 2. Sig's interaction with chemo in HER2 (within each HR subtype and in pooled stratified HR)
# 3. Sig's interaction with herceptin in HER2 (within each HR subtype and in pooled stratified HR)

# Note: By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to clim a significant treatment
# interaction.



# Limiting followup timeto 3 year
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clin_finher3yr <- survSplit(formula = Surv(DDFS_Time_Years, DDFS_Event) ~ .,
                            data = clin_finher  %>%
                              dplyr::filter(Subtype_IHC_2 == "HER2"),
                            cut = c(3),
                            start = "DDFS_Time_Start", # New variable introduced
                            end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                            zero = 0,
                            episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id)) %>%
  dplyr::filter(Interval_Id == 1)



#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_her2_3yr <- list()

# HER2 prognosis
# >>>>>>>>>>>>>>

finher_her2_3yr[["prog"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig prognosis + heterogenity
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),

          function(sig, event_type, xdata){ # subtype,

            get_adj_prog(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_End"),
              biomarker = sig,
              xdata = xdata
            )

            # # Prognostic model
            # paste("Surv(", time_variable, ",", event_variable, ") ~",
            #           biomarker, "+ strata(Hormone, Herceptin, Chemo)")

            # # Heterogenity/interaction model
            # paste("Surv(", time_variable, ",", event_variable, ") ~",
            #           biomarker,"* interaction(Hormone, Herceptin, Chemo)")

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")
  # %>%
  #   dplyr::filter(Arm != "DTX.TRA.HR")

)

names(finher_her2_3yr[["prog"]]) <- "DDFS"

# There were 50 or more warnings (use warnings() to see the first 50)
# warnings()
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#    Loglik converged before variable  2 ; coefficient may be infinite.

# Note !!!!!!!!!
# Warnings arise in DDFS and OS due to the no-event strata, DTX.TRA.HR. (Tested manually)
# No warnings from RFS (all strata has some events)!!!!



# HER2 Chemo interaction; all HER2 + strata(Hormone, Herceptin)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
finher_her2_3yr[["inter.chemo"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))



        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Chemo + strata(Hormone, Herceptin)")
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Chemo + strata(Hormone, Herceptin)")


            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")

)

names(finher_her2_3yr[["inter.chemo"]]) <- "DDFS"




# HER2 Chemo interaction; notra HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_3yr[["inter.chemo.notra"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant

            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "noTRA")

)

names(finher_her2_3yr[["inter.chemo.notra"]]) <- "DDFS"





# HER2 Chemo interaction; tra HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_3yr[["inter.chemo.tra"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant

            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "TRA")

)

names(finher_her2_3yr[["inter.chemo.tra"]]) <- "DDFS"

# There were 38 warnings (use warnings() to see them)
#
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#    Loglik converged before variable  2 ; coefficient may be infinite.



# HER2 TRA interaction; all HER2 + strata(Hormone, Chemo)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
finher_her2_3yr[["inter.tra"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))



        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Herceptin + strata(Hormone, Chemo)")
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin + strata(Hormone, Chemo)")


            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")

)

names(finher_her2_3yr[["inter.tra"]]) <- "DDFS"




# HER2 TRA interaction; dtx HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_3yr[["inter.tra.dtx"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Herceptin + strata(Hormone, Chemo)") # Chemo is redundant
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin + strata(Hormone, Chemo)") # Chemo is redundant

            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "DTX")

)

names(finher_her2_3yr[["inter.tra.dtx"]]) <- "DDFS"

# There were 38 warnings (use warnings() to see them)
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#     Loglik converged before variable  2 ; coefficient may be infinite.



# HER2 TRA interaction; nvb HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_3yr[["inter.tra.nvb"]] <- purrr::map(

  "DDFS",

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      c("HR-HER2+", "HR+HER2+", "HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", "TILsig", "TILsig_imm", "TILsig_fib",
            "Immune1", "Immune2", "Immune3",
            "Interferon1", "Interferon2", "Interferon3",
            # "Cholesterol1",
            "Cholesterol2",
            # "Cholesterol3",
            "Fibrosis1", "Fibrosis2", "Fibrosis3",
            "Proliferation1", "Proliferation2", "Proliferation3",
            "Tcell", #"CLymphocyte", "Bcell",
            # "NKcell", "Monocyte", "MDendritic",
            "Fibroblast"),


          function(biomarker, event_type, xdata){

            event_variable = str_c(event_type, "_Event")
            time_variable = str_c(event_type, "_Time_End")


            chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "+ Herceptin + strata(Hormone, Chemo)") # Chemo is redundant
            chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin + strata(Hormone, Chemo)") # Chemo is redundant

            get_adj_inter(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              xdata = xdata,
              chr_null_formula = chr_null_formula,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect = purrr::map(yy,~(.x$main_effect)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin = clin_finher3yr %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "NVB")

)

names(finher_her2_3yr[["inter.tra.nvb"]]) <- "DDFS"

#
# ==============================================================================




# 3. Analysis summary
# ==============================================================================

summarize_her2 <- function(x, prog = F){
  # x <- finher_her2_3yr$prog$DDFS
  # prog = T

  x <- purrr::map(x,
                  function(xx, prog){

                    if(prog){

                      nme <- c("Variable", "HR", "Low95","Up95", "P", "P_adj", "Heterogeneity", "Ph_ok")
                      xx <- xx$main_effect[ , nme] %>%
                        dplyr::mutate( HR = str_c(round(HR, digits = 2),
                                                  " (",
                                                  round(Low95,digits=2), "-",
                                                  round(Up95,digits=2),
                                                  ")"),
                                       # Note, for prognostic model the CI is in HR
                                       Heterogeneity = str_c(round(Heterogeneity, digits=2))) %>%
                        dplyr::rename(Het = "Heterogeneity")

                    } else {

                      nme <- c("Variable", "HR", "Low95","Up95", "P", "P_adj", "Ph_ok")
                      xx <- xx$main_effect[ , nme] %>%
                        dplyr::mutate( HR = str_c(round(HR, digits = 2),
                                                  " (",
                                                  round(exp(Low95),digits=2), "-",
                                                  round(exp(Up95),digits=2),
                                                  ")"))
                      # Note, for interaction model the CI is in coef

                    }

                    xx <- xx %>%
                      dplyr::mutate(
                        P = str_c(round(P, digits=2)),
                        P_adj = str_c(round(P_adj, digits=2))
                      ) %>%
                      dplyr::select(-c(Low95,Up95))
                  },
                  prog)


  x <- x[c("HER2","HR-HER2+","HR+HER2+")]

  x <- purrr::map(names(x),
                  function(subtype,x){

                    prefix <- switch(subtype,
                                     "HR-HER2+" = "neg",
                                     "HR+HER2+" = "pos",
                                     "HER2" = "all")
                    names(x[[subtype]]) <- str_c(prefix, "_", names(x[[subtype]]))
                    x[[subtype]]
                  },
                  x)


  bind_cols(x) %>%
    dplyr::rename(Vaiable = "all_Variable") %>%
    dplyr::select(-c(neg_Variable,pos_Variable))

}

ddfs <- list()

ddfs[["prog"]] <- summarize_her2(x = finher_her2_3yr$prog$DDFS, prog = T)
ddfs[["inter.chemo"]] <- summarize_her2(x = finher_her2_3yr$inter.chemo$DDFS, prog = F)
ddfs[["inter.chemo.notra"]] <- summarize_her2(x = finher_her2_3yr$inter.chemo.notra$DDFS, prog = F)
ddfs[["inter.chemo.tra"]] <- summarize_her2(x = finher_her2_3yr$inter.chemo.tra$DDFS, prog = F)
ddfs[["inter.tra"]] <- summarize_her2(x = finher_her2_3yr$inter.tra$DDFS, prog = F)
ddfs[["inter.tra.dtx"]] <- summarize_her2(x = finher_her2_3yr$inter.tra.dtx$DDFS, prog = F)
ddfs[["inter.tra.nvb"]] <- summarize_her2(x = finher_her2_3yr$inter.tra.nvb$DDFS, prog = F)

write_xlsx(ddfs, path = str_c(out_tables,"finher_her2_3yr_ddfs_summary.xlsx"))


#
# ==============================================================================




# 4. Save Robjects
# ==============================================================================

save(finher_her2_3yr, file = str_c(out_data,"finher_her2_3yr.RData"))

#
# ==============================================================================




