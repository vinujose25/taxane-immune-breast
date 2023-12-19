# finher_her2_split_survsplit.R


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
# load("results/data/expr_finher.RData") # for TILsig computation
# load("results/data/tilsig_clean.RData") # for TILsig computation
#
# # Update clin_finher with TILsig score
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sig <- tilsig_clean$ALL %>%
#   dplyr::select(Direction, Ncbi_gene_id1)
#
# # Compute module score and update clin_finher
# score <- get_module_score(
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


# SurvSplit
# >>>>>>>>>
clin_finher_split_ddfs <- survSplit(formula = Surv(DDFS_Time_Years, DDFS_Event) ~ .,
                               data = clin_finher  %>%
                                 dplyr::filter(Subtype_IHC_2 == "HER2"),
                               cut = c(3),
                               start = "DDFS_Time_Start", # New variable introduced
                               end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                               zero = 0,
                               episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))

clin_finher_split_os <- survSplit(formula = Surv(OS_Time_Years, OS_Event) ~ .,
                                    data = clin_finher  %>%
                                      dplyr::filter(Subtype_IHC_2 == "HER2"),
                                    cut = c(3),
                                    start = "OS_Time_Start", # New variable introduced
                                    end = "OS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                                    zero = 0,
                                    episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))

clin_finher_split_rfs <- survSplit(formula = Surv(RFS_Time_Years, RFS_Event) ~ .,
                                  data = clin_finher  %>%
                                    dplyr::filter(Subtype_IHC_2 == "HER2"),
                                  cut = c(3),
                                  start = "RFS_Time_Start", # New variable introduced
                                  end = "RFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                                  zero = 0,
                                  episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))

#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_her2_split <- list()

# HER2 prognosis
# >>>>>>>>>>>>>>

finher_her2_split[["prog"]] <- purrr::map(

  c("DDFS","OS", "RFS"),

  function(event_type, clin_ddfs, clin_os, clin_rfs){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os,
                   "RFS" = clin_rfs)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig prognosis + heterogenity
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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

            get_adj_prog2(
              event_variable = str_c(event_type, "_Event"),
              time_variable0 = str_c(event_type, "_Time_Start"),
              time_variable1 = str_c(event_type, "_Time_End"),
              time_split_id = "Interval_Id",
              biomarker = sig,
              xdata = xdata
            )

          },

          event_type,
          xdata
        )


        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2"),

  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2"),

  clin_rfs = clin_finher_split_rfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")

)

names(finher_her2_split[["prog"]]) <- c("DDFS","OS","RFS")

# There were 50 or more warnings (use warnings() to see the first 50)
# warnings()
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#    Loglik converged before variable  2 ; beta may be infinite.
#   4: In agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  ... :
#                   Ran out of iterations and did not converge



# HER2 Chemo interaction; all HER2 + strata(Hormone, Herceptin)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
finher_her2_split[["inter.chemo"]] <- purrr::map(

  c("DDFS","OS", "RFS"),

  function(event_type, clin_ddfs, clin_os, clin_rfs){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os,
                   "RFS" = clin_rfs)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))



        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"

            # chr_null_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
            #                           biomarker, "+ Chemo + strata(Hormone, Herceptin)")
            # chr_full_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
            #                           biomarker, "* Chemo + strata(Hormone, Herceptin)")
            chr_full_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Chemo *", time_split_id,
                                      "-", time_split_id, "+ strata(Hormone, Herceptin)")


            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )


          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2"),

  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2"),

  clin_rfs = clin_finher_split_rfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")

)

names(finher_her2_split[["inter.chemo"]]) <- c("DDFS","OS","RFS")
# Warning messages:
#   1: In sqrt(sum(x[c(idx2[1], i), "se(coef)"]^2) + (2 * vcov(m1)[idx2[1],  :
#       NaNs produced
#   2: In sqrt(sum(x[c(idx2[1], i), "se(coef)"]^2) + (2 * vcov(m1)[idx2[1],  :
#       NaNs produced
#   3: In sqrt(sum(x[c(idx2[1], i), "se(coef)"]^2) + (2 * vcov(m1)[idx2[1],  :
#       NaNs produced


# Just consider relevant data for mta til paper
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

finher_her2_split_relevant = finher_her2_split

save(finher_her2_split_relevant, file = str_c(out_data,"finher_her2_split_relevant.RData"))


ddfs <- list()

ddfs[["prog"]] <- summarize_her2_2(x = finher_her2_split_relevant$prog$DDFS$HER2, prog = T)
ddfs[["inter.chemo"]] <- summarize_her2_2(x = finher_her2_split_relevant$inter.chemo$DDFS$HER2, prog = F)

write_xlsx(ddfs, path = str_c(out_tables,"finher_her2_split_relevant_ddfs_summary.xlsx"))

os <- list()

os[["prog"]] <- summarize_her2_2(x = finher_her2_split_relevant$prog$OS$HER2, prog = T)
os[["inter.chemo"]] <- summarize_her2_2(x = finher_her2_split_relevant$inter.chemo$OS$HER2, prog = F)

write_xlsx(os, path = str_c(out_tables,"finher_her2_split_relevant_os_summary.xlsx"))

rfs <- list()

rfs[["prog"]] <- summarize_her2_2(x = finher_her2_split_relevant$prog$RFS$HER2, prog = T)
rfs[["inter.chemo"]] <- summarize_her2_2(x = finher_her2_split_relevant$inter.chemo$RFS$HER2, prog = F)

write_xlsx(rfs, path = str_c(out_tables,"finher_her2_split_relevant_rfs_summary.xlsx"))

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





# HER2 Chemo interaction; notra HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_split[["inter.chemo.notra"]] <- purrr::map(

  c("DDFS","OS"),

  function(event_type, clin_ddfs, clin_os){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"

            # chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
            #                           biomarker, "+ Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant
            chr_full_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Chemo *", time_split_id,
                                      "-",  time_split_id, "+ strata(Hormone, Herceptin)") # Herceptin is redundant


            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "noTRA"),
  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "noTRA")

)

names(finher_her2_split[["inter.chemo.notra"]]) <- c("DDFS","OS")
# Warning message:
# In sqrt(sum(x[c(idx2[1], i), "se(coef)"]^2) + (2 * vcov(m1)[idx2[1],  :
#  NaNs produced




# HER2 Chemo interaction; tra HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_split[["inter.chemo.tra"]] <- purrr::map(

  c("DDFS","OS"),

  function(event_type, clin_ddfs, clin_os){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"


            # chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
            #                           biomarker, "+ Chemo + strata(Hormone, Herceptin)") # Herceptin is redundant
            chr_full_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Chemo *", time_split_id,
                                      "-",  time_split_id, "+ strata(Hormone, Herceptin)") # Herceptin is redundant

            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Chemo",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "TRA"),
  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Herceptin == "TRA")

)

names(finher_her2_split[["inter.chemo.tra"]]) <- c("DDFS","OS")

# There were 50 or more warnings (use warnings() to see the first 50)
#
# Warning messages:
#   1: In fitter(X, Y, istrat, offset, init, control, weights = weights,  ... :
#    Loglik converged before variable  2 ; coefficient may be infinite.
#   2: In sqrt(sum(x[c(idx2[1], i), "se(coef)"]^2) + (2 * vcov(m1)[idx2[1],  ... :
#     NaNs produced





# HER2 TRA interaction; all HER2 + strata(Hormone, Chemo)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
finher_her2_split[["inter.tra"]] <- purrr::map(

  c("DDFS","OS"),

  function(event_type, clin_ddfs, clin_os){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))



        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"


            # chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
            #                           biomarker, "+ Herceptin + strata(Hormone, Chemo)")
            chr_full_formula <- paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin *", time_split_id,
                                      "-", time_split_id,"+ strata(Hormone, Chemo)")


            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2"),
  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")

)

names(finher_her2_split[["inter.tra"]]) <- c("DDFS","OS")




# HER2 TRA interaction; dtx HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_split[["inter.tra.dtx"]] <- purrr::map(

  c("DDFS","OS"),

  function(event_type, clin_ddfs, clin_os){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"

            # chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
            #                           biomarker, "+ Herceptin + strata(Hormone, Chemo)") # Chemo is redundant
            chr_full_formula <- paste("Surv(", time_variable0, ",",  time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin *", time_split_id,
                                      "-", time_split_id, "+ strata(Hormone, Chemo)") # Chemo is redundant

            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "DTX"),
  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "DTX")

)

names(finher_her2_split[["inter.tra.dtx"]]) <- c("DDFS","OS")

# There were 24 warnings (use warnings() to see them)
# Warning message:
#  In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#   Loglik converged before variable  5 ; beta may be infinite. te.



# HER2 TRA interaction; nvb HER2 + strata(Hormone)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_her2_split[["inter.tra.nvb"]] <- purrr::map(

  c("DDFS","OS"),

  function(event_type, clin_ddfs, clin_os){

    print(event_type)

    clin <- switch(event_type,
                   "DDFS" = clin_ddfs,
                   "OS" = clin_os)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      "HER2",

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig mta interaction
        yy <- purrr::map(

          c("TIL", #"TILsig", "TILsig_imm", "TILsig_fib",
            "TILsig_scaled",
            "TILsig_APP_Fc", "TILsig_Immune", "TILsig_IFNg", "TILsig_Innate",
            "TILsig_ECM", "TILsig_Adhesion",
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
            time_variable0 = str_c(event_type, "_Time_Start")
            time_variable1 = str_c(event_type, "_Time_End")
            time_split_id = "Interval_Id"


            # chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
            #                           biomarker, "+ Herceptin + strata(Hormone, Chemo)") # Chemo is redundant
            chr_full_formula <- paste("Surv(", time_variable0, ",",  time_variable1, ",", event_variable, ") ~",
                                      biomarker, "* Herceptin *", time_split_id,
                                      "-", time_split_id, "+ strata(Hormone, Chemo)") # Chemo is redundant

            get_adj_inter2(
              event_variable = event_variable, # to extract event summary
              biomarker = biomarker,
              interaction_variable = "Herceptin",
              time_split_id = time_split_id,
              xdata = xdata,
              chr_full_formula = chr_full_formula
            )

          },
          event_type,
          xdata
        )



        # Consolidate per sig stat and format
        list(
          main_effect1 = purrr::map(yy,~(.x$main_effect1)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = subtype
            ),

          main_effect2 = purrr::map(yy,~(.x$main_effect2)) %>%
            bind_rows() %>%
            dplyr::mutate(
              # P_adj = p.adjust(p = P, method = "BH"),
              # Discarding TIL from FDR coorection
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
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

    names(xx) <- "HER2" #c("HR-HER2+", "HR+HER2+", "HER2")

    xx
  },

  clin_ddfs = clin_finher_split_ddfs %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "NVB"),
  clin_os = clin_finher_split_os %>%
    dplyr::filter(Subtype_IHC_2 == "HER2" & Chemo == "NVB")

)

names(finher_her2_split[["inter.tra.nvb"]]) <- c("DDFS","OS")

#
# ==============================================================================




# 3. Analysis summary
# ==============================================================================

summarize_her2_2 <- function(x, prog = F){
  # x <- finher_her2_split$prog$DDFS$HER2
  # prog = T

  x <- x [names(x) != "ph_test"]
  x <- purrr::map(x,
                  function(xx, prog){

                    if(prog){

                      nme <- c("Variable", "HR", "Low95","Up95", "P", "P_adj", "Heterogeneity", "Ph_ok")
                      xx <- xx[ , nme] %>%
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
                      xx <- xx[ , nme] %>%
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


  names(x) <- str_replace(names(x),"main_effect","split")

  x <- purrr::map(names(x),
                  function(xname, x){

                    names(x[[xname]]) <- str_c(xname, "_", names(x[[xname]]))
                    x[[xname]]
                  },
                  x)


  bind_cols(x) %>%
    dplyr::rename(Vaiable = "split1_Variable") %>%
    dplyr::select(-c(split2_Variable))

}

ddfs <- list()

ddfs[["prog"]] <- summarize_her2_2(x = finher_her2_split$prog$DDFS$HER2, prog = T)
ddfs[["inter.chemo"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo$DDFS$HER2, prog = F)
ddfs[["inter.chemo.notra"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo.notra$DDFS$HER2, prog = F)
ddfs[["inter.chemo.tra"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo.tra$DDFS$HER2, prog = F)
ddfs[["inter.tra"]] <- summarize_her2_2(x = finher_her2_split$inter.tra$DDFS$HER2, prog = F)
ddfs[["inter.tra.dtx"]] <- summarize_her2_2(x = finher_her2_split$inter.tra.dtx$DDFS$HER2, prog = F)
ddfs[["inter.tra.nvb"]] <- summarize_her2_2(x = finher_her2_split$inter.tra.nvb$DDFS$HER2, prog = F)

write_xlsx(ddfs, path = str_c(out_tables,"finher_her2_split_ddfs_summary.xlsx"))

os <- list()

os[["prog"]] <- summarize_her2_2(x = finher_her2_split$prog$OS$HER2, prog = T)
os[["inter.chemo"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo$OS$HER2, prog = F)
os[["inter.chemo.notra"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo.notra$OS$HER2, prog = F)
os[["inter.chemo.tra"]] <- summarize_her2_2(x = finher_her2_split$inter.chemo.tra$OS$HER2, prog = F)
os[["inter.tra"]] <- summarize_her2_2(x = finher_her2_split$inter.tra$OS$HER2, prog = F)
os[["inter.tra.dtx"]] <- summarize_her2_2(x = finher_her2_split$inter.tra.dtx$OS$HER2, prog = F)
os[["inter.tra.nvb"]] <- summarize_her2_2(x = finher_her2_split$inter.tra.nvb$OS$HER2, prog = F)

write_xlsx(os, path = str_c(out_tables,"finher_her2_split_os_summary.xlsx"))


#
# ==============================================================================




# 4. Save Robjects
# ==============================================================================

save(finher_her2_split, file = str_c(out_data,"finher_her2_split.RData"))

#
# ==============================================================================




