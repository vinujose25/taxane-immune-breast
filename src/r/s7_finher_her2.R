# finher_her2.R


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
# 3. Analysis summary
# 4. Save Robjects


# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")



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
# 7 HR+HER2+     DTX.TRA.HR 16   2    0  0 12.50000   0.00000  0.000000 # no ddfs/os event
# 8 HR+HER2+     NVB.TRA.HR 30  10    9  4 33.33333  30.00000 13.333333

# Average size of strata, c(27,19,26,17,18,27,16,30) %>% mean() = 22.5

# Analysis
# Note that there are *eight* stratas in HER2 BC !!!!!!!!!!!!!!!
# 1. Sig's prognostic value in HER2 (with arm heterogenity; interaction with *eight* stratas)
# 2. Sig's interaction with chemo in HER2 in pooled data with HR/TRA stratification;
#     no within HR subtype or TRA treatment subgroup analysis due to low sample size.

# Note: By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to claim a significant treatment
# interaction.

#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_her2 <- list()

# HER2 prognosis
# >>>>>>>>>>>>>>

finher_her2[["prog"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      c("HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig prognosis + heterogenity
        yy <- purrr::map(

          c(
            "TIL",
            "scaled_Denovo_TILsig",
            "scaled_Denovo_Immune",
            "scaled_Denovo_ECM",

            "scaled_Gruosso2019_Immune",
            "scaled_Gruosso2019_Interferon",
            "scaled_Gruosso2019_Cholesterol",
            "scaled_Gruosso2019_Fibrosis",

            "scaled_General_Immune",
            "scaled_General_Interferon",
            "scaled_General_Cholesterol",
            "scaled_General_ECM",

            "scaled_MCPcounter_Tcell",
            "scaled_MCPcounter_Monocytic.Lineage",
            "scaled_MCPcounter_Fibroblasts",

            "scaled_Control_Proliferation",
            "scaled_Control_Non.breast.tissue",
            "scaled_Control_Human_Behaviour"
          ),

          function(sig, event_type, xdata){ # subtype,

            get_adj_prog(
              event_variable = str_c(event_type, "_Event"),
              time_variable = str_c(event_type, "_Time_Years"),
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
              # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
              P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
                p.adjust(method = "BH"),
              Subtype = "HER2"
              # Subtype = subtype
            ),

          ph_test = purrr::map(yy,~(.x$ph_test)) %>%
            bind_rows()
        )


      },
      event_type,
      clin
    )

    # names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")
    names(xx) <- c("HER2")

    xx
  },

  clin = clin_finher %>%
    dplyr::filter(Subtype_IHC_2 == "HER2")
  # %>%
  #   dplyr::filter(Arm != "DTX.TRA.HR")

)

names(finher_her2[["prog"]]) <- c("RFS", "DDFS", "OS")

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
finher_her2[["inter.chemo"]] <- purrr::map(

    c("RFS", "DDFS", "OS"),

    function(event_type, clin){

      print(event_type)

      xx <- purrr::map(

        # c("HR-HER2+", "HR+HER2+", "HER2"),
        c("HER2"),

        function(subtype, event_type, clin){

          print(subtype)

          subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
          xdata <- subset(x = clin,
                          subset = (clin[, subtype_varible] == subtype))



          # Per sig mta interaction
          yy <- purrr::map(

            c(
              "TIL",
              "scaled_Denovo_TILsig",
              "scaled_Denovo_Immune",
              "scaled_Denovo_ECM",

              "scaled_Gruosso2019_Immune",
              "scaled_Gruosso2019_Interferon",
              "scaled_Gruosso2019_Cholesterol",
              "scaled_Gruosso2019_Fibrosis",

              "scaled_General_Immune",
              "scaled_General_Interferon",
              "scaled_General_Cholesterol",
              "scaled_General_ECM",

              "scaled_MCPcounter_Tcell",
              "scaled_MCPcounter_Monocytic.Lineage",
              "scaled_MCPcounter_Fibroblasts",

              "scaled_Control_Proliferation",
              "scaled_Control_Non.breast.tissue",
              "scaled_Control_Human_Behaviour"
            ),


            function(biomarker, event_type, xdata){

              event_variable = str_c(event_type, "_Event")
              time_variable = str_c(event_type, "_Time_Years")


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

      # names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")
      names(xx) <- c("HER2")

      xx
    },

    clin = clin_finher %>%
      dplyr::filter(Subtype_IHC_2 == "HER2")

  )

names(finher_her2[["inter.chemo"]]) <- c("RFS", "DDFS", "OS")

#
# ==============================================================================




# 3. Analysis summary
# ==============================================================================

summarize_her2 <- function(x, prog = F){
  # x <- finher_her2$prog$DDFS
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

  return(x[[1]])

}

xout <- list()

xout[["ddfs"]] <- summarize_her2(x = finher_her2$prog$DDFS, prog = T)
xout[["os"]] <- summarize_her2(x = finher_her2$prog$OS, prog = T)
xout[["rfs"]] <- summarize_her2(x = finher_her2$prog$RFS, prog = T)

write_xlsx(xout, path = str_c(out_tables,"finher_her2_prognosis_summary.xlsx"))


xout <- list()

xout[["ddfs"]] <-  summarize_her2(x = finher_her2$inter.chemo$DDFS, prog = F)
xout[["os"]] <- summarize_her2(x = finher_her2$inter.chemo$OS, prog = F)
xout[["rfs"]] <- summarize_her2(x = finher_her2$inter.chemo$RFS, prog = F)

write_xlsx(xout, path = str_c(out_tables,"finher_her2_chemo_interaction_summary.xlsx"))

#
# ==============================================================================




# 4. Save Robjects
# ==============================================================================

save(finher_her2, file = str_c(out_data,"finher_her2.RData"))

#
# ==============================================================================




