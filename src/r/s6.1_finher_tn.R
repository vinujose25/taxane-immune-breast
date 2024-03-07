# s6.1_finher_tn.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Prognostic an Chemo interaction with TIL, TIL associated signatures and Immune celltypes.
# Account for arm heterogeneity.
# Prepare data for plotting.
# Generate sample plot for TN


# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis and chemo interaction in TN.
# 3. Analysis summary
# 4. Save Robjects



# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher_finneo.RData")


# Format clinical data
# >>>>>>>>>>>>>>>>>>>>

clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10) %>% # in units of 10% increments
  dplyr::filter(Subtype_IHC == "TN")


# Verify
clin_finher_finneo %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())
#   Subtype_IHC Subtype_IHC_2     N
# 1 TN          TN              120



# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

clin_finher_finneo %>%
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

#   Subtype            Arm  N RFS DDFS OS RFS_perc DDFS_perc  OS_perc
# 1      TN DTX.noTRA.noHR 60  14   13  9 23.33333  21.66667 15.00000
# 2      TN NVB.noTRA.noHR 60  21   19 14 35.00000  31.66667 23.33333


# Analysis
# 1. Sig's prognostic value in TN (with arm heterogeneity (DTX/NVB interaction))
# 2. Sig's interaction with chemo in TN

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

finher_tn <- list()

# TN prognosis
# >>>>>>>>>>>>


# Individual sig

finher_tn[["prog_individual.sig"]] <- purrr::map(

    c("RFS", "DDFS", "OS"),

    function(event_type, clin){

      print(event_type)

          # Per sig prognosis + heterogenity
          yy <- purrr::map(

            c(
              "TIL",

              "scaled_Denovo_TILsig",
              "scaled_Denovo_Immune",
              "scaled_Denovo_ECM",

              "scaled_Hamy2016_Immune",
              "scaled_Yang2018_Immune",

              "scaled_Gruosso2019_Interferon",
              "scaled_Farmer2009_MX1",
              "scaled_Hamy2016_Interferon",
              "scaled_Nirmal2018_Interferon",

              "scaled_Hamy2016_Ecm",
              "scaled_Naba2014_Ecmcore",
              "scaled_Triulzi2013_Ecm",

              "scaled_Sorrentino2014_Chol",

              "MCPcounter_Fibroblasts"
              ),

            function(sig, event_type, clin){ # subtype,

              get_adj_prog(
                event_variable = str_c(event_type, "_Event"),
                time_variable = str_c(event_type, "_Time_Years"),
                biomarker = sig,
                xdata = clin
              )

              # # Prognostic model
              # paste("Surv(", time_variable, ",", event_variable, ") ~",
              #           biomarker, "+ strata(Hormone, Herceptin, Chemo)")

              # # Heterogenity/interaction model
              # paste("Surv(", time_variable, ",", event_variable, ") ~",
              #           biomarker,"* interaction(Hormone, Herceptin, Chemo)")

            },

            event_type,
            clin
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
                Subtype = "TN"
              ),

            ph_test = purrr::map(yy,~(.x$ph_test)) %>%
              bind_rows()
          )

    },

    clin = clin_finher_finneo
  )

names(finher_tn[["prog_individual.sig"]]) <- c("RFS", "DDFS", "OS")



# Pooled sig

finher_tn[["prog_pooled.sig"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    # Per sig prognosis + heterogenity
    yy <- purrr::map(

      c(
        "TIL",

        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Fibroblasts"
      ),

      function(sig, event_type, clin){ # subtype,

        get_adj_prog(
          event_variable = str_c(event_type, "_Event"),
          time_variable = str_c(event_type, "_Time_Years"),
          biomarker = sig,
          xdata = clin
        )

        # # Prognostic model
        # paste("Surv(", time_variable, ",", event_variable, ") ~",
        #           biomarker, "+ strata(Hormone, Herceptin, Chemo)")

        # # Heterogenity/interaction model
        # paste("Surv(", time_variable, ",", event_variable, ") ~",
        #           biomarker,"* interaction(Hormone, Herceptin, Chemo)")

      },

      event_type,
      clin
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
          Subtype = "TN"
        ),

      ph_test = purrr::map(yy,~(.x$ph_test)) %>%
        bind_rows()
    )

  },

  clin = clin_finher_finneo
)

names(finher_tn[["prog_pooled.sig"]]) <- c("RFS", "DDFS", "OS")





# TN Chemo interaction
# >>>>>>>>>>>>>>>>>>>>


# individual sig

finher_tn[["inter.chemo_individual.sig"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    # Per sig mta interaction
    yy <- purrr::map(

      c(
        "TIL",

        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Hamy2016_Immune",
        "scaled_Yang2018_Immune",

        "scaled_Gruosso2019_Interferon",
        "scaled_Farmer2009_MX1",
        "scaled_Hamy2016_Interferon",
        "scaled_Nirmal2018_Interferon",

        "scaled_Hamy2016_Ecm",
        "scaled_Naba2014_Ecmcore",
        "scaled_Triulzi2013_Ecm",

        "scaled_Sorrentino2014_Chol",

        "MCPcounter_Fibroblasts"
        ),


      function(biomarker, event_type, clin){

        event_variable = str_c(event_type, "_Event")
        time_variable = str_c(event_type, "_Time_Years")


        chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                  biomarker, "+ Chemo")
        chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                  biomarker, "* Chemo")

        get_adj_inter(
          event_variable = event_variable, # to extract event summary
          biomarker = biomarker,
          interaction_variable = "Chemo",
          xdata = clin,
          chr_null_formula = chr_null_formula,
          chr_full_formula = chr_full_formula
        )

      },
      event_type,
      clin
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
          Subtype = "TN"
        ),

      ph_test = purrr::map(yy,~(.x$ph_test)) %>%
        bind_rows()
    )


  },

  clin = clin_finher_finneo

)

names(finher_tn[["inter.chemo_individual.sig"]]) <- c("RFS", "DDFS", "OS")



# pooled sig

finher_tn[["inter.chemo_pooled.sig"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    # Per sig mta interaction
    yy <- purrr::map(

      c(
        "TIL",

        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Fibroblasts"
      ),


      function(biomarker, event_type, clin){

        event_variable = str_c(event_type, "_Event")
        time_variable = str_c(event_type, "_Time_Years")


        chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                  biomarker, "+ Chemo")
        chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
                                  biomarker, "* Chemo")

        get_adj_inter(
          event_variable = event_variable, # to extract event summary
          biomarker = biomarker,
          interaction_variable = "Chemo",
          xdata = clin,
          chr_null_formula = chr_null_formula,
          chr_full_formula = chr_full_formula
        )

      },
      event_type,
      clin
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
          Subtype = "TN"
        ),

      ph_test = purrr::map(yy,~(.x$ph_test)) %>%
        bind_rows()
    )


  },

  clin = clin_finher_finneo

)

names(finher_tn[["inter.chemo_pooled.sig"]]) <- c("RFS", "DDFS", "OS")


#
# ==============================================================================



#
# # 3. Analysis summary
# # ==============================================================================
#
# summarize_tn <- function(x, prog = F){
#   # xx <- finher_tn$prog$DDFS
#   # prog = T
#
#   xx <- x
#
#   if(prog){
#
#     nme <- c("Variable", "HR", "Low95","Up95", "P", "P_adj", "Heterogeneity", "Ph_ok")
#     xx <- xx$main_effect[ , nme] %>%
#       dplyr::mutate( HR = str_c(round(HR, digits = 2),
#                                 " (",
#                                 round(Low95,digits=2), "-",
#                                 round(Up95,digits=2),
#                                 ")"),
#                      # Note, for prognostic model the CI is in HR
#                      Heterogeneity = str_c(round(Heterogeneity, digits=2))) %>%
#       dplyr::rename(Het = "Heterogeneity")
#
#   } else {
#
#     nme <- c("Variable", "HR", "Low95","Up95", "P", "P_adj", "Ph_ok")
#     xx <- xx$main_effect[ , nme] %>%
#       dplyr::mutate( HR = str_c(round(HR, digits = 2),
#                                 " (",
#                                 round(exp(Low95),digits=2), "-",
#                                 round(exp(Up95),digits=2),
#                                 ")"))
#     # Note, for interaction model the CI is in coef
#
#   }
#
#   xx <- xx %>%
#     dplyr::mutate(
#       P = str_c(round(P, digits=2)),
#       P_adj = str_c(round(P_adj, digits=2))
#     ) %>%
#     dplyr::select(-c(Low95,Up95))
#
#
#   return(xx)
#
# }
#
# xout <- list()
#
# xout[["ddfs"]] <- summarize_tn(x = finher_tn$prog_individual.sig$DDFS, prog = T)
# xout[["os"]] <- summarize_tn(x = finher_tn$prog_individual.sig$OS, prog = T)
# xout[["rfs"]] <- summarize_tn(x = finher_tn$prog_individual.sig$RFS, prog = T)
#
# write_xlsx(xout, path = str_c(out_tables,"finher_tn_prognosis_individual.sig_summary.xlsx"))
#
#
#
# xout <- list()
#
# xout[["ddfs"]] <- summarize_tn(x = finher_tn$prog_pooled.sig$DDFS, prog = T)
# xout[["os"]] <- summarize_tn(x = finher_tn$prog_pooled.sig$OS, prog = T)
# xout[["rfs"]] <- summarize_tn(x = finher_tn$prog_pooled.sig$RFS, prog = T)
#
# write_xlsx(xout, path = str_c(out_tables,"finher_tn_prognosis_pooled.sig_summary.xlsx"))
#
#
#
# xout <- list()
#
# xout[["ddfs"]] <-  summarize_tn(x = finher_tn$inter.chemo_individual.sig$DDFS, prog = F)
# xout[["os"]] <- summarize_tn(x = finher_tn$inter.chemo_individual.sig$OS, prog = F)
# xout[["rfs"]] <- summarize_tn(x = finher_tn$inter.chemo_individual.sig$RFS, prog = F)
#
# write_xlsx(xout, path = str_c(out_tables,"finher_tn_chemo_interaction_individual.sig_summary.xlsx"))
#
#
# xout <- list()
#
# xout[["ddfs"]] <-  summarize_tn(x = finher_tn$inter.chemo_pooled.sig$DDFS, prog = F)
# xout[["os"]] <- summarize_tn(x = finher_tn$inter.chemo_pooled.sig$OS, prog = F)
# xout[["rfs"]] <- summarize_tn(x = finher_tn$inter.chemo_pooled.sig$RFS, prog = F)
#
# write_xlsx(xout, path = str_c(out_tables,"finher_tn_chemo_interaction_pooled.sig_summary.xlsx"))
#
#
# #
# # ==============================================================================
#



# 4. Save Robjects
# ==============================================================================

save(finher_tn, file = str_c(out_data,"finher_tn.RData"))

#
# ==============================================================================


# Clear memory
# ==============================================================================

rm(clin_finher_finneo, finher_tn)

#
# ==============================================================================
