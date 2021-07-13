# finher_tn.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Prognostic an Chemo interaction with TIL, TIL associated signatures and Immune celltypes.
# Account for arm heterogenity.
# Prepare data for plotting.
# Generate sample plot for TN



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis and chemo interaction in TN.
# 3. Analysis summary
# 4. Save Robjects



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
                TIL = StrLy_Mean/10) %>%
  dplyr::filter(Subtype_IHC == "TN")

clin_finher %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 TN          TN              120



# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

clin_finher %>%
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
# 1. Sig's prognostic value in TN (with arm heterogenity (DTX/NVB interaction))
# 2. Sig's interaction with chemo in TN

# Note: By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to clim a significant treatment
# interaction.

#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_tn <- list()

# TN prognosis
# >>>>>>>>>>>>

finher_tn[["prog"]] <- purrr::map(

    c("RFS", "DDFS", "OS"),

    function(event_type, clin){

      print(event_type)

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

    clin = clin_finher
  )

names(finher_tn[["prog"]]) <- c("RFS", "DDFS", "OS")



# TN Chemo interaction
# >>>>>>>>>>>>>>>>>>>>
finher_tn[["inter.chemo"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

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

  clin = clin_finher

)


names(finher_tn[["inter.chemo"]]) <- c("RFS", "DDFS", "OS")

#
# ==============================================================================






# 3. Analysis summary
# ==============================================================================

# Analsysis summary: TN prognosis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# PH
finher_tn$prog$RFS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing
finher_tn$prog$DDFS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing
finher_tn$prog$OS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing


# prognosis
finher_tn$prog$RFS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
#         Coef        HR          P     Low95     Up95 Heterogeneity Variable  Endpoint Event   N P_adj Subtype
# 1 -0.1990437 0.8195141 0.05954192 0.6662475 1.008039     0.1819069      TIL RFS_Event    35 120    NA      TN
finher_tn$prog$DDFS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
#         Coef        HR          P      Low95      Up95 Heterogeneity Variable   Endpoint Event   N     P_adj Subtype
# 1 -0.2384062 0.7878826 0.03695126 0.62978589 0.9856668     0.1422721      TIL DDFS_Event    32 120        NA      TN
# 2 -1.2110289 0.2978906 0.09248506 0.07266994 1.2211215     0.7371714    Tcell DDFS_Event    32 120 0.6438521      TN
finher_tn$prog$OS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
# nothing

# !!!!!!!!!!!!!!!!!!
# Only TIL is prognostic with DDFS


# Analsysis summary: TN chemo interaction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# PH
finher_tn$inter.chemo$RFS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing
finher_tn$inter.chemo$DDFS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing
finher_tn$inter.chemo$OS$ph_test %>% dplyr::filter(p<.1) %>% as.data.frame()
# nothing


# interaction
finher_tn$inter.chemo$RFS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
#        Coef       HR          P        Low95     Up95                Therapy      Variable  Endpoint Event   N Ph_ok     P_adj Subtype
# 1 1.0233039 2.782372 0.06682110 -1.093378824 3.139987 TILsig_APP_Fc:ChemoDTX TILsig_APP_Fc RFS_Event    35 120  TRUE 0.2100092      TN
# 2 1.1969647 3.310055 0.03880712 -0.829763547 3.223693 TILsig_Immune:ChemoDTX TILsig_Immune RFS_Event    35 120  TRUE 0.2054705      TN
# 3 1.0048688 2.731549 0.05652231 -0.977443223 2.987181   TILsig_IFNg:ChemoDTX   TILsig_IFNg RFS_Event    35 120  TRUE 0.2072485      TN
# 4 1.3304993 3.782932 0.03807721 -0.817407709 3.478406 TILsig_Innate:ChemoDTX TILsig_Innate RFS_Event    35 120  TRUE 0.2054705      TN
# 5 0.7707341 2.161352 0.08571518 -1.164351226 2.705819       Immune2:ChemoDTX       Immune2 RFS_Event    35 120  TRUE 0.2357168      TN
# 6 1.3348210 3.799316 0.04669785 -0.757021025 3.426663   Interferon1:ChemoDTX   Interferon1 RFS_Event    35 120  TRUE 0.2054705      TN
# 7 1.9017051 6.697304 0.01100269  0.005308771 3.798101   Interferon2:ChemoDTX   Interferon2 RFS_Event    35 120  TRUE 0.2054705      TN
# 8 1.5766663 4.838798 0.02750725 -0.386488773 3.539821   Interferon3:ChemoDTX   Interferon3 RFS_Event    35 120  TRUE 0.2054705      TN
finher_tn$inter.chemo$DDFS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
#        Coef       HR          P        Low95     Up95                Therapy      Variable   Endpoint Event   N Ph_ok      P_adj Subtype
# 1  1.269340 3.558504 0.04116259 -0.961635210 3.500316 TILsig_APP_Fc:ChemoDTX TILsig_APP_Fc DDFS_Event    32 120  TRUE 0.10061967      TN
# 2  1.499065 4.477502 0.01609870 -0.648363225 3.646494 TILsig_Immune:ChemoDTX TILsig_Immune DDFS_Event    32 120  TRUE 0.08331451      TN
# 3  1.234452 3.436496 0.03029619 -0.849982288 3.318887   TILsig_IFNg:ChemoDTX   TILsig_IFNg DDFS_Event    32 120  TRUE 0.08331451      TN
# 4  1.409685 4.094666 0.02742959 -0.837773384 3.657143 TILsig_Innate:ChemoDTX TILsig_Innate DDFS_Event    32 120  TRUE 0.08331451      TN
# 5  1.558373 4.751084 0.05222046 -0.774326841 3.891073       Immune1:ChemoDTX       Immune1 DDFS_Event    32 120  TRUE 0.11488502      TN
# 6  1.222927 3.397117 0.02133537 -0.814162887 3.260017       Immune2:ChemoDTX       Immune2 DDFS_Event    32 120  TRUE 0.08331451      TN
# 7  1.264160 3.540118 0.03013899 -0.804739464 3.333060       Immune3:ChemoDTX       Immune3 DDFS_Event    32 120  TRUE 0.08331451      TN
# 8  1.826694 6.213314 0.01710431 -0.386190972 4.039580   Interferon1:ChemoDTX   Interferon1 DDFS_Event    32 120  TRUE 0.08331451      TN
# 9  2.126626 8.386520 0.01335468  0.187488561 4.065763   Interferon2:ChemoDTX   Interferon2 DDFS_Event    32 120  TRUE 0.08331451      TN
# 10 2.006531 7.437469 0.01317550 -0.005940085 4.019001   Interferon3:ChemoDTX   Interferon3 DDFS_Event    32 120  TRUE 0.08331451      TN
# 11 1.584030 4.874559 0.05932885 -0.559400215 3.727460  Cholesterol2:ChemoDTX  Cholesterol2 DDFS_Event    32 120  TRUE 0.11865770      TN
finher_tn$inter.chemo$OS$main_effect %>% dplyr::filter(P<.1) %>% as.data.frame()
#        Coef        HR           P      Low95     Up95                Therapy      Variable Endpoint Event   N Ph_ok      P_adj Subtype
# 1  2.840038 17.116418 0.006792984 -0.1880311 5.868107 TILsig_APP_Fc:ChemoDTX TILsig_APP_Fc OS_Event    23 120  TRUE 0.04682614      TN
# 2  2.424768 11.299604 0.005758224 -0.3396796 5.189215 TILsig_Immune:ChemoDTX TILsig_Immune OS_Event    23 120  TRUE 0.04682614      TN
# 3  1.956559  7.074940 0.010295993 -0.6621021 4.575220   TILsig_IFNg:ChemoDTX   TILsig_IFNg OS_Event    23 120  TRUE 0.04682614      TN
# 4  2.563581 12.982221 0.008808195 -0.3472161 5.474378 TILsig_Innate:ChemoDTX TILsig_Innate OS_Event    23 120  TRUE 0.04682614      TN
# 5  2.042410  7.709167 0.031075286 -0.7806577 4.865478       Immune1:ChemoDTX       Immune1 OS_Event    23 120  TRUE 0.06836563      TN
# 6  1.531483  4.625032 0.020990737 -0.9539417 4.016908       Immune2:ChemoDTX       Immune2 OS_Event    23 120  TRUE 0.05596067      TN
# 7  1.880162  6.554565 0.016117129 -0.6715944 4.431918       Immune3:ChemoDTX       Immune3 OS_Event    23 120  TRUE 0.05065384      TN
# 8  2.243918  9.430203 0.011198390 -0.4433457 4.931181   Interferon1:ChemoDTX   Interferon1 OS_Event    23 120  TRUE 0.04682614      TN
# 9  2.029300  7.608760 0.022893000 -0.2278734 4.286474   Interferon2:ChemoDTX   Interferon2 OS_Event    23 120  TRUE 0.05596067      TN
# 10 2.111066  8.257037 0.012770766 -0.2443011 4.466433   Interferon3:ChemoDTX   Interferon3 OS_Event    23 120  TRUE 0.04682614      TN


# !!!!!!!!!!!!!!!!!!
# TIL is not interacting with Chemo, but TILsig_imm
# Immune and interferon is significantly interacting with chemo DDFS, OS (RFS trend)


#
# ==============================================================================




# 4. Save Robjects
# ==============================================================================

save(finher_tn, file = str_c(out_data,"finher_tn.RData"))

#
# ==============================================================================
