# s4_exploring_interaction.R



# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the interaction between
#   1) Arm, Signature, and pCR within each Subtype, and
#   2) Subtype, Signature and pCR within each Arm.
# Account for study heterogenity.
# Prepare data necessary to plot prognostic plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load data.
# 2. Structure analysis
# 3. Interaction test
# 4. Analysis summary tables
# 5. Save Robjects odf interaction summaries




# 1. Load data
# ==============================================================================

load("results/data/clin_neoadj.RData")

dim(clin_neoadj)
# 1499  164


# clean Arm and factor levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
clin_neoadj$Arm_consolidated %>% factor() %>% levels()
# [1] "A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [2] "AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [3] "AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [4] "AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy"

clin_neoadj <- clin_neoadj %>%
  dplyr::mutate(
    Arm_consolidated =  Arm_consolidated %>%
      str_replace("\\+noTaxane", "") %>%
      str_replace("\\+Taxane", "+T") %>%
      str_replace("///Trastuzumab", "+TRA") %>%
      str_replace("///No_her2_agent", "") %>%
      str_replace("///No_hormone_therapy", "") %>%
      str_replace("///No_other_therapy", "") %>%
      factor(levels = c("AAA+T", "A0A+T", "AAA", "AAA+T+TRA"))
  )

clin_neoadj$Arm_consolidated %>% factor() %>% levels()
# [1] "AAA+T"     "A0A+T"     "AAA"       "AAA+T+TRA"



# set subtype ihc levels
# >>>>>>>>>>>>>>>>>>>>>>
clin_neoadj$Subtype_ihc %>% factor() %>% levels()
# [1] "HER2" "HR"   "TN"
clin_neoadj$Subtype_ihc2 %>% factor() %>% levels()
# [1] "HR"       "HR-HER2+" "HR+HER2+" "TN"

clin_neoadj <- clin_neoadj %>%
  dplyr::mutate(
    Subtype_ihc =  Subtype_ihc %>%
      factor(levels = c("TN", "HER2", "HR")),
    Subtype_ihc2 =  Subtype_ihc2 %>%
      factor(levels = c("TN", "HR-HER2+", "HR+HER2+", "HR"))
  )

clin_neoadj$Subtype_ihc %>% levels()
# [1] "TN" "HER2" "HR"
clin_neoadj$Subtype_ihc2 %>% levels()
# [1] "TN" "HR-HER2+" "HR+HER2+" "HR"
#
# ==============================================================================




# 2. Structure analysis
# ==============================================================================


# Summarize clin_neoadj
# >>>>>>>>>>>>>>>>>>>>>

# All neoadj regimen
clin_neoadj %>%
  dplyr::group_by(
    Subtype = Subtype_ihc,
    Hr,
    Arm = Arm_consolidated
  ) %>%
  dplyr::summarise(N = n()) %>%
  as.data.frame()
#    Subtype  Hr       Arm   N
# 1       TN neg     AAA+T 293
# 2       TN neg     A0A+T 150
# 3       TN neg       AAA  83

# 4     HER2 neg     AAA+T  80
# 5     HER2 neg AAA+T+TRA  58

# 6     HER2 pos     AAA+T  55
# 7     HER2 pos AAA+T+TRA  29

# 8       HR pos     AAA+T 521
# 9       HR pos     A0A+T 139
# 10      HR pos       AAA  91


# TN analysis
# 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 3. Per subtype global interaction: all-Arm * Sig * pCR


# HER2 analysis
# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR


# HR analysis
# 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 3. Per subtype global interaction: all-Arm * Sig * pCR


# Pan-Subtype (ALL) analysis
# 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
# 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
# 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR


# Analysis note: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# 1. Focus the analysis on AAA containing regimens due to large sample size.
# 2. Make plots accordingly.
# 3. Keep the other analysis as supplimentary.
# 4. Script structure:
#     Create geo_inter (interaction summary) list object.
#     Each elemet represent all anlysis done on each subtype.
#



# All neoadj regimen with Strata info
clin_neoadj %>%
  dplyr::group_by(
    Subtype = Subtype_ihc,
    Hr,
    Arm = Arm_consolidated,
    Strata # identical results when Strata =  Series_matrix_accession
  ) %>%
  dplyr::summarise(N = n(), Event = which(Response == 1) %>% length()) %>%
  dplyr::mutate(Event_percet = (Event/N) * 100) %>%
  as.data.frame()
#    Subtype  Hr       Arm             Strata   N Event Event_percet
# 1       TN neg     AAA+T     GSE20194:AAA+T  64    25    39.062500
# 2       TN neg     AAA+T     GSE20271:AAA+T  31     9    29.032258
# 3       TN neg     AAA+T     GSE23988:AAA+T  29    13    44.827586
# 4       TN neg     AAA+T     GSE25066:AAA+T 119    43    36.134454
# 5       TN neg     AAA+T     GSE32646:AAA+T  26    10    38.461538
# 6       TN neg     AAA+T     GSE42822:AAA+T  24    12    50.000000
# 7       TN neg     A0A+T     GSE25066:A0A+T  25     6    24.000000
# 8       TN neg     A0A+T     GSE41998:A0A+T 125    51    40.800000
# 9       TN neg       AAA       GSE20271:AAA  28     4    14.285714
# 10      TN neg       AAA       GSE22093:AAA  55    18    32.727273

# 11    HER2 neg     AAA+T     GSE20194:AAA+T  25    13    52.000000
# 12    HER2 neg     AAA+T     GSE32646:AAA+T  18     9    50.000000
# 13    HER2 neg     AAA+T     GSE50948:AAA+T  37     9    24.324324
# 14    HER2 neg AAA+T+TRA GSE42822:AAA+T+TRA  15    10    66.666667
# 15    HER2 neg AAA+T+TRA GSE50948:AAA+T+TRA  43    25    58.139535

# 16    HER2 pos     AAA+T     GSE20194:AAA+T  25     4    16.000000
# 17    HER2 pos     AAA+T     GSE32646:AAA+T  16     3    18.750000
# 18    HER2 pos     AAA+T     GSE50948:AAA+T  14     4    28.571429
# 19    HER2 pos AAA+T+TRA GSE42822:AAA+T+TRA   9     2    22.222222
# 20    HER2 pos AAA+T+TRA GSE50948:AAA+T+TRA  20     6    30.000000

# 21      HR pos     AAA+T     GSE20194:AAA+T 141     9     6.382979
# 22      HR pos     AAA+T     GSE20271:AAA+T  44     3     6.818182
# 23      HR pos     AAA+T     GSE23988:AAA+T  32     7    21.875000
# 24      HR pos     AAA+T     GSE25066:AAA+T 194    22    11.340206
# 25      HR pos     AAA+T     GSE32646:AAA+T  55     5     9.090909
# 26      HR pos     AAA+T     GSE42822:AAA+T  29     8    27.586207
# 27      HR pos     AAA+T     GSE50948:AAA+T  26     4    15.384615
# 28      HR pos     A0A+T     GSE25066:A0A+T  35     3     8.571429
# 29      HR pos     A0A+T     GSE41998:A0A+T 104    12    11.538462
# 30      HR pos       AAA       GSE20271:AAA  49     3     6.122449
# 31      HR pos       AAA       GSE22093:AAA  42    10    23.809524


# Note !!!!!
# The event rate in HR subtype is low possibily due to lack of hormone therapy

#
# ==============================================================================




# 3. Interaction tests
# ==============================================================================

# geo_inter <- vector(mode = "list", length = 4)
# names(geo_inter) <- c("TN", "HER2", "HR", "ALL")

geo_inter <- vector(mode = "list", length = 6)
names(geo_inter) <- c("TN", "HR-HER2+", "HR+HER2+","HER2", "HR", "ALL")


# Per-subtype: Arm, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# TN analysis
# @@@@@@@@@@@

# 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 3. Per subtype global interaction: all-Arm * Sig * pCR

# Naming convention:
# Name analysis based on backbone(non-interacting) therapy and interaction tested therpies.
# "." seperates back bone therpy from interaction tested therpaies/subtypes.
# The part left to the "." represents backbone therapy and
#   the part to the right represents interacion tested therapies or subtypes.

# Recommanded R naming convention:
# Name should only contains letters, numbers underscre and dot.
# Ref:


# Analysis structure
geo_inter$TN[["AAA.T_or_NoT"]] <- NA # main analysis
geo_inter$TN[["T.AAA_or_A0A"]] <- NA # supplemenatry analysis
geo_inter$TN[["Global"]] <- NA # supplemenatry analysis



# Interaction summary
geo_inter$TN <- purrr::map(

  #   Subtype       Arm   N
  # 1      TN     AAA+T 293
  # 2      TN     A0A+T 150
  # 3      TN       AAA  83


  c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
    "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
    "Global"), # 3. Per subtype global interaction: all-Arm * Sig * pCR

  function(analysis, clin){

    print(analysis)

    # Test data
    # clin = clin_neoadj
    # analysis = "AAA.T_or_NoT"
    # sig = "Immune1"

    xclin <- switch (
      analysis,

      # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
      AAA.T_or_NoT = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "AAA")),

      # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
      T.AAA_or_A0A = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "\\+T")),

      # 3. Per subtype global interaction: all-Arm * Sig * pCR
      Global = clin %>%
        dplyr::filter(Subtype_ihc == "TN")
    )

    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )


    # sig = "Immune1"
    # analysis = "AAA.T_or_NoT"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      # c("TILsig", "TILsig_imm", "TILsig_fib",
      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Arm_consolidated"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                              test_var = sig, inter_var = inter_var)



        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Arm_consolidated) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Arm_consolidated = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Arm_consolidated"),
                    by = "Interaction_var_levels")


        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){

            print(l)


            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm

              # Note in HR+HER2+ after filtering out strata with <10  samples,
              # only one dataset has TRA containing regime. Hence dataset heterogenity
              # assesment can't be done for this regimen. This situation will be
              # handled in estimate_prog_heterogenity()

              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Arm_consolidated == l),
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$TN) <- c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
                         "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
                         "Global") # 3. Per subtype global interaction: all-Arm * Sig * pCR



# HR-HER2+ analysis
# @@@@@@@@@@@@@@@@@

# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

# Analysis structure
geo_inter$`HR-HER2+`[["AAA_plus_T.TRA_or_NoTRA"]] <- NA # main analysis


# Interaction summary
geo_inter$`HR-HER2+` <- purrr::map(

  #   Subtype       Arm   N
  # 4    HER2     AAA+T 135
  # 5    HER2 AAA+T+TRA  87


  c("AAA_plus_T.TRA_or_NoTRA"), # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

  function(analysis, clin){

    print(analysis)

    # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR
    xclin <- clin %>%
      dplyr::filter(Subtype_ihc2 == "HR-HER2+" &
                      str_detect(Arm_consolidated, "AAA"))

    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )


    # sig = "Immune2"
    # analysis = "AAA_plus_T.TRA_or_NoTRA"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Arm_consolidated"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                                        test_var = sig, inter_var = inter_var)





        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Arm_consolidated) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Arm_consolidated = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Arm_consolidated"),
                    by = "Interaction_var_levels")

        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){

            print(l)

            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm
              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Arm_consolidated == l),
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$`HR-HER2+`) <- "AAA_plus_T.TRA_or_NoTRA"






# HR+HER2+ analysis
# @@@@@@@@@@@@@@@@@

# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

# Analysis structure
geo_inter$`HR+HER2+`[["AAA_plus_T.TRA_or_NoTRA"]] <- NA # main analysis


# Interaction summary
geo_inter$`HR+HER2+` <- purrr::map(

  #   Subtype       Arm   N
  # 4    HER2     AAA+T 135
  # 5    HER2 AAA+T+TRA  87


  c("AAA_plus_T.TRA_or_NoTRA"), # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

  function(analysis, clin){

    print(analysis)

    # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR
    xclin <- clin %>%
      dplyr::filter(Subtype_ihc2 == "HR+HER2+" &
                      str_detect(Arm_consolidated, "AAA")) %>%
      # discarding strata with <10 samples(n=9 and event=2)
      dplyr::filter(Series_matrix_accession != "GSE42822")

    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )


    # sig = "Immune2"
    # analysis = "AAA_plus_T.TRA_or_NoTRA"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Arm_consolidated"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                                        test_var = sig, inter_var = inter_var)





        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Arm_consolidated) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Arm_consolidated = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Arm_consolidated"),
                    by = "Interaction_var_levels")

        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){

            print(l)

            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm
              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Arm_consolidated == l),
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$`HR+HER2+`) <- "AAA_plus_T.TRA_or_NoTRA"



# HER2 analysis
# @@@@@@@@@@@@@

# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

# Analysis structure
geo_inter$HER2[["AAA_plus_T.TRA_or_NoTRA"]] <- NA # main analysis


# Interaction summary
geo_inter$HER2 <- purrr::map(

  #   Subtype       Arm   N
  # 4    HER2     AAA+T 135
  # 5    HER2 AAA+T+TRA  87


  c("AAA_plus_T.TRA_or_NoTRA"), # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

  function(analysis, clin){

    print(analysis)

    # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR
    xclin <- clin %>%
      dplyr::filter(Subtype_ihc == "HER2" &
                      str_detect(Arm_consolidated, "AAA"))

    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )


    # sig = "Immune2"
    # analysis = "AAA_plus_T.TRA_or_NoTRA"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession")),
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Hr + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Series_matrix_accession")),
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Hr + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Arm_consolidated"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                                        test_var = sig, inter_var = inter_var)





        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Arm_consolidated) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Arm_consolidated = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Arm_consolidated"),
                    by = "Interaction_var_levels")

        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){

            print(l)

            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm
              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Arm_consolidated == l),
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$HER2) <- "AAA_plus_T.TRA_or_NoTRA"




# HR analysis
# @@@@@@@@@@@

# 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 3. Per subtype global interaction: all-Arm * Sig * pCR


# Analysis structure
geo_inter$HR[["AAA.T_or_NoT"]] <- NA # main analysis
geo_inter$HR[["T.AAA_or_A0A"]] <- NA # supplemenatry analysis
geo_inter$HR[["Global"]] <- NA # supplemenatry analysis


# Interaction summary
geo_inter$HR <- purrr::map(

  #   Subtype       Arm   N
  # 6      HR     AAA+T 521
  # 7      HR     A0A+T 139
  # 8      HR       AAA  91


  c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
    "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
    "Global"), # 3. Per subtype global interaction: all-Arm * Sig * pCR

  function(analysis, clin){

    print(analysis)

    xclin <- switch (
      analysis,

      # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
      AAA.T_or_NoT = clin %>%
        dplyr::filter(Subtype_ihc == "HR"  &
                        str_detect(Arm_consolidated, "AAA")),

      # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
      T.AAA_or_A0A = clin %>%
        dplyr::filter(Subtype_ihc == "HR" &
                        str_detect(Arm_consolidated, "\\+T")),

      # 3. Per subtype global interaction: all-Arm * Sig * pCR
      Global = clin %>%
        dplyr::filter(Subtype_ihc == "HR")
    )


    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )


    # sig = "Cholesterol3"
    # analysis = "T.AAA_or_A0A"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Arm_consolidated"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                                        test_var = sig, inter_var = inter_var)





        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Arm_consolidated) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Arm_consolidated = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Arm_consolidated"),
                    by = "Interaction_var_levels")

        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){
            print(l)

            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm
              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Arm_consolidated == l),
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$HR) <- c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
                         "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
                         "Global") # 3. Per subtype global interaction: all-Arm * Sig * pCR






# Pan-subtype(ALL) per-arm: Subtype, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Pan-Subtype (ALL) analysis >>>>>
# 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
# 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
# 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR

# Analysis structure
geo_inter$ALL[["AAA_plus_T.TN_or_HER2_or_HR"]] <- NA # main analysis
geo_inter$ALL[["AAA_plus_T.TN_or_HRnHER2_or_HRpHER2_or_HR"]] <- NA # main analysis
geo_inter$ALL[["AAA.TN_or_HR"]] <- NA # main analysis
geo_inter$ALL[["AOA_plus_T.TN_or_HR"]] <- NA # supplementary analysis



# Interaction summary
geo_inter$ALL <- purrr::map(

  #   Subtype       Arm   N
  # 1      TN     AAA+T 293
  # 2      TN     A0A+T 150
  # 3      TN       AAA  83

  # 4    HER2     AAA+T 135
  # 5    HER2 AAA+T+TRA  87

  # 6      HR     AAA+T 521
  # 7      HR     A0A+T 139
  # 8      HR       AAA  91


  c("AAA_plus_T.TN_or_HER2_or_HR", # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
    "AAA_plus_T.TN_or_HRnHER2_or_HRpHER2_or_HR", # 2. Subtype interaction in AAA+Taxane: TN/HR-HER2+/HR-HER2+/HR * Sig * pCR
    # Note on the analysis: "AAA_plus_T.TN_or_HRnHER2_or_HRpHER2_or_HR"
    # No filtering of TRA regimen containg datset with n=9 in HR+HER2+ required.
    # As only non-TRA regimens were considered in this analysis.
    "AAA.TN_or_HR", # 3. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
    "AOA_plus_T.TN_or_HR"), # 4. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR

  function(analysis, clin){

    print(analysis)

    xclin <- switch (
      analysis,

      # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
      AAA_plus_T.TN_or_HER2_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "AAA\\+T") &
                        str_detect(Arm_consolidated, "TRA", negate = TRUE)),

      # 2. Subtype interaction in AAA+Taxane: TN/HR-HER2+/HR-HER2+/HR * Sig * pCR
      AAA_plus_T.TN_or_HRnHER2_or_HRpHER2_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "AAA\\+T") &
                        str_detect(Arm_consolidated, "TRA", negate = TRUE)) %>%
        dplyr::mutate(Subtype_ihc = Subtype_ihc2), # to incorporate HR based HER2 subtypes

      # 3. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
      AAA.TN_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "AAA") & str_detect(Arm_consolidated, "AAA\\+T", negate = TRUE)),

      # 4. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR
      AOA_plus_T.TN_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "A0A\\+T"))
    )

    xclin <- xclin %>%
      dplyr::mutate(
        # Droping levels as some part of the scripts depends on relevant factor levels
        Subtype_ihc = Subtype_ihc %>% droplevels(),
        Arm_consolidated = Arm_consolidated %>% droplevels()
      )

    print(table(xclin$Subtype_ihc))

    # sig = "Immune1"
    # analysis = "AAA_plus_T.TN_or_HER2_or_HR"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c(
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

        "MCPcounter_Tcell",
        "MCPcounter_Cyto.Lymphocyte",
        "MCPcounter_B.Lineage",
        "MCPcounter_NK.Cells",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Neutrophils",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts",

        "scaled_Control_Proliferation",
        "scaled_Control_Non.breast.tissue",
        "scaled_Control_Human_Behaviour"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Subtype_ihc")),
          formula =as.formula(paste("Response ~", sig, "+ Subtype_ihc + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          # formula =as.formula(paste("Response ~", sig, "* Subtype_ihc")),
          formula =as.formula(paste("Response ~", sig, "* Subtype_ihc + Series_matrix_accession")),
          data = xclin,
          family = "binomial"
        )

        inter_var <- "Subtype_ihc"

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                                        test_var = sig, inter_var = inter_var)





        # Estimate per arm and pooled data heterogenity

        inter_var_levels <- xinter %>%
          dplyr::select(Module_name) %>%
          dplyr::mutate(
            Interaction_var_levels = str_split_fixed(string = Module_name, pattern = ":", n = 2)[, 2] %>%
              str_replace(pattern = inter_var, "")
          )


        # Update with per arm sample size and event/pcr rate to 'inter_var_levels'
        response_sum  <- xclin %>%
          dplyr::group_by(Subtype_ihc) %>% # interaction variable
          dplyr::summarise(N_response = which(Response == 1) %>% length(),
                           N_patients = n(),
                           Percent_repsonse = (N_response/N_patients) * 100 %>% round(digits = 2),
                           .groups	= "keep")

        # adding row for All samples
        response_sum  <- bind_rows(
          response_sum,
          tibble(
            Subtype_ihc = "",
            N_response = response_sum$N_response %>% sum(),
            N_patients = response_sum$N_patients %>% sum(),
            Percent_repsonse = (N_response/N_patients) * 100
          )
        )

        inter_var_levels <- inter_var_levels %>%
          left_join(response_sum %>%
                      dplyr::rename(Interaction_var_levels = "Subtype_ihc"),
                    by = "Interaction_var_levels")


        # Heterogenity assessment
        xhet <- purrr::map_dfr(
          inter_var_levels$Interaction_var_levels,
          function(l, xclin, sig){
            print(l)

            if (l == "") {
              # estmate heterogenity for the pooled dataset
              estimate_prog_heterogenity(
                xclin = xclin,
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            } else {
              # estmate heterogenity per arm
              estimate_prog_heterogenity(
                xclin = xclin %>%
                  dplyr::filter(Subtype_ihc == l), # inter_var is "Subtype_ihc"
                sig = sig
              ) %>%
                dplyr::mutate(
                  Interaction_var_levels = l
                )
            }

          },
          xclin,
          sig
        )

        xhet <- xhet %>%
          dplyr::select(-Module_name) %>%
          left_join(inter_var_levels, by = "Interaction_var_levels")



        # Interaction + Heterogenity
        xinter <- xinter %>%
          left_join(xhet, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::mutate(
        P_inter_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), P, NA_real_),
          method = "BH"),
        P_prog_adj = p.adjust(
          p = if_else(str_detect(Module_name, ":"), NA_real_, P),
          method = "BH")
      )

  },

  clin = clin_neoadj
)

names(geo_inter$ALL) <-  c("AAA_plus_T.TN_or_HER2_or_HR", # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
                           "AAA_plus_T.TN_or_HRnHER2_or_HRpHER2_or_HR", # 2. Subtype interaction in AAA+Taxane: TN/HR-HER2+/HR-HER2+/HR * Sig * pCR
                           "AAA.TN_or_HR", # 3. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
                           "AOA_plus_T.TN_or_HR") # 4. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR




#
# ==============================================================================




# 4. Analysis summary tables
# ==============================================================================

# geo_inter

summarize_geo_inter <- function(x){

  # x <- geo_inter$TN$AAA.T_or_NoT




  x <- purrr::map(x,
                  function(xx, prog){

                    xx <- xx %>%
                      dplyr::filter(is.na(P_prog_adj))


                      nme <- c("Module_name", "Interaction_var_levels", "N_response", "N_patients", "Percent_repsonse",
                               "Estimate", "l95", "u95", "P", "P_inter_adj", "Q", "Q_p", "I2")
                      xx <- xx[ , nme] %>%
                        dplyr::mutate(
                          Percent_repsonse = round(Percent_repsonse, digits = 1) %>% str_c("%"),
                          Log_OR = str_c(round(Estimate, digits = 1),
                                         " (",
                                         round(l95,digits=1), "-",
                                         round(u95,digits=1),
                                         ")"),
                          P = round(P, digits = 3),
                          P_inter_adj = round(P_inter_adj, digits = 3),
                          Q = round(Q, digits = 3),
                          Q_p = round(Q_p, digits = 3),
                          I2 = round(I2, digits = 3)
                          )

                      xx %>%
                        dplyr::select(c("Module_name", "Interaction_var_levels",
                                        "N_response", "N_patients", "Percent_repsonse",
                                        "Log_OR", "P", "P_inter_adj", "Q", "Q_p", "I2"))


                    })

  return(x)

}

xout <- list()

xout[["TN"]] <- summarize_geo_inter(x = geo_inter$TN)$AAA.T_or_NoT
xout[["HR"]] <- summarize_geo_inter(x = geo_inter$HR)$AAA.T_or_NoT

write_xlsx(xout, path = str_c(out_tables,"geo_taxane_interaction_summary.xlsx"))



#
# ==============================================================================



# 5. Save Robjects odf interaction summaries
# ==============================================================================


# geo_inter : Interaction summary
save(geo_inter, file = str_c(out_data,"geo_inter.RData")) # old name: inter_sum


#
# ==============================================================================





