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
# 4. Save Robjects odf interaction summaries




# 1. Load data
# ==============================================================================

load("results/data/clin_neoadj.RData")
dim(clin_neoadj)
# 1499  152


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
    Arm = Arm_consolidated
  ) %>%
  dplyr::summarise(N = n()) %>%
  as.data.frame()

#   Subtype                                                                  Arm   N
# 1    HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 135
# 2    HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87

# 3      HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 139
# 4      HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
# 5      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521

# 6      TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
# 7      TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
# 8      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 293


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
# 1. Focus the analysis on AAA containing regimens due to latge sample size.
# 2. Make plots accordingly.
# 3. Keep the other analysis as supplimentary.
# 4. Script structure:
#     Create inter_sum (interaction summary) list object.
#     Each elemet represent all anlysis done on each subtype.
#



# All neoadj regimen with Strata info
clin_neoadj %>%
  dplyr::group_by(
    Subtype = Subtype_ihc,
    Arm = Arm_consolidated,
    Strata # identical results when Strata =  Series_matrix_accession
  ) %>%
  dplyr::summarise(N = n(), Event = which(Response == 1) %>% length()) %>%
  dplyr::mutate(Event_percet = (Event/N) * 100) %>%
  as.data.frame()
#    Subtype                                                                  Arm             Strata   N Event Event_percet
# 1     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20194:AAA+T  50    17    34.000000
# 2     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE32646:AAA+T  34    12    35.294118
# 3     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE50948:AAA+T  51    13    25.490196
# 4     HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE42822:AAA+T+TRA  24    12    50.000000
# 5     HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE50948:AAA+T+TRA  63    31    49.206349
# 6       HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE25066:A0A+T  35     3     8.571429
# 7       HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE41998:A0A+T 104    12    11.538462
# 8       HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy       GSE20271:AAA  49     3     6.122449
# 9       HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy       GSE22093:AAA  42    10    23.809524
# 10      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20194:AAA+T 141     9     6.382979
# 11      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20271:AAA+T  44     3     6.818182
# 12      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE23988:AAA+T  32     7    21.875000
# 13      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE25066:AAA+T 194    22    11.340206
# 14      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE32646:AAA+T  55     5     9.090909
# 15      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE42822:AAA+T  29     8    27.586207
# 16      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE50948:AAA+T  26     4    15.384615
# 17      TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE25066:A0A+T  25     6    24.000000
# 18      TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE41998:A0A+T 125    51    40.800000
# 19      TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy       GSE20271:AAA  28     4    14.285714
# 20      TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy       GSE22093:AAA  55    18    32.727273
# 21      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20194:AAA+T  64    25    39.062500
# 22      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20271:AAA+T  31     9    29.032258
# 23      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE23988:AAA+T  29    13    44.827586
# 24      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE25066:AAA+T 119    43    36.134454
# 25      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE32646:AAA+T  26    10    38.461538
# 26      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE42822:AAA+T  24    12    50.000000

# Note !!!!!
# The event rate in HR subtype is low possibily due to lack of hormone therapy


#
# ==============================================================================




# 3. Interaction tests
# ==============================================================================

inter_sum <- vector(mode = "list", length = 4)
names(inter_sum) <- c("TN", "HER2", "HR", "ALL")


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
inter_sum$TN[["AAA.T_or_NoT"]] <- NA # main analysis
inter_sum$TN[["T.AAA_or_A0A"]] <- NA # supplemenatry analysis
inter_sum$TN[["Global"]] <- NA # supplemenatry analysis


# Interaction summary
inter_sum$TN <- purrr::map(

  #   Subtype                                                                  Arm   N
  # 6      TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
  # 7      TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
  # 8      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 293

  c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
    "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
    "Global"), # 3. Per subtype global interaction: all-Arm * Sig * pCR

  function(analysis, clin){

    print(analysis)

    xclin <- switch (
      analysis,

      # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
      AAA.T_or_NoT = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "AAA")),

      # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
      T.AAA_or_A0A = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "\\+Taxane")),

      # 3. Per subtype global interaction: all-Arm * Sig * pCR
      Global = clin %>%
        dplyr::filter(Subtype_ihc == "TN")
    )

    # sig = "Immune2"
    # analysis = "Global"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        "Cholesterol1", "Cholesterol2", "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", "CLymphocyte", "Bcell",
        "NKcell", "Monocyte", "MDendritic", "Fibroblast"),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
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

names(inter_sum$TN) <- c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
                         "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
                         "Global") # 3. Per subtype global interaction: all-Arm * Sig * pCR



# HER2 analysis
# @@@@@@@@@@@@@

# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

# Analysis structure
inter_sum$HER2[["AAA_plus_T.TRA_or_NoTRA"]] <- NA # main analysis


# Interaction summary
inter_sum$HER2 <- purrr::map(

  #   Subtype                                                                  Arm   N
  # 1    HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 135
  # 2    HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87

  c("AAA_plus_T.TRA_or_NoTRA"), # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

  function(analysis, clin){

    print(analysis)

    # 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR
    xclin <- clin %>%
      dplyr::filter(Subtype_ihc == "HER2" &
                      str_detect(Arm_consolidated, "AAA"))

    # sig = "Immune2"
    # analysis = "AAA_plus_T.TRA_or_NoTRA"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        "Cholesterol1", "Cholesterol2", "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", "CLymphocyte", "Bcell",
        "NKcell", "Monocyte", "MDendritic", "Fibroblast"),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
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

names(inter_sum$HER2) <- "AAA_plus_T.TRA_or_NoTRA"





# HR analysis
# @@@@@@@@@@@

# 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 3. Per subtype global interaction: all-Arm * Sig * pCR


# Analysis structure
inter_sum$HR[["AAA.T_or_NoT"]] <- NA # main analysis
inter_sum$HR[["T.AAA_or_A0A"]] <- NA # supplemenatry analysis
inter_sum$HR[["Global"]] <- NA # supplemenatry analysis


# Interaction summary
inter_sum$HR <- purrr::map(

  #   Subtype                                                                  Arm   N
  # 3      HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 139
  # 4      HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
  # 5      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521


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
                        str_detect(Arm_consolidated, "\\+Taxane")),

      # 3. Per subtype global interaction: all-Arm * Sig * pCR
      Global = clin %>%
        dplyr::filter(Subtype_ihc == "HR")
    )

    # sig = "Cholesterol3"
    # analysis = "T.AAA_or_A0A"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        "Cholesterol1", "Cholesterol2", "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", "CLymphocyte", "Bcell",
        "NKcell", "Monocyte", "MDendritic", "Fibroblast"),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Arm_consolidated")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "* Arm_consolidated")),
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

names(inter_sum$HR) <- c("AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
                         "T.AAA_or_A0A", # 2. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
                         "Global") # 3. Per subtype global interaction: all-Arm * Sig * pCR






# Pan-subtype(ALL) per-arm: Subtype, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Pan-Subtype (ALL) analysis >>>>>
# 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
# 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
# 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR

# Analysis structure
inter_sum$ALL[["AAA_plus_T.TN_or_HER2_or_HR"]] <- NA # main analysis
inter_sum$ALL[["AAA.TN_or_HR"]] <- NA # main analysis
inter_sum$ALL[["AOA_plus_T.TN_or_HR"]] <- NA # supplementary analysis



# Interaction summary
inter_sum$ALL <- purrr::map(

  #   Subtype                                                                  Arm   N
  # 1    HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 135
  # 2    HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87

  # 3      HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 139
  # 4      HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
  # 5      HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521

  # 6      TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
  # 7      TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
  # 8      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 293


  c("AAA_plus_T.TN_or_HER2_or_HR", # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
    "AAA.TN_or_HR", # 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
    "AOA_plus_T.TN_or_HR"), # 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR

  function(analysis, clin){

    print(analysis)

    xclin <- switch (
      analysis,

      # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
      AAA_plus_T.TN_or_HER2_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "AAA\\+Taxane") &
                        str_detect(Arm_consolidated, "Trastuzumab", negate = TRUE)),

      # 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
      AAA.TN_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "AAA\\+noTaxane")),

      # 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR
      AOA_plus_T.TN_or_HR = clin %>%
        dplyr::filter(str_detect(Arm_consolidated, "A0A\\+Taxane"))
    )

    print(table(xclin$Subtype_ihc))

    # sig = "Cholesterol3"
    # analysis = "AAA_plus_T.TN_or_HER2_or_HR"
    # clin <- clin_neoadj

    # Per sig prognosis + heterogenity
    xx <- purrr::map(

      c("Immune1", "Immune2", "Immune3",
        "Interferon1", "Interferon2", "Interferon3",
        "Cholesterol1", "Cholesterol2", "Cholesterol3",
        "Fibrosis1", "Fibrosis2", "Fibrosis3",
        "Proliferation1", "Proliferation2", "Proliferation3",
        "Tcell", "CLymphocyte", "Bcell",
        "NKcell", "Monocyte", "MDendritic", "Fibroblast"),

      function(sig, xclin){

        print(sig)

        # Estimate interaction and prognosis
        m0_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Subtype_ihc")),
          data = xclin,
          family = "binomial"
        )

        m1_inter <- glm(
          formula =as.formula(paste("Response ~", sig, "* Subtype_ihc")),
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

names(inter_sum$ALL) <-   c("AAA_plus_T.TN_or_HER2_or_HR", # 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
                            "AAA.TN_or_HR", # 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
                            "AOA_plus_T.TN_or_HR") # 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR




#
# ==============================================================================



# 4. Save Robjects odf interaction summaries
# ==============================================================================


# inter_sum : Interaction summary
save(inter_sum, file = str_c(out_data,"inter_sum.RData"))


#
# ==============================================================================





