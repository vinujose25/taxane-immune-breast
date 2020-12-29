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

#    Subtype                                                                  Arm   N
# 1     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 161
# 2     HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87

# 3       HR   0A0+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy  60
# 4       HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 167
# 5       HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
# 6       HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521

# 7       TN A00+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  59
# 8       TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
# 9       TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
# 10      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 309


# TN analysis
# 1. Per subtype global interaction: all-Arm * Sig * pCR
# 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 3. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 4. Per subtype chemo-class-combination interaction: (AAA/A00)+noTaxane * Sig * pCR


# HER2 analysis
# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR


# HR analysis
# 1. Per subtype global interaction: all-Arm * Sig * pCR
# 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 3. Per subtype chemo-class-combination interaction: (AAA/A0A/0A0)+Taxane * Sig * pCR


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
  as.data.frame()

#
# ==============================================================================




# 3. Interaction tests
# ==============================================================================

inter_sum <- vector(mode = "list", length = 4)
names(inter_sum) <- c("TN", "HER2", "HR", "ALL")


# Per-subtype: Arm, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# TN analysis >>>>>
# 1. Per subtype global interaction: all-Arm * Sig * pCR
# 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 3. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
# 4. Per subtype chemo-class-combination interaction: (AAA/A00)+noTaxane * Sig * pCR

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
inter_sum$TN[["NoT.AAA_or_A00"]] <- NA # supplemenatry analysis
inter_sum$TN[["Global"]] <- NA # supplemenatry analysis


# Interaction summary
inter_sum$TN <- purrr::map(

  # 7       TN A00+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  59
  # 8       TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
  # 9       TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
  # 10      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 309

  c("AAA.T_or_NoT", # 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
    "T.AAA_or_A0A", # 3. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
    "NoT.AAA_or_A00", # 4. Per subtype chemo-class-combination interaction: (AAA/A00)+noTaxane * Sig * pCR
    "Global"), # 1. Per subtype global interaction: all-Arm * Sig * pCR

  function(analysis, clin){


    print(analysis)

    xclin <- switch (
      analysis,

      # 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
      AAA.T_or_NoT = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "AAA")),

      # 3. Per subtype chemo-class-combination interaction: (AAA/A0A)+Taxane * Sig * pCR
      T.AAA_or_A0A = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "\\+Taxane")),

      # 4. Per subtype chemo-class-combination interaction: (AAA/A00)+noTaxane * Sig * pCR
      NoT.AAA_or_A00 = clin %>%
        dplyr::filter(Subtype_ihc == "TN" &
                        str_detect(Arm_consolidated, "\\+noTaxane")),

      # 1. Per subtype global interaction: all-Arm * Sig * pCR
      Global = clin
    )

    # sig = "Immune1"
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

        xinter <- summarise_interaction(m1 = m1_inter, m0 = m0_inter,
                              test_var = sig, inter_var = "Arm_consolidated")





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

# There were 24 warnings (use warnings() to see them)
# Warning messages:
# 1: glm.fit: algorithm did not converge
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
# 3: glm.fit: algorithm did not converge

names(prog_sum) <- c("TN", "HER2", "HR", "ALL")


# !!!! Save prog_sum later in the script









# HER2 analysis >>>>>
# 1. Per subtype trastuzumab interaction: AAA+Taxane+/-Trastuzumab * Sig * pCR

inter_sum$HER2[["AAA_plus_T.TRA_or_NoTRA"]] <- NA # main analysis

#    Subtype                                                                  Arm   N
# 1     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 161
# 2     HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87





# HR analysis >>>>>
# 1. Per subtype global interaction: all-Arm * Sig * pCR
# 2. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR
# 3. Per subtype chemo-class-combination interaction: (AAA/A0A/0A0)+Taxane * Sig * pCR


inter_sum$HR[["AAA.T_or_NoT"]] <- NA # main analysis
inter_sum$HR[["T.AAA_or_A0A_or_0A0"]] <- NA # supoplementary analysis
inter_sum$HR[["Global"]] <- NA # supplementary analysis


# 3       HR   0A0+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy  60
# 4       HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 167
# 5       HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
# 6       HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521




# Pan-subtype(ALL) per-arm: Subtype, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Pan-Subtype (ALL) analysis >>>>>
# 1. Subtype interaction in AAA+Taxane: TN/HER2/HR * Sig * pCR
# 2. Subtype interaction in AAA+noTaxane: TN/HR * Sig * pCR
# 3. Subtype interaction in A0A+Taxane: TN/HR * Sig * pCR

inter_sum$ALL[["AAA_plus_T.TN_or_HER2_or_HR"]] <- NA # main analysis
inter_sum$ALL[["AAA.TN_or_HR"]] <- NA # main analysis
inter_sum$ALL[["AOA_plus_T.TN_or_HR"]] <- NA # supplementary analysis


#    Subtype                                                                  Arm   N
# 1     HER2   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 161
# 2     HER2     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  87

# 3       HR   0A0+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy  60
# 4       HR   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 167
# 5       HR AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
# 6       HR   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521

# 7       TN A00+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  59
# 8       TN   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
# 9       TN AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
# 10      TN   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 309


#
# ==============================================================================



# 4. Save Robjects odf interaction summaries
# ==============================================================================


#
# ==============================================================================





