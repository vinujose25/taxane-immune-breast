# s9.3_supp_firth_geo_exploring_interaction_finneo.R



# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Explore whether firth penalization reduces CI width



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load data.
# 2. Structure analysis
# 3. Explore chemo interaction
#   (standard logistic and firth penalized logistic).
# 4. Analysis summary figures




# 1. Load data
# ==============================================================================

load("results/data/clin_neoadj_finneo.RData")

dim(clin_neoadj_finneo)
# 1499  134


# clean Arm and factor levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
clin_neoadj_finneo$Arm_consolidated %>% factor() %>% levels()
# [1] "A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [2] "AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [3] "AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy"
# [4] "AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy"

clin_neoadj_finneo <- clin_neoadj_finneo %>%
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

clin_neoadj_finneo$Arm_consolidated %>% factor() %>% levels()
# [1] "AAA+T"     "A0A+T"     "AAA"       "AAA+T+TRA"



# set subtype ihc levels
# >>>>>>>>>>>>>>>>>>>>>>
clin_neoadj_finneo$Subtype_ihc %>% factor() %>% levels()
# [1] "HER2" "HR"   "TN"
clin_neoadj_finneo$Subtype_ihc2 %>% factor() %>% levels()
# [1] "HR"       "HR-HER2+" "HR+HER2+" "TN"

clin_neoadj_finneo <- clin_neoadj_finneo %>%
  dplyr::mutate(
    Subtype_ihc =  Subtype_ihc %>%
      factor(levels = c("TN", "HER2", "HR")),
    Subtype_ihc2 =  Subtype_ihc2 %>%
      factor(levels = c("TN", "HR-HER2+", "HR+HER2+", "HR"))
  )

clin_neoadj_finneo$Subtype_ihc %>% levels()
# [1] "TN" "HER2" "HR"
clin_neoadj_finneo$Subtype_ihc2 %>% levels()
# [1] "TN" "HR-HER2+" "HR+HER2+" "HR"
#
# ==============================================================================




# 2. Structure analysis
# ==============================================================================


# Summarize clin_neoadj
# >>>>>>>>>>>>>>>>>>>>>

# All neoadj regimen
clin_neoadj_finneo %>%
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




# Analysis note: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# 1. Focus the analysis on AAA containing regimens due to large sample size.
# 2. Discard AOA regimen
# 3. Discard HER2 subtype (no AAA+/-T interaction analysis possible)
# 4. Script structure:
#     Create geo_inter (interaction summary) list object.
#     Each element represent all analysis done on each subtype.
#




# All neoadj regimen with Strata info
clin_neoadj_finneo %>%
  dplyr::filter(Subtype_ihc != "HER2") %>%
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
# 7       TN neg     A0A+T     GSE25066:A0A+T  25     6    24.000000 # not used in analysis
# 8       TN neg     A0A+T     GSE41998:A0A+T 125    51    40.800000 # not used in analysis
# 9       TN neg       AAA       GSE20271:AAA  28     4    14.285714
# 10      TN neg       AAA       GSE22093:AAA  55    18    32.727273

# 11      HR pos     AAA+T     GSE20194:AAA+T 141     9     6.382979
# 12      HR pos     AAA+T     GSE20271:AAA+T  44     3     6.818182
# 13      HR pos     AAA+T     GSE23988:AAA+T  32     7    21.875000
# 14      HR pos     AAA+T     GSE25066:AAA+T 194    22    11.340206
# 15      HR pos     AAA+T     GSE32646:AAA+T  55     5     9.090909
# 16      HR pos     AAA+T     GSE42822:AAA+T  29     8    27.586207
# 17      HR pos     AAA+T     GSE50948:AAA+T  26     4    15.384615
# 18      HR pos     A0A+T     GSE25066:A0A+T  35     3     8.571429 # not used in analysis
# 19      HR pos     A0A+T     GSE41998:A0A+T 104    12    11.538462 # not used in analysis
# 20      HR pos       AAA       GSE20271:AAA  49     3     6.122449
# 21      HR pos       AAA       GSE22093:AAA  42    10    23.809524


# Note !!!!!
# The event rate in HR subtype is low possibily due to lack of hormone therapy

#
# ==============================================================================




# 3.Explore chemo interaction (standard logistic and firth penalized logistic).
# ==============================================================================


geo_inter <- vector(mode = "list", length = 2)
names(geo_inter) <-c("TN", "HR")


# Per-subtype: Arm, Signature, and pCR interaction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# TN analysis: individual.sig
# @@@@@@@@@@@@@@@@@@@@@@@@@@@



# Analysis structure
geo_inter$TN[["AAA.T_or_NoT"]] <- NA


geo_inter$TN <- purrr::map(

    "AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR

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

    xx <- purrr::map(

      c(
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

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Fibroblasts"
      ),

      function(sig, xclin){

        print(sig)


        get_std_firth_logist_out(
          formula_chr = str_c("Response ~ ", sig,
                              " * Arm_consolidated + Series_matrix_accession"),
          biomarker = sig,
          xdata = xclin
        )


      },

      xclin
    )

    xx %>% bind_rows()

  },

  clin = clin_neoadj_finneo
)


names(geo_inter$TN) <- "AAA.T_or_NoT"



# HR analysis:individual.sig
# @@@@@@@@@@@@@@@@@@@@@@@@@@


# Analysis structure
geo_inter$HR[["AAA.T_or_NoT"]] <- NA # main analysis


# Interaction summary
geo_inter$HR <- purrr::map(

  "AAA.T_or_NoT", # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR

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

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Fibroblasts"
      ),

      function(sig, xclin){

        print(sig)


        get_std_firth_logist_out(
          formula_chr = str_c("Response ~ ", sig,
                              " * Arm_consolidated + Series_matrix_accession"),
          biomarker = sig,
          xdata = xclin
        )


      },

      xclin
    )


    xx %>% bind_rows()

  },

  clin = clin_neoadj_finneo
)

names(geo_inter$HR) <- "AAA.T_or_NoT" # 1. Per subtype taxane interaction: AAA+/-Taxane * Sig * pCR


#
# ==============================================================================



# 4. Analysis summary plots
# ==============================================================================


# interaction model comparison

ggdf <- bind_rows(
  geo_inter$HR$AAA.T_or_NoT %>% dplyr::mutate(subtype = "HR"),
  geo_inter$TN$AAA.T_or_NoT %>% dplyr::mutate(subtype = "TN")
) %>%
  tidyr::gather("key", "value",
                "or", "ci.lo", "ci.up",
                "firth_or", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = geo_inter$HR$AAA.T_or_NoT$variable),
                group_var = factor(group_var,
                                   levels = c("ci.lo-Std", "ci.lo-Firth",
                                              "or-Std", "or-Firth",
                                              "ci.up-Std", "ci.up-Firth")))


p <- ggplot(ggdf, aes(x = variable, y = value, group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Logistic model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Gene-modules", y = "Estimates (OR and 95% CI)",
       title = "Firth vs Standard Logistic model estimates",
       subtitle = "GEO interaction models") +
  facet_wrap(~subtype, ncol = 1)

pdf(file = str_c(out_figures,"GEO_Interaction_Firth_vs_Std_Logistic.pdf"))
print(p)
dev.off()


p <- ggplot(ggdf, aes(x = variable, y = log(value), group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Logistic model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Gene-modules", y = "Log - Estimates (OR and 95% CI)",
       title = "Firth vs Standard Logistic model estimates",
       subtitle = "GEO interaction models") +
  facet_wrap(~subtype, ncol = 1)

pdf(file = str_c(out_figures,"GEO_Interaction_log_Firth_vs_Std_Logistic.pdf"))
print(p)
dev.off()

#
# ==============================================================================



# Clear memory
# ==============================================================================

rm(clin_neoadj_finneo, geo_inter)

#
# ==============================================================================
