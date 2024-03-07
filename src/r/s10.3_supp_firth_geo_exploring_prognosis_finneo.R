# s10.3_supp_firth_geo_exploring_prognosis_finneo.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Explore whether firth penalization reduces CI width



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and summarize clin_neoadj_finneo.
# 2. Explore prognosis.
# (standard logistic and firth penalized logistic).
# 3. Analysis summary plots



# 1. Load and format clin_neoadj_finneo.
# ==============================================================================

load("results/data/clin_neoadj_finneo.RData")




# Summarize clin_neoadj_finneo
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

#    Subtype  Hr                                                                  Arm   N
# 1     HER2 neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy  80
# 2     HER2 neg     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  58

# 3     HER2 pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy  55
# 4     HER2 pos     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy  29

# 5       HR pos   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 139
# 6       HR pos AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  91
# 7       HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 521

# 8       TN neg   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 150
# 9       TN neg AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy  83
# 10      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy 293


# Analysis
# 1. Per subtype prognosis (hr-arm-study heterogeneity)
# 2. Focus the analysis on AAA containing regimens due to large sample size.
# 2. Discard AOA regimen
# 3. Discard HER2 subtype (no AAA+/-T interaction analysis is performed in this subtype)
# 4. Script structure:
#     Create geo_inter (interaction summary) list object.
#     Each element represent all analysis done on each subtype.
#

# Analysis note 2:
# By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to clim a significant treatment
# interaction.


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

#    Subtype  Hr                                                                  Arm                 Strata   N Event Event_percet
# 1      HR pos   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:A0A+T  35     3     8.571429
# 2      HR pos   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE41998:A0A+T 104    12    11.538462
# 3      HR pos AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE20271:AAA  49     3     6.122449
# 4      HR pos AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE22093:AAA  42    10    23.809524
# 5      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20194:AAA+T 141     9     6.382979
# 6      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20271:AAA+T  44     3     6.818182
# 7      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE23988:AAA+T  32     7    21.875000
# 8      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:AAA+T 194    22    11.340206
# 9      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE32646:AAA+T  55     5     9.090909
# 10      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE42822:AAA+T  29     8    27.586207
# 11      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE50948:AAA+T  26     4    15.384615

# 12      TN neg   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:A0A+T  25     6    24.000000
# 13      TN neg   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE41998:A0A+T 125    51    40.800000
# 14      TN neg AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE20271:AAA  28     4    14.285714
# 15      TN neg AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE22093:AAA  55    18    32.727273
# 16      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20194:AAA+T  64    25    39.062500
# 17      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20271:AAA+T  31     9    29.032258
# 18      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE23988:AAA+T  29    13    44.827586
# 19      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:AAA+T 119    43    36.134454
# 20      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE32646:AAA+T  26    10    38.461538
# 21      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE42822:AAA+T  24    12    50.000000

# Note !!!!!
# The event rate in HR subtype is low possibily due to lack of hormone therapy


#
# ==============================================================================





# 2. Explore prognosis.
# ==============================================================================

# Per subtype prognosis (Series matrix + Arm (~ Strata) heterogeneity)


geo_prog <- purrr::map(

  c("TN", "HR"),


  function(subtype, clin){

    print(subtype)

    xclin <- clin %>%
        dplyr::filter(Subtype_ihc2 == subtype)

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
                              " + Arm_consolidated + Series_matrix_accession"),
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

names(geo_prog) <- c("TN", "HR")



#
# ==============================================================================



# 4. Analysis summary tables
# ==============================================================================


# prognostic model comparison

ggdf <- bind_rows(
  geo_prog$HR %>% dplyr::mutate(subtype = "HR"),
  geo_prog$TN %>% dplyr::mutate(subtype = "TN")
) %>%
  tidyr::gather("key", "value",
                "or", "ci.lo", "ci.up",
                "firth_or", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = geo_prog$HR$variable),
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
       subtitle = "GEO prognostic models") +
  facet_wrap(~subtype, ncol = 1)

pdf(file = str_c(out_figures,"GEO_Prognosis_Firth_vs_Std_Logistic.pdf"))
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
       subtitle = "GEO prognostic models") +
  facet_wrap(~subtype, ncol = 1)

pdf(file = str_c(out_figures,"GEO_Prognosis_log_Firth_vs_Std_Logistic.pdf"))
print(p)
dev.off()



#
# ==============================================================================





# Clear memory
# ==============================================================================

rm(clin_neoadj_finneo, geo_prog)

#
# ==============================================================================


