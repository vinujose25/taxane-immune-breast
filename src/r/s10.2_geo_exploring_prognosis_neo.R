# s10.2_geo_exploring_prognosis_neo.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the prognostic value of TIL associated signatures and Immune celltypes.
# Account for study heterogeneity.
# Prepare data necessary to plot prognostic plots.



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and summarize clin_neoadj_neo.
# 2. Explore prognosis.
# 3. Explore prognosis heterogeneity.
# 4. Analysis summary
# 5. Save Robjects



# 1. Load and format clin_neoadj_neo.
# ==============================================================================

load("results/data/clin_neoadj_neo.RData")




# Summarize clin_neoadj_neo
# >>>>>>>>>>>>>>>>>>>>>


# All neoadj regimen
clin_neoadj_neo %>%
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
clin_neoadj_neo %>%
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
# The event rate in HR subtype is low possibly due to lack of hormone therapy


#
# ==============================================================================





# 2. Explore prognosis.
# ==============================================================================

# Per subtype prognosis (Series matrix + Arm (~ Strata) heterogeneity)


# Reference
# >>>>>>>>>

# altmeta::metahet(y, s2, n.resam = 1000)
# y	: a numeric vector indicating the observed effect sizes in the collected studies; they are assumed to be normally distributed.
# s2	: a numeric vector indicating the within-study variances.#
# n.resam	: a positive integer indicating the number of resampling iterations for calculating p-values of test statistics and 95% confidence interval of heterogeneity measures.

# Reference: Negative value of I2 (I^2) set as zero
# https://www.researchgate.net/post/How-should-I-consider-a-negative-value-of-the-I-squared-during-the-heterogeneity-assessment-in-meta-analyses




# Analysis/plot structure
# >>>>>>>>>>>>>>>>>>>>>>>

# Basic prognosis plot structure ( per subtype - per sig - per arm - study heterogeneity)
# (!!! Note : This plot will be generated in the NEXT SECTION of the script)
#
# TN - sig1 - arm1 - prognosis - study heterogeneity
# TN - sig1 - arm2 - prognosis - study heterogeneity
# TN - sig1 - arm3 - prognosis - study heterogeneity
# TN - sig1 - arm4 - prognosis - study heterogeneity
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# TN - sig1 - armX - prognosis - study and arm heterogeneity (global)
# Note:
# In prognosis forest plot include the global prognosis (pooled arm) for sig1 per subtype.
# The global analysis account for study and arm heterogeneity,
# and it also assess the heterogeneity level.


# Higher level prognosis plot structure
# TN/HER2/HR/ALL (event summary for each subtype, event/size: e/n)
# (!!! Note : This plot will be generated in the THIS SECTION of the script)
#
# immune sig1 - armX - prognosis - study and arm heterogeneity (global)
# immune sig2 - armX - prognosis - study and arm heterogeneity (global)
# immune sig3 - armX - prognosis - study and arm heterogeneity (global)
# interferon sig1 - armX - prognosis - study and arm heterogeneity (global)
# interferon sig2 - armX - prognosis - study and arm heterogeneity (global)
# interferon sig3 - armX - prognosis - study and arm heterogeneity (global)
# fibrosis sig1 - armX - prognosis - study and arm heterogeneity (global)
# fibrosis sig2 - armX - prognosis - study and arm heterogeneity (global)
# fibrosis sig3 - armX - prognosis - study and arm heterogeneity (global)
# cholesterol sig1 - armX - prognosis - study and arm heterogeneity (global)
# cholesterol sig2 - armX - prognosis - study and arm heterogeneity (global)
# cholesterol sig3 - armX - prognosis - study and arm heterogeneity (global)
# proliferation sig1 - armX - prognosis - study and arm heterogeneity (global)
# proliferation sig2 - armX - prognosis - study and arm heterogeneity (global)
# proliferation sig3 - armX - prognosis - study and arm heterogeneity (global)


geo_prog <- vector(mode = "list", length = 2)
names(geo_prog) <- c("individual.sig", "pooled.sig")

geo_prog <- purrr::map(
  geo_prog,
  ~{
    x <- vector(mode = "list", length = 2)
    names(x) <-  c("TN", "HR")
    x
  }
)



# Get prognostic effect and heterogenity measurement
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Prognosis: individual sig
geo_prog$individual.sig <- purrr::map(

  c("TN", "HR"),


  function(subtype, clin){

    print(subtype)

    xclin <- switch (
      subtype,

      # ALL = clin,
      #
      # HER2 = clin %>%
      #   dplyr::filter(Subtype_ihc == subtype),
      #
      # "HR+HER2+" = clin %>%
      #   dplyr::filter(Subtype_ihc2 == subtype & Strata != "GSE42822:AAA+T+TRA"),
      # # The above strata has only 9 samples with 2 events
      # # If included in analysis, the above strata will throw the following warnings.
      # # Warning messages:
      # # 1: glm.fit: algorithm did not converge
      # # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred

      clin %>%
        dplyr::filter(Subtype_ihc2 == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
    # clin <- clin_neoadj_neo

    # Per sig prognosis + heterogenity
    xx <- purrr::map(


      c(
        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Gruosso2019_Immune",
        "scaled_Hamy2016_Immune",
        "scaled_Teschendorff2007_Immune",
        "scaled_Yang2018_Immune",
        "scaled_Desmedt2008_STAT1",

        "scaled_Gruosso2019_Interferon",
        "scaled_Farmer2009_MX1",
        "scaled_Hamy2016_Interferon",
        "scaled_Nirmal2018_Interferon",

        "scaled_Gruosso2019_Cholesterol",
        "scaled_Ehmsen2019_Chol",
        "scaled_Simigdala2016_Chol",
        "scaled_Sorrentino2014_Chol",

        "scaled_Gruosso2019_Fibrosis",
        "scaled_Hamy2016_Ecm",
        "scaled_Naba2014_Ecmcore",
        "scaled_Triulzi2013_Ecm",

        "MCPcounter_Tcell",
        "MCPcounter_B.Lineage",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts"
        ),

      function(sig, xclin){

        print(sig)

        # Estimate prognosis
        m1_prog <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Strata")),
          # Note that in the formula the term "Strata" combines dataset and arm effect.
          # Although, using arm and dataset term seperatly would distinguishes the arm amd dataset effect, we
          # are not interesed in it at this stage.
          # This model tests the prognostic value of sig irrespective of arm and dataset
          # Alternative formula:
          formula = as.formula(
            switch(subtype,
                   # For HER2 adjust for HR status
                   "HER2" = paste("Response ~", sig, "+ Hr + Arm_consolidated + Series_matrix_accession"),
                   "ALL" = paste("Response ~", sig, "+ Subtype_ihc2 + Arm_consolidated + Series_matrix_accession"),
                   paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession"))
          ),
          # formula =as.formula(paste("Response ~", sig, "+ strata(Strata)")),
          # Both model will give identical p value

          data = xclin,
          family = "binomial"
        )
        # summary(m1_prog)$coefficients

        m1_prog_clean <- cbind(
          summary(m1_prog)$coefficients[sig, , drop = F],
          confint.default(m1_prog)[sig, , drop = F]
        ) %>%
          as_tibble(rownames = "Module_name")


        # Estimate heterogeneity
        # Note that heterogeneity assessed w.r.t each strata
        # For each strata, the prognostic value is assessed and the heterogeneity is measured
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig)


        # Prognosis + Heterogeneity
        m1_prog_clean %>%
          left_join(m1_het_clean, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::rename(
        Std_error = "Std. Error",
        Z_value = "z value",
        P = "Pr(>|z|)",
        Lower_ci = "2.5 %",
        Upper_ci = "97.5 %"
      ) %>%
      dplyr::mutate(
        P_adj = p.adjust(p = P, method = "BH")
      )

  },

  clin = clin_neoadj_neo
)

names(geo_prog$individual.sig) <- c("TN", "HR")



# Prognosis: pooled sig

geo_prog$pooled.sig <- purrr::map(

  c("TN", "HR"),


  function(subtype, clin){

    print(subtype)

    xclin <- switch (
      subtype,

      # ALL = clin,
      #
      # HER2 = clin %>%
      #   dplyr::filter(Subtype_ihc == subtype),
      #
      # "HR+HER2+" = clin %>%
      #   dplyr::filter(Subtype_ihc2 == subtype & Strata != "GSE42822:AAA+T+TRA"),
      # # The above strata has only 9 samples with 2 events
      # # If included in analysis, the above strata will throw the following warnings.
      # # Warning messages:
      # # 1: glm.fit: algorithm did not converge
      # # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred

      clin %>%
        dplyr::filter(Subtype_ihc2 == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
    # clin <- clin_neoadj_neo

    # Per sig prognosis + heterogenity
    xx <- purrr::map(


      c(
        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Tcell",
        "MCPcounter_B.Lineage",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts"
      ),

      function(sig, xclin){

        print(sig)

        # Estimate prognosis
        m1_prog <- glm(
          # formula =as.formula(paste("Response ~", sig, "+ Strata")),
          # Note that in the formula the term "Strata" combines dataset and arm effect.
          # Although, using arm and dataset term seperatly would distinguishes the arm amd dataset effect, we
          # are not interesed in it at this stage.
          # This model tests the prognostic value of sig irrespective of arm and dataset
          # Alternative formula:
          formula = as.formula(
            switch(subtype,
                   # For HER2 adjust for HR status
                   "HER2" = paste("Response ~", sig, "+ Hr + Arm_consolidated + Series_matrix_accession"),
                   "ALL" = paste("Response ~", sig, "+ Subtype_ihc2 + Arm_consolidated + Series_matrix_accession"),
                   paste("Response ~", sig, "+ Arm_consolidated + Series_matrix_accession"))
          ),
          # formula =as.formula(paste("Response ~", sig, "+ strata(Strata)")),
          # Both model will give identical p value

          data = xclin,
          family = "binomial"
        )
        # summary(m1_prog)$coefficients

        m1_prog_clean <- cbind(
          summary(m1_prog)$coefficients[sig, , drop = F],
          confint.default(m1_prog)[sig, , drop = F]
        ) %>%
          as_tibble(rownames = "Module_name")


        # Estimate heterogeneity
        # Note that heterogeneity assessed w.r.t each strata
        # For each strata, the prognostic value is assessed and the heterogeneity is measured
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig)


        # Prognosis + Heterogeneity
        m1_prog_clean %>%
          left_join(m1_het_clean, by = "Module_name")

      },

      xclin
    )


    # Consolidate per sig stat and format
    bind_rows(xx) %>%
      dplyr::rename(
        Std_error = "Std. Error",
        Z_value = "z value",
        P = "Pr(>|z|)",
        Lower_ci = "2.5 %",
        Upper_ci = "97.5 %"
      ) %>%
      dplyr::mutate(
        P_adj = p.adjust(p = P, method = "BH")
      )

  },

  clin = clin_neoadj_neo
)

names(geo_prog$pooled.sig) <- c("TN", "HR")


# !!!! Save geo_prog late r in the script

#
# ==============================================================================





# 3. Explore prognosis heterogeneity.
# ==============================================================================


geo_prog_het <- vector(mode = "list", length = 2)
names(geo_prog_het) <- c("individual.sig", "pooled.sig")

geo_prog_het <- purrr::map(
  geo_prog_het,
  ~{
    x <- vector(mode = "list", length = 2)
    names(x) <-  c("TN", "HR")
    x
  }
)




# Explore per arm prognostic effect and study heterogeneity
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Per-subtype, per-sig strata (study-arm) heterogeneity assessment


# Individual sig >>>>>>>>>>>

geo_prog_het$individual.sig <- purrr::map(

  c("TN", "HR"),

  function(subtype, clin){

    print(subtype)

    xclin <- switch (
      subtype,

      ALL = clin,

      HER2 = clin %>%
        dplyr::filter(Subtype_ihc == subtype),

      "HR+HER2+" = clin %>%
        dplyr::filter(Subtype_ihc2 == subtype & Strata != "GSE42822:AAA+T+TRA"),
      # The above strata has only 9 samples with 2 events
      # If included in analysis, the above strata will throw the following warnings.
      # Warning messages:
      # 1: glm.fit: algorithm did not converge
      # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred

      clin %>%
        dplyr::filter(Subtype_ihc2 == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
    # clin <- clin_neoadj_neo

    # Per-subtype, per-sig heterogenity assesment
    xx <- purrr::map(

      c(
        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Gruosso2019_Immune",
        "scaled_Hamy2016_Immune",
        "scaled_Teschendorff2007_Immune",
        "scaled_Yang2018_Immune",
        "scaled_Desmedt2008_STAT1",

        "scaled_Gruosso2019_Interferon",
        "scaled_Farmer2009_MX1",
        "scaled_Hamy2016_Interferon",
        "scaled_Nirmal2018_Interferon",

        "scaled_Gruosso2019_Cholesterol",
        "scaled_Ehmsen2019_Chol",
        "scaled_Simigdala2016_Chol",
        "scaled_Sorrentino2014_Chol",

        "scaled_Gruosso2019_Fibrosis",
        "scaled_Hamy2016_Ecm",
        "scaled_Naba2014_Ecmcore",
        "scaled_Triulzi2013_Ecm",

        "MCPcounter_Tcell",
        "MCPcounter_B.Lineage",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts"
        ),

      function(sig, xclin){

        print(sig)


        # Per-strata level !!!!
        # Per-subtype, per-sig, per-strata prognosis assesment
        xxx <- purrr::map(
          xclin$Strata %>% unique(),
          function(xstrata, sig, xclin){

            # Estimate prognosis for individual strata
            m1_prog <- glm(
              formula =as.formula(paste("Response ~", sig)),
              data = xclin %>% dplyr::filter(Strata == xstrata),
              family = "binomial"
            )

            m1_prog_clean <- cbind(
              summary(m1_prog)$coefficients[sig, , drop = F],
              confint.default(m1_prog)[sig, , drop = F]
            ) %>%
              as_tibble(rownames = "Module_name") %>%
              dplyr::mutate(Strata = xstrata)


            # Pseudo code/analysis/data !!!
            # Estimate heterogenity
            # Note that m1_het_clean has structure idenical to output from
            # estimate_prog_heterogenity(). This is to make make tibbles congruent
            # when heterogenity stat is added for pooled data.
            m1_het_clean <-
              tibble(
                Module_name = sig,
                Q = NA,
                Q_p = NA,
                I2 = NA,
                I2_lower_ci = NA,
                I2_upper_ci = NA
              ) %>%
              # Sample size informaion is appended to use in plot
              dplyr::mutate(
                N_response = xclin %>%
                  dplyr::filter(Strata == xstrata & Response ==1) %>%
                  nrow(),

                N_patients = xclin %>%
                  dplyr::filter(Strata == xstrata) %>%
                  nrow()
              )


            # Prognosis + Heterogenity
            m1_prog_clean %>%
              left_join(m1_het_clean, by = "Module_name")

            # Note that joining by Module_name is okay as within each unique strata
            # only one module is present.
            # In other words, the logic is per-sig - per-strata, hence joining by
            # either Module_name or Strata is feasible.

          },

          sig,
          xclin
        )

        names(xxx) <- xclin$Strata %>% unique()



        # Pooled strata level !!!!
        # Per-subtype, per-sig, pooled-strata prognosis assesment
        m1_prog <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Strata")),
          data = xclin,
          family = "binomial"
        )

        m1_prog_clean <- cbind(
          summary(m1_prog)$coefficients[sig, , drop = F],
          confint.default(m1_prog)[sig, , drop = F]
        ) %>%
          as_tibble(rownames = "Module_name") %>%
          dplyr::mutate(Strata = "All")

        # Estimate heterogenity due to strata
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig) %>%
          # Sample size informaion is appended to use in plot
          dplyr::mutate(
            N_response = xclin %>%
              dplyr::filter(Response ==1) %>%
              nrow(),

            N_patients = xclin %>%
              nrow()
          )

        # Prognosis + Heterogenity of pooled strata appended to
        # -individual strata level analysis summary
        xxx[["All"]] <-  m1_prog_clean %>%
          left_join(m1_het_clean, by = "Module_name")



        # Consolidate strata level prognosis results of each sig
        bind_rows(xxx) %>%
          dplyr::rename(
            Std_error = "Std. Error",
            Z_value = "z value",
            P = "Pr(>|z|)",
            Lower_ci = "2.5 %",
            Upper_ci = "97.5 %"
          ) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH")
          )


      },
      xclin
    )

    # names(xx) <- c("TILsig", "TILsig_imm", "TILsig_fib",
    names(xx) <-   c(
      "scaled_Denovo_TILsig",
      "scaled_Denovo_Immune",
      "scaled_Denovo_ECM",

      "scaled_Gruosso2019_Immune",
      "scaled_Hamy2016_Immune",
      "scaled_Teschendorff2007_Immune",
      "scaled_Yang2018_Immune",
      "scaled_Desmedt2008_STAT1",

      "scaled_Gruosso2019_Interferon",
      "scaled_Farmer2009_MX1",
      "scaled_Hamy2016_Interferon",
      "scaled_Nirmal2018_Interferon",

      "scaled_Gruosso2019_Cholesterol",
      "scaled_Ehmsen2019_Chol",
      "scaled_Simigdala2016_Chol",
      "scaled_Sorrentino2014_Chol",

      "scaled_Gruosso2019_Fibrosis",
      "scaled_Hamy2016_Ecm",
      "scaled_Naba2014_Ecmcore",
      "scaled_Triulzi2013_Ecm",

      "MCPcounter_Tcell",
      "MCPcounter_B.Lineage",
      "MCPcounter_Monocytic.Lineage",
      "MCPcounter_Myeloid.Dendritic",
      "MCPcounter_Endothelial",
      "MCPcounter_Fibroblasts"
      )
    xx

  },

  clin = clin_neoadj_neo
)

names(geo_prog_het$individual.sig) <- c("TN", "HR")



# pooled sig >>>>>>>>>>>

geo_prog_het$pooled.sig <- purrr::map(

  c("TN", "HR"),

  function(subtype, clin){

    print(subtype)

    xclin <- switch (
      subtype,

      ALL = clin,

      HER2 = clin %>%
        dplyr::filter(Subtype_ihc == subtype),

      "HR+HER2+" = clin %>%
        dplyr::filter(Subtype_ihc2 == subtype & Strata != "GSE42822:AAA+T+TRA"),
      # The above strata has only 9 samples with 2 events
      # If included in analysis, the above strata will throw the following warnings.
      # Warning messages:
      # 1: glm.fit: algorithm did not converge
      # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred

      clin %>%
        dplyr::filter(Subtype_ihc2 == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
    # clin <- clin_neoadj_neo

    # Per-subtype, per-sig heterogenity assesment
    xx <- purrr::map(

      c(
        "scaled_Denovo_TILsig",
        "scaled_Denovo_Immune",
        "scaled_Denovo_ECM",

        "scaled_Pooled_Immune",
        "scaled_Pooled_Interferon",
        "scaled_Pooled_Fibrosis",
        "scaled_Pooled_Cholesterol",

        "MCPcounter_Tcell",
        "MCPcounter_B.Lineage",
        "MCPcounter_Monocytic.Lineage",
        "MCPcounter_Myeloid.Dendritic",
        "MCPcounter_Endothelial",
        "MCPcounter_Fibroblasts"
      ),

      function(sig, xclin){

        print(sig)


        # Per-strata level !!!!
        # Per-subtype, per-sig, per-strata prognosis assesment
        xxx <- purrr::map(
          xclin$Strata %>% unique(),
          function(xstrata, sig, xclin){

            # Estimate prognosis for individual strata
            m1_prog <- glm(
              formula =as.formula(paste("Response ~", sig)),
              data = xclin %>% dplyr::filter(Strata == xstrata),
              family = "binomial"
            )

            m1_prog_clean <- cbind(
              summary(m1_prog)$coefficients[sig, , drop = F],
              confint.default(m1_prog)[sig, , drop = F]
            ) %>%
              as_tibble(rownames = "Module_name") %>%
              dplyr::mutate(Strata = xstrata)


            # Pseudo code/analysis/data !!!
            # Estimate heterogenity
            # Note that m1_het_clean has structure idenical to output from
            # estimate_prog_heterogenity(). This is to make make tibbles congruent
            # when heterogenity stat is added for pooled data.
            m1_het_clean <-
              tibble(
                Module_name = sig,
                Q = NA,
                Q_p = NA,
                I2 = NA,
                I2_lower_ci = NA,
                I2_upper_ci = NA
              ) %>%
              # Sample size informaion is appended to use in plot
              dplyr::mutate(
                N_response = xclin %>%
                  dplyr::filter(Strata == xstrata & Response ==1) %>%
                  nrow(),

                N_patients = xclin %>%
                  dplyr::filter(Strata == xstrata) %>%
                  nrow()
              )


            # Prognosis + Heterogeneity
            m1_prog_clean %>%
              left_join(m1_het_clean, by = "Module_name")

            # Note that joining by Module_name is okay as within each unique strata
            # only one module is present.
            # In other words, the logic is per-sig - per-strata, hence joining by
            # either Module_name or Strata is feasible.

          },

          sig,
          xclin
        )

        names(xxx) <- xclin$Strata %>% unique()



        # Pooled strata level !!!!
        # Per-subtype, per-sig, pooled-strata prognosis assessment
        m1_prog <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Strata")),
          data = xclin,
          family = "binomial"
        )

        m1_prog_clean <- cbind(
          summary(m1_prog)$coefficients[sig, , drop = F],
          confint.default(m1_prog)[sig, , drop = F]
        ) %>%
          as_tibble(rownames = "Module_name") %>%
          dplyr::mutate(Strata = "All")

        # Estimate heterogeneity due to strata
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig) %>%
          # Sample size information is appended to use in plot
          dplyr::mutate(
            N_response = xclin %>%
              dplyr::filter(Response ==1) %>%
              nrow(),

            N_patients = xclin %>%
              nrow()
          )

        # Prognosis + Heterogenity of pooled strata appended to
        # -individual strata level analysis summary
        xxx[["All"]] <-  m1_prog_clean %>%
          left_join(m1_het_clean, by = "Module_name")



        # Consolidate strata level prognosis results of each sig
        bind_rows(xxx) %>%
          dplyr::rename(
            Std_error = "Std. Error",
            Z_value = "z value",
            P = "Pr(>|z|)",
            Lower_ci = "2.5 %",
            Upper_ci = "97.5 %"
          ) %>%
          dplyr::mutate(
            P_adj = p.adjust(p = P, method = "BH")
          )


      },
      xclin
    )

    # names(xx) <- c("TILsig", "TILsig_imm", "TILsig_fib",
    names(xx) <-   c(
      "scaled_Denovo_TILsig",
      "scaled_Denovo_Immune",
      "scaled_Denovo_ECM",

      "scaled_Pooled_Immune",
      "scaled_Pooled_Interferon",
      "scaled_Pooled_Fibrosis",
      "scaled_Pooled_Cholesterol",

      "MCPcounter_Tcell",
      "MCPcounter_B.Lineage",
      "MCPcounter_Monocytic.Lineage",
      "MCPcounter_Myeloid.Dendritic",
      "MCPcounter_Endothelial",
      "MCPcounter_Fibroblasts"
    )
    xx

  },

  clin = clin_neoadj_neo
)

names(geo_prog_het$pooled.sig) <- c("TN", "HR")



#
# ==============================================================================






# # 4. Analysis summary tables
# # ==============================================================================
#
# # geo_inter
#
# summarize_geo_prog <- function(x){
#
#   # x <- geo_prog
#
#
#   x <- purrr::map(x,
#                   function(xx, prog){
#
#                     nme <- c("Module_name", "Estimate", "Lower_ci", "Upper_ci", "P", "P_adj", "Q", "Q_p", "I2")
#                     xx <- xx[ , nme] %>%
#                       dplyr::mutate(
#                         Log_OR = str_c(round(Estimate, digits = 1),
#                                        " (",
#                                        round(Lower_ci,digits=1), "-",
#                                        round(Upper_ci,digits=1),
#                                        ")"),
#                         P = round(P, digits = 3),
#                         P_adj = round(P_adj, digits = 3),
#                         Q = round(Q, digits = 3),
#                         Q_p = round(Q_p, digits = 3),
#                         I2 = round(I2, digits = 3)
#                       )
#
#                     xx %>%
#                       dplyr::select(c("Module_name",
#                                       "Log_OR", "P", "P_adj", "Q", "Q_p", "I2"))
#
#
#                   })
#
#   return(x)
#
# }
#
# xout <- summarize_geo_prog(x = geo_prog$individual.sig)
# write_xlsx(xout, path = str_c(out_tables,"geo_prognosis_individual.sig_summary_v2.2.xlsx"))
#
# xout <- summarize_geo_prog(x = geo_prog$pooled.sig)
# write_xlsx(xout, path = str_c(out_tables,"geo_prognosis_pooled.sig_summary_v2.2.xlsx"))
#
#
#
# #
# # ==============================================================================





# 5. Save Robjects
# ==============================================================================

# geo_prog : Prognosis summary
geo_prog_neo <- geo_prog
save(geo_prog_neo, file = str_c(out_data,"geo_prog_neo.RData")) # old name: prog_sum

# geo_prog_het: Prognosis heterogeneity summary
geo_prog_het_neo <- geo_prog_het
save(geo_prog_het_neo, file = str_c(out_data,"geo_prog_het_neo.RData")) # old name: prog_sum_het

#
# ==============================================================================



# Clear memory
# ==============================================================================

rm(clin_neoadj_neo, geo_prog_neo, geo_prog_het_neo, geo_prog, geo_prog_het)

#
# ==============================================================================
