# s3_exploring_prognosis.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the prognostic value of TIL associated signatures and Immune celltypes.
# Account for study heterogenity.
# Prepare data necessary to plot prognostic plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and summarize clin_neoadj_finneo.
# 2. Explore prognosis.
# 3. Explore prognosis heterogenity.
# 4. Analysis summary
# 5. Save Robjects



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
# 1. Per subtype prognosis (hr-arm-study heterogenity)
# 2. Pan-subtype prognosis (subtype-hr-arm-study heterogenity)

# Analysis note 1:
# For the main analysis consider HER2 as a whole (i.e. both HR- and HR+ HER2) with
# stratification based on hr-arm-study.
# Further as a supplementary analysis, repeat the same analysis (as in HER2 whole)
# in each HR-HER2+ and HR+HER2+ independetly.

# Analysis note 2:
# By definition, prognosis is the natural course of disease independent of
# treatment regimen. Hence limiting analysis to specific treatment regimen is
# meaning less. Further, it is wrong to interpret interaction from prognostic model
# even if prognostic effect is limited to a single treatment regimen and not to others.
# An interaction test is the minimum requirement to clim a significant treatment
# interaction.


# All neoadj regimen with Strata info
clin_neoadj_finneo %>%
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
# 1     HER2 neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20194:AAA+T:HR-  25    13    52.000000
# 2     HER2 neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE32646:AAA+T:HR-  18     9    50.000000
# 3     HER2 neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE50948:AAA+T:HR-  37     9    24.324324
# 4     HER2 neg     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE42822:AAA+T+TRA:HR-  15    10    66.666667
# 5     HER2 neg     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE50948:AAA+T+TRA:HR-  43    25    58.139535

# 6     HER2 pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE20194:AAA+T:HR+  25     4    16.000000
# 7     HER2 pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE32646:AAA+T:HR+  16     3    18.750000
# 8     HER2 pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy     GSE50948:AAA+T:HR+  14     4    28.571429
# 9     HER2 pos     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE42822:AAA+T+TRA:HR+   9     2    22.222222
# Note that the above strata with 9 samples create convergense issue !!!!
# 10    HER2 pos     AAA+Taxane///Trastuzumab///No_hormone_therapy///No_other_therapy GSE50948:AAA+T+TRA:HR+  20     6    30.000000

# 11      HR pos   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:A0A+T  35     3     8.571429
# 12      HR pos   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE41998:A0A+T 104    12    11.538462
# 13      HR pos AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE20271:AAA  49     3     6.122449
# 14      HR pos AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE22093:AAA  42    10    23.809524
# 15      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20194:AAA+T 141     9     6.382979
# 16      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20271:AAA+T  44     3     6.818182
# 17      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE23988:AAA+T  32     7    21.875000
# 18      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:AAA+T 194    22    11.340206
# 19      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE32646:AAA+T  55     5     9.090909
# 20      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE42822:AAA+T  29     8    27.586207
# 21      HR pos   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE50948:AAA+T  26     4    15.384615

# 22      TN neg   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:A0A+T  25     6    24.000000
# 23      TN neg   A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE41998:A0A+T 125    51    40.800000
# 24      TN neg AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE20271:AAA  28     4    14.285714
# 25      TN neg AAA+noTaxane///No_her2_agent///No_hormone_therapy///No_other_therapy           GSE22093:AAA  55    18    32.727273
# 26      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20194:AAA+T  64    25    39.062500
# 27      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE20271:AAA+T  31     9    29.032258
# 28      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE23988:AAA+T  29    13    44.827586
# 29      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE25066:AAA+T 119    43    36.134454
# 30      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE32646:AAA+T  26    10    38.461538
# 31      TN neg   AAA+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy         GSE42822:AAA+T  24    12    50.000000

# Note !!!!!
# The event rate in HR subtype is low possibily due to lack of hormone therapy


#
# ==============================================================================





# 2. Explore prognosis.
# ==============================================================================

# Per subtype prognosis (Sereis matrix + Arm (~ Strata) heterogenity)


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

# Basic prognosis plot strucutre ( per subtype - per sig - per arm - study heterogenity)
# (!!! Note : This plot will be generated in the NEXT SECTION of the script)
#
# TN - sig1 - arm1 - prognosis - study heterogenity
# TN - sig1 - arm2 - prognosis - study heterogenity
# TN - sig1 - arm3 - prognosis - study heterogenity
# TN - sig1 - arm4 - prognosis - study heterogenity
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# TN - sig1 - armX - prognosis - study and arm heterogenity (global)
# Note:
# In prognosis forest plot include the global prognosis (pooled arm) for sig1 per subtype.
# The global analysis account for study and arm heterogenity,
# and it also asseses the heterogenity level.


# Higher level prognosis plot strucutre
# TN/HER2/HR/ALL (event summary for each subtype, event/size: e/n)
# (!!! Note : This plot will be generated in the THIS SECTION of the script)
#
# immune sig1 - armX - prognosis - study and arm heterogenity (global)
# immune sig2 - armX - prognosis - study and arm heterogenity (global)
# immune sig3 - armX - prognosis - study and arm heterogenity (global)
# interferon sig1 - armX - prognosis - study and arm heterogenity (global)
# interferon sig2 - armX - prognosis - study and arm heterogenity (global)
# interferon sig3 - armX - prognosis - study and arm heterogenity (global)
# fibrosis sig1 - armX - prognosis - study and arm heterogenity (global)
# fibrosis sig2 - armX - prognosis - study and arm heterogenity (global)
# fibrosis sig3 - armX - prognosis - study and arm heterogenity (global)
# cholesterol sig1 - armX - prognosis - study and arm heterogenity (global)
# cholesterol sig2 - armX - prognosis - study and arm heterogenity (global)
# cholesterol sig3 - armX - prognosis - study and arm heterogenity (global)
# proliferation sig1 - armX - prognosis - study and arm heterogenity (global)
# proliferation sig2 - armX - prognosis - study and arm heterogenity (global)
# proliferation sig3 - armX - prognosis - study and arm heterogenity (global)


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
    # clin <- clin_neoadj_finneo

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


        # Estimate heterogenity
        # Note that heterogenity assessed w.r.t each strata
        # For each strata, the prognostic value is assessed and the heterogenity is measured
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig)


        # Prognosis + Heterogenity
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

  clin = clin_neoadj_finneo
)

names(geo_prog$individual.sig) <- c("TN", "HR")



# Prognosis: pooled sig

geo_prog$pooled.sig <- purrr::map(

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
    # clin <- clin_neoadj_finneo

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


        # Estimate heterogenity
        # Note that heterogenity assessed w.r.t each strata
        # For each strata, the prognostic value is assessed and the heterogenity is measured
        m1_het_clean <- estimate_prog_heterogenity(xclin = xclin, sig = sig)


        # Prognosis + Heterogenity
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

  clin = clin_neoadj_finneo
)

names(geo_prog$pooled.sig) <- c("TN", "HR")


# !!!! Save geo_prog late r in the script

#
# ==============================================================================





# 3. Explore prognosis heterogenity.
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




# Explore per arm prognostic effect and study heterogenity
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Per-subtype, per-sig strata (study-arm) heterogenity assesment


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
    # clin <- clin_neoadj_finneo

    # Per-subtype, per-sig heterogenity assesment
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
      )
    xx

  },

  clin = clin_neoadj_finneo
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
    # clin <- clin_neoadj_finneo

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

      "scaled_Pooled_Immune",
      "scaled_Pooled_Interferon",
      "scaled_Pooled_Fibrosis",
      "scaled_Pooled_Cholesterol",

      "MCPcounter_Fibroblasts"
    )
    xx

  },

  clin = clin_neoadj_finneo
)

names(geo_prog_het$pooled.sig) <- c("TN", "HR")



#
# ==============================================================================






# 4. Analysis summary tables
# ==============================================================================

# geo_inter

summarize_geo_prog <- function(x){

  # x <- geo_prog


  x <- purrr::map(x,
                  function(xx, prog){

                    nme <- c("Module_name", "Estimate", "Lower_ci", "Upper_ci", "P", "P_adj", "Q", "Q_p", "I2")
                    xx <- xx[ , nme] %>%
                      dplyr::mutate(
                        Log_OR = str_c(round(Estimate, digits = 1),
                                       " (",
                                       round(Lower_ci,digits=1), "-",
                                       round(Upper_ci,digits=1),
                                       ")"),
                        P = round(P, digits = 3),
                        P_adj = round(P_adj, digits = 3),
                        Q = round(Q, digits = 3),
                        Q_p = round(Q_p, digits = 3),
                        I2 = round(I2, digits = 3)
                      )

                    xx %>%
                      dplyr::select(c("Module_name",
                                      "Log_OR", "P", "P_adj", "Q", "Q_p", "I2"))


                  })

  return(x)

}

xout <- summarize_geo_prog(x = geo_prog$individual.sig)
write_xlsx(xout, path = str_c(out_tables,"geo_prognosis_individual.sig_summary_v2.1.xlsx"))

xout <- summarize_geo_prog(x = geo_prog$pooled.sig)
write_xlsx(xout, path = str_c(out_tables,"geo_prognosis_pooled.sig_summary_v2.1.xlsx"))



#
# ==============================================================================





# 5. Save Robjects
# ==============================================================================

# geo_prog : Prognosis summary
save(geo_prog, file = str_c(out_data,"geo_prog_v2.1.RData")) # old name: prog_sum

# geo_prog_het: Prognosis heterogenity summary
save(geo_prog_het, file = str_c(out_data,"geo_prog_het_v2.1.RData")) # old name: prog_sum_het

#
# ==============================================================================



