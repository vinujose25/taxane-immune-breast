# s3_exploring_prognosis.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Explore the prognostic value of TIL associated signatures and Immune celltypes.
# Account for study heterogenity.
# Prepare data necessary to plot prognostic plots.



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and summarize clin_neoadj.
# 2. Explore prognosis.
# 3. Explore prognosis heterogenity.
# 4. Save Robjects



# 1. Load and format clin_neoadj.
# ==============================================================================

load("results/data/clin_neoadj.RData")


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


# Analysis
# 1. Per subtype prognosis (arm-study heterogenity)
# 2. Pan-subtype prognosis (subtype-arm-study heterogenity)


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



# Get prognostic effect and heterogenity measurement
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Prognosis summary
prog_sum <- purrr::map(

  c("TN", "HER2", "HR", "ALL"),

  function(subtype, clin){

    print(subtype)

    xclin <- switch (
      subtype,
      ALL = clin,
      clin %>%
        dplyr::filter(Subtype_ihc == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
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

        # Estimate prognosis
        m1_prog <- glm(
          formula =as.formula(paste("Response ~", sig, "+ Strata")),
          data = xclin,
          family = "binomial"
        )

        m1_prog_clean <- cbind(
          summary(m1_prog)$coefficients[sig, , drop = F],
          confint.default(m1_prog)[sig, , drop = F]
        ) %>%
          as_tibble(rownames = "Module_name")


        # Estimate heterogenity
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

  clin = clin_neoadj
)

# There were 24 warnings (use warnings() to see them)
# Warning messages:
# 1: glm.fit: algorithm did not converge
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
# 3: glm.fit: algorithm did not converge

names(prog_sum) <- c("TN", "HER2", "HR", "ALL")


# !!!! Save prog_sum later in the script

#
# ==============================================================================





# 3. Explore prognosis heterogenity.
# ==============================================================================


# Explore per arm prognostic effect and study heterogenity
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Per-subtype, per-sig study-arm heterogenity assesment

prog_het_sum <- purrr::map(

  c("TN", "HER2", "HR", "ALL"),

  function(subtype, clin){

    print(subtype)

    xclin <- switch(
      subtype,
      ALL = clin,
      clin %>%
        dplyr::filter(Subtype_ihc == subtype)
    )

    # sig = "Immune1"
    # subtype = "TN"
    # clin <- clin_neoadj

    # Per-subtype, per-sig heterogenity assesment
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

    names(xx) <- c("Immune1", "Immune2", "Immune3",
                   "Interferon1", "Interferon2", "Interferon3",
                   "Cholesterol1", "Cholesterol2", "Cholesterol3",
                   "Fibrosis1", "Fibrosis2", "Fibrosis3",
                   "Proliferation1", "Proliferation2", "Proliferation3",
                   "Tcell", "CLymphocyte", "Bcell",
                   "NKcell", "Monocyte", "MDendritic", "Fibroblast")
    xx

  },

  clin = clin_neoadj
)

# There were 48 warnings (use warnings() to see them)
# Warning messages:
# 1: glm.fit: algorithm did not converge
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
# 3: glm.fit: algorithm did not converge

names(prog_het_sum) <- c("TN", "HER2", "HR", "ALL")

#
# ==============================================================================




# 4. Save Robjects
# ==============================================================================

# prog_sum : Prognosis summary
save(prog_sum, file = str_c(out_data,"prog_sum.RData"))

# prog_het_sum: Prognosis heterogenity summary
save(prog_het_sum, file = str_c(out_data,"prog_het_sum.RData"))

#
# ==============================================================================



