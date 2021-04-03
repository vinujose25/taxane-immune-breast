# functions.R


filter_data <- function(x, type = "filter_data"){

  # What the function does?
  # Filters data accrording to the criteria mentioned in the manuscript.
  # Depends on the clinical tibble (internal) variable name and values.
  # Strinctly an internal function, or not a general purpose function.
  #
  # Study/sample filtering criteria:
  # The pooled dataset is filtered using the following criteria;
  # a) Discard samples with missing IHC status for HR and HER2
  #
  # For each IHC subtype (HR,HER2,TN):
  # b) Discard samples with missing treatment data
  # c) Discard samples with missing clinical endpoint (pCR/DFS) data
  # d) Discard matched samples from longitudnal studies (series matrices)
  #   (i.e. consider only pre-treatment expression profiles),
  # e) Discard studies from a treatment regimen with less than 20 samples.
  # f) Discard treatment regimens with sample size less than 50 per subtype, and
  # g) Discard treatment regimens with samples only from a single study per subtype,


  #
  # The present function performs steps b - g.



  # Script structure:
  # For each IHC subtype (HR,HER2,TN):
  # b) Discard samples with missing treatment data
  # c) Discard samples with missing clinical endpoint (pCR/DFS) data
  # d) Discard matched samples from longitudnal studies (series matrices)
  #   (i.e. consider only pre-treatment expression profiles),
  # e) Discard studies from a treatment regimen with less than 20 samples per subtype.
  # f) Discard treatment regimens with sample size less than 50 per subtype, and
  # g) Discard treatment regimens with samples only from a single study per subtype,

# Note:
# Less than 20 samples will crreate convergense issue and
# large variance in estimates which are problematic in prognosis, interaction
# and heterogenisty assesments


  # x: clinical data tibble
  # type: desription of x, used for printing informative messages.

  # x = clin %>%
  #   dplyr::filter(Regimen_updated == "neoadj" & Subtype_ihc == "HR")
  # type = "neoadj"

  cat(
    type,
    ": Started with",
    x %>% nrow(),
    "no.of samples.\n"
  )



  # b) Discard samples with missing treatment data
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  cat(
    type,
    ": Discarding",
    x %>%
      dplyr::filter(is.na(Arm_consolidated)) %>%
      nrow(),
    "samples with missing treatment data.\n"
  )

  x <- x %>%
    dplyr::filter(! is.na(Arm_consolidated))



  # c) Discard samples with missing clinical endpoint (pCR/DFS) data
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if( all(x$Regimen_updated == "neoadj" )){

    cat(
      type,
      ": Discarding",
      x %>%
        dplyr::filter(is.na(Response)) %>%
        nrow(),
      "samples with missing clinical response data.\n"
    )

    x <- x %>%
      dplyr::filter(! is.na(Response))

  } else if( all(x$Regimen_updated == "adj" )){

    cat(
      type,
      ": Discarding",
      x %>%
        dplyr::filter(is.na(Event_dfs)) %>%
        nrow(),
      "samples with missing clinical follow-up data.\n"
    )

    x <- x %>%
      dplyr::filter(! is.na(Event_dfs))
  }





  # d) Discard matched samples from longitudnal studies (series matrices)
  #   (i.e. consider only pre-treatment expression profiles),
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  cat(
    type,
    ": Discarding",
    x %>%
      dplyr::filter(Timepoint == "Mid_treatment" | Timepoint == "Post_treatment") %>%
      nrow(),
    "matched samples from longitudnal studies.\n"
  )

  x <- x %>%
    dplyr::filter(is.na(Timepoint) | Timepoint == "Pre_treatment")
  # Timepoint == NA implies pre-treatment sample by default.
  # Timepoint is set to pre/mid/post-treatment only for longitudnal studies.




  # e) Discard studies from a treatment regimen with less than 20 samples.
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  xsum <- x %>%
    dplyr::group_by(Arm_consolidated, Series_matrix_accession) %>%
    dplyr::summarise(n_samp = n(), .groups = "keep") %>%
    # dplyr::filter(n_samp < 10) # old cutoff
  dplyr::filter(n_samp < 20) # recommanded
    # dplyr::filter(n_samp < 15) # test

  cat(
    type,
    ": Discarding",
    xsum$n_samp %>% sum(),
    "samples from",
    xsum %>% nrow(),
    "series matrices with <20 samples per treatment regimen.\n" # recommanded
    # "series matrices with <15 samples per treatment regimen.\n" # test
  )

  x <- x %>%
    dplyr::mutate(
      # index is introduced to aid filtering
      index = str_c(Arm_consolidated, Series_matrix_accession)
    ) %>%
    dplyr::filter(
      ! (index %in% all_of(str_c(xsum$Arm_consolidated, xsum$Series_matrix_accession)))
    ) %>%
    dplyr::select(-index)



  # f) Discard treatment regimens with sample size less than 50 per subtype, and
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  xsum <- x %>%
    dplyr::group_by(Arm_consolidated ) %>%
    dplyr::summarise(n_dataset = unique(Series_matrix_accession) %>% length(),
                     n_samp = n(), .groups = "keep") %>%
    dplyr::filter(n_samp < 50)

  cat(
    type,
    ": Discarding",
    xsum$n_samp %>% sum(),
    "samples from",
    xsum %>% nrow(),
    "treatment regimens with <50 samples.\n"
  )

  x <- x %>%
    dplyr::filter( ! (Arm_consolidated %in% all_of(xsum$Arm_consolidated)) )




  # g) Discard treatment regimens with samples only from a single study per subtype,
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  xsum <- x %>%
    dplyr::group_by(Arm_consolidated ) %>%
    dplyr::summarise(n_dataset = unique(Series_matrix_accession) %>% length(),
                     n_samp = n(), .groups = "keep") %>%
    dplyr::filter(n_dataset == 1)

  cat(
    type,
    ": Discarding",
    xsum$n_samp %>% sum(),
    "samples from",
    xsum %>% nrow(),
    "treatment regimens unique to a single study.\n"
  )

  x <- x %>%
    dplyr::filter( ! (Arm_consolidated %in% all_of(xsum$Arm_consolidated)) )



  # Output
  # >>>>>>
  cat(
    type,
    ": End with",
    x %>% nrow(),
    "no.of samples.\n"
  )

  return(x)

}



# Note: Copied from treatment-curated-breast-dataset/src/r/functions.R
# updated t_tibble()
t_tibble <- function(x, names_x_desc = NULL){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  # Tranposes a gene expression tibble.
  # By definition, row names are not recommanded for tibble, instead any
  # rowname must be included as a column in the tibble. For instance the
  # as_tibble() function has "rownames" option to give names for existing rownames
  # so as to include the rownames as the first colum of a tibble. However,
  # rownames are conventional and usually included in a tibble as the first column.
  # In this case, the default transpose function may not work properly as
  # transposing a tibble may cause data type conflicts and may convert the
  # entrire colum of a transposed tibble as the with the type of first column of
  # non-transposed tibble. This will be case with gene expression tibbles with
  # first column as gene symbols, and column names are sample names.
  # Note that this function is defined to use with gene expression tibbles.
  # Gene expression tibbles have the identical data type for all columns except
  # the first column representing rownames.

  # input
  # >>>>>
  # x: A tibble to transpose. The 1st column is considered as rownames.
  # names_x_desc: A character string that describes what names(x) represents.
  #               This character string will be used as the name of the first
  #               column of the trassposed tibble.

  # output
  # >>>>>>
  # A transposed tibble


  # Note: t(x) will set all data to single type

  # If names_x_desc is NULL, it is set as the name of the 1st column of x
  if (is.null(names_x_desc)) {
    names_x_desc <- names(x)[1]
  }

  # Extract rownames of x. This will later used to set column names of t(x)
  rownames_x <- x %>% dplyr::select(1) %>% tibble::deframe()


  tx <- t(x %>% dplyr::select(-1))
  colnames(tx) <- rownames_x

  return(as_tibble(tx, rownames = names_x_desc))
}


# Note: Copied from treatment-curated-breast-dataset/src/r/functions.R
get_module_score <- function(x, module_list, by = "Entrez_Gene"){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  # Compute gene module score as a weighted average

  # input
  # >>>>>
  # x: A tibble with genes on column and samples on rows.
  #     1st column is considered as sample names.
  # module_list: List of gene modules. Each gene module must contain at least
  #               two columns; a gene id column with column name exactly as the
  #               value of the "by" input and a direction column with
  #               column name as "Direction". Further, the gene id
  #               in the gene module must match the gene ids of x, ie. column
  #               names of x. The Direction in the gene module represents the
  #               gene's association with the phenotype that the
  #               gene module represents. For instance, for a gene module
  #               representing TP53 mutation srtatus a direction of +1 suggests
  #               up regulated in mutant phenotype and -1 represents
  #               down regulated in mutant phenotype.
  # by: The column name of gene module which contain gene id. The type gene id
  #     in the gene module must match that of column names of x.

  # output
  # >>>>>>
  # A tibble of module score. Column names represnts module names and the 1st
  # column contains sample ids.


  sig_score <- purrr::map_dfr( # map_dfr
    module_list,
    function(sig, xdata, by){

      sig <- sig %>%
        dplyr::filter_at(
          by,
          dplyr::all_vars(. %in% names(xdata))
        )


      xdata <- xdata %>%
        dplyr::select(
          c(names(xdata)[1],
            sig %>% dplyr::select(by) %>% tibble::deframe())
        )


      score <- (xdata %>% dplyr::select(-1) %>% as.matrix()) %*% (sig %>% dplyr::select(Direction) %>% as.matrix())

      # score/nrow(sig) # this old code generates matrix
      as.numeric(score/nrow(sig))
    },
    xdata = x,
    by = by
  )

  sig_score <- dplyr::bind_cols(x %>% dplyr::select(1),
                                sig_score)

}



# Validation of module subset in independent dataset
validate_gene_modules <- function(module_full_list, module_subset_list, validation_data){

  # module_full_list: list of full modules; each module dataframe must contain
  #                   minimum two columns "Ncbi_gene_id" and "Direction"
  # module_subset_list: list of module subsets; each module dataframe must contain
  #                     minimum two columns "Ncbi_gene_id" and "Direction"
  # validation_data: dataset for validation; gene on columns and samples on rows

  # # Test data
  # module_full_list = module_list
  # module_subset_list = module_list_subset
  # validation_data = tcga # 1073 x 19406

  score_full <- get_module_score(x = validation_data,
                                 module_list = module_full_list,
                                 by = "Ncbi_gene_id")
  score_subset <- get_module_score(x = validation_data,
                                   module_list = module_subset_list,
                                   by = "Ncbi_gene_id")
  pearson <- cor(score_full[,-1], score_subset[,-1], method = "pearson") %>% diag()
  spearman <- cor(score_full[,-1], score_subset[,-1], method = "spearman") %>% diag()

  tibble(
    Module_id = names(module_full_list),
    Module_full_size = purrr::map_int(module_full_list, nrow),
    Module_subset_size = purrr::map_int(module_subset_list, nrow),
    Module_subset_size_percent = (Module_subset_size/Module_full_size) * 100,
    Pearson = pearson,
    Spearman = spearman
  )
}



filter_no_event_strata <- function(x, strata_variable, event_variable){

  # strata_variable = "Arm"
  # event_variable = "DDFS_Event"
  # x = clin_finher

  # summarize x
  xsum <- x %>%
    dplyr::group_by(Strata = get(strata_variable)) %>%
    dplyr::summarise(Event = (get(event_variable) == 1) %>% sum(),
                     N = n(), .groups = "keep") %>%
    dplyr::filter(Event == 0)


  if(nrow(xsum) == 0) {

    str_c(
      "No arms with zero events. No arms were filtered.") %>%
      print()

    return(x)

  } else {

    str_c(
      "Filtering arms with zero events. Samples per arm filtered: ",
      str_c(xsum$N, " (", xsum$Strata, ")", collapse = ", ")
    ) %>%
      print()

    return(x %>% dplyr::filter( !(get(strata_variable) %in% xsum$Strata)) )

  }

}


# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_prog <- function(
  event_variable,
  time_variable,
  biomarker,
  strata_variable,
  xdata
){

  # # Test data
  # event_variable = "DDFS_Event"
  # time_variable = "DDFS_Time_Years"
  # biomarker = "Immune1"
  # strata_variable = "Arm"
  # # xdata = clin_finher
  #
  # subtype = "HR-HER2+"
  # print(subtype)
  # subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
  #
  # # xdata <- filter_no_event_strata(
  # #   x = subset(x = clin,
  # #              subset = (clin[, subtype_varible] == subtype)),
  # #   strata_variable = "Arm",
  # #   event_variable = str_c(event_type, "_Event"))
  #
  # xdata = subset(x = clin_finher,
  #                subset = (clin_finher[, subtype_varible] == subtype))

  m1 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "+ strata(", strata_variable, ")")
    ),
    data = xdata
    )

  # testing heterogenity w.r.t strata
  m1_het0 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "+", strata_variable)
    ),
    data = xdata
    )

  m1_het1 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "*", strata_variable)
    ),
    data = xdata
    )

  # summary(m1)$coefficients
  # summary(m1_het0)$coefficients
  # summary(m1_het1)$coefficients
  # lrtest(m1_het0, m1_het1)[, "Pr(>Chisq)"] %>% na.omit()

  # Summary dataframe
  x <- cbind(
    summary(m1)$coefficients,
    summary(m1)$conf.int[, c("lower .95", "upper .95"), drop = F],
    lrtest(m1_het0, m1_het1)[2, "Pr(>Chisq)", drop = F]
  )

  x <- data.frame(
    x,
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    stringsAsFactors = F,
    check.names = F
  )

  x <- x[, c("exp(coef)", "Pr(>|z|)", "lower .95", "upper .95", "Pr(>Chisq)", "Variable", "Endpoint", "Event", "N")]
  names(x) <- c("HR", "P", "Low95", "Up95", "Heterogeneity", "Variable", "Endpoint", "Event", "N")

  x
}



# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_inter <- function(
  event_variable,
  time_variable,
  biomarker,
  interaction_variable,
  strata_variable,
  xdata
){

  # event_variable: variable (column) name present in xdata. e.g. "DDFS_Event"
  # time_variable: variable (column) name present in xdata. e.g. "DDFS_Time_Years"
  # biomarker: variable (column) name present in xdata. e.g. "Immune1"
  # interaction_variable: variable (column) name present in xdata. e.g. "Chemo"
  # strata_variable: variable (column) name present in xdata. e.g. "Stratum_hr_tra"
  # xdata: a tibble/data frame


  # Test data
  # event_variable = "DDFS_Event"
  # time_variable = "DDFS_Time_Years"
  # biomarker = "Immune1"
  # interaction_variable = "Chemo"
  # strata_variable = "Stratum_hr_tra"
  # xdata = subset(x = clin_finher,
  #                subset = (clin_finher[, "Subtype_IHC"] == "HR+HER2+"))


  # The original idea is to filtering out entire starat with no events, to eliminate warnings !!
  # e.g. FEC+DTX+Hormone strata in HR+HER2+ subtype
  # However the coxph documntation stated that coeffcient of only the no-event strata
  # will get affected (estimates are not usable), the estimates of the
  # remaining strata with some events are all reliable.
  # Hence, all strata (all data) is used.
  # Ref: https://rdrr.io/cran/survival/man/coxph.html (See convergense section)

  # Filtering no-event-strata
  # xdata <- filter_no_event_strata(
  #   x = subset(x = clin_finher,
  #              subset = (clin_finher[, "Subtype_IHC_2"] == "HER2")),
  #   strata_variable = "Arm",
  #   event_variable = "DDFS_Event")


  # interaction null model
  m0 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "+", interaction_variable, "+ strata(", strata_variable,")")
    ),
    data = xdata,
    model = TRUE, # used to identify the reference/non-reference(factor levels) regimes
    x = TRUE # used to identify the reference/non-reference(factor levels) regimes
  )

  # interaction full model
  m1 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "*", interaction_variable, "+ strata(", strata_variable, ")")
    ),
    data = xdata,
    model = TRUE, # used to identify the reference/non-reference(factor levels) regimes
    x = TRUE # used to identify the reference/non-reference(factor levels) regimes
  )


  # summary(m1)$coefficients
  # summary(m1_het0)$coefficients
  # summary(m1_het1)$coefficients
  # lrtest(m1_het0, m1_het1)[, "Pr(>Chisq)"] %>% na.omit()

  x <- data.frame(
    summarise_interaction_cox(m1, m0, test_var = biomarker,
                              inter_var = interaction_variable) %>%
      dplyr::mutate(Therapy = Module_name),
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    stringsAsFactors = F,
    check.names = F
  )

  x <- x[, c("Exp_coef", "P", "L95", "U95", "Therapy", "Variable", "Endpoint", "Event", "N")]
  names(x) <- c("HR", "P", "Low95", "Up95", "Therapy", "Variable", "Endpoint", "Event", "N")

  x
}


# clean_names <- function(x){
#
#   # Ref:http://edrub.in/CheatSheets/cheatSheetStringr.pdf
#
#   colnames(x) <- colnames(x) %>%
#     stringr::str_replace_all(pattern = "([:blank:]|[:punct:]|[<>|])+", replacement = "_") %>%
#     # [<>|] is included as [:punct:] did not recognize > & |
#     stringr::str_replace_all(pattern = "_$", replacement = "")
#   # trimming tail "_"
#
#   x
# }

estimate_prog_heterogenity <- function(xclin, sig, perm = 100){

  # xclin <- clin_neoadj %>%
  #   dplyr::filter(Subtype_ihc == "HR" &
  #                   Strata != "GSE22226_GPL4133:A0A+T" &
  #                   Strata != "GSE21974:A0A+T" &
  #                   Arm_consolidated == "A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy")
  # sig = "Cholesterol3"
  # perm = 100

  xeffect <-  purrr::map(
    xclin$Strata %>% unique(),

    # Strata =  Series_matrix_accession * Arm_consolidated
    function(study, xclin){

      m1_het <- glm(
        formula =as.formula(paste("Response ~", sig)),
        data = xclin %>%
          dplyr::filter(Strata == study),
        family = "binomial"
      )

      m1_het_sum <- summary(m1_het)

      list(
        Estimate = m1_het_sum$coefficients[sig, "Estimate"],
        Variance = m1_het_sum$cov.unscaled[sig,sig]
      )

    },
    xclin
  )

  xeffect <- bind_cols(
    tibble(Strata = xclin$Strata %>% unique()),
    bind_rows(xeffect)
  )


  xhet <- altmeta::metahet(y = xeffect$Estimate, s2 = xeffect$Variance, n.resam = perm)

  tibble(
    Module_name = sig,
    Q = xhet$Q,
    Q_p = xhet$p.Q,
    I2 = xhet$I2,
    I2_lower_ci = xhet$ci.I2[1],
    I2_upper_ci = xhet$ci.I2[2]
  )

}


summarise_interaction <- function(m1, m0, test_var, inter_var){

  # Funtion to extract and summarize interaction effect ready to plot

  # m1: full model with interaction term
  # m0: null model without interaction term
  # test_var: variable tested for interaction effect (signature)
  # inter_var: variable interacting with test var
  # Note that this script expects the "test_var" to be a continuous variable,
  # and the "inter_var" to be a factor or a character vector that can be
  # coerased to be a factor.

  # m0 = m0_inter
  # m1 = m1_inter
  # test_var = sig
  # inter_var = "Arm_consolidated"


  # interaction effect
  m1_coef <- summary(m1)$coefficients
  m1_ci <- confint.default(m1)
  # subsetting relevent terms
  # (test_var(sig) and interaction terms with test_var)
  idx <- str_detect(rownames(m1_ci), test_var) %>% which()

  if (identical(rownames(m1_ci)[idx], rownames(m1_coef)[idx])) {

    x <- cbind(m1_coef[idx, ], m1_ci[idx, ])

    # Update wald p with loglikelihood p
    p <- lrtest(m0, m1)$"Pr(>Chisq)" %>% na.omit() %>% as.numeric()

    # Obsolete code logic
    # Both pvalues in x is set to p, to aid in plotting
    # x[, "Pr(>|z|)"] <- p #

    # New code logic
    # The pvalues can be replicated at any point later in the script,
    # but deletion is difficult.
    x[, "Pr(>|z|)"] <- NA # Resetting original pvalues
    x[1, "Pr(>|z|)"] <- p # Setting likihood ratio test interaction pvalue


    # Obsolete logic !!!!!!!!!!!!
    # nrow(x): the last row will be of interaction term

    # New logic !!!!!!!!!!!!!!!!!
    # The first row is the effect from reference arm.
    # All other rows represents change with respect to the reference row
    # All other rows required to adjust it effect for independet comparison in plots
    # Adjustment 1: add reference effect to each other row's effect
    # Adjustment 2: adjust other row's standard error
    # Adjustment 3: adjust other row's 95% ci


    # Updating estimate, std.err and ci for all relevant rows except the 1st row (reference row)

    idx2 <- (1:nrow(x))[-1]

    for(i in idx2){

      # Update Estimate(Log-OR)
      # Interaction term contains change from reference
      # To get actual Estiate(Log-OR), add interaction Estiate(Log-OR)
      # to reference Estiate(Log-OR)
      x[i, "Estimate"] <- sum(x[c(1, i), "Estimate"])

      # Note in the following code descriptions, coef ~ Estimate(Log-OR).
      # To get CI of sum of coefs find Std.Error of sum of coefs
      # Std.Error of sum of coefs: sqrt( sum(coef variances) + 2 * their covariance)
      # Sum of variances reference from wikipedia, search term: "Variance"
      se <- sqrt(
        sum(x[c(1, i), "Std. Error"] ^ 2) +
          (2 * vcov(m1)[rownames(x)[1], rownames(x)[i]])
      )
      x[i, "Std. Error"] <- se
      x[i, "2.5 %"] <- x[i, "Estimate"] - (1.96 * se)
      x[i, "97.5 %"] <- x[i, "Estimate"] + (1.96 * se)

    }

  } else {
    # As an error catching logic
    x <- rep(NA, 6)
    names(x) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "2.5 %", "97.5 %")
    #vector names used in following code
  }

  # Appending reference arm detials to rownames(x)[1]
  if (class(m1$data[,inter_var] %>% tibble::deframe()) == "factor") {
    ref_arm <- (m1$data[,inter_var] %>% tibble::deframe() %>% levels())[1]
  } else {
    ref_arm <- (m1$data[,inter_var] %>% tibble::deframe() %>% factor() %>% levels())[1]
  }
  ref_arm <- str_c(inter_var, ref_arm)

  if (str_detect(rownames(x)[1], ":")) {

    # Making sure the row is reference row.
    # For reference row ":" is not present.
    stop("Error in interaction summary preparation.")

  } else {

    rownames(x)[1] <- str_c(rownames(x)[1], ":", ref_arm)
  }



  # Adding prognostic effect
  m0_coef <- summary(m0)$coefficients
  m0_ci <- confint.default(m0)

  x <- rbind(x, cbind(m0_coef[test_var, , drop=F], m0_ci[test_var, , drop = F]))
  # Note the irrelevant pvalue of prognostic effect is also added to x.
  # Displaying this prognostic pvalue in interaction plot add confusion.
  # Ignore the prognostic pvalue.


  x %>%
    as_tibble(rownames = "Module_name") %>%
    dplyr::rename(
      Std_error = "Std. Error",
      Z_value = "z value",
      P = "Pr(>|z|)",
      l95 = "2.5 %",
      u95 = "97.5 %"
    )

} # end of fun


summarise_interaction_cox <- function(m1, m0, test_var, inter_var){

  # Funtion to extract and summarize interaction effect ready to plot

  # m1: full model with interaction term
  # m0: null model without interaction term
  # test_var: variable tested for interaction effect (signature)
  # inter_var: variable interacting with test var
  # Note that this script expects the "test_var" to be a continuous variable,
  # and the "inter_var" to be a factor or a character vector that can be
  # coerased to be a factor.

  # Test data
  # m0 = m0_inter
  # m1 = m1_inter
  # test_var = sig
  # inter_var = "Arm_consolidated"

  # Test data
  # m0 = m0
  # m1 = m1
  # test_var = "Immune1"
  # inter_var = "Chemo"


  # interaction effect summary
  # >>>>>>>>>>>>>>>>>>>>>>>>>>
  m1_coef <- summary(m1)$coefficients
  m1_ci <- confint.default(m1)



  # subsetting relevent terms
  # >>>>>>>>>>>>>>>>>>>>>>>>>
  # (i.e. test_var (sig) and terms interacting with test_var)

  # Obsolete, error prone code
  # idx <- str_detect(rownames(m1_ci), test_var) %>% which()
  # Convert ineger index based code to row/col name based one.
  # Indices are error prone !!!!!!!!!!!!!!!!

  idx <- rownames(m1_ci)[str_detect(rownames(m1_ci), test_var)]
  # Note
  # idx[1] : reference arm
  # The remaining terms are interaction terms



  # Updating estimate, std.err and ci for all relevant rows except the 1st row
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # if (identical(rownames(m1_ci)[idx], rownames(m1_coef)[idx])) {
  if (identical(rownames(m1_ci), rownames(m1_coef))) {

    x <- cbind(m1_coef[idx, ], m1_ci[idx, ])

    # Update wald p with loglikelihood p
    p <- lrtest(m0, m1)$"Pr(>Chisq)" %>% na.omit() %>% as.numeric()

    # Obsolete code logic
    # Both pvalues in x is set to p, to aid in plotting
    # x[, "Pr(>|z|)"] <- p #

    # New code logic
    # The pvalues can be replicated at any point later in the script,
    # but deletion is difficult.
    x[, "Pr(>|z|)"] <- NA # Resetting original pvalues
    x[idx[1], "Pr(>|z|)"] <- p # Setting likihood ratio test interaction pvalue


    # Obsolete logic !!!!!!!!!!!!
    # nrow(x): the last row will be of interaction term

    # New logic !!!!!!!!!!!!!!!!!
    # The first row is the effect from reference arm.
    # All other rows represents change with respect to the reference row
    # All other rows required to adjust it effect for independet comparison in plots
    # Adjustment 1: add reference effect to each other row's effect
    # Adjustment 2: adjust other row's standard error
    # Adjustment 3: adjust other row's 95% ci


    # Updating estimate, std.err and ci for all relevant rows except the 1st row (reference row)


    for(i in idx[-1]){

      # idx[-1] suggest all terms excluding reference term/arm

      # Update coef(Log-HR)
      # Interaction term contains change from reference
      # To get actual coef(Log-HR), add interaction coef(Log-HR)
      # to reference Estiate(Log-OR)
      x[i, "coef"] <- sum(x[c(idx[1], i), "coef"]) # idx[1] represents reference arm

      # Note in the following code descriptions, coef ~ Estimate(Log-OR).
      # To get CI of sum of coefs find Std.Error of sum of coefs
      # Std.Error of sum of coefs: sqrt( sum(coef variances) + 2 * their covariance)
      # Sum of variances reference from wikipedia, search term: "Variance"
      se <- sqrt(
        sum(x[c(idx[1], i), "se(coef)"] ^ 2) +
          (2 * vcov(m1)[idx[1], i])
      )
      x[i, "se(coef)"] <- se
      x[i, "2.5 %"] <- x[i, "coef"] - (1.96 * se)
      x[i, "97.5 %"] <- x[i, "coef"] + (1.96 * se)

    }

  } else {
    # As an error catching logic
    x <- rep(NA, 6)
    names(x) = c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "2.5 %", "97.5 %")
    #vector names used in following code
  }


  # Appending reference arm detials to rownames(x)[1]
  non_ref_arm <- grep(x = colnames(m1$x), pattern = ":", value = T)
  non_ref_arm <- str_split_fixed(string = non_ref_arm, pattern = ":", n = 2)[,2]


  ref_arm <- setdiff(
    unique(m1$model[, inter_var]),
    str_replace(string = non_ref_arm, pattern = inter_var, replacement = "")
  )
  ref_arm <- str_c(inter_var, ref_arm)


  # if (class(m1$data[,inter_var] %>% tibble::deframe()) == "factor") {
  #   ref_arm <- (m1$data[,inter_var] %>% tibble::deframe() %>% levels())[1]
  # } else {
  #   ref_arm <- (m1$data[,inter_var] %>% tibble::deframe() %>% factor() %>% levels())[1]
  # }
  # ref_arm <- str_c(inter_var, ref_arm)

  if (str_detect(rownames(x)[1], ":")) {

    # Making sure the 1st row is reference row.
    # For reference row ":" is not present.
    stop("Error in interaction summary preparation.")

  } else {

    rownames(x)[1] <- str_c(rownames(x)[1], ":", ref_arm)
  }


  # Discard adding prognostic effect in interaction test

  # # Adding prognostic effect
  # m0_coef <- summary(m0)$coefficients
  # m0_ci <- confint.default(m0)
  #
  # x <- rbind(x, cbind(m0_coef[test_var, , drop=F], m0_ci[test_var, , drop = F]))
  # # Note the irrelevant pvalue of prognostic effect is also added to x.
  # # Displaying this prognostic pvalue in interaction plot add confusion.
  # # Ignore the prognostic pvalue.


  x %>%
    as_tibble(rownames = "Module_name") %>%
    dplyr::rename(
      Exp_coef = "exp(coef)",
      Std_error = "se(coef)",
      Z_value = "z",
      P = "Pr(>|z|)",
      l95 = "2.5 %",
      u95 = "97.5 %"
    ) %>%
    dplyr::rename_all(~str_to_title(.x))

} # end of fun
