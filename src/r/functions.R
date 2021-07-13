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

# Note: "averaged the expression of positively and negatively associated genes
# independently and computed their differences"
get_module_score_2 <- function(x, module_list, by = "Entrez_Gene"){

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

      sig_up <- sig %>%
        dplyr::filter(Direction == 1)

      sig_dn <- sig %>%
        dplyr::filter(Direction == -1)

      if(nrow(sig_up) >= 1){
        xdata_up <- xdata %>%
          dplyr::select(
            c(names(xdata)[1],
              sig_up %>% dplyr::select(by) %>% tibble::deframe())
          )
        score_up <- (xdata_up %>% dplyr::select(-1) %>% as.matrix()) %*% (sig_up %>% dplyr::select(Direction) %>% as.matrix())
        # score/nrow(sig) # this old code generates matrix
        score_up <- as.numeric(score_up/nrow(sig_up))
      } else {
        score_up <- 0
      }

      if(nrow(sig_dn) >= 1){
      xdata_dn <- xdata %>%
        dplyr::select(
          c(names(xdata)[1],
            sig_dn %>% dplyr::select(by) %>% tibble::deframe())
        )
      score_dn <- (xdata_dn %>% dplyr::select(-1) %>% as.matrix()) %*% (sig_dn %>% dplyr::select(Direction) %>% as.matrix())
      # score/nrow(sig) # this old code generates matrix
      score_dn <- as.numeric(score_dn/nrow(sig_dn))
      } else {
        score_dn <- 0
      }


      score_up + score_dn

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

  score_full <- get_module_score_2(x = validation_data,
                                 module_list = module_full_list,
                                 by = "Ncbi_gene_id")
  score_subset <- get_module_score_2(x = validation_data,
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


tertile <- function(x){
  q1 <- quantile(x, probs = 0.33, na.rm=T)
  q2 <- quantile(x, probs = 0.66, na.rm=T)
  dplyr::if_else( x <= q1, "low",
                  if_else(x <= q2, "mid", "high")) %>%
    factor(levels = c("low", "mid", "high"))
}

binirize <- function(x){
  q1 <- quantile(x, probs = 0.25, na.rm=T)
  dplyr::if_else( x <= q1, "low", "high") %>%
    factor(levels = c("low", "high"))
}

binirize_q2 <- function(x){
  q2 <- quantile(x, probs = 0.5, na.rm=T)
  dplyr::if_else( x <= q2, "low", "high") %>%
    factor(levels = c("low", "high"))
}

# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_prog <- function(
  event_variable,
  time_variable,
  biomarker,
  xdata
){

  # # Test data
  # event_variable = "DDFS_Event"
  # time_variable = "DDFS_Time_Years"
  # biomarker = "Immune1"
  # xdata = clin_finher %>%
  #   dplyr::mutate( Proliferation1 = binirize(Proliferation1) %>%
  #                    factor(levels = c("low", "high"))) %>%
  #   dplyr::filter(Subtype_IHC == "TN")


  m1 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "+ strata(Hormone, Herceptin, Chemo)")
    ),
    data = xdata
  )


 # For testing heterogenity(interaction) w.r.t strata
 m1_het0 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker, "+ interaction(Hormone, Herceptin, Chemo)")
    ),
    data = xdata
    )

  m1_het1 <- coxph( # m1_het1
    formula = as.formula(
      paste("Surv(", time_variable, ",", event_variable, ") ~",
            biomarker,"* interaction(Hormone, Herceptin, Chemo)")
    ),
    data = xdata
    )


  # Summary dataframe
  # if( class( xdata[,biomarker] %>% deframe() ) == "factor" ){
  #   # To accomodate factor Proliferation
  #   #
  #   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #   # This if condition code can be removed after careful checking according to new analysis
  #   rowname = str_c(biomarker, (xdata[,biomarker] %>% deframe() %>% levels())[2])
  #   x <- cbind(
  #     summary(m1)$coefficients[rowname, , drop=F],
  #     summary(m1)$conf.int[rowname, c("lower .95", "upper .95"), drop = F],
  #     lrtest(m1_het0, m1_het1)[2, "Pr(>Chisq)", drop = F]
  #   )
  # } else {
    x <- cbind(
      summary(m1)$coefficients[biomarker, , drop=F],
      summary(m1)$conf.int[biomarker, c("lower .95", "upper .95"), drop = F],
      lrtest(m1_het0, m1_het1)[2, "Pr(>Chisq)", drop = F]
    )
  # }



  ph_test <- cox.zph(m1)$table %>%
    tibble::as_tibble(rownames = "Test") %>%
    dplyr::mutate(Model = biomarker)

  x <- data.frame(
    x,
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    Ph_ok = if_else(any(ph_test$p < 0.05), F, T),
    stringsAsFactors = F,
    check.names = F
  )


  x <- x[, c("coef", "exp(coef)", "Pr(>|z|)", "lower .95", "upper .95", "Pr(>Chisq)", "Variable", "Endpoint", "Event", "N", "Ph_ok")]
  names(x) <- c("Coef", "HR", "P", "Low95", "Up95", "Heterogeneity", "Variable", "Endpoint", "Event", "N", "Ph_ok")

  list(
    main_effect =x,
    ph_test = ph_test
  )

}



# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_inter <- function(
  event_variable, # to extract event summary
  biomarker, # to extract model estimates
  interaction_variable, # to extract model estimates
  xdata,
  chr_null_formula, # confounding variable will change according to subtype
  chr_full_formula # confounding variable will change according to subtype
){

  # event_variable: variable (column) name present in xdata. e.g. "DDFS_Event"
  # biomarker: variable (column) name present in xdata. e.g. "Immune1"
  # interaction_variable: variable (column) name present in xdata. e.g. "Chemo"
  # xdata: a tibble/data frame
  # chr_null_formula: formula for null model(no interaction)
  # chr_full_formula: formula for null model(interaction)

  # # Test data
  # event_variable = "OS_Event"
  # biomarker = "Immune1"
  # interaction_variable = "Chemo"
  # xdata = clin_finher %>%
  #   # dplyr::mutate( Proliferation1 = binirize(Proliferation1) %>%
  #   #                  factor(levels = c("low", "high"))) %>%
  #   dplyr::filter(Subtype_IHC_2 == "HER2")
  # chr_null_formula <- paste("Surv( OS_Time_Years, OS_Event ) ~",
  #                           biomarker, "+ Chemo + strata(Hormone, Herceptin)")
  # chr_full_formula <- paste("Surv( OS_Time_Years, OS_Event ) ~",
  #                           biomarker, "* Chemo + strata(Hormone, Herceptin)")

  # The original idea is to filtering out entire starat with no events, to eliminate warnings !!
  # e.g. FEC+DTX+Hormone strata in HR+HER2+ subtype
  # However the coxph documntation stated that coeffcient of only the no-event strata
  # will get affected (estimates are not usable), the estimates of the
  # remaining strata with some events are all reliable.
  # Hence, all strata (all data) is used.
  # Ref: https://rdrr.io/cran/survival/man/coxph.html (See convergense section)


  # interaction null model
  m0 <- coxph(
    formula = as.formula(chr_null_formula),
    data = xdata,
    model = TRUE, # used to identify the reference/non-reference(factor levels) regimes
    x = TRUE # used to identify the reference/non-reference(factor levels) regimes
  )

  # interaction full model
  m1 <- coxph(
    formula = as.formula(chr_full_formula),
    data = xdata,
    model = TRUE, # used to identify the reference/non-reference(factor levels) regimes
    x = TRUE # used to identify the reference/non-reference(factor levels) regimes
  )



  # if( class( xdata[,biomarker] %>% deframe() ) == "factor" ){
  #
  #   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #   # This if condition code can be discarded according to new analysis
  #
  #   # To accomodate factor Proliferation    #
  #   # Note: for interaction this code will generate interaction extimates only for one level of the biomarker.
  #   # i.e For Proliferation1high:DTX and Proliferation1high:NVB,
  #   # but not for Proliferation1low:DTX and Proliferation1low:NVB
  #
  #   rowname = str_c(biomarker, (xdata[,biomarker] %>% deframe() %>% levels())[2])
  #   x <- summarise_interaction_cox(
  #     m1, m0, test_var = rowname, inter_var = interaction_variable
  #   ) %>%
  #     dplyr::mutate(Therapy = Module_name)
  #
  # } else {

    x <- summarise_interaction_cox(
      m1, m0, test_var = biomarker, inter_var = interaction_variable
    ) %>%
      dplyr::mutate(Therapy = Module_name)

  # }

  ph_test <- cox.zph(m1)$table %>%
    tibble::as_tibble(rownames = "Test") %>%
    dplyr::mutate(Model = biomarker)

  x <- data.frame(
    x,
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    Ph_ok = if_else(any(ph_test$p < 0.05), F, T),
    stringsAsFactors = F,
    check.names = F
  )

  x <- x[, c("Coef", "Exp_coef", "P", "L95", "U95", "Therapy", "Variable", "Endpoint", "Event", "N", "Ph_ok")]
  names(x) <- c("Coef", "HR", "P", "Low95", "Up95", "Therapy", "Variable", "Endpoint", "Event", "N", "Ph_ok")

  list(
    main_effect =x,
    ph_test = ph_test
  )

}


# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_prog2 <- function(
  event_variable,
  time_variable0,
  time_variable1,
  time_split_id,
  biomarker,
  xdata
){

  # # # Test data
  # event_variable = "DDFS_Event"
  # time_variable0 = "DDFS_Time_Start"
  # time_variable1 = "DDFS_Time_End"
  # time_split_id = "Interval_Id"
  # biomarker = "Immune1"
  # xdata = clin_finher_split


  m1 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
            biomarker, "*", time_split_id, "-", time_split_id,  "+ strata(Hormone, Herceptin, Chemo)")
    ),
    data = xdata
  )


  # For testing heterogenity(interaction) w.r.t strata
  m1_het0 <- coxph(
    formula = as.formula(
      paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
            biomarker, "*", time_split_id, "+ interaction(Hormone, Herceptin, Chemo)",
            "-", time_split_id)
    ),
    data = xdata
  )

  m1_het1 <- coxph( # m1_het1
    formula = as.formula(
      paste("Surv(", time_variable0, ",", time_variable1, ",", event_variable, ") ~",
            biomarker, "*", time_split_id, "* interaction(Hormone, Herceptin, Chemo)",
            "-", time_split_id)
    ),
    data = xdata
  )


  xterm <- names(m1$coefficients)[str_detect(names(m1$coefficients),time_split_id)]

  x <- cbind(
    summary(m1)$coefficients[, , drop=F],
    summary(m1)$conf.int[, c("lower .95", "upper .95"), drop = F],
    lrtest(m1_het0, m1_het1)[2, "Pr(>Chisq)", drop = F]
  )

  # updating interaction term
  x[xterm, "coef"] <- x[ , "coef"] %>% sum()
  x[, "exp(coef)"] <- x[ , "coef"] %>% exp()

  # updating interaction SE and CI
  se <- sqrt(
      sum(x[, "se(coef)"] ^ 2) +
        (2 * vcov(m1)[biomarker, xterm])
    )
  x[xterm, "se(coef)"] <- se
  x[xterm, "lower .95"] <- (x[xterm, "coef"] - (1.96 * se)) %>% exp()
  x[xterm, "upper .95"] <- (x[xterm, "coef"] + (1.96 * se)) %>% exp()


  ph_test <- cox.zph(m1)$table %>%
    tibble::as_tibble(rownames = "Test") %>%
    dplyr::mutate(Model = biomarker)

  x <- data.frame(
    x,
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    Ph_ok = if_else(any(ph_test$p < 0.05), F, T),
    stringsAsFactors = F,
    check.names = F
  )


  x <- x[, c("coef", "exp(coef)", "Pr(>|z|)", "lower .95", "upper .95", "Pr(>Chisq)", "Variable", "Endpoint", "Event", "N", "Ph_ok")]
  names(x) <- c("Coef", "HR", "P", "Low95", "Up95", "Heterogeneity", "Variable", "Endpoint", "Event", "N", "Ph_ok")


  list(
    main_effect1 = x[biomarker, , drop=F],
    main_effect2 = x[xterm, , drop=F],
    ph_test = ph_test
  )

}



# Test finher adjuvant prognosis with heterogeneity assesed using cox interaction test.
get_adj_inter2 <- function(
  event_variable, # to extract event summary
  biomarker, # to extract model estimates
  interaction_variable, # to extract model estimates
  time_split_id, # to extract model estimates
  xdata,
  # chr_null_formula, # confounding variable will change according to subtype
  chr_full_formula # confounding variable will change according to subtype
){

  # event_variable: variable (column) name present in xdata. e.g. "DDFS_Event"
  # biomarker: variable (column) name present in xdata. e.g. "Immune1"
  # interaction_variable: variable (column) name present in xdata. e.g. "Chemo"
  # time_split_id: variable (column) name present in xdata. e.g. "Interval_Id
  # xdata: a tibble/data frame
  # chr_null_formula: formula for null model(no interaction)
  # chr_full_formula: formula for null model(interaction)

  # # Test data
  # event_variable = "DDFS_Event"
  # biomarker = "Immune1"
  # interaction_variable = "Chemo"
  # time_split_id = "Interval_Id"
  # xdata = clin_finher_split %>%
  #   dplyr::filter(Subtype_IHC_2 == "HER2")
  # # chr_null_formula <- paste("Surv( OS_Time_Years, OS_Event ) ~",
  # #                           biomarker, "+ Chemo + strata(Hormone, Herceptin)")
  # chr_full_formula <- paste("Surv( DDFS_Time_Start, DDFS_Time_End, DDFS_Event ) ~",
  #                           biomarker, "* Chemo *", time_split_id,
  #                           "-", time_split_id, "+ strata(Hormone, Herceptin)")


  # The original idea is to filtering out entire starat with no events, to eliminate warnings !!
  # e.g. FEC+DTX+Hormone strata in HR+HER2+ subtype
  # However the coxph documntation stated that coeffcient of only the no-event strata
  # will get affected (estimates are not usable), the estimates of the
  # remaining strata with some events are all reliable.
  # Hence, all strata (all data) is used.
  # Ref: https://rdrr.io/cran/survival/man/coxph.html (See convergense section)


  # interaction full model
  m1 <- coxph(
    formula = as.formula(chr_full_formula),
    data = xdata,
    model = TRUE, # used to identify the reference/non-reference(factor levels) regimes
    x = TRUE # used to identify the reference/non-reference(factor levels) regimes
  )


  x <- summarise_interaction_cox2(
    m1, test_var = biomarker, inter_var = interaction_variable, time_split_id = time_split_id
  ) %>%
    dplyr::mutate(Therapy = Module_name)

  ph_test <- cox.zph(m1)$table %>%
    tibble::as_tibble(rownames = "Test") %>%
    dplyr::mutate(Model = biomarker)

  x <- data.frame(
    x,
    Variable = biomarker,
    Endpoint = event_variable,
    Event = xdata[,event_variable] %>% sum(na.rm = T),
    N = nrow(xdata),
    Ph_ok = if_else(any(ph_test$p < 0.05), F, T),
    stringsAsFactors = F,
    check.names = F
  )

  x <- x[, c("Coef", "Exp_coef", "P", "L95", "U95", "Therapy", "Variable", "Endpoint", "Event", "N", "Ph_ok")]
  names(x) <- c("Coef", "HR", "P", "Low95", "Up95", "Therapy", "Variable", "Endpoint", "Event", "N", "Ph_ok")

  list(
    main_effect1 = x[str_detect(x$Therapy,str_c(time_split_id,"1")), ],
    main_effect2 = x[str_detect(x$Therapy,str_c(time_split_id,"2")), ],
    ph_test = ph_test
  )

}


estimate_prog_heterogenity <- function(xclin, sig, perm = 100){

  # xclin <- clin_neoadj %>%
  #   dplyr::filter(Subtype_ihc == "HR" &
  #                   Strata != "GSE22226_GPL4133:A0A+T" &
  #                   Strata != "GSE21974:A0A+T" &
  #                   Arm_consolidated == "A0A+Taxane///No_her2_agent///No_hormone_therapy///No_other_therapy")
  # sig = "Cholesterol3"
  # perm = 100


  # Minimum two datasets are required for dataset hetrerogenity assessment
  if( length(xclin$Strata %>% unique()) == 1 ){

     return(

       tibble(
        Module_name = sig,
        Q = NA,
        Q_p = NA,
        I2 = NA,
        I2_lower_ci = NA,
        I2_upper_ci = NA
      )

    )

  }


  xeffect <-  purrr::map(
    xclin$Strata %>% unique(),

    # Strata =  Series_matrix_accession * Arm_consolidated
    function(study, xclin){

      m1_het <- glm(
        formula = as.formula( paste("Response ~", sig) ),
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

  # Test data
  # m0 = m0_inter
  # m1 = m1_inter
  # test_var = sig
  # inter_var = "Arm_consolidated"


  # interaction effect
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

      # i =  "Immune1:Arm_consolidatedAAA"

      # idx[-1] suggest all terms excluding reference term/arm

      # Update Estimate(Log-OR)
      # Interaction term contains change from reference
      # To get actual Estiate(Log-OR), add interaction Estiate(Log-OR)
      # to reference Estiate(Log-OR)
      x[i, "Estimate"] <- sum(x[c(idx[1], i), "Estimate"]) # idx[1] represents reference arm

      # Note in the following code descriptions, coef ~ Estimate(Log-OR).
      # To get CI of sum of coefs find Std.Error of sum of coefs
      # Std.Error of sum of coefs: sqrt( sum(coef variances) + 2 * their covariance)
      # Sum of variances reference from wikipedia, search term: "Variance"
      se <- sqrt(
        sum(x[c(idx[1], i), "Std. Error"] ^ 2) +
        (2 * vcov(m1)[[idx[1], i]])
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
  # Note that this function is identical to summarise_interaction(), except
  # for columnnames of cox model output.

  # m1: full model with interaction term
  # m0: null model without interaction term
  # test_var: variable tested for interaction effect (signature)
  # inter_var: variable interacting with test var
  # Note that this script expects the "test_var" to be a continuous variable,
  # and the "inter_var" to be a factor or a character vector that can be
  # coerased to be a factor.

  # # Test data
  # m0 = m0
  # m1 = m1
  # test_var = biomarker
  # inter_var = interaction_variable


  # interaction effect summary
  # >>>>>>>>>>>>>>>>>>>>>>>>>>
  m1_coef <- summary(m1)$coefficients
  m1_ci <- confint.default(m1)



  # subsetting relevent terms
  # >>>>>>>>>>>>>>>>>>>>>>>>>
  # (i.e. test_var (sig) and terms interacting with test_var)

  idx <- rownames(m1_ci)[str_detect(rownames(m1_ci), test_var)] # this will select Sig and Sig:Arm2
  # Note
  # idx[1] : reference arm
  # The remaining terms are interaction terms


  # Updating estimate, std.err and ci for all relevant rows except the 1st row
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if (identical(rownames(m1_ci), rownames(m1_coef))) {

    x <- cbind(m1_coef[idx, ], m1_ci[idx, ])

    # Update wald p with loglikelihood p
    p <- lrtest(m0, m1)$"Pr(>Chisq)" %>% na.omit() %>% as.numeric()

    # Obsolete code logic
    # Both pvalues in x is set to p, to aid in plotting
    # x[, "Pr(>|z|)"] <- p #

    # New code logic
    # The pvalues can be replicated at any point later in the script,
    # but deletion is difficult. Hence likilihoog pvalue is set for reference term
    x[, "Pr(>|z|)"] <- NA # Resetting original pvalues
    x[idx[1], "Pr(>|z|)"] <- p # Setting likihood ratio test interaction pvalue


    # Updating estimate, std.err and ci for all relevant rows except the 1st row (reference row)

    # The first row is the effect from reference arm.
    # All other rows represents change with respect to the reference row
    # All other rows required to adjust it effect for independet comparison in plots
    # Adjustment 1: add reference effect to each other row's effect
    # Adjustment 2: adjust other row's standard error
    # Adjustment 3: adjust other row's 95% ci

    for(i in idx[-1]){

      # Note: idx[-1] suggest all terms excluding reference term/arm

      # Update coef(Log-HR)
      # Interaction term contains change from reference
      # To get actual coef(Log-HR), add interaction coef(Log-HR)
      # to reference Estiate(Log-OR)
      x[i, "coef"] <- sum(x[c(idx[1], i), "coef"]) # idx[1] represents reference arm
      # update exp(coef) ~ HR
      x[i, "exp(coef)"] <- exp(x[i, "coef"])

      # Note in the following code descriptions, coef ~ Estimate(Log-OR).
      # To get CI of sum of coefs find Std.Error of sum of coefs
      # Std.Error of sum of coefs: sqrt( sum(coef variances) + 2 * their covariance)
      # Sum of variances reference from wikipedia, search term: "Variance"
      # Also see:
      # https://stats.stackexchange.com/questions/401439/how-to-calculate-treatment-effect-and-its-confidence-interval-for-subgroups-in-a
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
    # vector names used in following code
  }


  # Appending reference arm detials to rownames(x)[1]
  non_ref_arm <- grep(x = colnames(m1$x), pattern = ":", value = T)
  non_ref_arm <- str_split_fixed(string = non_ref_arm, pattern = ":", n = 2)[,2]

  ref_arm <- setdiff(
    unique(m1$model[, inter_var]),
    str_replace(string = non_ref_arm, pattern = inter_var, replacement = "")
  )
  ref_arm <- str_c(inter_var, ref_arm)

  if (str_detect(rownames(x)[1], ":")) {
    # Making sure the 1st row is reference row.
    # For reference row ":" is not present.
    stop("Error in interaction summary preparation.")
  } else {
    rownames(x)[1] <- str_c(rownames(x)[1], ":", ref_arm)
  }



  x %>%
    as_tibble(rownames = "Module_name") %>%
    dplyr::rename(
      Exp_coef = "exp(coef)",
      Std_error = "se(coef)",
      Z_value = "z",
      P = "Pr(>|z|)",
      L95 = "2.5 %", # coef 95% CI
      U95 = "97.5 %" # coef 95% CI
    ) %>%
    dplyr::rename_all(~str_to_title(.x))

} # end of fun


summarise_interaction_cox2 <- function(m1, test_var, inter_var, time_split_id){


  # Funtion to extract and summarize interaction effect ready to plot
  # Note that this function is identical to summarise_interaction(), except
  # for columnnames of cox model output.

  # m1: full model with interaction term
  # m0: null model without interaction term
  # test_var: variable tested for interaction effect (signature)
  # inter_var: variable interacting with test var
  # time_split_id:
  # Note that this script expects the "test_var" to be a continuous variable,
  # and the "inter_var" to be a factor or a character vector that can be
  # coerased to be a factor.

  # # Test data
  # m1 = m1
  # test_var = biomarker
  # inter_var = interaction_variable
  # time_split_id = time_split_id

  # interaction effect summary
  # >>>>>>>>>>>>>>>>>>>>>>>>>>
  m1_coef <- summary(m1)$coefficients
  m1_ci <- confint.default(m1)



  # subsetting relevent terms
  # >>>>>>>>>>>>>>>>>>>>>>>>>
  # (i.e. test_var (sig) and terms interacting with test_var)

  idx <- rownames(m1_ci)[str_detect(rownames(m1_ci), test_var)] # this will select Sig and Sig:Arm2
  # Note
  # idx[1] : reference arm
  # The remaining terms are interaction terms

  # updating estimates in two sets(the difference in two sets is the reference term used)
  # in idx1 reference term "biomarker"
  idx1 <- idx[!(str_detect(idx,inter_var) & str_detect(idx,time_split_id))]
  # in idx1 reference term "biomarker:time_split_id" (Interval_Id2)
  idx2 <- idx[str_detect(idx,time_split_id)]


  # Updating estimate, std.err and ci for all relevant rows except the 1st row
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if (identical(rownames(m1_ci), rownames(m1_coef))) {

    x <- cbind(m1_coef[idx, ], m1_ci[idx, ])

    # Wald test pvalue of individual terms are used.
    # No loglikelihood p vlaue
    # Setting non-interaction term pval to NA.
    # The pvalues can be replicated at any point later in the script,
    # but deletion is difficult.
    i <- str_detect(rownames(x),inter_var, negate = T)
    x[i, "Pr(>|z|)"] <- NA


    # Updating estimate, std.err and ci for all relevant rows except the 1st row (reference row)

    # The first row is the effect from reference arm.
    # All other rows represents change with respect to the reference row
    # All other rows required to adjust it effect for independet comparison in plots
    # Adjustment 1: add reference effect to each other row's effect
    # Adjustment 2: adjust other row's standard error
    # Adjustment 3: adjust other row's 95% ci

    # The above logic is repeated in wo sets(idx1 and idx2)
    # In the 1st set (idx1), the reference is "biomarker"
    # In the 2nd set (idx2), the reference is "biomarker:time_split_id"

    # updating set1
    for(i in idx1[-1]){

      # Note: idx1[-1] suggest all terms excluding reference term/arm

      # Update coef(Log-HR)
      # Interaction term contains change from reference
      # To get actual coef(Log-HR), add interaction coef(Log-HR)
      # to reference Estiate(Log-OR)
      x[i, "coef"] <- sum(x[c(idx1[1], i), "coef"]) # idx1[1] represents reference arm
      # update exp(coef) ~ HR
      x[i, "exp(coef)"] <- exp(x[i, "coef"])

      # Note in the following code descriptions, coef ~ Estimate(Log-OR).
      # To get CI of sum of coefs find Std.Error of sum of coefs
      # Std.Error of sum of coefs: sqrt( sum(coef variances) + 2 * their covariance)
      # Sum of variances reference from wikipedia, search term: "Variance"
      # Also see:
      # https://stats.stackexchange.com/questions/401439/how-to-calculate-treatment-effect-and-its-confidence-interval-for-subgroups-in-a
      se <- sqrt(
        sum(x[c(idx1[1], i), "se(coef)"] ^ 2) +
          (2 * vcov(m1)[idx1[1], i])
      )
      x[i, "se(coef)"] <- se
      x[i, "2.5 %"] <- x[i, "coef"] - (1.96 * se) # ci of coef
      x[i, "97.5 %"] <- x[i, "coef"] + (1.96 * se) # ci of coef

    }

    # original
    #                                     coef   exp(coef) se(coef)          z  Pr(>|z|)     2.5 %    97.5 %
    # Immune1                       -1.3827801  0.25088011 1.211099 -1.1417563        NA -3.756491 0.9909306
    # Immune1:ChemoNVB               2.3332841 10.31175064 1.544564  1.5106429 0.1308795 -0.694005 5.3605731
    # Immune1:Interval_Id2           0.9591351  2.60943862 2.408564  0.3982186        NA -3.761564 5.6798341
    # Immune1:ChemoNVB:Interval_Id2 -3.6450245  0.02612077 2.841388 -1.2828326 0.1995507 -9.214042 1.9239928

    # after set1
    #                                    coef  exp(coef)  se(coef)          z  Pr(>|z|)      2.5 %    97.5 %
    # Immune1                       -1.382780 0.25088011 1.2110991 -1.1417563        NA -3.7564908 0.9909306
    # Immune1:ChemoNVB               0.950504 2.58701316 0.9811567  1.5106429 0.1308795 -0.9725632 2.8735711
    # Immune1:Interval_Id2          -0.423645 0.65465625 2.0819270  0.3982186        NA -4.5042220 3.6569320
    # Immune1:ChemoNVB:Interval_Id2 -3.645024 0.02612077 2.8413875 -1.2828326 0.1995507 -9.2140417 1.9239928

    # after set2
    #                                    coef  exp(coef)  se(coef)          z  Pr(>|z|)      2.5 %     97.5 %
    # Immune1                       -1.382780 0.25088011 1.2110991 -1.1417563        NA -3.7564908  0.9909306
    # Immune1:ChemoNVB               0.950504 2.58701316 0.9811567  1.5106429 0.1308795 -0.9725632  2.8735711
    # Immune1:Interval_Id2          -0.423645 0.65465625 2.0819270  0.3982186        NA -4.5042220  3.6569320
    # Immune1:ChemoNVB:Interval_Id2 -4.068669 0.01710013 0.9370571 -1.2828326 0.1995507 -5.9053013 -2.2320375

    # after term expansion
    #                                    coef  exp(coef)  se(coef)          z  Pr(>|z|)      2.5 %     97.5 %
    # Immune1:ChemoDTX:Interval_Id1 -1.382780 0.25088011 1.2110991 -1.1417563        NA -3.7564908  0.9909306
    # Immune1:ChemoNVB:Interval_Id1  0.950504 2.58701316 0.9811567  1.5106429 0.1308795 -0.9725632  2.8735711
    # Immune1:ChemoDTX:Interval_Id2 -0.423645 0.65465625 2.0819270  0.3982186        NA -4.5042220  3.6569320
    # Immune1:ChemoNVB:Interval_Id2 -4.068669 0.01710013 0.9370571 -1.2828326 0.1995507 -5.9053013 -2.2320375


    # updating set2 ()
    for(i in idx2[-1]){

      # Note: idx2[-1] is updated in set1

      x[i, "coef"] <- sum(x[c(idx2[1], i), "coef"]) # idx2[1] represents reference arm
      # update exp(coef) ~ HR
      x[i, "exp(coef)"] <- exp(x[i, "coef"])


      se <- sqrt(
        sum(x[c(idx2[1], i), "se(coef)"] ^ 2) +
          (2 * vcov(m1)[idx2[1], i])
      )


      x[i, "se(coef)"] <- se
      x[i, "2.5 %"] <- x[i, "coef"] - (1.96 * se)
      x[i, "97.5 %"] <- x[i, "coef"] + (1.96 * se)

    }

  } else {
    # As an error catching logic
    x <- rep(NA, 6)
    names(x) = c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "2.5 %", "97.5 %")
    # vector names used in following code
  }


  # Expanding individual term names to reflect full term meaning
  # Note: only implemented for two factor levels !!!!!!!!
  inter_var_val <- str_c(inter_var, unique(m1$model[, inter_var]))
  time_split_id_val <- str_c(time_split_id, unique(m1$model[, time_split_id]))

  if(any(str_detect(rownames(x), inter_var_val[1])))
    inter_var_val[1] <- NA
  if(any(str_detect(rownames(x), inter_var_val[2])))
    inter_var_val[2] <- NA


  if(any(str_detect(rownames(x), time_split_id_val[1])))
    time_split_id_val[1] <- NA
  if(any(str_detect(rownames(x), time_split_id_val[2])))
    time_split_id_val[2] <- NA

  inter_var_val <- na.omit(inter_var_val) %>% as.character()
  time_split_id_val <- na.omit(time_split_id_val) %>% as.character()

  # expanding each rownames, using unique patterns
  if (str_detect(rownames(x)[1], ":", negate = T)) {
    # Making sure the 1st row is reference row.
    # For reference row ":" is not present.
    rownames(x)[1] <- str_c(rownames(x)[1], ":", inter_var_val, ":", time_split_id_val)
  }
  if (str_detect(rownames(x)[2], time_split_id, negate = T)) {
    rownames(x)[2] <- str_c(rownames(x)[2], ":", time_split_id_val)
  }
  if (str_detect(rownames(x)[3], inter_var, negate = T)) {
    xrowname <- str_split(rownames(x)[3],":") %>% unlist()
    rownames(x)[3] <- str_c(xrowname[1], ":", inter_var_val, ":", xrowname[2])
  }




  x %>%
    as_tibble(rownames = "Module_name") %>%
    dplyr::rename(
      Exp_coef = "exp(coef)",
      Std_error = "se(coef)",
      Z_value = "z",
      P = "Pr(>|z|)",
      L95 = "2.5 %", # coef 95% CI
      U95 = "97.5 %" # coef 95% CI
    ) %>%
    dplyr::rename_all(~str_to_title(.x))

} # end of fun
