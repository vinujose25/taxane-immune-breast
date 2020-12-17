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
  # e) Discard treatment regimens with samples only from a single study,
  # f) Discard treatment regimens with sample size less than 50, and
  # g) Discard studies from a treatment regimen with less than 10 samples.
  #
  # The present function performs steps b - g.



  # Script structure:
  # For each IHC subtype (HR,HER2,TN):
  # b) Discard samples with missing treatment data
  # c) Discard samples with missing clinical endpoint (pCR/DFS) data
  # d) Discard matched samples from longitudnal studies (series matrices)
  #   (i.e. consider only pre-treatment expression profiles),
  # e) Discard treatment regimens with samples only from a single study,
  # f) Discard treatment regimens with sample size less than 50, and
  # g) Discard studies from a treatment regimen with less than 10 samples.



  # x: clinical data tibble
  # type: desription of x, used for printing informative messages.

  # x = clin %>%
  #   dplyr::filter(Regimen_updated == "neoadj")
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



  # e) Discard treatment regimens with samples only from a single study
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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



  # f) Discard studies from a treatment regimen with less than 10 samples
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  xsum <- x %>%
    dplyr::group_by(Arm_consolidated, Series_matrix_accession) %>%
    dplyr::summarise(n_samp = n(), .groups = "keep") %>%
    dplyr::filter(n_samp < 10)

  cat(
    type,
    ": Discarding",
    xsum$n_samp %>% sum(),
    "samples from",
    xsum %>% nrow(),
    "series matrices with <10 samples per treatment regimen.\n"
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




  # g) Discard treatment regimens with sample size less than 50
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
