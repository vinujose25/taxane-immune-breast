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
