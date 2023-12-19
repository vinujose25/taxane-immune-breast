# s7_supp_firth_finher_her2_v2.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Explore whether firth penalization reduces CI width




# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis, chemo interaction, and herceptin interaction in HER2.
# (standard cox and firth penalized cox).
# 3. Analysis summary figures



# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher_finneo.RData")



# Format clinical data
# >>>>>>>>>>>>>>>>>>>>
clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10)

clin_finher_finneo  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 HR-HER2+    HER2             89
# 2 HR+HER2+    HER2             91



# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

clin_finher_finneo %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(
    Subtype = Subtype_IHC, # Subtype_IHC_2
    Arm
  ) %>%
  dplyr::summarise(N = n(),
                   RFS = which(RFS_Event == 1) %>% length(),
                   DDFS = which(DDFS_Event == 1) %>% length(),
                   OS = which(OS_Event == 1) %>% length()) %>%
  dplyr::mutate(
    RFS_perc = (RFS/N)*100,
    DDFS_perc = (DDFS/N)*100,
    OS_perc = (OS/N)*100
  ) %>%
  as.data.frame()

#   Subtype            Arm  N RFS DDFS OS RFS_perc DDFS_perc   OS_perc
# 1 HR-HER2+ DTX.noTRA.noHR 27   7    7  5 25.92593  25.92593 18.518519
# 2 HR-HER2+ NVB.noTRA.noHR 19   6    6  6 31.57895  31.57895 31.578947
# 3 HR-HER2+   DTX.TRA.noHR 26   3    3  2 11.53846  11.53846  7.692308
# 4 HR-HER2+   NVB.TRA.noHR 17   6    5  3 35.29412  29.41176 17.647059
# 5 HR+HER2+   DTX.noTRA.HR 18   4    4  2 22.22222  22.22222 11.111111
# 6 HR+HER2+   NVB.noTRA.HR 27   7    6  3 25.92593  22.22222 11.111111
# 7 HR+HER2+     DTX.TRA.HR 16   2    0  0 12.50000   0.00000  0.000000 # no ddfs/os event
# 8 HR+HER2+     NVB.TRA.HR 30  10    9  4 33.33333  30.00000 13.333333

# Average size of strata, c(27,19,26,17,18,27,16,30) %>% mean() = 22.5


# Remove NAs from model variables as coxphf() don't accept NAs
clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::filter_at(c( "TIL",

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

                      "MCPcounter_Fibroblasts",

                      "scaled_Denovo_TILsig",
                      "scaled_Denovo_Immune",
                      "scaled_Denovo_ECM",

                      "scaled_Pooled_Immune",
                      "scaled_Pooled_Interferon",
                      "scaled_Pooled_Fibrosis",
                      "scaled_Pooled_Cholesterol",

                      "RFS_Event",
                      "RFS_Time_Years",
                      "DDFS_Event",
                      "DDFS_Time_Years",
                      "OS_Event",
                      "OS_Time_Years",

                      "Hormone",
                      "Herceptin",
                      "Chemo"
  ),
  ~(!is.na(.x)))

# N = 170

#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_her2 <- list()

# HER2 prognosis
# >>>>>>>>>>>>>>


# >>>>>>>>>>> Individual.sig

finher_her2[["prog"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

    xx <- purrr::map(

      # c("HR-HER2+", "HR+HER2+", "HER2"),
      c("HER2"),

      function(subtype, event_type, clin){

        print(subtype)

        subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
        xdata <- subset(x = clin,
                        subset = (clin[, subtype_varible] == subtype))


        # Per sig prognosis + heterogenity
        yy <- purrr::map(

          c(
            "TIL",

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

          function(sig, event_type, xdata){ # subtype,


            get_std_firth_cox_out(
              formula_chr = str_c("Surv(", event_type, "_Time_Years, ",
                                  event_type,"_Event ) ~ ", sig,
                                  " + strata(Hormone, Herceptin, Chemo)"),
              biomarker = sig,
              xdata = clin
            )

            # # Prognostic model
            # paste("Surv(", time_variable, ",", event_variable, ") ~",
            #           biomarker, "+ strata(Hormone, Herceptin, Chemo)")

            # # Heterogenity/interaction model
            # paste("Surv(", time_variable, ",", event_variable, ") ~",
            #           biomarker,"* interaction(Hormone, Herceptin, Chemo)")

          },

          event_type,
          xdata
        )

        yy %>% bind_rows()

      },
      event_type,
      clin
    )

    # names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")
    names(xx) <- c("HER2")

    xx
  },

  clin = clin_finher_finneo
)

names(finher_her2[["prog"]]) <- c("RFS", "DDFS", "OS")



# HER2 Chemo interaction; all HER2 + strata(Hormone, Herceptin)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# >>>>>>>>>>> Individual.sig

finher_her2[["inter"]] <- purrr::map(

    c("RFS", "DDFS", "OS"),

    function(event_type, clin){

      print(event_type)

      xx <- purrr::map(

        # c("HR-HER2+", "HR+HER2+", "HER2"),
        c("HER2"),

        function(subtype, event_type, clin){

          print(subtype)

          subtype_varible = if_else(subtype == "HER2", "Subtype_IHC_2", "Subtype_IHC")
          xdata <- subset(x = clin,
                          subset = (clin[, subtype_varible] == subtype))



          # Per sig mta interaction
          yy <- purrr::map(

            c(
              "TIL",

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


            function(sig, event_type, xdata){

              # event_variable = str_c(event_type, "_Event")
              # time_variable = str_c(event_type, "_Time_Years")
              #
              #
              #   chr_null_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
              #                             biomarker, "+ Chemo + strata(Hormone, Herceptin)")
              #   chr_full_formula <- paste("Surv(", time_variable, ",", event_variable, ") ~",
              #                             biomarker, "* Chemo + strata(Hormone, Herceptin)")
              #
              #
              # get_adj_inter(
              #   event_variable = event_variable, # to extract event summary
              #   biomarker = biomarker,
              #   interaction_variable = "Chemo",
              #   xdata = xdata,
              #   chr_null_formula = chr_null_formula,
              #   chr_full_formula = chr_full_formula
              # )


              get_std_firth_cox_out(
                formula_chr = str_c("Surv(", event_type, "_Time_Years, ",
                                    event_type,"_Event ) ~ ", sig,
                                    " * Chemo + strata(Hormone, Herceptin)"),
                biomarker = sig,
                xdata = xdata
              )
            },
            event_type,
            xdata
          )



          # # Consolidate per sig stat and format
          # list(
          #   main_effect = purrr::map(yy,~(.x$main_effect)) %>%
          #     bind_rows() %>%
          #     dplyr::mutate(
          #       # P_adj = p.adjust(p = P, method = "BH"),
          #       # Discarding TIL from FDR coorection
          #       # P_adj = if_else(str_detect(Variable, "TIL"), NA_real_, P) %>%
          #       P_adj = if_else(Variable == "TIL", NA_real_, P) %>%
          #         p.adjust(method = "BH"),
          #       Subtype = subtype
          #     ),
          #
          #   ph_test = purrr::map(yy,~(.x$ph_test)) %>%
          #     bind_rows()
          # )

          yy %>% bind_rows()


        },
        event_type,
        clin
      )

      # names(xx) <- c("HR-HER2+", "HR+HER2+", "HER2")
      names(xx) <- c("HER2")

      xx
    },

    clin = clin_finher_finneo

  )

names(finher_her2[["inter"]]) <- c("RFS", "DDFS", "OS")
# Warning message:
#   In coxphf(formula = as.formula(formula_chr), data = xdata, maxit = 1000) :
#   Convergence in estimating profile likelihood CI or p-values not attained for all variables.
# Consider re-run with smaller maxstep and larger maxit.


#
# ==============================================================================



# 3. Analysis summary
# ==============================================================================


# prognostic model comparison

ggdf <- bind_rows(
  finher_her2$prog$RFS$HER2 %>% dplyr::mutate(event = "RFS"),
  finher_her2$prog$DDFS$HER2 %>% dplyr::mutate(event = "DDFS"),
  finher_her2$prog$OS$HER2 %>% dplyr::mutate(event = "OS")
) %>%
  tidyr::gather("key", "value",
                "hr", "ci.lo", "ci.up",
                "firth_hr", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = finher_her2$prog$RFS$HER2$variable),
                group_var = factor(group_var,
                                   levels = c("ci.lo-Std", "ci.lo-Firth",
                                              "hr-Std", "hr-Firth",
                                              "ci.up-Std", "ci.up-Firth")))


p <- ggplot(ggdf, aes(x = variable, y = value, group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Cox model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "TIL-H&E / Gene-modules", y = "Estimates (HR and 95% CI)",
       title = "Firth vs Standard Cox estimates",
       subtitle = "FinHER HER2 prognostic models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"Firth_vs_Std_COX-FinHER_HER2_Prognosis.pdf"))
print(p)
dev.off()


p <- ggplot(ggdf, aes(x = variable, y = log(value), group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Cox model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "TIL-H&E / Gene-modules", y = "Log - Estimates (HR and 95% CI)",
       title = "Firth vs Standard Cox estimates",
       subtitle = "FinHER HER2 prognostic models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"Firth_vs_Std_COX-FinHER_HER2_Prognosis_log.pdf"))
print(p)
dev.off()


# interaction model comparison

ggdf <- bind_rows(
  finher_her2$inter$RFS$HER2 %>% dplyr::mutate(event = "RFS"),
  finher_her2$inter$DDFS$HER2 %>% dplyr::mutate(event = "DDFS"),
  finher_her2$inter$OS$HER2 %>% dplyr::mutate(event = "OS")
) %>%
  tidyr::gather("key", "value",
                "hr", "ci.lo", "ci.up",
                "firth_hr", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = finher_her2$inter$RFS$HER2$variable),
                group_var = factor(group_var,
                                   levels = c("ci.lo-Std", "ci.lo-Firth",
                                              "hr-Std", "hr-Firth",
                                              "ci.up-Std", "ci.up-Firth")))


p <- ggplot(ggdf, aes(x = variable, y = value, group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Cox model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "TIL-H&E / Gene-modules", y = "Estimates (HR and 95% CI)",
       title = "Firth vs Standard Cox estimates",
       subtitle = "FinHER HER2 interaction models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"Firth_vs_Std_COX-FinHER_HER2_Interaction.pdf"))
print(p)
dev.off()


p <- ggplot(ggdf, aes(x = variable, y = log(value), group = group_var)) +
  geom_point(aes(color = group_var, shape = model)) +
  geom_line(aes(color = group_var)) +
  scale_shape_manual(values = c(1,3)) +
  guides(shape = guide_legend(title = "Cox model"),
         color = guide_legend(title = "Estimate")) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "TIL-H&E / Gene-modules", y = "Log - Estimates (HR and 95% CI)",
       title = "Firth vs Standard Cox estimates",
       subtitle = "FinHER HER2 interaction models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"Firth_vs_Std_COX-FinHER_HER2_Interaction_log.pdf"))
print(p)
dev.off()


#
# ==============================================================================




