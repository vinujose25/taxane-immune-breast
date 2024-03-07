# s6.2_supp_firth_finher_tn.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Explore whether "coxphf: Cox Regression with Firth's Penalized Likelihood" reduces CI width.
# Firth's Penalization did not reduce CI width. !!!!!!!!!!!




# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore prognosis and chemo interaction in TN
#   (standard cox and firth penalized cox).
# 3. Analysis summary figures




# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher_finneo.RData")


# Format clinical data
# >>>>>>>>>>>>>>>>>>>>

clin_finher_finneo <- clin_finher_finneo %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10) %>%
  dplyr::filter(Subtype_IHC == "TN")


# Verify
clin_finher_finneo %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())
#   Subtype_IHC Subtype_IHC_2     N
# 1 TN          TN              120



# Summarize clinical data
# >>>>>>>>>>>>>>>>>>>>>>>

clin_finher_finneo %>%
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

#   Subtype            Arm  N RFS DDFS OS RFS_perc DDFS_perc  OS_perc
# 1      TN DTX.noTRA.noHR 60  14   13  9 23.33333  21.66667 15.00000
# 2      TN NVB.noTRA.noHR 60  21   19 14 35.00000  31.66667 23.33333


# Remove NAs from model variables as coxphf() don't accept NAs
clin_finher_finneo <- clin_finher_finneo %>%
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

dim(clin_finher_finneo) #[1] 118 119 # 2 NAs

#
# ==============================================================================




# 2. Explore prognosis and chemo interaction in TN.
# ==============================================================================

finher_tn <- list()

# TN prognosis
# >>>>>>>>>>>>


# Individual sig

finher_tn[["prog"]] <- purrr::map(

    c("RFS", "DDFS", "OS"),

    function(event_type, clin){

      print(event_type)

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

            function(sig, event_type, clin){ # subtype,

              # get_adj_prog(
              #   event_variable = str_c(event_type, "_Event"),
              #   time_variable = str_c(event_type, "_Time_Years"),
              #   biomarker = sig,
              #   xdata = clin
              # )
              get_std_firth_cox_out(
                formula_chr = str_c("Surv(", event_type, "_Time_Years, ",
                                    event_type,"_Event ) ~ ", sig,
                                    " + strata(Hormone, Herceptin, Chemo)"),
                biomarker = sig,
                xdata = clin
              )

            },

            event_type,
            clin
          )

          # merge per sig output from get_std_firth_cox_out()
          yy %>% bind_rows()

    },

    clin = clin_finher_finneo
  )

names(finher_tn[["prog"]]) <- c("RFS", "DDFS", "OS")






# TN Chemo interaction
# >>>>>>>>>>>>>>>>>>>>


# individual sig

finher_tn[["inter"]] <- purrr::map(

  c("RFS", "DDFS", "OS"),

  function(event_type, clin){

    print(event_type)

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


      function(sig, event_type, clin){

       get_std_firth_cox_out(
          formula_chr = str_c("Surv(", event_type, "_Time_Years, ",
                              event_type,"_Event ) ~ ", sig,
                              " * Chemo"),
          biomarker = sig,
          xdata = clin
        )

      },
      event_type,
      clin
    )



    # merge per sig output from get_std_firth_cox_out()
    yy %>% bind_rows()

  },

  clin = clin_finher_finneo

)
# "OS" interaction model with "MCPcounter_Fibroblasts" raises the following warning
# with maxstep = 0.5 and maxit = 1000 !!!!!
# Warning message:
#   In coxphf(formula = as.formula(formula_chr), data = xdata, maxit = 1000) :
#   Convergence in estimating profile likelihood CI or p-values not attained for all variables.
# Consider re-run with smaller maxstep and larger maxit.

names(finher_tn[["inter"]]) <- c("RFS", "DDFS", "OS")


#
# ==============================================================================



# 3. Analysis summary
# ==============================================================================

# prognostic model comparison

ggdf <- bind_rows(
  finher_tn$prog$RFS %>% dplyr::mutate(event = "RFS"),
  finher_tn$prog$DDFS %>% dplyr::mutate(event = "DDFS"),
  finher_tn$prog$OS %>% dplyr::mutate(event = "OS")
) %>%
  tidyr::gather("key", "value",
                "hr", "ci.lo", "ci.up",
                "firth_hr", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = finher_tn$prog$RFS$variable),
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
       subtitle = "FinHER TN prognostic models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"FinHER_TN_Prognosis_Firth_vs_Std_COX.pdf"))
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
       subtitle = "FinHER TN prognostic models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"FinHER_TN_Prognosis_log_Firth_vs_Std_COX.pdf"))
print(p)
dev.off()


# interaction model comparison

ggdf <- bind_rows(
  finher_tn$inter$RFS %>% dplyr::mutate(event = "RFS"),
  finher_tn$inter$DDFS %>% dplyr::mutate(event = "DDFS"),
  finher_tn$inter$OS %>% dplyr::mutate(event = "OS")
) %>%
  tidyr::gather("key", "value",
                "hr", "ci.lo", "ci.up",
                "firth_hr", "firth_ci.lo", "firth_ci.up") %>%
  dplyr::mutate(model = if_else(str_detect(key, "firth"),"Firth","Std"),
                key2 = str_replace(key, "firth_", ""),
                group_var = str_c(key2, "-", model),
                variable = factor(variable, levels = finher_tn$inter$RFS$variable),
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
       subtitle = "FinHER TN interaction models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"FinHER_TN_Interaction_Firth_vs_Std_COX.pdf"))
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
       subtitle = "FinHER TN interaction models") +
  facet_wrap(~event, ncol = 1)

pdf(file = str_c(out_figures,"FinHER_TN_Interaction_Firth_vs_Std_COX_log.pdf"))
print(p)
dev.off()


#
# ==============================================================================




# Clear memory
# ==============================================================================

rm(clin_finher_finneo, finher_tn)

#
# ==============================================================================

