# finher_her2_ph_fix.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Acounting for non-PH in TRA interaction models in HER2
# The following models failed PH asumptions
# TIL, Fibrosis1-3, and Fibroblast.

# !!!!!!!!!!!! Implement as an R markdown report file.!!!!!!!!!!!!!!!!!!!!!!!!!!



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Explore non-PH in TRA interaction models (TIL, Fibrosis1-3, and Fibroblast)


# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")
load("results/data/expr_finher.RData") # for TILsig computation
load("results/data/tilsig_clean.RData") # for TILsig computation

# Update clin_finher with TILsig score
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sig <- tilsig_clean$ALL %>%
  dplyr::select(Direction, Ncbi_gene_id1)

# Compute module score and update clin_finher
score <- get_module_score(
  x = expr_finher %>% t_tibble(names_x_desc = "CEL_filename"),
  module_list = list(TILsig_imm = sig %>%
                       dplyr::filter(Direction == 1) %>%
                       dplyr::mutate(Direction = 1), # average imm
                     TILsig_fib = sig %>%
                       dplyr::filter(Direction == -1) %>%
                       dplyr::mutate(Direction = 1), # average fib
                     TILsig = sig), # weighted average
  by = "Ncbi_gene_id1"
) %>%
  dplyr::mutate(TILsig_imm = TILsig_imm %>% genefu::rescale(q=0.05),
                TILsig_fib = TILsig_fib %>% genefu::rescale(q=0.05),
                TILsig = TILsig %>% genefu::rescale(q=0.05))

clin_finher <- clin_finher %>% left_join(score, by = "CEL_filename")



# Format clinical data
# >>>>>>>>>>>>>>>>>>>>

clin_finher <- clin_finher %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10)

clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2") %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 HR-HER2+    HER2             89
# 2 HR+HER2+    HER2             91

#
# ==============================================================================



# 2. Explore non-PH in ddfs
# ==============================================================================


xdata <- clin_finher  %>%
  dplyr::filter(Subtype_IHC_2 == "HER2")

xformula <- str_c("Surv(DDFS_Time_Years, DDFS_Event) ~ TIL + TILsig + TILsig_fib",
                  "+ Fibrosis1 + Fibrosis2 + Fibrosis3 + Fibroblast ",
                  "+ Herceptin + Hormone + Chemo")

xdata_split <- survSplit(formula = as.formula(xformula),
                         data = clin_finher  %>%
                           dplyr::filter(Subtype_IHC_2 == "HER2"),
                         cut = c(3),
                         # cut = c(2,4),
                         start = "DDFS_Time_Start", # New variable introduced
                         end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                         zero = 0,
                         episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))

xdata_split2 <- survSplit(formula = as.formula(xformula),
                          data = clin_finher  %>%
                            dplyr::filter(Subtype_IHC_2 == "HER2"),
                          # cut = c(3),
                          cut = c(2,4),
                          start = "DDFS_Time_Start", # New variable introduced
                          end = "DDFS_Time_End", # Renaming of DDFS_Time_Years to avoid confusion
                          zero = 0,
                          episode = "Interval_Id") %>%
  dplyr::mutate(Interval_Id = factor(Interval_Id))


# Regular models with PH test
# >>>>>>>>>>>>>>>>>>>>>>>>>>>

pdf("results/figures/TRA_interaction_regular_model_ph_diagnosis.pdf", width = 5, height = 5, onefile = T)
for(i in c("TIL", "TILsig", "TILsig_fib", "Fibrosis1", "Fibrosis2", "Fibrosis3", "Fibroblast")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i, "* Herceptin + strata(Hormone, Chemo)")),
    data = xdata
  )

  print(str_c(i,": summary"))
  print(summary(m1)$coefficients)
  print(str_c(i,": ph test"))
  print(cox.zph(m1))
  print("")

  plot(cox.zph(m1), main = i)
}

dev.off()

# [1] "TIL: summary"
# coef exp(coef)  se(coef)          z   Pr(>|z|)
# TIL               0.1545398 1.1671207 0.1087553  1.4209865 0.15532068
# HerceptinTRA      0.3068122 1.3590857 0.5208642  0.5890445 0.55583143
# TIL:HerceptinTRA -0.3851188 0.6803698 0.2068467 -1.8618563 0.06262335
# [1] "TIL: ph test"
# chisq df     p
# TIL           0.0299  1 0.863
# Herceptin     1.5919  1 0.207
# TIL:Herceptin 4.4797  1 0.034 !!!
# GLOBAL        5.2280  3 0.156

# [1] "TILsig: summary"
# coef exp(coef)  se(coef)         z   Pr(>|z|)
# TILsig               1.4719092 4.3575466 0.8800898  1.672454 0.09443493
# HerceptinTRA         0.2616801 1.2991109 0.6416496  0.407824 0.68340288
# TILsig:HerceptinTRA -1.5333492 0.2158117 1.2676984 -1.209554 0.22645024
# [1] "TILsig: ph test"
# chisq df      p
# TILsig            2.80  1 0.0945 !!!
# Herceptin         2.25  1 0.1338
# TILsig:Herceptin  9.96  1 0.0016 !!!
# GLOBAL           12.52  3 0.0058 !!!

# [1] "TILsig_fib: summary"
# coef exp(coef)  se(coef)          z  Pr(>|z|)
# TILsig_fib              -1.1355348 0.3212503 0.8066747 -1.4076737 0.1592277
# HerceptinTRA            -0.9935912 0.3702447 0.7098919 -1.3996374 0.1616219
# TILsig_fib:HerceptinTRA  1.1014700 3.0085853 1.1782165  0.9348621 0.3498593
# [1] "TILsig_fib: ph test"
# chisq df      p
# TILsig_fib            1.971  1 0.1604
# Herceptin             2.341  1 0.1260
# TILsig_fib:Herceptin  0.147  1 0.7018
# GLOBAL               11.522  3 0.0092 !!!

# [1] "Fibrosis1: summary"
# coef exp(coef)  se(coef)         z  Pr(>|z|)
# Fibrosis1              -1.117045 0.3272454 0.7797113 -1.432639 0.1519610
# HerceptinTRA           -1.096172 0.3341476 0.6918568 -1.584392 0.1131045
# Fibrosis1:HerceptinTRA  1.330474 3.7828355 1.1561780  1.150752 0.2498344
# [1] "Fibrosis1: ph test"
# chisq df     p
# Fibrosis1            0.9688  1 0.325
# Herceptin            2.3133  1 0.128
# Fibrosis1:Herceptin  0.0742  1 0.785
# GLOBAL              10.5627  3 0.014 !!!

# [1] "Fibrosis2: summary"
# coef exp(coef) se(coef)         z   Pr(>|z|)
# Fibrosis2              -1.141659 0.3192888 0.844550 -1.351796 0.17644059
# HerceptinTRA           -1.257261 0.2844321 0.700838 -1.793939 0.07282289
# Fibrosis2:HerceptinTRA  1.705078 5.5018152 1.222635  1.394593 0.16313872
# [1] "Fibrosis2: ph test"
# chisq df      p
# Fibrosis2            1.787  1 0.1813
# Herceptin            2.176  1 0.1401
# Fibrosis2:Herceptin  0.279  1 0.5971
# GLOBAL              12.725  3 0.0053 !!!

# [1] "Fibrosis3: summary"
# coef exp(coef)  se(coef)         z   Pr(>|z|)
# Fibrosis3              -1.280833 0.2778058 0.8385525 -1.527433 0.12665326
# HerceptinTRA           -1.393960 0.2480909 0.7368035 -1.891902 0.05850404
# Fibrosis3:HerceptinTRA  1.878457 6.5434030 1.2368129  1.518789 0.12881571
# [1] "Fibrosis3: ph test"
# chisq df      p
# Fibrosis3            1.454  1 0.2278
# Herceptin            2.198  1 0.1382
# Fibrosis3:Herceptin  0.137  1 0.7112
# GLOBAL              12.489  3 0.0059 !!!

# [1] "Fibroblast: summary"
# coef exp(coef)  se(coef)         z   Pr(>|z|)
# Fibroblast              -1.614822 0.1989261 0.8593839 -1.879046 0.06023822
# HerceptinTRA            -1.651459 0.1917699 0.7818421 -2.112267 0.03466357
# Fibroblast:HerceptinTRA  2.240205 9.3952614 1.2673178  1.767675 0.07711533
# [1] "Fibroblast: ph test"
# chisq df      p
# Fibroblast            1.025  1 0.3113
# Herceptin             2.105  1 0.1468
# Fibroblast:Herceptin  0.051  1 0.8214
# GLOBAL               12.065  3 0.0072 !!!




# Fixing PH survSplit
# >>>>>>>>>>>>>>>>>>>


pdf("results/figures/TRA_interaction_ph_fixed_model_diagnosis.pdf", width = 5, height = 5, onefile = T)

# TIL (survSplit)
m1 <- coxph(
  # formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ TIL * Herceptin + TIL : Herceptin : Interval_Id + strata(Hormone, Chemo),
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ TIL * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
#                                      coef exp(coef)  se(coef)           z   Pr(>|z|)
# TIL                            0.16931765 1.1844963 0.1135306  1.49138297 0.13586098
# HerceptinTRA                   0.31156086 1.3655549 0.6690365  0.46568589 0.64144035
# TIL:HerceptinTRA              -0.58153373 0.5590403 0.3139666 -1.85221510 0.06399494
# TIL:Interval_Id2              -0.31993794 0.7261941 0.4539518 -0.70478400 0.48094467
# HerceptinTRA:Interval_Id2     -0.05875974 0.9429333 1.2109153 -0.04852506 0.96129779
# TIL:HerceptinTRA:Interval_Id2  0.67289995 1.9599127 0.5823259  1.15553836 0.24787003
cox.zph(m1)
#                             chisq df     p
# TIL                       0.00776  1 0.930
# Herceptin                 0.00594  1 0.939
# TIL:Herceptin             1.74000  1 0.187
# TIL:Interval_Id           2.71919  1 0.099
# Herceptin:Interval_Id     0.29648  1 0.586
# TIL:Herceptin:Interval_Id 1.65435  1 0.198
# GLOBAL                    5.53875  6 0.477
plot(cox.zph(m1), main = "TIL - ph fixed")


# TILsig (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ TILsig * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
#                                        coef    exp(coef)  se(coef)           z   Pr(>|z|)
# TILsig                            1.3846334   3.99336183 0.9599120  1.44245879 0.14917299
# HerceptinTRA                      0.6500074   1.91555505 0.7061366  0.92051229 0.35730511
# TILsig:HerceptinTRA              -3.9086791   0.02006699 1.7181157 -2.27498022 0.02290712
# TILsig:Interval_Id2               0.1283604   1.13696265 2.3849402  0.05382121 0.95707761
# HerceptinTRA:Interval_Id2        -0.7812637   0.45782708 1.5944138 -0.49000062 0.62413346
# TILsig:HerceptinTRA:Interval_Id2  4.8837166 132.12079713 3.1341309  1.55823632 0.11917724
cox.zph(m1)
#                               chisq df    p
# TILsig                       0.0887  1 0.77
# Herceptin                    0.0135  1 0.91
# TILsig:Herceptin             2.2642  1 0.13
# TILsig:Interval_Id           0.3392  1 0.56
# Herceptin:Interval_Id        0.1774  1 0.67
# TILsig:Herceptin:Interval_Id 0.4441  1 0.51
# GLOBAL                       8.6612  6 0.19
plot(cox.zph(m1), main = "TILsig - ph fixed")


# TILsig_fib (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ TILsig_fib * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
#                                            coef    exp(coef)  se(coef)          z   Pr(>|z|)
# TILsig_fib                           -1.1579000  0.314145204 0.8929219 -1.2967539 0.19471588
# HerceptinTRA                         -2.9328556  0.053244776 1.1512808 -2.5474719 0.01085066
# TILsig_fib:HerceptinTRA               3.4811668 32.497617208 1.6136777  2.1572875 0.03098327
# TILsig_fib:Interval_Id2               0.2650276  1.303466940 2.1114156  0.1255213 0.90011088
# HerceptinTRA:Interval_Id2             4.0113795 55.223000082 1.7403207  2.3049657 0.02116849
# TILsig_fib:HerceptinTRA:Interval_Id2 -4.9153235  0.007333345 2.8406231 -1.7303681 0.08356454
cox.zph(m1)
#                                    chisq df    p
# TILsig_fib                       0.00166  1 0.97
# Herceptin                        0.01425  1 0.90
# TILsig_fib:Herceptin             0.31624  1 0.57
# TILsig_fib:Interval_Id           0.55378  1 0.46
# Herceptin:Interval_Id            0.17441  1 0.68
# TILsig_fib:Herceptin:Interval_Id 0.00107  1 0.97
# GLOBAL                           7.36955  6 0.29
plot(cox.zph(m1), main = "TILsig_fib - ph fixed")


# Fibrosis1 (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ Fibrosis1 * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
#                                          coef    exp(coef)  se(coef)          z    Pr(>|z|)
# Fibrosis1                           -1.338717  0.262181741 0.8790062 -1.5229896 0.127761300
# HerceptinTRA                        -2.869235  0.056742329 1.0819241 -2.6519743 0.008002264
# Fibrosis1:HerceptinTRA               3.543737 34.595975988 1.5666682  2.2619578 0.023700012
# Fibrosis1:Interval_Id2               1.383584  3.989172645 2.0628187  0.6707249 0.502395808
# HerceptinTRA:Interval_Id2            4.200672 66.731148328 1.7455476  2.4065066 0.016105911
# Fibrosis1:HerceptinTRA:Interval_Id2 -5.218660  0.005414578 2.7393331 -1.9050842 0.056769123
cox.zph(m1)
#                                  chisq df    p
# Fibrosis1                       0.0876  1 0.77
# Herceptin                       0.0122  1 0.91
# Fibrosis1:Herceptin             0.4490  1 0.50
# Fibrosis1:Interval_Id           0.0103  1 0.92
# Herceptin:Interval_Id           0.1679  1 0.68
# Fibrosis1:Herceptin:Interval_Id 0.0766  1 0.78
# GLOBAL                          5.3388  6 0.50
plot(cox.zph(m1), main = "Fibrosis1 - ph fixed")


# Fibrosis2 (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ Fibrosis2 * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
# coef    exp(coef)  se(coef)         z    Pr(>|z|)
# Fibrosis2                           -1.266554  0.281801165 0.9427606 -1.343452 0.179125692
# HerceptinTRA                        -2.902828  0.054867822 1.0337529 -2.808048 0.004984272
# Fibrosis2:HerceptinTRA               3.817985 45.512389601 1.5701798  2.431559 0.015034005
# Fibrosis2:Interval_Id2               1.309280  3.703505980 2.2079517  0.592984 0.553191885
# HerceptinTRA:Interval_Id2            4.246337 69.849061244 1.7735071  2.394316 0.016651388
# Fibrosis2:HerceptinTRA:Interval_Id2 -5.618790  0.003629028 2.9193171 -1.924693 0.054267732
cox.zph(m1)
# chisq df    p
# Fibrosis2                       0.2852  1 0.59
# Herceptin                       0.0126  1 0.91
# Fibrosis2:Herceptin             0.8532  1 0.36
# Fibrosis2:Interval_Id           0.0257  1 0.87
# Herceptin:Interval_Id           0.1752  1 0.68
# Fibrosis2:Herceptin:Interval_Id 0.2229  1 0.64
# GLOBAL                          7.2566  6 0.30
plot(cox.zph(m1), main = "Fibrosis2 - ph fixed")



# Fibrosis3 (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ Fibrosis3 * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
# coef   exp(coef) se(coef)          z    Pr(>|z|)
# Fibrosis3                           -1.477846  0.22812862 0.943107 -1.5669969 0.117115422
# HerceptinTRA                        -3.236579  0.03929812 1.119675 -2.8906416 0.003844563
# Fibrosis3:HerceptinTRA               4.221007 68.10204163 1.644057  2.5674334 0.010245448
# Fibrosis3:Interval_Id2               1.523076  4.58630963 2.206437  0.6902875 0.490013411
# HerceptinTRA:Interval_Id2            4.549572 94.59189738 1.852457  2.4559668 0.014050616
# Fibrosis3:HerceptinTRA:Interval_Id2 -5.809586  0.00299867 2.906902 -1.9985493 0.045657145
cox.zph(m1)
# chisq df    p
# Fibrosis3                       0.3201  1 0.57
# Herceptin                       0.0132  1 0.91
# Fibrosis3:Herceptin             0.7066  1 0.40
# Fibrosis3:Interval_Id           0.1808  1 0.67
# Herceptin:Interval_Id           0.1730  1 0.68
# Fibrosis3:Herceptin:Interval_Id 0.1644  1 0.69
# GLOBAL                          5.9376  6 0.43
plot(cox.zph(m1), main = "Fibrosis3 - ph fixed")


# Fibroblast (survSplit)
m1 <- coxph(
  formula = Surv(DDFS_Time_Start, DDFS_Time_End, DDFS_Event) ~ Fibroblast * Herceptin * Interval_Id - Interval_Id + strata(Hormone, Chemo),
  data = xdata_split
)

summary(m1)$coefficients
# coef    exp(coef) se(coef)          z    Pr(>|z|)
# Fibroblast                           -1.848251 1.575124e-01 0.968727 -1.9079173 0.056401906
# HerceptinTRA                         -3.531242 2.926854e-02 1.196416 -2.9515170 0.003162171
# Fibroblast:HerceptinTRA               4.576281 9.715237e+01 1.721651  2.6580763 0.007858810
# Fibroblast:Interval_Id2               1.807002 6.092157e+00 2.268850  0.7964398 0.425776434
# HerceptinTRA:Interval_Id2             4.738707 1.142863e+02 1.980074  2.3931969 0.016702275
# Fibroblast:HerceptinTRA:Interval_Id2 -5.853329 2.870327e-03 2.969553 -1.9711145 0.048710785
cox.zph(m1)
# chisq df    p
# Fibroblast                       0.3020  1 0.58
# Herceptin                        0.0143  1 0.90
# Fibroblast:Herceptin             0.6942  1 0.40
# Fibroblast:Interval_Id           0.2141  1 0.64
# Herceptin:Interval_Id            0.1828  1 0.67
# Fibroblast:Herceptin:Interval_Id 0.1835  1 0.67
# GLOBAL                           6.8638  6 0.33
plot(cox.zph(m1), main = "Fibroblast - ph fixed")

dev.off()




# Checking non-linerar terms
# >>>>>>>>>>>>>>>>>>>>>>>>>>

for(i in c("TIL", "TILsig", "TILsig_fib", "Fibrosis1", "Fibrosis2", "Fibrosis3", "Fibroblast")){

  m1 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~", i, "+ Herceptin + strata(Hormone, Chemo)")),
    data = xdata
  )

  m2 <- coxph(
    formula = as.formula(str_c("Surv(DDFS_Time_Years, DDFS_Event) ~ ns(", i, ", df = 3) + Herceptin + strata(Hormone, Chemo)")),
    data = xdata
  )

  print(str_c(i,"; nonlin test p: ", lrtest(m1,m2)$'Pr(>Chisq)' %>% na.omit() %>% round(digits=2)))

}
# No- nonlin term is present
# [1] "TIL; nonlin test p: 0.6"
# [1] "TILsig; nonlin test p: 0.15"
# [1] "TILsig_fib; nonlin test p: 0.15"
# [1] "Fibrosis1; nonlin test p: 0.93"
# [1] "Fibrosis2; nonlin test p: 0.48"
# [1] "Fibrosis3; nonlin test p: 0.7"
# [1] "Fibroblast; nonlin test p: 0.66"


#
# ==============================================================================

