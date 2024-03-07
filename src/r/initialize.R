# initialize.R

# backup

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# The idea of this script is to initialize folder structure for the project and
# to load necessary R packages.


# Script structure
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/figures")){
  dir.create(path = "results/figures", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/data")){
  dir.create(path = "results/data", recursive = TRUE, mode = "0700")
}

out_figures <- "results/figures/"
out_data <- "results/data/"


#
# ==============================================================================




# 2. Library: Loading of depended R packages and local functions
# ==============================================================================

# Public
library(tidyverse)
library(MCPcounter) # MCPcounter.estimate: produce cell type estimates
library(genefu) # rescale() based on quantiles
library(writexl) # write_xlsx: wrirte xlsx files in multiple sheets
library(readxl) # read_excel: Read xls and xlsx files
library(altmeta) # metahet: Meta-Analysis Heterogeneity Measures
library(ggpubr) # ggarrange: Arrange nultiple ggplots
                # annotate_figure: Annotate figures including: i) ggplots, ii) arranged ggplots from ggarrange(), grid.arrange() and plot_grid().
library(lmtest) # lrtest: likilihood ratio test
library(survival) # pspline(), coxph()
library(survminer) # ggsurvplot()
# NOTE: surminer depends on  ggpubr package which has ggarrange() function
# Try using ggarrange to arrange ggaurvplot()
library(ggh4x) # guide_stringlegend(),
library(coxphf) # suggested by Stefan
library(logistf) # suggested by Stefan


# library(ggstance) # position_dodgev() # vertical dodging
# library(hablar) # retype: Transforms all elements into simple classes
# library(rmarkdown) # render: To render the Rmarkdown file into specific output format
# library(openxlsx) #
# library(meta) # metabin: Meta-analysis of binary outcome data
# library(MetaGxBreast)
# library(pdftools) # qpdf::pdf_subset: Split, Combine and Compress PDF Files
# library(splines) # ns()
# library(coxme)
# library(lme4)
# library(grid) # unit()
# library(survcomp) # combine.test()
# library(pROC) # auc
# library(GOfuncR) # get_anno_genes()
# library(Homo.sapiens) # to support GOfuncR::get_anno_genes()
# library(reporter) # supsc() https://cran.r-project.org/web/packages/reporter/vignettes/reporter-super.html
# library(ggsci) # scale_color_jco() Ref : https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#:~:text=In%20this%20case%2C%20we'll,that%20are%20color%2Dblind%20friendly.



# Private
source("src/r/functions.R")

#
# ==============================================================================



