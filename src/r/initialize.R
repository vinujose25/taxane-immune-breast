# main.R

# backup

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# The idea of this script is to initialize folder structure for the project and
# to load necessary R packages.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/tables")){
  dir.create(path = "results/tables", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/figures")){
  dir.create(path = "results/figures", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/data")){
  dir.create(path = "results/data", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/documents")){
  dir.create(path = "results/documents", recursive = TRUE, mode = "0700")
}


out_tables <- "results/tables/"
out_figures <- "results/figures/"
out_data <- "results/data/"
out_documents <- "results/documents/"


#
# ==============================================================================




# 2. Library: Loading of depended R packages and local functions
# ==============================================================================

# Public
library(tidyverse)
library(ggstance) # position_dodgev() # vertical dodging
library(MCPcounter) # MCPcounter.estimate: produce cell type estimates
# library(hablar) # retype: Transforms all elements into simple classes
library(rmarkdown) # render: To render the Rmarkdown file into specific output format
# library(openxlsx) #

library(writexl) # write_xlsx: wrirte xlsx files in multiple sheets
# library(readxl) # read_excel: Read xls and xlsx files
library(altmeta) # metahet: Meta-Analysis Heterogeneity Measures
# library(meta) # metabin: Meta-analysis of binary outcome data
# library(MetaGxBreast)
library(ggpubr) # ggarrange: Arrange nultiple ggplots
                # annotate_figure: Annotate figures including: i) ggplots, ii) arranged ggplots from ggarrange(), grid.arrange() and plot_grid().
library(pdftools) # qpdf::pdf_subset: Split, Combine and Compress PDF Files
library(genefu) # rescale() based on quantiles
library(lmtest) # lrtest: likilihood ratio test
library(survival) # pspline(), coxph()
library(splines) # ns()
# library(coxme)
# library(lme4)
library(survminer) # ggsurvplot()
# NOTE: surminer depends on  ggpubr package which has ggarrange() function
# Try using ggarrange to arrange ggaurvplot()
library(grid) # unit()
library(survcomp) # combine.test()
library(pROC) # auc
library(ggh4x) # guide_stringlegend(),
library(GOfuncR) # get_anno_genes()
library(Homo.sapiens) # to support GOfuncR::get_anno_genes()

library(coxphf) # suggested by Stefan
library(logistf) # suggested by Stefan


# library(ggsci) # scale_color_jco() Ref : https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#:~:text=In%20this%20case%2C%20we'll,that%20are%20color%2Dblind%20friendly.



# Private
source("src/r/functions.R")

#
# ==============================================================================



