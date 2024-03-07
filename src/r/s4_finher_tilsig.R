# finher_tilsig.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Generate TILsig from log2 TIL counts from TN, HER2 and TN+HER2
# Decipher TILsig based on its associated biological processes !!! (Deciphering will
# remove any genes associated with confounding factors involved
# in the TIL's semi quantitative measurement process) !!!!!



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Generate TILsig and compare TILsig based on geneset and enriched GO BP terms
# 3. Decipher TILsig based on associated biological processes
# 4. Save Robjects



# 1. Load and format clinical data.
# ==============================================================================

load("results/data/clin_finher.RData")
load("results/data/expr_finher.RData")

clin_finher <- clin_finher %>%
  dplyr::mutate(Arm = interaction(Chemo, Herceptin, Hormone),
                TIL = StrLy_Mean/10)

clin_finher %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarize(N = n())

#   Subtype_IHC Subtype_IHC_2     N
# 1 TN          TN              120
# 2 HR-HER2+    HER2             89
# 3 HR+HER2+    HER2             91


dim(expr_finher) # 3350 301
expr_finher[1:3,1:3]
#   Ncbi_gene_id `A1764-001.CEL.gz` `A1764-002.CEL.gz`
# 1 ncbi_8355                  5.25               5.92
# 2 ncbi_126282                8.37               8.72
# 3 ncbi_55199                 7.02               5.81

#
# ==============================================================================



# 2. Generate TILsig and compare TILsig based on geneset and enriched GO BP terms
# ==============================================================================

tilsig <- purrr::map(

  c("TN", "HER2", "ALL"),

  function(subtype, clin, expr){

    print(subtype)

    xclin <- switch(subtype,
                    TN = clin %>% dplyr::filter(Subtype_IHC_2 == "TN"),
                    HER2 = clin %>% dplyr::filter(Subtype_IHC_2 == "HER2"),
                    ALL = clin)

    xexpr <- expr %>%
      dplyr::select(Ncbi_gene_id, all_of(xclin$CEL_filename)) %>%
      t_tibble(names_x_desc = "CEL_filename")

    xformula <- switch(subtype,
                       TN = "log2(StrLy_Mean + 1) ~ gene",
                       HER2 = "log2(StrLy_Mean + 1) ~ gene + Subtype_IHC", #HR+HER2+ format
                       ALL = "log2(StrLy_Mean + 1) ~ gene + Subtype_IHC") #HR+HER2+ format

    sig <- purrr::map_df(
      xexpr[,-1],
      function(gene, xclin, xformula){
        x <- lm(formula = as.formula(xformula), data = xclin) %>% summary()
        x$coefficients["gene", , drop = T] # matrix will give unexpected results
      },
      xclin,
      xformula
    )

    sig %>%
      dplyr::rename(Std_error = "Std. Error",
                    T_value = "t value",
                    P = "Pr(>|t|)") %>%
      dplyr::mutate(P_adj = p.adjust(p = P, method = "BH"),
                    Direction = if_else(Estimate < 0, -1, 1),
                    Ncbi_gene_id1 = names(xexpr)[-1],
                    Ncbi_gene_id2 = str_replace(Ncbi_gene_id1,"ncbi_",replacement = ""))

  },
  clin = clin_finher,
  expr = expr_finher
)

names(tilsig) <- c("TN", "HER2", "ALL")

tilsig_clean <- purrr::map(
  tilsig,
  ~(.x %>% dplyr::filter(P_adj<.01 & str_detect(Ncbi_gene_id2, "///", negate = T)))
)

purrr::map(tilsig_clean,dim)
# $TN
# [1] 484   8
# $HER2
# [1] 171   8
# $ALL
# [1] 994   8



# gene set agreement
# >>>>>>>>>>>>>>>>>>

tilsig_agreement(x = tilsig_clean$TN, y = tilsig_clean$ALL)
# [1] "% of x(n=484) in y(n=994): 95 (n=460); % common genes with same direction: 100"
tilsig_agreement(x = tilsig_clean$HER2, y = tilsig_clean$ALL)
# [1] "% of x(n=171) in y(n=994): 98 (n=168); % common genes with same direction: 100"
tilsig_agreement(x = tilsig_clean$HER2, y = tilsig_clean$TN)
# [1] "% of x(n=171) in y(n=484): 78 (n=134); % common genes with same direction: 100"
tilsig_agreement(x = bind_rows(tilsig_clean$TN, tilsig_clean$HER2), y = tilsig_clean$ALL)
# [1] "% of x(n=655) in y(n=994): 96 (n=628); % common genes with same direction: 100"

# There is >78% agreement between different versions of TILsig.
# Further, among common genes the direction is identical.


# biological agreement
# >>>>>>>>>>>>>>>>>>>>

tilsig_print(x = tilsig_clean$TN, label = "tilsig_tn", path = out_data)
tilsig_print(x = tilsig_clean$HER2, label = "tilsig_her2", path = out_data)
tilsig_print(x = tilsig_clean$ALL, label = "tilsig_all", path = out_data)


# background for david analysis
x = expr_finher %>%
  dplyr::select(Ncbi_gene_id) %>%
  dplyr::filter(str_detect(Ncbi_gene_id,"///",negate = T)) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id,"ncbi_",""))

write.table(x = x,
            sep = "\t",
            row.names = F,
            col.names = T,
            file = str_c(out_data, "tilsig_background_genes_for_overrepresentation_analysis.txt"))


# load DAVID results
x <- list.files("results/data/tilsig_david_results/")
david <- vector(mode = "list", length = length(x))
names(david) = x

for(i in names(david)){
  david[[i]] <- read_delim(file = str_c("results/data/tilsig_david_results/", i),delim = "\t")
}
names(david) <- str_replace(names(david),"\\.txt","") %>%
  str_replace("dn","down") %>%
  str_to_upper()

david_clean <- purrr::map(names(david),
                          function(i, david){
                            x <- david[[i]] %>%
                              dplyr::filter( FDR < .05) %>%
                              dplyr::select(Term)
                            names(x) <- i
                            x
                          },
                          david)
names(david_clean) <- names(david)

# padding "" terms to make equal no. of rows in each list to apply bind_cols()
max_terms <- purrr::map(david_clean, nrow) %>% unlist() %>% max()
david_clean <- purrr::map(david_clean,
                          function(x, max_terms){

                            if(nrow(x) < max_terms){
                              xx <- tibble(nme = rep("", times = max_terms - nrow(x)))
                              names(xx) <- names(x)

                              return(bind_rows(x, xx))

                              } else {

                              return(x)

                            }

                          },max_terms)

x <- list()
x[["TN+HER2"]] <- bind_cols(david_clean[str_detect(names(david_clean), "ALL")])
x[["TN"]] <- bind_cols(david_clean[str_detect(names(david_clean), "TN")])
x[["HER2"]] <- bind_cols(david_clean[str_detect(names(david_clean), "HER2")])

write_xlsx(x,path = str_c(out_data, "tilsig_david_results_fdr.05.xlsx"))

# Similar biological process across up and down regulated gene independently

#
# ==============================================================================



# 3. Decipher TILsig based on associated biological processes
# ==============================================================================

# common enriched terms in up geneset in all versions of tilsig
up_terms <- david_clean[str_detect(names(david_clean),"UP")] %>%
  unlist() %>% as.character() %>% table()
up_terms <- names(up_terms)[up_terms == 3]
# [1] "GO:0002479~antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent"
# [2] "GO:0006955~immune response"
# [3] "GO:0006956~complement activation"
# [4] "GO:0006958~complement activation, classical pathway"
# [5] "GO:0038095~Fc-epsilon receptor signaling pathway"
# [6] "GO:0045087~innate immune response"
# [7] "GO:0050776~regulation of immune response"
# [8] "GO:0060333~interferon-gamma-mediated signaling pathway"


# common enriched terms in up geneset in all versions of tilsig
dn_terms <- david_clean[str_detect(names(david_clean),"DOWN")] %>%
  unlist() %>% as.character() %>% table()
dn_terms <- names(dn_terms)[dn_terms == 3]
# [1] "GO:0007155~cell adhesion"
# [2] "GO:0030198~extracellular matrix organization"


# Extracts common up and down enriched geneset from all versions of TILsig
#
# up sig
up_sig <- purrr::map(
  david[str_detect(names(david),"UP")],
  function(x,up_terms){
    x <- x %>%
      dplyr:::filter(Term %in% all_of(up_terms))

    y <- x$Genes
    names(y) <- x$Term

    purrr::map(y,function(x){
      tibble(Ncbi_gene_id = str_split(x, ", ") %>% unlist(),
             Direction = 1) %>%
        dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))
    })
  },
  up_terms
)

# down sig
dn_sig <- purrr::map(
  david[str_detect(names(david),"DOWN")],
  function(x,dn_terms){
    x <- x %>%
      dplyr:::filter(Term %in% all_of(dn_terms))

    y <- x$Genes
    names(y) <- x$Term

    purrr::map(y,function(x){
      tibble(Ncbi_gene_id = str_split(x, ", ") %>% unlist(),
             Direction = 1) %>%
        dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))
    })
  },
  dn_terms
)


# up gene agremment
purrr::map_df(up_terms, function(nme,up_sig){
  list(
    term = nme %>% str_sub(end = 40),
    size = str_c("all:", up_sig$ALL_UP[[nme]] %>% nrow(), ", ",
             "her2:", up_sig$HER2_UP[[nme]] %>% nrow(), ", ",
             "tn:", up_sig$TN_UP[[nme]] %>% nrow()),
    her2_in_all_prop = sum(up_sig$HER2_UP[[nme]]$Ncbi_gene_id %in% up_sig$ALL_UP[[nme]]$Ncbi_gene_id) / nrow(up_sig$HER2_UP[[nme]]),
    tn_in_all_prop = sum(up_sig$TN_UP[[nme]]$Ncbi_gene_id %in% up_sig$ALL_UP[[nme]]$Ncbi_gene_id) / nrow(up_sig$TN_UP[[nme]])
  )
},up_sig)
#   term                                     size                   her2_in_all_prop tn_in_all_prop
# 1 GO:0002479~antigen processing and presen all:22, her2:9, tn:12                 1              1
# 2 GO:0006955~immune response               all:29, her2:19, tn:20                1              1
# 3 GO:0006956~complement activation         all:8, her2:7, tn:7                   1              1
# 4 GO:0006958~complement activation, classi all:9, her2:8, tn:8                   1              1
# 5 GO:0038095~Fc-epsilon receptor signaling all:21, her2:9, tn:12                 1              1
# 6 GO:0045087~innate immune response        all:22, her2:12, tn:15                1              1
# 7 GO:0050776~regulation of immune response all:15, her2:11, tn:11                1              1
# 8 GO:0060333~interferon-gamma-mediated sig all:14, her2:8, tn:10                 1              1

# ALL version of siglist contains all genes from her2 and tn
# Hence consider siglist from ALL version


# down gene agremment
purrr::map_df(dn_terms, function(nme,dn_sig){
  list(
    term = nme %>% str_sub(end = 40),
    size = str_c("all:", dn_sig$ALL_DOWN[[nme]] %>% nrow(), ", ",
                 "her2:", dn_sig$HER2_DOWN[[nme]] %>% nrow(), ", ",
                 "tn:", dn_sig$TN_DOWN[[nme]] %>% nrow()),
    her2_in_all_prop = sum(dn_sig$HER2_DOWN[[nme]]$Ncbi_gene_id %in% dn_sig$ALL_DOWN[[nme]]$Ncbi_gene_id) / nrow(dn_sig$HER2_DOWN[[nme]]),
    tn_in_all_prop = sum(dn_sig$TN_DOWN[[nme]]$Ncbi_gene_id %in% dn_sig$ALL_DOWN[[nme]]$Ncbi_gene_id) / nrow(dn_sig$TN_DOWN[[nme]])
  )
},dn_sig)
#   term                                     size                   her2_in_all_prop tn_in_all_prop
# 1 GO:0007155~cell adhesion                 all:50, her2:14, tn:29                1          0.966
# 2 GO:0030198~extracellular matrix organiza all:47, her2:15, tn:29                1          1


# TILsig_bp gene agreement summary
#
# ALL version of siglist contains all genes from her2 and tn (one gene from tn is not present)
# Hence consider siglist from ALL version

# Further, due to similarity in genesets and biological signals between different
# versions of TILsigup/TILsig_dn (from TN, HER2, TN+HER2),
# use TILsig generated from TN+HER2 samples for downstream analysis.



# tilsig_bp from ALL dataset
# >>>>>>>>>>>>>>>>>>>>>>>>>>

tilsig_bp <- c(up_sig$ALL_UP, dn_sig$ALL_DOWN)

tilsig_bp_legend <- tibble(
  Term = c(up_sig$ALL_UP %>% names(), dn_sig$ALL_DOWN %>% names()),
  Direction = c(rep("up_regulated", length(up_sig$ALL_UP)),
                rep("down_regulated", length(dn_sig$ALL_DOWN)))
) %>%
  dplyr::mutate(
    Signame = purrr::map_chr(Term,
                             ~switch(.x,
                                     "GO:0002479~antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent" = "APP",
                                     "GO:0006955~immune response" = "Immune",
                                     "GO:0050776~regulation of immune response" = "Immune_reg.",
                                     "GO:0038095~Fc-epsilon receptor signaling pathway" = "Fc_epsilon",
                                     "GO:0060333~interferon-gamma-mediated signaling pathway" = "IFN_gamma",
                                     "GO:0045087~innate immune response" = "Immune_innate",
                                     "GO:0006958~complement activation, classical pathway" = "Complement1",
                                     "GO:0006956~complement activation" ="Complement2",
                                     "GO:0030198~extracellular matrix organization" = "ECM",
                                     "GO:0007155~cell adhesion" = "Adhesion"
                                     ))
  )

if(identical(names(tilsig_bp), tilsig_bp_legend$Term)){
  names(tilsig_bp) <- tilsig_bp_legend$Signame
}
names(tilsig_bp)
# [1] "APP"           "Immune"        "Immune_reg."   "Fc_epsilon"
# [5] "IFN_gamma"     "Immune_innate" "Complement1"   "Complement2"
# [9] "ECM"           "Adhesion"


# ### added on 22 mar 2022


# tilsig_bp geneset agreement plot

module_overlap <- module_gene_overlap(lst = tilsig_bp)

module_annot = module_overlap %>%
  dplyr::select(Module_name, Module_name2) %>%
  dplyr::mutate(
    Module_group = purrr::map_chr(
      Module_name,
      ~case_when(
        .x %in% c("APP", "Immune", "Immune_reg.", "Fc_epsilon", "IFN_gamma",
                  "Immune_innate", "Complement1", "Complement2") ~ "Immune",
        .x %in% c("ECM", "Adhesion") ~ "ECM"
      )),
    Module_group_col = purrr::map_chr(
      Module_group,
      ~case_when(
        .x == "Immune" ~ "tomato",
        .x == "ECM" ~ "steelblue"
      ))
  )


p_agreement <- module_gene_overlap_plot(

  module_overlap = module_overlap,
  overlap_col = c(low = "gold", mid = "white", high = "blue"),
  module_label = module_annot$Module_name2,
  module_group = module_annot$Module_group,
  module_group_col = module_annot$Module_group_col,
  module_group_vline = c(0.5, 8.5,10.5),
  module_group_hline = c(0.5, 8.5,10.5),
  line_col = "gray30",
  ticksize.x = 5,
  ticksize.y = 4.5,
  angle = 45,
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2

)



# Module-score
x <- expr_finher %>% t_tibble(names_x_desc = "Sample_id") # Transpose dataset
score <- get_module_score_2(
  x = x,
  module_list = tilsig_bp,
  by = "Ncbi_gene_id"
)
# rm(x)

p_correlation <- module_correlation_plot(
  module_score = score, # module score output from get_module_score_2()
  score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
  module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
  module_group_col = module_annot$Module_group_col, # char vector of module group color, preserve order
  module_group_vline = c(0.5, 8.5,10.5), # xintercept for vlines to group heatmap
  module_group_hline = c(0.5, 8.5,10.5), # yintercept for hlines to group heatmap
  line_col = "gray30",
  ticksize.x = 5,
  ticksize.y = 4.5,
  angle = 45, # x axis angle
  axis.text.size = 7.5,
  axis.text.x.vjust = 1,
  plotarea_text_size = 2
)




# plot printin
pdf(file = str_c(out_figures,"Module_tilsig_bp_overlap_and_cocorrelation.pdf"),
    width = 6, height = 7.5)
print(
  ggarrange( p_agreement, p_correlation,
             ncol = 1,
             nrow = 2,
             labels = c("A", "B"),
             align = "hv",
             widths = 1,
             # heights = c(.61,.39),
             legend = "right",
             common.legend = F
  ),
  newpage = F
)
dev.off()


# # write_down tilsig_bp gene modules
#
# write_xlsx(
#   x = tilsig_bp %>%
#     purrr::map(
#       ~(.x %>%
#           dplyr::mutate(
#             Ncbi_gene_id = str_replace(
#               string = Ncbi_gene_id,
#               pattern = "ncbi_",
#               replacement = ""
#             )
#           )
#       )
#     ),
#   path = str_c(out_data, "tilsig_bp_modules.xlsx")
# )



# Summary: tilsig_bp geneset agreement and co-correlation plots
# The up regulated BPs were all related to immune signalling and are highly correlated
# - with low to high overlap in genes.
# Similar effect is observed with down regulated BPs, all related to ECM and cell adhesion to ECM.
# Hence all unregulated/down regulated BP genesets were pooled independently.


# Merge up/down regulated gene sets independently
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


tilsig_bp_merged <- list(
  Immune = bind_rows(tilsig_bp[c("APP", "Immune", "Immune_reg.",
                                 "Fc_epsilon", "IFN_gamma",
                                 "Immune_innate", "Complement1", "Complement2")]),
  ECM = bind_rows(tilsig_bp[c("ECM", "Adhesion")])
)
purrr::map(tilsig_bp_merged,nrow) %>% unlist()
# Immune    ECM
# 140     97


tilsig_bp_merged <- purrr::map(tilsig_bp_merged,
                               ~(.x %>% dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE)))
purrr::map(tilsig_bp_merged,nrow) %>% unlist()
# Immune    ECM
# 71     78

# ### end edit on 22 mar 2022



#
# ==============================================================================



# 4. Save Robjects
# ==============================================================================

save(tilsig, file = str_c(out_data,"tilsig.RData")) # Original TILsig no filtering
save(tilsig_clean, file = str_c(out_data,"tilsig_clean.RData")) # P_adj<.01
save(tilsig_bp, file = str_c(out_data, "tilsig_bp.RData"))
save(tilsig_bp_merged, file = str_c(out_data, "tilsig_bp_merged.RData"))

#
# ==============================================================================


# Clear memory
# ==============================================================================

rm(clin_finher, expr_finher,
   tilsig, tilsig_clean,
   david, david_clean,
   up_sig,up_terms,
   dn_sig, dn_terms,
   tilsig_bp, tilsig_bp_merged, tilsig_bp_legend)

#
# ==============================================================================
