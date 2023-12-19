# finher_tilsig.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# Generate TILsig from log2 TIL counts from TN, HER2 and TN+HER2
# Compare TILsig from TN, HER2 and TN+HER2 w.r.t gene set, and GO BP enrichment
# Decipher TILsig based on asociated biological processes !!! (Deciphering will
# remove any genes accociated with confounding factors involved
# in the TIL's semi quantitative measurement process) !!!!!



# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Load and format clinical data.
# 2. Generate TILsig and compare TILsig based on geneset and enriched GO BP terms
# 3. Decipher TILsig based on asociated biological processes
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

# purrr::map(tilsig_clean,dim)
# $TN
# [1] 484   8
# $HER2
# [1] 171   8
# $ALL
# [1] 994   8



# gene set agreement
# >>>>>>>>>>>>>>>>>>

sig_agreement <- function(x,y){
  xx = x %>% inner_join(y %>% dplyr::select(Ncbi_gene_id1, Direction), by = "Ncbi_gene_id1")
  gper = (nrow(xx)/nrow(x)) * 100
  dper = ((purrr::map2_lgl(xx$Direction.x, xx$Direction.y,~(.x == .y)) %>% sum()) / nrow(xx)) * 100

  str_c("% of x(n=", nrow(x),") in y(n=", nrow(y),"): ", round(gper), " (n=",nrow(xx),"); ",
        "% common genes with same direction: ", round(dper))
}

sig_agreement(x = tilsig_clean$TN, y = tilsig_clean$ALL)
# [1] "% of x(n=484) in y(n=994): 95 (n=460); % common genes with same direction: 100"
sig_agreement(x = tilsig_clean$HER2, y = tilsig_clean$ALL)
# [1] "% of x(n=171) in y(n=994): 98 (n=168); % common genes with same direction: 100"
sig_agreement(x = tilsig_clean$HER2, y = tilsig_clean$TN)
# [1] "% of x(n=171) in y(n=484): 78 (n=134); % common genes with same direction: 100"
sig_agreement(x = bind_rows(tilsig_clean$TN, tilsig_clean$HER2), y = tilsig_clean$ALL)
# [1] "% of x(n=655) in y(n=994): 96 (n=628); % common genes with same direction: 100"

# There is >78% agreement between different versions of TILsig.
# Further, among common genes the direction is identical.


# biological agreement
# >>>>>>>>>>>>>>>>>>>>

print_sig <- function(x, label = "sig"){
  up <- x %>% dplyr::filter(Direction == 1)
  dn <- x %>% dplyr::filter(Direction == -1)
  write.table(x = up %>% dplyr::select(Ncbi_gene_id2),
              sep = "\t",
              row.names = F,
              col.names = F,
              file = str_c("results/tables/",label,"_up.txt"))
  write.table(x = dn %>% dplyr::select(Ncbi_gene_id2),
              sep = "\t",
              row.names = F,
              col.names = F,
              file = str_c("results/tables/",label,"_dn.txt"))
  write.table(x = bind_rows(up,dn) %>% dplyr::select(Ncbi_gene_id2),
              sep = "\t",
              row.names = F,
              col.names = F,
              file = str_c("results/tables/",label,".txt"))
}


print_sig(x = tilsig_clean$TN, label = "tilsig_tn")
print_sig(x = tilsig_clean$HER2, label = "tilsig_her2")
print_sig(x = tilsig_clean$ALL, label = "tilsig_all")


# background for david analysis
x = expr_finher %>%
  dplyr::select(Ncbi_gene_id) %>%
  dplyr::filter(str_detect(Ncbi_gene_id,"///",negate = T)) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id,"ncbi_",""))

write.table(x = x,
            sep = "\t",
            row.names = F,
            col.names = T,
            file = str_c("results/tables/tilsig_background_genes.txt"))


# load david results
x <- list.files("results/tables/david_results/")
david <- vector(mode = "list", length = length(x))
names(david) = x

for(i in names(david)){
  david[[i]] <- read_delim(file = str_c("results/tables/david_results/", i),delim = "\t")
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

write_xlsx(x,path = "results/tables/david_results_fdr.05.xlsx")

# Similar bioloical proceses across up and down regulated gene independently

#
# ==============================================================================



# 3. Decipher TILsig based on asociated biological processes
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

# ALL version of siglist contins all genes fropm her2 and tn
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
# ALL version of siglist contins all genes fropm her2 and tn (one gene from tn is not present)
# Hence consider siglist from ALL version

# Further, due to similarity in genesets and biological signals betwern different
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
# names(tilsig_bp)
# [1] "APP"           "Immune"        "Immune_reg."   "Fc_epsilon"
# [5] "IFN_gamma"     "Immune_innate" "Complement1"   "Complement2"
# [9] "ECM"           "Adhesion"




# ### added on 22 mar 2022



# tilsig_bp geneset agreement and co-correlation plots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# # Functions
#
# module_gene_overlap <- function(
#
#   lst # list of modules for which gene overlap is to be measured
#
#   ){
#
#   # lst = tilsig_bp
#
#   # Overlap data frame
#   odf <- matrix(data = NA, nrow = length(lst), ncol = length(lst)) %>%
#     as.data.frame()
#   rownames(odf) = names(lst)
#   colnames(odf) = names(lst)
#
#   for(irow in names(lst)){
#     for(icol in names(lst)){
#       odf[irow, icol] =
#         (intersect(lst[[irow]]$Ncbi_gene_id, lst[[icol]]$Ncbi_gene_id) %>% length()) /
#         (nrow(lst[[irow]]))
#       # odf[irow, icol] contains proportion of "irow" sig genes present in "icol" sig
#     }
#   }
#
#   odf <- odf %>%
#     tibble::as_tibble(rownames = "Module_name") %>%
#     dplyr::mutate(Module_name2 = str_c(names(lst)," (",purrr::map_int(lst, nrow),")") ) %>%
#     dplyr::select(c(Module_name, Module_name2), !c(Module_name, Module_name2))
#
#   return(odf)
#
#   # odf interpretation:
#   # each row contains proportion of row signature genes present in column signature
# }
#
#
# module_gene_overlap_plot <- function(
#
#   module_overlap, # output from module_gene_overlap()
#   overlap_col,    # char vector of heatmap color in the order of low, mid, high
#   module_label,   # alternative names to use in axis annotation
#   module_group,   # char vector of module group name, preserve module order as of module_overlap df
#   module_group_col, # char vector of module group color, preserve order
#   module_group_vline = NULL, # xintercept for vlines to group heatmap
#   module_group_hline = NULL, # yintercept for hlines to group heatmap
#   line_col = "gray30",
#   ticksize.x = 5,
#   ticksize.y = 4.5,
#   angle = 45, # x axis angle
#   plotarea_text_size = 2
#
# ){
#
#   # module_overlap = module_overlap
#   # overlap_col = c(low = "gold", mid = "white", high = "blue")
#   # module_label = module_annot$Module_name2
#   # module_group = module_annot$Module_group
#   # module_group_col = module_annot$Module_group_col
#   # module_group_vline = c(0.5, 8.5,10.5)
#   # module_group_hline = c(0.5, 8.5,10.5)
#   # line_col = "gray30"
#   # ticksize.x = 5
#   # ticksize.y = 4.5
#   # angle = 45
#   # plotarea_text_size = 2
#
#
#   if(!identical(module_overlap$Module_name2, module_label)){
#     stop("module_label and module_overlap order not identical.")
#   }
#
#
#   ggdf <- module_overlap %>%
#     tidyr::gather(key = "Column", value = "Proportion", -Module_name, -Module_name2) %>%
#     dplyr::rename(Row = "Module_name") %>%
#     dplyr::mutate(Proportion = round(Proportion, digits = 1),
#                   Row = factor(Row, levels = unique(Row)),
#                   Column = factor(Column, levels = levels(Row)))
#
#
#   p <- ggplot(data = ggdf, aes(x = Column, y = Row)) +
#     geom_raster(aes(fill = Proportion)) +
#     scale_fill_gradient2(low = overlap_col[1],
#                          mid = overlap_col[2],
#                          high = overlap_col[3],
#                          midpoint = 0.5)+
#     geom_text(aes(label = Proportion), size = plotarea_text_size) +
#     # Hack to include module group names
#     geom_text(data = tibble(Module_group = module_group),
#               aes(x = Inf, y = Inf, label = "", color = Module_group))+
#     scale_color_manual(values = unique(module_group_col),
#                        labels = unique(module_group))+
#     scale_x_discrete(expand = c(0, 0), labels = unique(ggdf$Module_name2)) + # expand() remove extra space
#     scale_y_discrete(expand = c(0, 0), labels = unique(ggdf$Module_name2)) + # expand() remove extra space
#     # End hack
#     guides(fill = guide_colorbar(title.position = "top"),
#            color = ggh4x::guide_stringlegend(title.position = "top",
#                                       nrow = 4,
#                                       face = "bold")) +
#     theme(axis.text.x.bottom = element_text(angle = angle, hjust = 1, vjust = 1),
#           legend.position = "bottom",
#           axis.text = element_text(colour = module_group_col),
#           axis.ticks = element_line(colour = module_group_col),
#           axis.ticks.x = element_line(size = ticksize.x),
#           axis.ticks.y = element_line(size = ticksize.y)) +
#     labs(x = "Reference modules", y = "Test modules",
#          title = "Proportion of test modules present in reference module") +
#
#     geom_hline(yintercept = module_group_hline, color = line_col) +
#     geom_vline(xintercept = module_group_vline, color = line_col)
#
#   return(p)
# }
#
#
#
# module_correlation_plot <- function(
#   module_score = score, # module score output from get_module_score_2()
#   score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan"),
#   module_group = module_annot$Module_group, # char vector of module group name, preserve module order as of module score df
#   module_group_col = module_annot$Module_group_col, # char vector of module group color, preserve order
#   module_group_vline = c(0.5, 8.5,10.5), # xintercept for vlines to group heatmap
#   module_group_hline = c(0.5, 8.5,10.5), # yintercept for hlines to group heatmap
#   line_col = "gray30",
#   ticksize.x = 5,
#   ticksize.y = 4.5,
#   angle = 45, # x axis angle
#   plotarea_text_size = 2
# ){
#
#   # module_score = score
#   # score_col = c(low = "darkgoldenrod", mid = "white", high = "darkcyan")
#   # module_group = module_annot$Module_group
#   # module_group_col = module_annot$Module_group_col
#   # module_group_vline = c(0.5, 8.5,10.5)
#   # module_group_hline = c(0.5, 8.5,10.5)
#   # line_col = "gray30"
#   # ticksize.x = 5
#   # ticksize.y = 4.5
#   # angle = 45
#   # plotarea_text_size = 2
#
#
#
#   # Co-correlaion
#   ggdf <- module_score %>%
#     dplyr::select(-1) %>%
#     cor(method = "pearson") %>%
#     as_tibble(rownames = "Row") %>%
#     tidyr::gather(key = "Column", value = "Pearson", -1) %>%
#     dplyr::mutate(Pearson = round(Pearson, digits = 1),
#                   Row = factor(Row, levels = unique(Row)),
#                   Column = factor(Column, levels = levels(Row)))
#
#
#   p <- ggplot(data = ggdf, aes(x = Column, y = Row)) +
#     geom_raster(aes(fill = Pearson)) +
#     scale_fill_gradient2(low = score_col[1],
#                          mid = score_col[2],
#                          high = score_col[3],
#                          midpoint = 0,
#                          limits = c(-1,1))+
#     geom_text(aes(label = Pearson), size = plotarea_text_size) +
#
#     scale_x_discrete(expand = c(0,0))+
#     scale_y_discrete(expand = c(0,0))+
#
#     # ! Hack to set legend for module class
#     geom_text(data = tibble(Module_group = module_group),
#               aes(x = Inf, y= -Inf, label = "", color = Module_group)) +
#     scale_color_manual(values = unique(module_group_col),
#                        labels = unique(module_group)) +
#     # ! end hack
#
#     guides(fill = guide_colorbar(title.position = "top"),
#            color = guide_stringlegend(title = "Module group",
#                                       face = "bold",
#                                       nrow = 4,
#                                       title.position = "top") ) +
#
#     theme(axis.text.x.bottom = element_text(angle = angle, hjust = 1, vjust = 1),
#           legend.position = "bottom",
#           axis.text = element_text(colour = module_group_col),
#           axis.ticks = element_line(colour = module_group_col),
#           axis.ticks.x = element_line(size = ticksize.x),
#           axis.ticks.y = element_line(size = ticksize.y),
#           panel.grid = element_blank()) + #axis.title = element_blank()
#     labs(x = "Gene modules", y = "Gene modules",
#          title = "Pearson correlation between gene modules") +
#
#
#     geom_hline(yintercept = module_group_hline, color = line_col) +
#     geom_vline(xintercept = module_group_vline, color = line_col)
#
#
#   return(p)
#
# }



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


# write_down tilsig_bp gene modules

write_xlsx(
  x = tilsig_bp %>%
    purrr::map(
      ~(.x %>%
          dplyr::mutate(
            Ncbi_gene_id = str_replace(
              string = Ncbi_gene_id,
              pattern = "ncbi_",
              replacement = ""
            )
          )
      )
    ),
  path = "results/tables/tilsig_bp_modules.xlsx"
)



# Summary: tilsig_bp geneset agreement and co-correlation plots
# The up regulated BPs were all related to immune signalling and are highly correlated
# - with low to high overlap in genes.
# Similar effect is oberved with down regulated BPs, all related to ECM and cell adhision to ECM.
# Hence all upregulated/downregulated BP genesets were pooled independently.


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


# Commented on 22 mar 2022
# The below code become obsolete, as in the previous code its is decided to pool
# all upreguled/downregulated BP genesets togother independetly.
# !!!!!! Note that the R object "tilsig_bp_merged" now contains TILsig, Immune and ECM signatures

#
# # tilsig_bp gene agrremnt martix
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# x <- matrix(data=NA, ncol = length(tilsig_bp), nrow = length(tilsig_bp))
# colnames(x) <- names(tilsig_bp)
# rownames(x) <- names(tilsig_bp)
# for(row in rownames(x)){
#   for(col in colnames(x)){
#   x[row, col] = sum(tilsig_bp[[row]]$Ncbi_gene_id %in% tilsig_bp[[col]]$Ncbi_gene_id) / nrow(tilsig_bp[[row]])
#   x[row, col] = round(x[row, col], digits = 2)
#   }
# }
# x
# #               APP Immune Immune_reg. Fc_epsilon IFN_gamma Immune_innate Complement1 Complement2  ECM Adhesion
# # APP           1.00   0.23        0.23       0.68      0.23          0.23        0.00        0.00 0.00      0.0
# # Immune        0.17   1.00        0.31       0.14      0.28          0.17        0.17        0.17 0.00      0.0
# # Immune_reg.   0.33   0.60        1.00       0.27      0.53          0.27        0.27        0.27 0.00      0.0
# # Fc_epsilon    0.71   0.19        0.19       1.00      0.00          0.10        0.19        0.19 0.00      0.0
# # IFN_gamma     0.36   0.57        0.57       0.00      1.00          0.21        0.00        0.00 0.00      0.0
# # Immune_innate 0.23   0.23        0.18       0.09      0.14          1.00        0.27        0.23 0.00      0.0
# # Complement1   0.00   0.56        0.44       0.44      0.00          0.67        1.00        0.89 0.00      0.0
# # Complement2   0.00   0.62        0.50       0.50      0.00          0.62        1.00        1.00 0.00      0.0
# # ECM           0.00   0.00        0.00       0.00      0.00          0.00        0.00        0.00 1.00      0.4
# # Adhesion      0.00   0.00        0.00       0.00      0.00          0.00        0.00        0.00 0.38      1.0
#
# purrr::map(tilsig_bp, nrow) %>% unlist()
# # APP        Immune   Immune_reg.    Fc_epsilon
# # 22            29            15            21
# # IFN_gamma Immune_innate   Complement1   Complement2
# # 14            22             9             8
# # ECM      Adhesion
# # 47            50
#
#
# # Merge genesets with > 60% agreement
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# # APP-Fc # Fc receptor binding process part of APP, https://www.frontiersin.org/articles/10.3389/fimmu.2020.01393/full , https://pubmed.ncbi.nlm.nih.gov/10631953/ , https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4145246/
# # Imm-Imm.reg
# # IFN-gamma
# # Innate-Complement1-2 # complement part of innate, https://pubmed.ncbi.nlm.nih.gov/16234578/
# # ECM
# # Adhesion
#
# tilsig_bp_merged <- list(
#   APP_Fc = bind_rows(tilsig_bp[c("APP","Fc_epsilon")]),
#   Immune = bind_rows(tilsig_bp[c("Immune", "Immune_reg.")]),
#   IFN_gamma = tilsig_bp$IFN_gamma,
#   Innate = bind_rows(tilsig_bp[c("Immune_innate", "Complement1", "Complement2")]),
#   ECM = tilsig_bp$ECM,
#   Adhesion = tilsig_bp$Adhesion
# )
# purrr::map(tilsig_bp_merged,nrow) %>% unlist()
# # APP_Fc    Immune IFN_gamma    Innate       ECM  Adhesion
# # 43        44        14        39        47        50
#
#
# tilsig_bp_merged <- purrr::map(tilsig_bp_merged,~(.x %>% dplyr::distinct(Ncbi_gene_id)))
# purrr::map(tilsig_bp_merged,nrow) %>% unlist()
# # APP_Fc    Immune IFN_gamma    Innate       ECM  Adhesion
# # 28        35        14        25        47        50
#
#
#
# # tilsig_bp_merged gene agrremnt martix
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# x <- matrix(data=NA, ncol = length(tilsig_bp_merged), nrow = length(tilsig_bp_merged))
# colnames(x) <- names(tilsig_bp_merged)
# rownames(x) <- names(tilsig_bp_merged)
# for(row in rownames(x)){
#   for(col in colnames(x)){
#     x[row, col] = sum(tilsig_bp_merged[[row]]$Ncbi_gene_id %in% tilsig_bp_merged[[col]]$Ncbi_gene_id) / nrow(tilsig_bp_merged[[row]])
#     x[row, col] = round(x[row, col], digits = 2)
#   }
# }
# x
# #           APP_Fc Immune IFN_gamma Innate  ECM Adhesion
# # APP_Fc      1.00   0.32      0.18   0.36 0.00      0.0
# # Immune      0.26   1.00      0.31   0.23 0.00      0.0
# # IFN_gamma   0.36   0.79      1.00   0.21 0.00      0.0
# # Innate      0.40   0.32      0.12   1.00 0.00      0.0
# # ECM         0.00   0.00      0.00   0.00 1.00      0.4
# # Adhesion    0.00   0.00      0.00   0.00 0.38      1.0
#
#
# purrr::map(tilsig_bp_merged, nrow) %>% unlist()
# # APP_Fc    Immune IFN_gamma    Innate    ECM  Adhesion
# # 28        35        14        25        47        50


#
# ==============================================================================



# 3. Save Robjects
# ==============================================================================

save(tilsig, file = str_c(out_data,"tilsig.RData")) # Original TILsig no filtering
save(tilsig_clean, file = str_c(out_data,"tilsig_clean.RData")) # P_adj<.01
save(tilsig_bp, file = "results/data/tilsig_bp.RData")
save(tilsig_bp_merged, file = "results/data/tilsig_bp_merged.RData")

#
# ==============================================================================



