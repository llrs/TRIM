cd <- setwd("..")

# Load the helper file
source("helper_functions.R")
# Read files
otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)
setwd(cd)


# keep <- !grepl("28_T52_T_DM_CH", meta$Sample_Code) # Remove outlier See PCA
# meta <- meta[keep, ]
# otus_s <- otus_s[keep, ]
# otus_i <- otus_i[keep, ]

##### STATegRa #####
keepCol <- allComb(meta, c("CD_Aftected_area", "Time"))
keepCol <- cbind(keepCol, allComb(meta, c("CD_Aftected_area")))
# keepColon <- meta$CD_Aftected_area == "COLON"
# keepT0 <- meta$Time == "T106"
# keep <- keepColon & keepT0
pdf(paste0("Figures/", today, "_STATegRa_plots.pdf"))
library("gridExtra")
for (i in seq_len(ncol(keepCol))){
  keep <- keepCol[, i]
  message("\n", colnames(keepCol)[i])
  # Create ExpressionSet objects
  expr <- as.matrix(t(otus_i))[, keep]
  colnames(expr) <- rownames(meta[keep, ])
  
  eS_i <- createOmicsExpressionSet(expr, pData = meta[keep, ])
  
  expr <- as.matrix(t(otus_s))[, keep]
  colnames(expr) <- rownames(meta[keep, ])
  
  eS_s <- createOmicsExpressionSet(expr, pData = meta[keep, ])
  
  stopifnot(ncol(assay(eS_s)) == ncol(assay(eS_i)))
  
  # Selecting components
  tryCatch({
    cc <- selectCommonComps(t(otus_i[keep, ]), t(otus_s[keep, ]), Rmax = 3)
  PCA.selection(t(otus_i[keep, ]), fac.sel = "single%", varthreshold = 0.03)$numComps
  PCA.selection(t(otus_s[keep, ]), fac.sel = "single%", varthreshold = 0.03)$numComps
  (ms <- modelSelection(list(eS_i, eS_s), Rmax = 7, fac.sel = "single%",
                        varthreshold = 0.03))
  }, error = function(e){message(e)})

  if (ms$common != 0) {
    # grid.arrange(cc$pssq, cc$pratios, ncol=2)
    # Omics Integration
    discoRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
                                  Names = c("Intestinal", "Stools"),
                                  method = "DISCOSCA",
                                  Rcommon = ms$common,
                                  Rspecific = ms$dist,
                                  center = TRUE, scale = TRUE)
    if (is.null(unlist(discoRes@VAF$dist)) | 
        length(discoRes@VAF$dist$Block1) != length(1:discoRes@distComps[1])){
      warning("Malfunctioning plotVAF")
    } else {
      tryCatch({plotVAF(discoRes, main = colnames(keepCol)[i])}, 
               error = function(e){
                 print("Error:")
                 message(e)
               },
               warning = function(w){
                 print("Warning:")
                 message(w)
               })
    }
  
  # jiveRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
  #                              Names = c("Intestinal", "Stools"),
  #                              method = "JIVE",
  #                              Rcommon = ms$common,
  #                              Rspecific = ms$dist,
  #                              center=TRUE, scale=TRUE)
  # o2plsRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
  #                               Names = c("Intestinal", "Stools"),
  #                               method="O2PLS", 
  #                               Rcommon = ms$common,
  #                               Rspecific = ms$dist,
  #                               center=TRUE, scale=TRUE, weight=TRUE)
  
    tryCatch({
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "common",
                    combined = TRUE, block = "1", color = "HSCT_responder"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "1", color = "HSCT_responder"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "2", color = "HSCT_responder"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "common",
                    combined = TRUE, block = "1", color = "Endoscopic_Activity"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "1", color = "Endoscopic_Activity"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "2", color = "Endoscopic_Activity"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "common",
                    combined = TRUE, block = "1", color = "Involved_Healthy"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "1", color = "Involved_Healthy"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "scores", type = "both",
              combined = TRUE, block = "2", color = "Involved_Healthy"))
      },
      error = function(e){message(e)})
    
    tryCatch({
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "1", color = "HSCT_responder") + ggtitle(subtitle = colnames(keepCol)[i]))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "2", color = "HSCT_responder"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "1", color = "Endoscopic_Activity"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "2", color = "Endoscopic_Activity"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "1", color = "Involved_Healthy"))
      print(plotRes(object = discoRes, comps = c(1, 1), what = "loadings", type = "both",
              combined = TRUE, block = "2", color = "Involved_Healthy"))},
      error = function(e){print("Error:")
        message(e)}, 
      warning = function(w){
        print("Warning:")
        message(w)})
  } else {
    message("No common components")
  }
}

