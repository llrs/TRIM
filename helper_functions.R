# Helper functions and variables

library("ggplot2")
library("RGCCA")
library("STATegRa")
library("reshape2")
library("ggforce")
library("MultiAssayExperiment")

#' Create a circle
#' @param center The position where the center of the circle is
#' @param diameter Arbitrary units of diameter of the circle
#' @param npoints Number of points of the circle (aka: definition of the circle)
#' @return A data.frame with the position of the points
circleFun <- function(center = c(-1, 1), diameter = 1, npoints = 100) {
  r <- diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circle <- circleFun(c(0, 0), 2, npoints = 100)

#' Clean and prepare the data from IMNGS
#'
#' Divides the taxonomy into a new matrix for each otu
#'
#' @param taxonomy Last column of files a string ; separated with domain,
#' phylum, vlass, order, family, genus and species.
#' @param otus The name of the rows
#'
#' @return
#' A matrix with the taxonomic information ready for the package phylo
taxonomy <- function(taxonomy, otus) {
  taxonomy <- sapply(taxonomy, strsplit, split = ";")
  names(taxonomy) <- otus
  otus_tax <- t(sapply(taxonomy, "[", seq(max(sapply(taxonomy, length)))))
  if (ncol(otus_tax) == 7) {
    colnames(otus_tax) <- c(
      "Domain", "Phylum", "Class", "Order",
      "Family", "Genus", "Species")
  } else if (ncol(otus_tax) == 6 ) {
    colnames(otus_tax) <- c(
      "Domain", "Phylum", "Class", "Order",
      "Family", "Genus") 
  }
  # Remove spaces
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "\\s", replacement = "")
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "[;:]", replacement = "")
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "^([a-z]__)", replacement = "")
  otus_tax[otus_tax == ""] <- NA # Remove empty cells
  otus_tax
}


# Check the taxonomy
# https://stackoverflow.com/q/7943695/2886003
#' Check if a vector is in the matrix
#' @param x The matrix to check if it is in the matrix
#' @param matrix The matrix where x is looked up.
#'
#' @return A logical vector saying if the vector is in the matrix
fastercheck <- function(x, matrix) {
  nc <- ncol(matrix)
  rec.check <- function(r, i, id) {
    id[id] <- matrix[id, i] %in% r[i]
    if (i < nc & any(id)) {
      rec.check(r, i + 1, id)
    } else {
      any(id)
    }
  }
  apply(x, 1, rec.check, 1, rep(TRUE, nrow(matrix)))
}

tol21rainbow <- c(
  "#771155", "#AA4488", "#CC99BB",
  "#114477", "#4477AA", "#77AADD",
  "#117777", "#44AAAA", "#77CCCC",
  "#117744", "#44AA77", "#88CCAA",
  "#777711", "#AAAA44", "#DDDD77",
  "#774411", "#AA7744", "#DDAA77",
  "#771122", "#AA4455", "#DD7788"
)

colors <- c(
  "#a692d2", "#6de14d", "#5a3bcb", "#b2e145", "#b844dd", "#50a93e",
  "#c640b6", "#6de697", "#472383", "#ddce3e", "#8766d6", "#8a9a37",
  "#517cd1", "#e24428", "#7dddcf", "#da425f", "#56ad7b", "#df4e98",
  "#c9e393", "#6a2d6d", "#de8d2f", "#3d3b76", "#d9ba73", "#d57ed2",
  "#3a682c", "#93295e", "#8f9c73", "#3c2032", "#d5d0bf", "#782e27",
  "#7ab4d5", "#ba5132", "#4d847e", "#de8f7d", "#2c3a29", "#d6abca",
  "#615325", "#ab6583", "#a47735", "#506280", "#89685f"
)

#' Calculates the angle between to slopes
#'
#' @param x,y Slope of the lines
#' @note The default compares the first slope with the slope 1
#' @return The smaller angle between the slopes in degrees
angle <- function(x, y = 1) {
  atan(abs((x - y) / (1 + y * x))) * 180 / pi
}

#' Calculate the distance between a line and a point
#'
#' The line is defined by points b and d.
#' @param p Point c(x, y)
#' @param b,d Points c(x, y) defining the line to calculate the distance with.
#' @return The units of distance between the point and the line
#' @note Change the d point to change the direction of the diagnoal
dist2d <- function(p, b = c(0, 0), d = c(1, 1)) {
  v1 <- b - d
  v2 <- p - b
  m <- cbind(v1, v2)
  abs(det(m)) / sqrt(sum(v1 * v1))
}

today <- format(Sys.time(), "%Y%m%d")

#' Compares the taxonomy of the otus
#'
#' Given two taxonomy tables find which one is in which one
#' @param taxa_1, taxa_2 taxonomic table as described in taxonomy
#' @return A matrix of nrow(taxa_1) * nrow(taxa_2)
contingency_taxa <- function(taxa_1, taxa_2) {

  # Check if the input is factor or character
  txc1 <- apply(taxa_1, 2, is.character)
  txc2 <- apply(taxa_2, 2, is.character)
  txf1 <- apply(taxa_1, 2, is.factor)
  txf2 <- apply(taxa_2, 2, is.factor)

  if (!all(txc1 | txf1) | !all(txc2 | txf2)) {
    stop("Taxa should be factors or characters (NAs allowed) ")
  }
  l <- sapply(rownames(taxa_2), function(x) {
    nas_1 <- sum(is.na(taxa_1[x, ]))
    nas_2 <- sum(is.na(taxa_2[x, ]))
    (colSums(apply(taxa_1, 1, `==`, taxa_2[x, ]), na.rm = TRUE) +
      min(nas_2, nas_1)) / 7
  })
  l
}

#' Z-score to p-value
#'
#' Calculates the p-value of a z-score
#' @param z the normalized value
#' @param one.sided Either NULL, - or +
#' @return The p-value
convert.z.score <- function(z, one.sided = NULL) {
  # https://www.biostars.org/p/17227/#136907
  if (!is.null(one.sided)) {
    if (!one.sided %in% c("+", "-")) {
      stop("one.sided should be NULL or + or -")
    }
  }

  if (!is.numeric(z)) {
    stop("z should be numeric")
  }

  if (is.null(one.sided)) {
    pval <- pnorm(-abs(z))
    pval <- 2 * pval
  } else if (one.sided == "-") {
    pval <- pnorm(z)
  } else if (one.sided == "+") {
    pval <- pnorm(-z)
  }
  return(pval)
}

#' Calculates the z-score of two correlations
#'
#' @param r1,r2 The correlation coefficients
#' @param n1,n2 The number of samples used to calculate the correlations
#' @return The z score of the comparison of the coefficients.
compare.correlations <- function(r1, r2, n1, n2) {
  r.norm <- function(r) {
    log(abs((1 + r) / (1 - r))) / 2
  }
  r1n <- r.norm(r1)
  r2n <- r.norm(r2)
  (r1n - r2n) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
}

#' Logical vectors of meta data
#'
#' Given the names of the columns of the data calculates the logical vectors of
#' each subset of the data
#' @param columns names of the columns to be used
#' @param data A \code{data.frame} with factors
#' @return a matrix with the logical values of each combination of the levels of
#' the columns given for the data in each column.
#' @note If some rows are all FALSE it means some values are NA.
allComb <- function(data, columns) {
  if (is.null(dim(data))) {
    stop("data should be a data.frame or a matrix")
  }

  if (!is.character(columns)) {
    stop("columns should be a character")
  }

  # if (!length(columns) >= 2) {
  #   stop("Several columns should be used")
  # }

  keep <- columns %in% colnames(data)
  if (sum(keep) == 0) {
    stop("Names of columns not present on data")
  } else if (sum(keep) != length(keep)) {
    warning("Columns:", paste(columns[!keep]), "are not present on data")
  }

  data <- data[, columns[keep]]

  if (is.null(dim(data))) {
    data <- as.factor(data)
    out <- sapply(levels(data), function(x) {
      data == x
    })
  } else {
    data <- sapply(data, as.factor)
    lvl <- sapply(as.data.frame(data), levels)
    if (is.matrix(lvl)) {
      lvl <- as.data.frame(lvl)
    }

    comb <- expand.grid(sapply(lvl, as.factor))

    comb2 <- apply(comb, 1, paste0, collapse = "_")
    out <- apply(comb, 1, function(x) {
      # Repeat the terms as much as the data
      combT <- sapply(x, function(y) {
        rep(y, nrow(data))
      })
      # Compare the data with the levels
      check <- rowSums(data == combT)
      # Convert to logical by performing an &
      o <- check == 2
      o[is.na(o)] <- FALSE
      o
    })
    colnames(out) <- comb2
  }
  out
}

#' Make a rectangle at those cells
#'
#' @param tfMat matrix to take the dimensions from and with TRUE at the cells
#' to be higlighted
#' @param border Type of border
#' @return A side effect of the adding a rectangle to the heatmap
#' @references
#' https://stackoverflow.com/a/45891116/2886003
#' https://stackoverflow.com/a/7980349/2886003
makeRects <- function(tfMat, border) {
  nx <- nrow(tfMat)
  ny <- ncol(tfMat)
  cAbove <- expand.grid(seq_len(nx), seq_len(ny))[tfMat, ]
  cA_1 <- table(cAbove[, 1])

  # In case all the positions are consecutive
  # draw a big rectangle by keeping the x positions
  if (all(names(cA_1) %in% seq_len(max(cA_1)))) {
    xl <- rep(min(as.numeric(names(cA_1))) - 0.49, max(cA_1) * length(cA_1))
    xr <- rep(max(as.numeric(names(cA_1))) + 0.49, max(cA_1) * length(cA_1))
  } else {
    xl <- cAbove[, 1] - 0.49
    xr <- cAbove[, 1] + 0.49
  }
  # In case all the positions are consecutive
  # draw a big rectangle by keeping the y positions
  cA_2 <- table(cAbove[, 2])
  if (all(names(cA_2) %in% seq_len(max(cA_2)))) {
    yb <- rep(min(as.numeric(names(cA_2))) - 0.49, max(cA_2) * length(cA_2))
    yt <- rep(max(as.numeric(names(cA_2))) + 0.49, max(cA_2) * length(cA_2))
  } else {
    yb <- cAbove[, 2] - 0.49
    yt <- cAbove[, 2] + 0.49
  }
  rect(xl, yb, xr, yt, border = border, lwd = 3)
}

#' Select variable from bootstrapping
#'
#' @param x List of the summary statsitics of the bootstrapping
#' @return The names of the selected variables
selectVar <- function(x) {
  varNames <- rownames(x)
  x <- x[, "freq"]
  names(x) <- varNames
  if (length(unique(x)) == 1) {
    names(x)
  } else {
    q <- quantile(x, na.rm = TRUE)
    names(x)[x > q["75%"] & !is.na(x)]
  }
}


# Select the significant relationships ####
#' Two sided test
#'
#' Test in a vector from a permutation if there is a relationship or not.
#' Assumes that the distribution is symmetric around 0.
#' @param z Vector of of the permutations
#' @param y Value of the test
#' @return The p-value
two.sided <- function(y, z) {
  stopifnot(length(y) == 1)
  stopifnot(sum(!is.na(z)) > 2)
  greater <- sum(abs(z) >= abs(y), na.rm = TRUE)

  (1 + greater) / (1 + sum(!is.na(z)))
}


#' Normalize the metadata of the RNA
#'
#' Performs some modifications to the data.frame
#' @param meta Meta
#' @return A data.frame with some modifications
meta_r_norm <- function(meta) {

  # Duplicated label: we don't know which one is from where!!
  meta$CD_Aftected_area[meta$Sample_Code == "22_T52_T_DM_III"] <- NA
  meta$Exact_location[meta$Sample_Code == "22_T52_T_DM_III"] <- NA
  
  meta$`Sample Name_RNA` <- toupper(as.character(meta$`Sample Name_RNA`))
  meta$Sample_Code_uDNA <- toupper(as.character(meta$Sample_Code_uDNA))

  # Add people IDs (some people has several Patient ID)
  # This affects transplant too
  meta$ID <- meta$Patient_ID
  meta$ID[meta$Patient_ID %in% c("15", "23")] <- "15/23"
  meta$ID[meta$Patient_ID %in% c("33", "36")] <- "33/36"
  meta$ID[meta$Patient_ID %in% c("29", "35")] <- "29/35"
  meta$ID[meta$Patient_ID %in% c("17", "40")] <- "17/40"
  meta$ID <- as.factor(meta$ID)

  # We don't know yet if the newest samples are responders or not (yet)
  meta$HSCT_responder[meta$ID %in% c("38", "40", "41")] <- NA
  meta$HSCT_responder[meta$IBD == "CONTROL"] <- NA

  # Pre transplantament post and baseline three phases of the treatment
  meta$Transplant <- "Post"
  meta$Transplant[meta$Patient_ID %in% c("15", "33", "29")] <- "Pre"
  meta$Transplant[meta$Time %in% c("T0", "S0")] <- "Baseline"
  meta$Transplant[meta$Time %in% c("C")] <- NA

  meta$Active_area[meta$Involved_Healthy == "HEALTHY"] <- "HEALTHY"

  meta$Treatment <- "NO"
  meta$Treatment[grep("W TMT| W Surgery", meta$Endoscopic_Activity)] <- "YES"
  meta$Treatment[meta$Treatment == "C"] <- NA

  # Surgery post transplant
  meta$Surgery[!is.na(meta$Treatment)] <- "NO"
  meta$Surgery[grep("Surgery", meta$Endoscopic_Activity)] <- "YES"

  # Correct clinical scores (mainly for controls)
  meta$SESCD_local[meta$IBD == "CONTROL"] <- 0
  meta$SESCD_global[meta$IBD == "CONTROL"] <- 0
  meta$CDAI[meta$IBD == "CONTROL"] <- 0
  meta$DiagDate[meta$IBD == "CONTROL"] <- 0
  meta$Transplant[meta$IBD == "CONTROL"] <- "Baseline"
  
  # Add the date of the diagnosis
  strDates <- c(
    "13" = "15-feb-02", "14" = "15-jul-97", "15/23" = "15-jul-06",
    "16" = "15-oct-06", "17/40" = "15-nov-02", "18" = "15-jul-99",
    "19" = "15-jul-99", "20" = "19-jun-09", "21" = "15-jun-06",
    "22" = "15-ago-09", "25" = "15-feb-99", "26" = "15-mar-00",
    "27" = "15-abr-07", "28" = "15-abr-01", "29/35" = "15-jun-04",
    "30" = "15-jun-89", "31" = "15-jun-03", "33/36" = "15-jun-05",
    "34" = "15-jun-03", "37" = "15-jun-06", "38" = "15-feb-12",
    "39" = "31-ago-11", "41" = "15-Jun-03"
  )
  dates <- as.Date(strDates, "%d-%b-%y")
  names(dates) <- names(strDates)
  o <- match(meta$ID, names(dates))
  meta <- cbind(meta, "DiagDate" = dates[o])
  diagTime <- as.Date(meta$DATE_SAMPLE, "%d/%m/%Y") - dates[o]
  diagTime <- as.numeric(diagTime / 365)
  diagTime[is.na(diagTime)] <- 0 # If no diagnosis (controls) set to 0
  meta <- cbind(meta, diagTime)
  AgeDiag <- as.numeric(dates[o] - 
                as.Date(meta$Birth_date, "%d/%m/%Y"))/365.25
  meta <- cbind(meta, "AgeDiag" = AgeDiag)

  return(meta)
}

#' Normalize 16S biopsies metadata
#'
#' @param meta The metadata
#' @return The dataframe with the metadata corrected
meta_i_norm <- function(meta) {
  meta <- meta[, apply(meta, 2, function(x) {
    length(unique(x)) != 1
  })]
  meta$ID <- meta$Patient_ID
  meta$ID[meta$Patient_ID %in% c("15", "23")] <- "15/23"
  meta$ID[meta$Patient_ID %in% c("33", "36")] <- "33/36"
  meta$ID[meta$Patient_ID %in% c("29", "35")] <- "29/35"
  meta$ID[meta$Patient_ID %in% c("17", "40")] <- "17/40"
  meta$ID <- as.factor(meta$ID)

  # There is a mislabeling on those tubes, we don't know which is which
  meta$CD_Aftected_area[meta$Sample_Code == "22_T52_T_DM_III"] <- NA

  # We don't know yet if the newest samples are responders or not
  meta$HSCT_responder[(meta$ID %in% c("38", "40", "41"))] <- NA

  return(meta)
}

#' Filter expressions
#'
#' @param expr Input RNAseq data (not microbiota)
#' @return A matrix without low expressed genes and with low variance
norm_RNAseq <- function(expr) {
  # Remove low expressed genes
  expr <- expr[rowSums(expr != 0) >= (0.25 * ncol(expr)), ]
  expr <- expr[rowMeans(expr) > quantile(rowMeans(expr), prob = 0.1), ]

  # Filter genes by variance
  SD <- apply(expr, 1, sd)
  CV <- sqrt(exp(SD ^ 2) - 1)
  expr[CV > quantile(CV, probs = 0.1), ]
}

#' Normalize 16S stools metadata
#'
#' @param meta The metadata
#' @return The dataframe with the metadata corrected
meta_s_norm <- function(meta) {
  # Clean the metadata
  meta <- meta[, apply(meta, 2, function(x) {
    length(unique(x)) != 1
  })]
  meta$ID <- meta$Patient_ID
  meta$ID[meta$Patient_ID %in% c("15", "23")] <- "15/23"
  meta$ID[meta$Patient_ID %in% c("33", "36")] <- "33/36"
  meta$ID[meta$Patient_ID %in% c("29", "35")] <- "29/35"
  meta$ID <- as.factor(meta$ID)
  meta
}
#' Filters a matrix by abundance
#'
#' @param x A matrix of the microbiota (absolute values)
#' @param perc A numeric value between 0 and 100 of a percentage
#' @return
filter_abundance <- function(x, perc = 0.5) {

}
