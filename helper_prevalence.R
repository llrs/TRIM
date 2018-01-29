# Code for testing the prevalence with fisher tables
source("helper_functions.R")

#' Test prevalence
#'
#' Iterate in the tables to reach
#' @param presence Matrix with presence of certain microorganism
#' @param absence Matrix with the abscence of certain microorganism
#'
prevalence <- function(presence, absence) {
  stopifnot(all(rownames(presence) == rownames(absence)))
  sapply(rownames(presence), function(i) {
    m <- rbind(
      P = presence[i, ],
      A = absence[i, ]
    )
    f <- fisher.test(m, workspace = 2e8)
    # ch <- chisq.test(m)
    # c(f$p.value, ch$p.value)
    f$p.value
  })
}

#' Calculates the presence or absence of a microorganism
#'
#' It takes into account only samples with more than 0.5% of relative presence
#' of the microorganism
#' @param table Input data with the samples in columns and microorganism in rows
#' @param meta Metadata associated with those samples (it assumes they are in
#' the same order as in table)
#' @param columns Column to make constrasts from
#' @return a list with the matrices of presence and absence
prevalence_tab <- function(table, meta, columns) {
  stopifnot(all(colnames(table) == rownames(meta)))
  prevalence <- prop.table(as.matrix(table), 2) > 0.005
  subSets <- allComb(meta, columns)
  subSets[is.na(subSets)] <- FALSE
  totalSamples <- colSums(subSets, na.rm = TRUE)
  subSets <- subSets[, totalSamples >= 1]
  totalSamples <- totalSamples[totalSamples >= 1]
  presence <- prevalence %*% subSets
  absence <- matrix(
    totalSamples, nrow = nrow(presence),
    byrow = TRUE, ncol = ncol(presence)
  ) - presence
  list(presence = presence, absence = absence)
}


prevalence_2factors <- function(table, meta, columns) {
  stopifnot(all(colnames(table) == rownames(meta)))
  stopifnot(length(columns) == 2)
  res <- prevalence_tab(table, meta, columns)
  levels <- lapply(meta[, columns], function(x) {
    levels(as.factor(x))
  })
  aux <- function(y) {
    # To convert a long line to a matrix
    apply(
      y, 1, matrix, ncol = length(levels[[1]]), nrow = length(levels[[2]]),
      dimnames = levels, byrow = TRUE
    )
  }
  out <- lapply(res, aux)
  sapply(out, `names<-`, names(res$presence))
}

# Calculates the ratio of p-values of the prevalence between samples
ratio <- function(columns, data, indices, meta) {
  a <- t(data[indices, ]) # allows boot to select sample
  bindices <- !seq_len(nrow(data)) %in% indices
  b <- t(data[bindices, ])
  Ameta <- meta[indices, ]
  Bmeta <- meta[bindices, ]

  AsubSets <- allComb(Ameta, columns)
  if (!is.matrix(AsubSets)) {
    return(NA)
  }
  AtotalSamples <- colSums(AsubSets)
  AsubSets <- AsubSets[, AtotalSamples >= 1]
  AtotalSamples <- AtotalSamples[AtotalSamples >= 1]

  Apresence <- a %*% AsubSets
  AtotalSamplesm <- matrix(
    AtotalSamples, nrow = nrow(Apresence),
    byrow = TRUE, ncol = ncol(Apresence)
  )
  Aabsence <- AtotalSamplesm - Apresence

  BsubSets <- allComb(Bmeta, columns)
  if (!is.matrix(BsubSets)) {
    return(NA)
  }
  BtotalSamples <- colSums(BsubSets)
  BsubSets <- BsubSets[, BtotalSamples >= 1]
  BtotalSamples <- BtotalSamples[BtotalSamples >= 1]

  Bpresence <- b %*% BsubSets
  BtotalSamplesm <- matrix(
    BtotalSamples, nrow = nrow(Bpresence),
    byrow = TRUE, ncol = ncol(Bpresence)
  )
  Babsence <- BtotalSamplesm - Bpresence



  # Fisher test and ratio calculation
  sapply(rownames(presence), function(i) {
    Am <- rbind(
      P = Apresence[i, ],
      A = Aabsence[i, ]
    )
    Af <- fisher.test(Am, workspace = 2e8)

    Bm <- rbind(
      P = Bpresence[i, ],
      A = Babsence[i, ]
    )
    Bf <- fisher.test(Bm, workspace = 2e8)

    # ch <- chisq.test(m)
    # c(f$p.value, ch$p.value)
    Af$p.value / Bf$p.value
  })
}

#' Analyze the data for the relationship in time
#'
#' It looks in general and in colon and ILEUM specifics
#' @param input the table of microorganisms
#' @param meta the metadata it has encoded things only for intestinal
#' @return A list of the p-values the first two are for all comparisons
#' the following are for ileum and colon.
comp <- function(input, meta) {
  out <- vector("list", 4)

  removeControls <- meta$IBD == "CD"
  keepResponders <- meta$HSCT_responder %in% "YES"
  keepNonResponders <- meta$HSCT_responder %in% "NO"
  keepTime <- meta$Time %in% c("T0", "T26", "T52")

  res <- prevalence_tab(
    input[, keepTime & removeControls & keepResponders],
    meta[keepTime & removeControls & keepResponders, ],
    "Time"
  )
  # Fisher test
  Responders_Time <- prevalence(res$presence, res$absence)
  Responders_Time <- p.adjust(Responders_Time, "BH")

  out[[1]] <- Responders_Time

  res <- prevalence_tab(
    input[, keepTime & removeControls & keepNonResponders],
    meta[keepTime & removeControls & keepNonResponders, ],
    "Time"
  )
  # Fisher test
  NonResponders_Time <- prevalence(res$presence, res$absence)
  NonResponders_Time <- p.adjust(NonResponders_Time, "BH")

  out[[2]] <- NonResponders_Time

  ## Ileum ####
  ### Respondres ####
  keepILEUM <- meta$CD_Aftected_area %in% "ILEUM"
  keepCOLON <- meta$CD_Aftected_area %in% "COLON"

  res <- prevalence_tab(
    input[, keepTime & removeControls & keepResponders & keepILEUM],
    meta[keepTime & removeControls & keepResponders & keepILEUM, ],
    "Time"
  )
  # Fisher test
  Responders_Time <- prevalence(res$presence, res$absence)
  Responders_Time <- p.adjust(Responders_Time, "BH")

  ### Non Responders
  res <- prevalence_tab(
    input[, keepTime & removeControls & keepNonResponders & keepILEUM],
    meta[keepTime & removeControls & keepNonResponders & keepILEUM, ],
    "Time"
  )
  # Fisher test
  NonResponders_Time <- prevalence(res$presence, res$absence)
  NonResponders_Time <- p.adjust(NonResponders_Time, "BH")

  out[[3]] <- list(
    Responder = Responders_Time,
    NonResponders = NonResponders_Time
  )

  ## COLON ####
  ### Respondres ####
  keepILEUM <- meta$CD_Aftected_area %in% "ILEUM"
  keepCOLON <- meta$CD_Aftected_area %in% "COLON"

  res <- prevalence_tab(
    input[, keepTime & removeControls & keepResponders & keepCOLON],
    meta[keepTime & removeControls & keepResponders & keepCOLON, ],
    "Time"
  )
  # Fisher test
  Responders_Time <- prevalence(res$presence, res$absence)
  Responders_Time <- p.adjust(Responders_Time, "BH")

  ### Non Responders
  res <- prevalence_tab(
    input[, keepTime & removeControls & keepNonResponders & keepCOLON],
    meta[keepTime & removeControls & keepNonResponders & keepCOLON, ],
    "Time"
  )
  # Fisher test
  NonResponders_Time <- prevalence(res$presence, res$absence)
  NonResponders_Time <- p.adjust(NonResponders_Time, "BH")

  out[[4]] <- list(
    Responder = Responders_Time,
    NonResponders = NonResponders_Time
  )

  names(out) <- c("Responders", "NonResponders", "Ileum", "Colon")
  as.data.frame(out)
}

response_time <- function(table_org, meta) {
  # Test the prevalence between controls and non controls ####
  res <- prevalence_tab(table_org, meta, "IBD")
  IBD <- prevalence(res$presence, res$absence)

  IBD[is.na(IBD)] <- 1
  IBD <- p.adjust(IBD, "BH")

  out <- data.frame("ControlsVsIBD" = IBD)

  # Differences between responders and non responders at time 0 ####
  removeControls <- meta$IBD != "CONTROL"
  removeTimes <- meta$Time == "T0"
  keep <- removeControls & removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T0 <- prevalence(res$presence, res$absence)
  T0[is.na(T0)] <- 1
  T0 <- p.adjust(T0, "BH")

  out <- cbind(out, "T0_RvsNR" = T0)

  # Differences between responders and non responders at time 0 ####
  removeControls <- meta$IBD == "CONTROL"
  removeTimes <- meta$Time == "T0"
  keep <- removeControls | removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T0 <- prevalence(res$presence, res$absence)
  T0[is.na(T0)] <- 1
  T0 <- p.adjust(T0, "BH")

  out <- cbind(out, "T0vsC" = T0)



  # Differences between responders and non responders at time T26 ####
  removeControls <- meta$IBD != "CONTROL"
  removeTimes <- meta$Time == "T26"
  keep <- removeControls & removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T26 <- prevalence(res$presence, res$absence)
  T26[is.na(T26)] <- 1
  T26 <- p.adjust(T26, "BH")

  out <- cbind(out, "T26_RvsNR" = T26)

  removeControls <- meta$IBD == "CONTROL"
  removeTimes <- meta$Time == "T26"
  keep <- removeControls | removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T26 <- prevalence(res$presence, res$absence)
  T26[is.na(T26)] <- 1
  T26 <- p.adjust(T26, "BH")

  out <- cbind(out, "T26vsC" = T26)

  # Differences between responders and non responders at time T52 ####
  removeControls <- meta$IBD != "CONTROL"
  removeTimes <- meta$Time == "T52"
  keep <- removeControls & removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T52 <- prevalence(res$presence, res$absence)
  T52[is.na(T52)] <- 1
  T52 <- p.adjust(T52, "BH")

  out <- cbind(out, "T52_RvsNR" = T52)

  removeControls <- meta$IBD == "CONTROL"
  removeTimes <- meta$Time == "T52"
  keep <- removeControls | removeTimes
  res <- prevalence_tab(table_org[, keep], meta[keep, ], "HSCT_responder")

  # Fisher test
  T52 <- prevalence(res$presence, res$absence)
  T52[is.na(T52)] <- 1
  T52 <- p.adjust(T52, "BH")

  cbind(out, "T52vsC" = T52)
}

#' Look for prevalene in combinations of two
#' 
#' @param table is the data
#' @param meta is the metadata
#' @param columns is the columns to make comparisons from
#' @return A table with all the pairwise comparisons
comb_prevalence <- function(table, meta, columns) {
  above <- rowSums(prop.table(as.matrix(table), 2) > 0.005)
  
  if (sum(above > 0) == 0) {
    return(data.frame(0))
  }
  table <- table[above > 0, ]
  
  out <- prevalence_tab(table, meta, columns)
  func <- function(x) {
    res <- prevalence(out$presence[, x], out$absence[, x])
    res[is.na(res)] <- 1
    res <- p.adjust(res, "BH")
    res
  }
  o <- combn(colnames(out$presence), 2, FUN = func)
  colnames(o) <- combn(colnames(out$presence), 2, paste, collapse = "_&_")
  rownames(o) <- rownames(out$presence)
  o
}
