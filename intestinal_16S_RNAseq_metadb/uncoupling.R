
# Prepare data #####
otus <- read.delim("../uncoupling_data/SupplementalTableS005.txt", 
                   check.names = FALSE, row.names = 1)
transcripts <- read.delim("../uncoupling_data/SupplementaryTableS002.txt", 
                   check.names = FALSE, row.names = 1)
# From the body of the second email from Rob Haesler
metadata <- read.table(text = "C1675   CD_SIG_b        F       39
C1492   CD_SIG_b        F       39
D0742   CD_SIG_b        F       32
D0743   CD_SIG_b        M       41
C1670   CD_SIG_nb       F       47
C1673   CD_SIG_nb       M       35
D0745   CD_SIG_nb       F       35
D0746   CD_SIG_nb       F       35
C1671   CD_TILE_b       M       35
D0755   CD_TILE_b       F       32
D0756   CD_TILE_b       F       41
D0757   CD_TILE_b       F       42
D0758   CD_TILE_b       M       28
C1674   CD_TILE_nb      F       39
C1686   CD_TILE_nb      F       18
C1672   CD_TILE_nb      F       39
C1669   CD_TILE_nb      F       47
D0759   CD_TILE_nb      F       35
D0760   CD_TILE_nb      F       25
C1680   Ctrl_SIG_nb     M       38
C1684   Ctrl_SIG_nb     M       60
C1656   Ctrl_SIG_nb     F       17
C1689   Ctrl_SIG_nb     F       72
C1663   Ctrl_SIG_nb     F       28
D0768   Ctrl_SIG_nb     F       52
C1679   Ctrl_TILE_nb    M       38
C1683   Ctrl_TILE_nb    M       60
C1655   Ctrl_TILE_nb    F       17
C1659   Ctrl_TILE_nb    M       51
C1662   Ctrl_TILE_nb    F       28
D0761   Ctrl_TILE_nb    F       52
C1668   DCtrl_SIG_b     F       63
C1682   DCtrl_SIG_b     F       44
C1493   DCtrl_SIG_b     F       43
D0748   DCtrl_SIG_b     M       68
C1676   DCtrl_SIG_nb    F       28
C1649   DCtrl_SIG_nb    F       33
D0749   DCtrl_SIG_nb    F       51
D0769   DCtrl_SIG_nb    F       45
C1677   DCtrl_TILE_b    F       28
D0762   DCtrl_TILE_b    M       68
C1648   DCtrl_TILE_nb   F       33
C1681   DCtrl_TILE_nb   F       44
C1664   DCtrl_TILE_nb   F       43
D0763   DCtrl_TILE_nb   F       24
D0764   DCtrl_TILE_nb   F       51
C1661   UC_SIG_b        M       40
C1666   UC_SIG_b        M       37
C2485   UC_SIG_b        F       25
C1658   UC_SIG_b        M       32
C1494   UC_SIG_nb       F       42
D0751   UC_SIG_nb       F       30
D0752   UC_SIG_nb       F       55
D0753   UC_SIG_nb       M       44
D0754   UC_SIG_nb       F       27
C1665   UC_TILE_b       M       37
D0765   UC_TILE_b       F       24
C1660   UC_TILE_nb      M       40
C1650   UC_TILE_nb      M       42
C1653   UC_TILE_nb      F       34
C1652   UC_TILE_nb      M       25
C2484   UC_TILE_nb      F       25
C1678   UC_TILE_nb      F       42", header = FALSE)

colnames(metadata) <- c("ID", "Label", "Gender", "AgeAtDateOfSampling")

# Fix ids of samples of the transcripts to match metadata
o <- strsplit(colnames(transcripts), "_")
ids_trans <- vapply(o, function(x){x[[1]]}, FUN.VALUE = character(1L))
stopifnot(length(ids_trans) == 63)
stopifnot(all(ids_trans %in% metadata$ID))
colnames(transcripts) <- ids_trans

# Fix ids of samples of the otus to match metadata
o <- strsplit(colnames(otus)[1:59], "_")
ids_trans <- vapply(o, function(x){x[[length(x)]]}, FUN.VALUE = character(1L))
stopifnot(length(ids_trans) == 63)
stopifnot(all(ids_trans %in% metadata$ID))
colnames(otus)[1:59] <- ids_trans

length(intersect(colnames(otus)[1:59], colnames(transcripts)))

paste(setdiff(metadata$ID, ids_trans), collapse = ", ")
Zeros <- rowSums(otus[, 1:59] == 0)
NonZeros <- rowSums(otus[, 1:59] != 0)

stopifnot(all(otus[, "CountZeros"] >= Zeros ))
stopifnot(all(otus[, "CountNonZeros"] >= NonZeros ))

