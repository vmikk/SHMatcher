#!/usr/bin/env Rscript

## Convert SH-matches into UC-like format

library(data.table)

## Load matches
datt <- fread(file = "matches_out_all.csv", sep = "\t", header = TRUE)

## Rename columns
setnames(
  x = datt,
  old = c("SH code (0.5)", "SH code (1.0)", "SH code (1.5)", "SH code (2.0)", "SH code (2.5)", "SH code (3.0)"),
  new = c("SH_005", "SH_010", "SH_015", "SH_020", "SH_025", "SH_030"))

## Subset data
clz <- c("seq_accno", "SH_005", "SH_010", "SH_015", "SH_020", "SH_025", "SH_030")
datt <- datt[, ..clz]

## Print summary
cat("Number of sequences in the dataset: ", length(unique(datt$seq_accno)), "\n")
cat("Number of 0.5% SHs: ", length(unique(datt$SH_005)), "\n")
cat("Number of 1.0% SHs: ", length(unique(datt$SH_010)), "\n")
cat("Number of 1.5% SHs: ", length(unique(datt$SH_015)), "\n")
cat("Number of 2.0% SHs: ", length(unique(datt$SH_020)), "\n")
cat("Number of 2.5% SHs: ", length(unique(datt$SH_025)), "\n")
cat("Number of 3.0% SHs: ", length(unique(datt$SH_030)), "\n")


## Function to add leading zeros
leading_zero <- function(x){
  nn <- paste("%0", nchar(max(x)), "d", sep = "")
  res <- sprintf(nn, x)
  return(res)
}

## Function to create UC table
create_uc <- function(x){
	# x = two-column data.table; e.g., x <- datt[, .(seq_accno, SH_020)]

  cl <- colnames(x)[2]
  tr <- sub(pattern = "SH_", replacement = "", x = cl)
  setnames(x, old = cl, new = "SH")

  ## Rename SHs (consecutively by the number of sequences per SH)
  SHH <- x[, .(N = .N), by = "SH"]
  setorder(x = SHH, -N, SH)
  SHH[ , Num  := .I ]
  SHH[ , Numz := leading_zero(Num) ]
  SHH[ , SHID := paste0("SH_", tr, "_", Numz) ]
  SHH[ , Numz := NULL ]

  ## Add new ID to the match table
  x <- merge(x = x, y = SHH, by = "SH", all.x = TRUE)
  setorder(x = x, SHID)

  x[, Record := "H" ]   # 1. Record type == hit
  x[, Num := Num - 1 ]  # 2. Cluster number (0-based)
  # x[, .(N) ]          # 3. Cluster size
  x[, Sim := (100-as.numeric(tr)/10) ] # 4. Percent identity
  x[, Strand := "+" ]   # 5. Strand
  x[, V6 := "." ]       # 6. Not used
  x[, V7 := "." ]       # 7. Not used
  x[, V8 := "." ]       # 8. Compressed alignment or the symbol '=' (equals sign)
  # x[, .(seq_accno) ]  # 9. Label of query sequence
  # x[, .(SHID)]        # 10. Label of target sequence

  clz <- c("Record", "Num", "N", "Sim", "Strand", "V6", "V7", "V8", "seq_accno", "SHID")
  x <- x[ , ..clz]

  return(x)
}


## Create UC files for all thresholds
UC_005 <- create_uc(x = datt[, .(seq_accno, SH_005)])
UC_010 <- create_uc(x = datt[, .(seq_accno, SH_010)])
UC_015 <- create_uc(x = datt[, .(seq_accno, SH_015)])
UC_020 <- create_uc(x = datt[, .(seq_accno, SH_020)])
UC_025 <- create_uc(x = datt[, .(seq_accno, SH_025)])
UC_030 <- create_uc(x = datt[, .(seq_accno, SH_030)])


## Export UC files
export_uc <- function(uc, trsh = "005"){
  fwrite(
    x = uc,
    file = paste0("UC_", trsh, ".uc.gz"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,
    compress = "gzip")
}

export_uc(uc = UC_005, trsh = "005")
export_uc(uc = UC_010, trsh = "010")
export_uc(uc = UC_015, trsh = "015")
export_uc(uc = UC_020, trsh = "020")
export_uc(uc = UC_025, trsh = "025")
export_uc(uc = UC_030, trsh = "030")

