#! /usr/bin/env Rscript 
#
# calculateJetset.R
#
# 2012-04-30
# Aron Eklund
#
# Calculate jetset scores
#  argument should be one of 'hgu95av2', 'hgu133a', 'hgu133plus2', 'u133x3p'

args <- commandArgs(TRUE)
stopifnot(length(args) == 1)

## get REFSEQ2EG version
refseq2eg.version <- readLines('REFSEQ2EG.version')

## get human.rna version
human.rna.version <- readLines('human.rna.version')

## refseq-entrez lookup
message('Reading REFSEQ2EG.txt')
refseq2eg.df <- read.delim('REFSEQ2EG.txt', colClasses = 'character')
refseq2eg <- refseq2eg.df$entrez
names(refseq2eg) <- refseq2eg.df$refseq
rm(refseq2eg.df)

## refseq length lookup
message('Reading human.rna.len')
refseq.length.df <- read.delim('human.rna.len',
  header = FALSE, stringsAsFactors = FALSE)
refseq.length <- refseq.length.df[, 2]
names(refseq.length) <- refseq.length.df[,1]
rm(refseq.length.df)

## calculate jetset scores from a data frame of blast results; 
## x = data frame with all hits for a single probe set
## n = number of probes in the probe set
## cuthi, cutlo = thresholds for defining strong and moderate hits
getScores <- function(x, n, cuthi = 48, cutlo = 32) {
  ## Filter out hits below the lower threshold, and hits on the wrong strand
  x <- x[x$bitScore >= cutlo & x$sStart < x$sEnd, ]
  ## Map Refseq hits to Entrez gene ID ("eg")
  x$refseqNoVersion <- sub('\\..*$', '', x$subjectID, perl = TRUE)
  x$eg <- refseq2eg[x$refseqNoVersion]
  ## Filter out hits that are not matched to a gene
  x <- x[!is.na(x$eg), ]   
  ## Calculate the target and specificity score 
  if(nrow(x) > 0) {
    # For each probe, identify the targeted genes (if unique)
    detectedGenes <- sapply(split(x, x$queryID), function(y) {
      eg.strong <- unique(y$eg[y$bitScore >= cuthi])
      eg.mod <- unique(y$eg)
      if(length(eg.strong) == 1 && length(eg.mod) == 1) { 
        eg.strong
      } else {
        NA
      }
    })
    detectedGenesCounts <- sort(table(detectedGenes), decreasing = TRUE)
    if(length(detectedGenesCounts) > 1 &&  detectedGenesCounts[1] == detectedGenesCounts[2]) {
      best <- NA  ## >1 target genes with same score
    } else {
      best <- detectedGenesCounts[1]
    }
  } else {  ## no blast hits
    best <- NA
  }
  if(!is.na(best)) {
    target <- names(best)
    specificity <- unname(best / n)
  } else {
    target <- NA
    specificity <- NA
  }
  ## Now calculate Coverage and Processivity requirement
  if(!is.na(target)) {
    ## identify the subset of x that is both Strong and On-target
    xso <- x[x$eg == target & x$bitScore > cuthi, ]  
    # coverage
    tgt.allrefseq <- names(refseq2eg)[refseq2eg == target]
    tgt.allrefseq.nm <- grep('^[NX][MR]', tgt.allrefseq, value = TRUE)  # limit to RNA
    xq.refseq <- lapply(split(xso$refseqNoVersion, xso$queryID), unique)
    xq.refseq.table <- table(unlist(xq.refseq))
    xq.refseq.table.detected <- names(xq.refseq.table)[xq.refseq.table > (n/2)]
    coverage <- mean(tgt.allrefseq.nm %in% xq.refseq.table.detected)
    # process
    dists <- refseq.length[xso$subjectID] - xso$sStart
    process = median(dists)
  } else {
    coverage = NA
    process = NA
  }
  list(target = target, process = process, specificity = specificity, coverage = coverage)  
}  

# apply "getScores" over each probe set in blast results
# NOte: this uses the unix command "cut" and thus may not work on Windows
blastapply <- function(filename, cdf) {
  nProbes <- sapply(as.list(cdf), nrow)
  nProbes <- nProbes[order(names(nProbes))]
  message('Prescanning ', filename, '...', appendLF = FALSE)
  p <- pipe(paste('cut -f 1', filename))
  all.queryID <- scan(p, what = 'character')
  close(p)
  all.probeset <- sub(':.*$', '', all.queryID, perl = TRUE) 
  idx <- factor(all.probeset)
  grouplength <- c(0, which(idx[-1] != idx[-length(idx)]), length(idx))
  groupsize <- diff(grouplength)
  message('found ', length(groupsize), ' probe sets')
  nms <- as.character(idx[grouplength[-1]])
  open(con <- file(filename))
  on.exit(close(con))
  n <- length(nProbes)
  out <- data.frame(nProbes = nProbes, 
    EntrezID = as.character(rep(NA, n)), process = as.integer(rep(NA, n)), 
    specificity = as.numeric(rep(NA, n)), coverage = as.numeric(rep(NA, n)), 
    row.names = names(nProbes), stringsAsFactors = FALSE) 
  message('Processing file, wait for ', length(groupsize) %/% 1000, ' dots:')
  for (i in 1:length(groupsize)) {
    if(i %% 1000 == 0) message('.', appendLF = FALSE)
    x <- read.delim(con, nrows = groupsize[i],
      header = FALSE, stringsAsFactors = FALSE,
      col.names = c("queryID", "subjectID", "percentId", "alignLength", "mismatches", 
            "gapOpenings", "qStart", "qEnd", "sStart", "sEnd", 
            "eVal", "bitScore")
    )
    intendedProbesetID <- nms[i]
    x.probesetIDs <- sub(':.*', '', x$queryID, perl = TRUE) 
    stopifnot(all(x$probeset == intendedProbesetID))
    sc <- getScores(x, n = nProbes[intendedProbesetID])
    out[intendedProbesetID, 'EntrezID'] <- sc$target
    out[intendedProbesetID, 'process'] <- sc$process
    out[intendedProbesetID, 'specificity'] <- sc$specificity
    out[intendedProbesetID, 'coverage'] <- sc$coverage
  }
  message('done')
  out
}

go <- function(chip) {
  cdfname <- paste(chip, 'cdf', sep = '')
  blastname <- paste('blastresult.', chip, '.refseq.txt', sep = '')
  objname <- paste('scores.', chip, sep = '')
  savename <- paste(objname, '.RData', sep = '')
  library(cdfname, character.only = TRUE)
  message('Calculating jetset scores for ', chip)
  tempScores <- blastapply(filename = blastname, cdf = get(cdfname))
  p <- ifelse(chip == 'u133x3p', 300, 600)
  temp.robust <- (1 - (1/p)) ^ tempScores$process
  temp.overall <- tempScores$specificity * tempScores$coverage * temp.robust
  sortedScores <- tempScores[order(temp.overall, decreasing = TRUE), ]
  attr(sortedScores, 'version.Refseq') <- human.rna.version
  attr(sortedScores, 'version.org.Hs.eg.db') <- refseq2eg.version
  assign(objname, sortedScores)
  save(list = objname, file = savename)
}

go(args)
message('All finished, no problem')
q(save = 'no')

