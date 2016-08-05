# makeREFSEQ2EG.R
# Aron
# 2011-12-08
#

## first update  org.Hs.eg.db
source("http://bioconductor.org/biocLite.R")
update.packages(repos=biocinstallRepos(), ask=FALSE, 
  checkBuilt=TRUE, oldPkgs = 'org.Hs.eg.db')

library(Biobase)
library(org.Hs.eg.db) 

refseq2eg.list <- as.list(org.Hs.egREFSEQ2EG)
refseq2eg <- sapply(refseq2eg.list, function(x) x[1])
out <- data.frame(refseq = names(refseq2eg), 
        entrez = refseq2eg, stringsAsFactors = FALSE)
write.table(out, file = 'REFSEQ2EG.txt', sep = '\t',
  row.names = FALSE, quote = FALSE)

info <- package.version('org.Hs.eg.db')
writeLines(info, 'REFSEQ2EG.version')

