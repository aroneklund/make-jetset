# 2012-05-02
# Aron Eklund
#
# makeCSVfiles.R
#
# make csv files for website from R packages


library(jetset)
library(org.Hs.eg.db)
ver <-  packageDescription("jetset", fields = "Version")

plat <- c('hgu95av2', 'hgu133a', 'hgu133plus2', 'u133x3p')

# scores
for (p in plat) {
    scores.name <- paste('jetset.scores', p, sep = '.')
    scores <- jscores(p)
    scores.csv.name <- paste(scores.name, '_', ver, '.csv', sep = '')
    scores.zip.name <- paste(scores.csv.name, '.zip', sep = '')
    isBest = !is.na(scores$EntrezID) & !duplicated(scores$EntrezID)
    new.scores <- data.frame(probeset = rownames(scores),
                       scores,
                       best = isBest,
                       stringsAsFactors = FALSE)
    write.csv(new.scores, file = scores.csv.name, 
      row.names = FALSE, na = '--')
    system(paste('zip', scores.zip.name, scores.csv.name))
    unlink(scores.csv.name)
}

