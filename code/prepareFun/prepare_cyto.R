###                    prepare_cyto.R                      ###
### ====================================================== ###
# This R script is to prepare gene-gene interaction files from gene network matrix



### 1. KEGG CGP information
source.gene <- c()
target.gene <- c()
interact <- c()
gene.name <- row.names(CGP.mat)

for (i in 1:nrow(CGP.mat)) {
      for (j in 1:ncol(CGP.mat)) {
            if (CGP.mat[i,j] != 0) {
                  source.gene <- c(source.gene, gene.name[i])
                  target.gene <- c(target.gene, gene.name[j])
                  interact <- c(interact, CGP.mat[i,j])
            }
      }
}

KEGG.CGP <- data.frame(source = source.gene, target = target.gene, action = interact)

# convert KEGG id to Gene Symbol
KEGG.CGP$source <- gsub("hsa:", "", KEGG.CGP$source)
KEGG.CGP$target <- gsub("hsa:", "", KEGG.CGP$target)
KEGG.CGP$source <- sapply(mget(KEGG.CGP$source, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
KEGG.CGP$target <- sapply(mget(KEGG.CGP$target, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
KEGG.CGP <- KEGG.CGP[order(KEGG.CGP$source),]
write.csv(KEGG.CGP, file = "KEGG.CGP.csv", row.names = FALSE)



### 2. KEGG GP information
source.gene <- c()
target.gene <- c()
interact <- c()
gene.name <- row.names(GP.mat)

for (i in 1:nrow(GP.mat)) {
      for (j in 1:ncol(GP.mat)) {
            if (GP.mat[i,j] != 0) {
                  source.gene <- c(source.gene, gene.name[i])
                  target.gene <- c(target.gene, gene.name[j])
                  interact <- c(interact, GP.mat[i,j])
            }
      }
}

KEGG.GP <- data.frame(source = source.gene, target = target.gene, action = interact)

# convert KEGG id to Gene Symbol
KEGG.GP$source <- gsub("hsa:", "", KEGG.GP$source)
KEGG.GP$target <- gsub("hsa:", "", KEGG.GP$target)
KEGG.GP$source <- sapply(mget(KEGG.GP$source, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
KEGG.GP$target <- sapply(mget(KEGG.GP$target, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
KEGG.GP <- KEGG.GP[order(KEGG.GP$source),]
write.csv(KEGG.GP, file = "KEGG.GP.csv", row.names = FALSE)



### 3. Marginal correlation GP
source.gene <- c()
target.gene <- c()
interact <- c()
gene.name <- row.names(mGP.mat)

for (i in 1:nrow(mGP.mat)) {
      for (j in 1:ncol(mGP.mat)) {
            if (mGP.mat[i,j] != 0) {
                  source.gene <- c(source.gene, gene.name[i])
                  target.gene <- c(target.gene, gene.name[j])
                  interact <- c(interact, mGP.mat[i,j])
            }
      }
}

Mar.mGP <- data.frame(source = source.gene, target = target.gene, action = interact)

# convert KEGG id to Gene Symbol
Mar.mGP$source <- gsub("hsa:", "", Mar.mGP$source)
Mar.mGP$target <- gsub("hsa:", "", Mar.mGP$target)
Mar.mGP$source <- sapply(mget(Mar.mGP$source, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
Mar.mGP$target <- sapply(mget(Mar.mGP$target, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
Mar.mGP <- Mar.mGP[order(Mar.mGP$source),]
write.csv(Mar.mGP, file = "Mar.GP.csv", row.names = FALSE)



### 4. Partial correlation GP
source.gene <- c()
target.gene <- c()
interact <- c()
gene.name <- row.names(mGP.mat)

for (i in 1:nrow(mGP.mat)) {
      for (j in 1:ncol(mGP.mat)) {
            if (mGP.mat[i,j] != 0) {
                  source.gene <- c(source.gene, gene.name[i])
                  target.gene <- c(target.gene, gene.name[j])
                  interact <- c(interact, mGP.mat[i,j])
            }
      }
}

Par.mGP <- data.frame(source = source.gene, target = target.gene, action = interact)

# convert KEGG id to Gene Symbol
Par.mGP$source <- gsub("hsa:", "", Par.mGP$source)
Par.mGP$target <- gsub("hsa:", "", Par.mGP$target)
Par.mGP$source <- sapply(mget(Par.mGP$source, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
Par.mGP$target <- sapply(mget(Par.mGP$target, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
Par.mGP <- Par.mGP[order(Par.mGP$source),]
write.csv(Par.mGP, file = "Par.GP.csv", row.names = FALSE)




###############################################
### clean gene network

raw.net <- read.csv("Par.GP.csv")
raw.net <- as.matrix(raw.net)

# remove duplicated interaction
par.geneNet <- c()
for (i in 1:nrow(raw.net)) {
      par.geneNet <- rbind(par.geneNet, sort(drop(raw.net[i,])))
}
par.geneNet <- par.geneNet[!(duplicated(par.geneNet)),]

# remove self interaction
idx <- c()
for (i in 1:nrow(par.geneNet)) {
      if (par.geneNet[i,1] == par.geneNet[i,2]) {
            idx <- c(idx, i)
      }
}
par.geneNet <- par.geneNet[-idx,]

write.table(par.geneNet, file = "Par.GP.txt", row.names = FALSE, sep = "\t")
