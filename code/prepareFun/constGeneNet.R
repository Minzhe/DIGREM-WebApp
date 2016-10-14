###                    constGeneNet.R                      ###
### ====================================================== ###
# This R script is to construct gene network matrix
library(org.Hs.eg.db)
library(KEGGgraph)


constGeneNet <- function(file) {
      rawInteract <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
      
      source.gene <- sapply(mget(rawInteract[,1], org.Hs.egSYMBOL2EG, ifnotfound = NA), "[[", 1)
      target.gene <- sapply(mget(rawInteract[,2], org.Hs.egSYMBOL2EG, ifnotfound = NA), "[[", 1)
      source.gene <- paste("hsa:", source.gene, sep = "")
      target.gene <- paste("hsa:", target.gene, sep = "")
      
      if (ncol(rawInteract) > 2) {
            interact <- as.numeric(rawInteract[,3])
            if (!all(interact %in% c(-1, 0, 1))) {
                  stop("Interaction type not recognized. Should be 1, 0 or -1 for activation, no interaction or inhibition respectively.")
            }
      } else {
            cat("No interaction type specified, default all for activation ...")
            interact <- rep(1, nrow(rawInteract))
      }
      
      # construct GP matrix
      gene.list <- union(source.gene, target.gene)
      GP.mat <- matrix(0, length(gene.list), length(gene.list), dimnames = list(gene.list, gene.list))
      
      for (i in 1:length(source.gene)) {
            GP.mat[source.gene[i],target.gene[i]] <- interact[i]
      }
      
      return(GP.mat)
}