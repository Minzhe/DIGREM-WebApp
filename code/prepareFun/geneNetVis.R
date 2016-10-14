interact.CGP.mat <- c("source", "target")

for (i in 1:nrow(CGP.mat)) {
      for (j in 1:ncol(CGP.mat)) {
            if (CGP.mat[i,j] != 0) {
                  
                  temp.interact <- c(row.names(CGP.mat)[i], colnames(CGP.mat)[j])
                  interact.CGP.mat <- rbind(interact.CGP.mat, temp.interact)
                  
            }
      }
}

colnames(interact.CGP.mat) <- c("source", "target")
interact.CGP.mat <- interact.CGP.mat[-1,]
row.names(interact.CGP.mat) <- NULL

write.table(interact.CGP.mat, "~/Project/Drug_combination/WebApp/data/KEGG.CGP.net.txt", sep = "\t", row.names = FALSE)
