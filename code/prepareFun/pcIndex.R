### PC-index
##############################
# function to calculate pc-index, written by Minzhe, modified from Sangin's code.


pcIndex <- function(gStd.data, predRank.data, num = 10000, seed = 122){
      
      ### read gold standard file
      drugPair <- paste(gStd.data$drugA, gStd.data$drugB, sep=" & ")
      
      rank.gStd <- sort(gStd.data$eob, index.return = TRUE)$ix
      names(rank.gStd) <- drugPair[rank.gStd]
      
      eob <- gStd.data$eob[rank.gStd]
      eobErr <- gStd.data$eobErr[rank.gStd]

      
      ### read the prediction rank
      rank.pred <- rank.gStd*0
      drugPair.pred <- paste(predRank.data$drugA, predRank.data$drugB, sep = " & ")
      drugPair.pred.rev <- paste(predRank.data$drugB, predRank.data$drugA, sep = " & ")
      
      for(i in 1:max(rank.gStd)) {
            if (names(rank.pred)[i] %in% drugPair.pred) {
                  rank.pred[i] <- predRank.data$Rank[which(drugPair.pred == names(rank.gStd)[i])]
            } else if (names(rank.pred)[i] %in% drugPair.pred.rev) {
                  rank.pred[i] <- predRank.data$Rank[which(drugPair.pred.rev == names(rank.gStd)[i])]
            }
      }

      ### compute weighted c-index, and p-value
      pMatrix <- prob.mat(eob,eobErr)
      c.index <- con.idx(x = rank.pred, y = rank.gStd, pMatrix = pMatrix)

      if (num > 999) {
            null.dist <- rep(0, num)
            for (i in 1:num) {
                  null.dist[i] <- con.idx(x = sample(rank.gStd), y = rank.gStd, pMatrix = pMatrix)
            }
            p.value <- sum(null.dist>=c.index)/num
      } else { p.value = NULL }
      
      return(list(c.index = c.index, p.value = p.value))       
  
}

### Function for making probility matrix
# ==========================================
prob.mat <- function(x, x.sd){
      
      n <- length(x)
      x.tmp <- matrix(x, nrow = n, ncol = n)
      x.mat <- x.tmp - t(x.tmp) # i-j
      
      sd.tmp <- matrix(x.sd, nrow = n, ncol = n)
      sd.mat <- sqrt(sd.tmp^2 + t(sd.tmp)^2)
      
      erf <- function(a) pchisq(2*a^2, 1) * sign(a)
      pMatrix <- 0.5 * (1 + erf(x.mat/sd.mat))
      
      return(pMatrix)
      
}


### Function for computing concordance index
# ============================================
con.idx <- function(x, y, pMatrix) {
      
      if(length(x) != length(y)) {
            cat("Error messange: the size of two vectors are different.\n")
            break
      }
      n <- length(x)
      x.mat <- matrix(x, nrow = n, ncol = n)
      y.mat <- matrix(y, nrow = n, ncol = n)
      
      Con  <- sign(x.mat-t(x.mat)) == sign(y.mat-t(y.mat))
      wCon <- Con*(1-t(pMatrix)) + (1-Con)*t(pMatrix)
      pIndex <- sum(wCon[lower.tri(wCon)])/choose(n,2)
      
      return(pIndex)
}



