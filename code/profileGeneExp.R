###                     profileGeneExp.R                   ###
### ====================================================== ###
# This R function is to profile drug treated gene expression data




profileGeneExp <- function(file) {
      
      library(preprocessCore)
            
	### Read drug treated gene expression data
	geneExp <- readGeneExp.csv(file)
	
	### Average duplicated drug data
	geneExp.mat <- geneExp.aveDrug(geneExp)
	
	### Normalize data
	geneExp.mat <- normGeneExp(geneExp.mat)
	
	### Subtract nagtive control effect
	geneExp.mat <- geneExp.diffCtl(geneExp.mat)
	
	### Convert probe level to gene level
	geneExp.profile <- prob2gene(geneExp.mat)
	
	return(geneExp.profile)
      
}


### Read drug treated gene expression data
readGeneExp.csv <- function(file) {
      geneExp <- read.csv(file = file, header = FALSE, stringsAsFactors = FALSE)
      return(geneExp)
}


### Average duplicated drug data in gene expression file
geneExp.aveDrug <- function(geneExp) {
      
      ### Check duplicated drugs
      if (any(duplicated(as.character(geneExp[1,-1])))) {
            cat("***Duplicated drug name found, average duplicated data.\n")
      }
      
      ### Print drug list
      drugName <- sort(unique(as.character(geneExp[1,-1])))
      cat("***Drug list you provided:\n")
      for (i in 1:length(drugName)) {
            if (drugName[i] != "Neg_control") cat(drugName[i], "\n")
      }
      
      ### Average dupliacte drug data (if any)
      geneExp.mat <- data.frame(matrix(NA, nrow = nrow(geneExp)-1, ncol = length(drugName)+1))
      colnames(geneExp.mat) <- c("Gene", drugName)
      geneExp.mat$Gene <- geneExp[-1,1]
      
      for (i in 1:length(drugName)) {
            tmpDrug.data <- geneExp[-1, geneExp[1,] == drugName[i]]
            tmpDrug.data <- apply(as.data.frame(tmpDrug.data), 2, as.numeric)
            tmpDrug.mean <- apply(tmpDrug.data, 1, mean)
            geneExp.mat[drugName[i]] <- tmpDrug.mean
      }
      
      return(geneExp.mat)
}


### Quantile normalize gene expression data
normGeneExp <- function(geneExp) {
      library(preprocessCore)
      
      drugName <- colnames(geneExp)[-1]; geneName <- geneExp[,1]
      geneExp.mat <- normalize.quantiles(as.matrix(geneExp[,-1]))
      colnames(geneExp.mat) <- drugName
      geneExp.mat <- data.frame(Gene = geneName, geneExp.mat, stringsAsFactors = FALSE, check.names = FALSE)
      
      return(geneExp.mat)
}


### Subtract nagtive control effect
geneExp.diffCtl <- function(geneExp) {
      
      neg_control <- geneExp[,"Neg_control"]
      for (i in 2:ncol(geneExp)) {
            geneExp[,i] <- geneExp[,i] - neg_control
      }
      
      ### delete DMSO column
      geneExp$Neg_control <- NULL

      return(geneExp)
}


### Conver probe level to gene level
prob2gene <- function(geneExp) {
      
      geneName <- sapply(strsplit(geneExp$Gene, " ", fixed = TRUE), "[[", 1)
      geneExp$Gene <- geneName
      
      geneExp.mean <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = mean)
      geneExp.max <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = max)
      geneExp.min <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = min)
      
      geneExp.profile <- list(geneExp.mean = geneExp.mean, geneExp.max = geneExp.max, geneExp.min = geneExp.min)
      return(geneExp.profile)
}



    
