###                    process_pipline.R                   ###
### ====================================================== ###
# This R script is a pipline to analyze drug pair synergy

library(argparse)
library(preprocessCore)
library(org.Hs.eg.db)
library(KEGGgraph)
library(ggplot2)
library(grid)


### 0. Parse comandline argument
parser <- ArgumentParser(description = "This pipline is to analyze drug pair syngergy.")
parser$add_argument("jobID", type = "character", help = "Passing the job ID")
parser$add_argument("-p", "--pathway", type = "integer", default = 1, help =  "Specify the pathway information want to use: 1 for KEGG pathway information, 2 for constructed lymphoma gene network(partial corelation) information, 3 for constructed lymphoma gene network(marginal corelation) information. Default is 1.")

pipArgs <- parser$parse_args()
jobID <- pipArgs$jobID
pathway <- pipArgs$pathway

logPath <- paste("report/log", jobID, ".txt", sep = "")
dosePath <- paste("data/doseRes", jobID, ".csv", sep = "")
geneExpPath <- paste("data/GeneExpr", jobID, ".csv", sep = "")
rankPath <- paste("report/pred.pairRank", jobID, ".csv", sep = "")
p.heat.path <- paste("report/score_heatmap", jobID, ".jpeg", sep = "")
p.bar.path <- paste("report/score_rank", jobID, ".jpeg", sep = "")
statPath <- paste("report/status", jobID, ".txt", sep = "")



### 1. Start tracking status
if (!dir.exists("report/")) {
      dir.create("report/")
      file.create(logPath)
} else if (!file.exists(logPath)) {
      file.create(logPath)
}
f <- file(logPath, open = "wt")
sink(f, type = "message")


### 2. Load functions
source("code/doseRes.R")
source("code/profileGeneExp.R")
source("code/scoring.R")
source("code/plotting.R")



### 3. Read data

# read dose response data
doseRes <- readDoseRes.csv(dosePath)

# prepare drug treated gene expression data
geneExpDiff <- profileGeneExp(geneExpPath)

# load pathway information



### 4. Analyze drug pair synergy
if (pathway == 1) {
      load("data/pathInfo.Rdata")
      res <- scoring(geneExpDiff = geneExpDiff$geneExp.max, doseRes = doseRes, CGP.mat = CGP.mat, GP.mat = GP.mat, fold = 0.6)
}
if (pathway == 2) {
      load("data/gLassoGeneNet.Rdata")
      res <- scoring(geneExpDiff = geneExpDiff$geneExp.max, doseRes = doseRes, CGP.mat = mCGP.mat, GP.mat = mGP.mat, fold = 0.6)
}
if (pathway == 3) {
      load("data/marGeneNet.Rdata")
      res <- scoring(geneExpDiff = geneExpDiff$geneExp.max, doseRes = doseRes, CGP.mat = mCGP.mat, GP.mat = mGP.mat, fold = 0.6)
}
score.rank <- res$scoreRank
write.csv(score.rank, rankPath, row.names = FALSE)




### 5. Make plots

# plot heatmap
p.heat <- pair.ggheat(pred.pair = score.rank)
ggsave(p.heat, filename =  p.heat.path, width = 8, height = 6)

# plot barplot
p.bar <- pair.ggbar(pred.pair = score.rank)
ggsave(p.bar, filename = p.bar.path, width = 8, height = 7)



### 6. Check error
sink(type = "message")
close(f)

Lines <- readLines(logPath)
Err <- any(grepl("Error", Lines))
if (Err) {
      write("Fail", file = statPath)
} else {
      write("Success", file = statPath)
}