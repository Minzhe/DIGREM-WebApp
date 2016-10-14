              process_pipline.R
=================================================

This R script is the pipline of DIGRE algorithm to analyze the drug pair synergistic score, using drug treated gene expression data and dose response curve data.

--------------   Dependency   ------------------
--- R ---

Following packages need to be installed:
1. argparse
2. proto
3. getopt
4. preprocessCore
5. org.Hs.eg.db
6. KEGGgraph
7. ggplot2
8. grid

--- python 2 ---
Following modules need to be installed:
1. argparse	 


--------------   Process   -----------------------
This script mainly do the following things:
1. Load functions
2. Read data
3. Analyze drug pair synergy
4. Make plots
5. Check errors


-----------------   Before start   ----------------------
To run this script, you need to:
1. Change the working directory of the "process_pipline.R", which is to change the path of the setwd("...") to your own path that contains this file.
2. Have the drug dose response data named as "doseRes.csv" in the data folder.
3. Have the drug treated gene expression data named as "data/GeneExpr.csv" in the data folder.
4. Have the pathway information named as "pathInfo.Rdata" in the data folder.
5. The result of the pipline would be generated in the report folder.

Detail function see in code folder.
Detail input data requirement see in the data folder.