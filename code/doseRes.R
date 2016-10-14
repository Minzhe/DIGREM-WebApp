###                        doseRes.R                       ###
### ====================================================== ###
# This R function is to read dose response curve data

readDoseRes.csv <- function(file) {
      
      doseRes <- read.csv(file = file, row.names = 1, check.names = FALSE)
      
      return(doseRes)
}