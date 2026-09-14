library("NCmisc")

knitr::purl("260422_Nocodozole_treatment_erm1_smFISH_wNT002.Rmd", output = "temp_script.R")
list.functions.in.file("temp_script.R")

getwd()
