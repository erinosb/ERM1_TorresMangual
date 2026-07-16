260527

Sam Zavislan-Pullaro

Project: Naly's RNAi experiments with LP306 membrane marker strain 

Purpose: Quantify signal across plasma membrane in 2-cell embryos and generate plots 

Goal: Run the R script and troubleshoot/adopt for RNAi experiments. This one will be for bmk-1 RNAi against the erm-1 set-3 control (no RNAi)

Steps: 

1) Copy datafiles from fiji macros membrane quantification script output

PATH to folder containing each rep for erm1 set3 controls (no RNAi): 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260521_RNAi_treatments/LP306_no_RNAi

PATH to folder containing each rep for bmk-1 RNAi: 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260521_RNAi_treatments/BMK-1_RNAi

Note: Each datafile is found within the 03_output directory for each rep

PATH to input folder for the R script: 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260525_R_analysis/bmk-1_RNAi/01_input

Names of input files: 

controls: 

2026-5-21_162_datafile_for_221011_LP306_no_treatment.txt
2026-5-21_1642_datafile_for_240610_LP306_erm-1.txt
2026-5-26_1047_datafile_for_240621_LP306_erm-1_set-3.txt
2026-5-26_1130_datafile_for_240701_LP306_erm-1_set-3.txt

bmk-1 RNAi: 

2026-5-26_1149_datafile_for_220719_LP306_bmk-1_RNAi.txt
2026-5-26_1221_datafile_for_230408_LP306_BMK-1_rep1.txt
2026-5-26_1242_datafile_for_231107_LP306_BMK-1_rep2.txt
2026-5-26_1449_datafile_for_231113_LP306_BMK-1_rep3.txt

Erin had a naming convention for each image name in the nocodazole experiments. I manually modified each image name in the txt files to follow this convention so that extracting salient information from each image file name is consistent for each file: 

date_strain_treatment_R3D.dv 









