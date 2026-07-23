260617

Sam Zavislan-Pullaro

Project: Naly's RNAi experiments with LP306 membrane marker strain 

Purpose: Quantify signal across plasma membrane in 2-cell embryos and generate plots 

Goal: Run the R script and troubleshoot/adopt for RNAi experiments. This one will be for bmk-1 RNAi against the erm-1 set-3 control (L4440). L4440 is an empty vector control for RNAi. 

Steps: 

1) Copy datafiles from fiji macros membrane quantification script output

PATH to folder containing each rep for erm1 set3 controls (L4440): 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260521_RNAi_treatments/L4440_RNAi_control

PATH to folder containing each rep for bmk-1 RNAi: 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260521_RNAi_treatments/BMK-1_RNAi

Note: Each datafile is found within the 03_output directory for each rep

PATH to input folder for the R script: 

/Users/samzavislan/Desktop/onish/people/naly/projects/erm-1_LP306_RNAi/SZP260525_R_analysis/bmk-1_L4440/01_input

Names of input files: 

controls: 

2026-6-17_1450_datafile_for_230713_L4440_RNAi_rep1.txt
2026-6-17_153_datafile_for_230828_L4440_RNAi_rep2.txt
2026-6-17_158_datafile_for_240122_LP306_L4440-rep3.txt

bmk-1 RNAi: 

2026-5-26_1149_datafile_for_220719_LP306_bmk-1_RNAi.txt
2026-5-26_1221_datafile_for_230408_LP306_BMK-1_rep1.txt
2026-5-26_1242_datafile_for_231107_LP306_BMK-1_rep2.txt
2026-5-26_1449_datafile_for_231113_LP306_BMK-1_rep3.txt

Erin had a naming convention for each image name in the nocodazole experiments. I manually modified each image name in the txt files to follow this convention so that extracting salient information from each image file name is consistent for each file: 

date_strain_treatment_R3D.dv 









