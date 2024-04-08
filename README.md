# Codes for symptom-based depression subtypes

This repo contains the codes and collect data for replicating the results in our paper about exploring the biological underpinnings for MDD subtypes.

## inputs

This folder contains the raw subject list table (subject_merge_MDD_HC.xlsx), and the generated files for each analysis pipeline (begins wiht "Analysis"). Moreover, the gene rank of each subtype is presented, which was used to run GO analysis.

## scripts

This folder contains all codes used to analyze the data.

In **/matlab** sub-folder, there are codes for running time delay estimation. 

In **/python** sub-folder, it contains python script to clean and extract the BOLD signal, and to run abagen as well. 

In **/R** sub-folder, there are R codes to run the main analysis, including LPA, PLS and statistics.

## outputs

This folder contains the files generated in R codes, including figures, tables and gene expression ranks.
