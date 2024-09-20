# tssH
This repository supports [Symbiotic T6SS affects horizontal transmission of Paraburkholderia bonniea among 
Dictyostelium discoideum amoeba hosts](). RNA-sequencing data (raw and processed) are available from NCBI GEO 
GSE276651. 

It contains two directories:

**phenotype**, which contains
*   **tssh_pheno_analysis.R** - R code used to run statistical analyses using the data files and generate figures
*   **tssH.pheno_fitness.tsv** - Data from host fitness experiment with columns in the following order
    *   date of experiment
    *   host identity
    *   symbiont identity
    *   MOI (multiplicity of infection) of host-symbiont pairing
    *   percent of host spores infected by symbiont in sample
    *   percent of spores produced by sample relative to mean of uninfected controls
    *   symbiont treatment (wildtype or mutant)
*   **tssH.pheno_transmission.tsv** - Data from symbiont transmission experiment with columns in the following 
order
    *   date of experiment
    *   host identity
    *   symbiont identity
    *   MOI of host-symbiont pairing
    *   percent of host spores infected by symbiont in sample
    *   percent of infected spores that were previously uninfected
    *   symbiont treatment (wildtype or mutant)
*   **tssH.pheno_rescue.tsv** - Data from mutant rescue experiment with columns in the following order
    *   date of experiment
    *   host identity
    *   symbiont identity
    *   MOI of host-symbiont pairing
    *   percent of host spores infected by symbiont in sample
    *   percent of infected spores that were previously uninfected
    *   symbiont treatment (wildtype, mutant, rescue)

**RNA-seq**, which contains
*   **alignment_to_count_summary.txt** - Summary of pipeline to generate count tables from raw illumina reads
*   **tssh_de_analysis.R** - R code used to run statistical analyses using count files and generate figures
*   **gene_association.dictybase.filter.gostats.txt** - Processed GO annotation file for _D. discoideum_ genome
*   **endo_phago_path.tsv** - Table of _D. discoideum_ genes with known roles in endocytosis and phagocytosis, 
with columns in the following order
    *   gene symbol
    *   [dictyBase](http://dictybase.org/) gene id
    *   decription of protein coded by gene
*   **tssH.de_*.txt** - Three DESeq2 summary tables for the time series contrasts performed in the study
*   **tssH.GO_*.txt** - Six GOstats summary tables for the time series contrasts performed in the study; tables 
for up- and down-regulated groups of genes are separated
*   **tssH.phago_*.txt** - Five tables of phagocytosis-related genes clustered by DEGreport based on expression 
change over time in the study, with columns in the following order
    *   normalized transcript abundance
    *   timepoint (negative, 090 minutes-post-infection, 360 minutes-post-infection, and positive approximately 
2250 minutes-post-infection)
    *   gene symbol

