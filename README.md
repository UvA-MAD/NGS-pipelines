# NGS-pipelines
MAD-RBAB Pipelines for analysis of NGS experiments. Combines en inherits material form basicQC and rnaxcount.  
**Usage**

 - Clone project in to Script directory of experiment 
 - Pick the workflow you want to run.  
 - Edit the snakefile to set relevant variables such as databases and species
 - Use snakemake -s [ workflow ]  to run the workflow.

Additional snakemake flags which might be useful:
-j the number of concurrent jobs

files 
-----
**Workflows**

- Snakefile.qc_basic  : workflow to generate the qc file
- Snakefile.microrna  : process a microRNA experiment
- Snakefile.mrna_euk  : process a eukaryotic mRNA experiment
- Snakefile.mrna_pro  : process a prokaryotic mRNA experiment

**Suporting files**

- LICENSE : GPL licence
- README.md : this file
- sRNA_tools.py : python fuction library use in mapping and counting
- basicQC.Rmd : Input file for basic QC markdown file
- get_experiment_fq.py : script to generate fastq files per sample based on design file and bam files from run
- .Rprofile,  packrat/init.R, packrat/packrat.lock: files in packratsetup
