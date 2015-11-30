import os

#domain for which pipeline is run pick use "prokaryote"  or "eukaryote"
DOMAIN= "xxx"

# species for which pipeline is run
SPECIES = "xxx"

# in an ideal world the database is found by species alone 
GENOME_DB="/xxx/xxx.dna.genome.fa"

# fasta file with transcriptome: mandatory for eukaryotes not used for prokaryotes
TRANSCRIPTOME_DB="/xxx/xxx.transcripts.fa"

#  genome annotation with transcriptome: mandatory for prokaryotes not used for eukaryotes
GENOME_ANNOT="/xxx/xxx.gff"

#location of the experiment (variable only used in this file) 
EXPERIMENT_DIR = "../.."

# directory with raw reads in fastq format
FQ_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/raw")

# samples to be processed if unchanged: all files with fastq extention in FQ_DIR
SAMPLES =  [s[:-6] for s in os.listdir(FQ_DIR) if s.endswith(".fastq")]

# map to use for temporary files. Removal is not automatic  
SCRATCH_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/mapped")
# map with scratch results 
RESULT_DIR = os.path.join(EXPERIMENT_DIR,"Results/mapping")

# reference spike-in sequences
SPIKES_REF = "/mad/MAD-RBAB/05_Reference-db/external/ERCC/ERCC92.fa"

