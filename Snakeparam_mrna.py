
import os

#domain for which pipeline is run pick use "prokaryote"  or "eukaryote"
DOMAIN= "xxx"

# species for which pipeline is run
SPECIES = "xxx"

# fasta file containing 
# in an ideal world the database is found by species alone 
GENOME_DB="/xxx/xxx/xxx.fasta" 

# url which points to genome (used if GENOME_DB is not found) 
GENOME_FASTA_URL="xxx://xxx.yyy.zzz/abc"


# fasta file with transcriptome: mandatory for eukaryotes not used for prokaryotes
TRANSCRIPTOME_DB="/xxx/xxx/xxx.fa"

#  genome annotation with transcriptome: mandatory for prokaryotes not used for eukaryotes
GENOME_ANNOT="/xxx/xxx.gff"

# url which points to genome annotation (used if GENOME_ANNOT is not found) 
GENOME_GFF_URL="xxx://xxx.yyy.zzz/abc"

#location of the experiment (variable only used in this file) 
EXPERIMENT_DIR = "../.."

# directory with raw reads in fastq format
FQ_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/raw")

# samples to be processed if unchanged: all files with fastq extention in FQ_DIR
SAMPLES =  [s[:-6] for s in os.listdir(FQ_DIR) if s.endswith(".fastq")]

# map to use for temporary files.
# Removal is not automatic  
# snakemake 
SCRATCH_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/mrna_map")

# map with scratch results 
RESULT_DIR = os.path.join(EXPERIMENT_DIR,"Results/mrna_map")

# reference spike-in sequences
SPIKES_REF = "/mad/MAD-RBAB/05_Reference-db/external/ERCC/ERCC92.fa"


# triming of reads. 
# reads shorter than MIN_READLEN wil be removed from analysis 
# reads longer than MAX_READLEN will be tail clipped to a length of MAX_READLEN
MIN_READLEN ="25"
MAX_READLEN ="50"
