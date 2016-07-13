import os

#####################################################################
# Filesystem parameters
#####################################################################

#location of the experiment (variable only used in this file) 
EXPERIMENT_DIR = "../.."

# map to use for temporary files.  Removal is not automatic  
SCRATCH_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/mrna_map")

# Location of final results 
RESULT_DIR = os.path.join(EXPERIMENT_DIR,"Results/mrna_map")

# directory with raw reads in fastq format
FQ_DIR = os.path.join(EXPERIMENT_DIR,"Scratch/raw")

# samples to be processed if unchanged: all files with fastq extention in FQ_DIR
SAMPLES =  [s[:-6] for s in os.listdir(FQ_DIR) if s.endswith(".fastq")]


#####################################################################
# Experiment parameters
#####################################################################

#domain for which pipeline is run pick use "prokaryote"  or "eukaryote"
DOMAIN= "xxx"

# species for which pipeline is run
SPECIES = "xxx"

# fasta file containing genome 
# in an ideal world the database is found by species alone 
GENOME_DB=os.path.join(SCRATCH_DIR,"genome.fa")
# fasta file with transcriptome: mandatory for eukaryotes not used for prokaryotes
TRANSCRIPTOME_DB=os.path.join(SCRATCH_DIR,"transcriptome.fa")

# url which points to gzipped fasta file with genome and transcriptome
# (used if GENOME_DB or TRANSCRIPTOME is not found) 
# example  mouse data from ensembl
#ENSEMBL_84_URL="ftp://ftp.ensembl.org/pub/release-84/fasta/"
#GENOME_FASTA_URL=ENSEMBL_84_URL+"mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
#CDNA_FASTA_URL=  ENSEMBL_84_URL+"mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
#NCRNA_FASTA_URL= ENSEMBL_84_URL+"mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"

GENOME_ZIP_FASTA_URL=""
CDNA_ZIP_FASTA_URL=""
NCRNA_ZIP_FASTA_URL=""

#  genome annotation with transcriptome: mandatory for prokaryotes not used for eukaryotes
GENOME_ANNOT="/xxx/xxx.gff"

# url which points to zipped gff genome annotation (used if GENOME_ANNOT is not found) 
GENOME_ZIP_GFF_URL="xxx://xxx.yyy.zzz/abc"

# reference spike-in sequences
SPIKES_REF = "/mad/MAD-RBAB/05_Reference-db/external/ERCC/ERCC92.fa"


# triming of reads. 
# reads shorter than MIN_READLEN wil be removed from analysis 
# reads longer than MAX_READLEN will be tail clipped to a length of MAX_READLEN
MIN_READLEN ="25"
# In earlier versions of the pipeline reads were "illuminised": truncated at 50 nt
# By setting the the MAX_READLEN to a high value no reads are truncated. 
MAX_READLEN ="99999"
