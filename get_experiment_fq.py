import subprocess
import argparse
import pandas as pd
import re
import os
import sys


proton_run_store   = "/mad/MAD-RBAB/04_MAD-RBAB-runs/data"
illumina_run_store = "/mad/MAD-RBAB/04_MAD-RBAB-runs/illumina/data"

def illumina_gather_sample(sample_info):
    """
    For each sampleid generate fastq file(s).
    """
    sample_data, output_dir = sample_info
    sampleid = list(sample_data['sampleid'])[0]
    for pend in ["1","2"] : 
      file_pattern = '%s_S[1-9][0-9]*_R%s_001.fastq.gz'
      sample_gz = []
      for index, row in sample_data.iterrows():
          iid = row['runid']
          fa_gz_regexp = re.compile(file_pattern % (sampleid,pend))
          subdir, dirs, files = next(os.walk(os.path.join(illumina_run_store, iid)))
          for f in files:
              if fa_gz_regexp.match(f):
                  gz_file = os.path.join(illumina_run_store, iid, f)
                  sample_gz.append(gz_file)
      sample_fastq = os.path.join(output_dir, str(sampleid) + '_R'+pend+'.fastq.gz')
      if (len(sample_gz) > 0) : 
         zcat_command = ['zcat'] + sample_gz
         zip_command = ['gzip','-c'] 
         ps = subprocess.Popen(zcat_command,stdout=subprocess.PIPE)
         with open(sample_fastq, 'w') as fout:
            subprocess.call(zip_command,stdin=ps.stdout,stdout=fout)
            ps.wait()


def proton_sample_bam2fq(sample_info):
    """
    For each sampleid generate fastq file(s) from bam.
    """
    # unpack arguments
    sample_data, output_dir = sample_info
    # patttern of bam file names
    bam_pattern = '^IonXpress(RNA)?_0%s_rawlib.(basecaller.)?bam$'
    sampleid = list(sample_data['sampleid'])[0]
    sample_bams = []
    for index, row in sample_data.iterrows():
        rid = row['runid']
        bc = row['barcode'][2:]
        subdir, dirs, files = next(os.walk(os.path.join(proton_run_store, rid)))

        bam_regexp = re.compile(bam_pattern % bc)
        for f in files:
            if bam_regexp.match(f):
                bam_file = os.path.join(proton_run_store, rid, f)
                sample_bams.append(bam_file)

    # generate fastq files from bam files for this sample
    sample_fastqs = []
    for i, bam_file in enumerate(sample_bams):
        fastq_file = os.path.join(
            ".",
            '_'.join((str(sampleid), str(i))) + '.fastq'
            )
        sample_fastqs.append(fastq_file)
        subprocess.call(["bamToFastq","-i", bam_file, "-fq", fastq_file])
    if len(sample_fastqs) == 0:
        print("Could not find data for "+rid+" Barcode number "+bc); 
        sys.exit(1)
     
    # concatenate fastq for same sample from multiple runs
    sample_fastq = os.path.join(output_dir, str(sampleid) + '.fastq')
    cat_command = ['cat'] + sample_fastqs
    with open(sample_fastq, 'w') as fout:
        subprocess.call(cat_command, stdout=fout)

    # remove fastq  files from before concatenation
    [os.remove(f) for f in sample_fastqs]



def main(argv):
    # setup argument parse
    parser = argparse.ArgumentParser(
        description="Generate per sample fastq files from sequencer output.",
        epilog="Finds files for samples described in design file."
               "If sample has been run in muliple runs, fq files are combine to one."
               "Generated fastq files are placed in 'output_dir'.")

    parser.add_argument('-d', '--design_file',
                        type=argparse.FileType('r'),
                        help="path to design file")
    parser.add_argument('-a', '--output_dir',
                        type=str,
                        help="path to target directory with analysis")

    # print help if no command line arguments have been used
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # parse arguments
    args = parser.parse_args(argv)

    # parse design csv, which is tab separated file
    df = pd.read_csv(args.design_file, sep="\t")

    # first three column names should be case insensitive
    # therefore convert to lower case
    colnames = list(df.columns)
    colnames = [s.lower() for s in colnames]
    df.columns = colnames

    # get unique sample names
    # sample name can repeat it design file
    samples = list(set(df['sampleid']))

    # create dir for fastq files
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

# Now find out the type of experiment Proton or illumina 
    idtypes=list(set([s[0:3] for s in df['runid']]))
    if len(idtypes) != 1: 
        print("ERROR: get_experiment.py currently can not handle multi-platform designs."); 
        sys.exit(1)
    if idtypes == ["RID"]: 
       for sampleid in samples:
           sample_data = df[df['sampleid'] == sampleid]
           proton_sample_bam2fq((sample_data, args.output_dir))
    elif idtypes == ["IID"]:
       for sampleid in samples:
           sample_data = df[df['sampleid'] == sampleid]
           illumina_gather_sample((sample_data, args.output_dir))
    else:
        print("ERROR: get_experiment.py unknown sequencing platform id:"+idtypes[0]); 
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
