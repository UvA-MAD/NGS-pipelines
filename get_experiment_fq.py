import subprocess
import argparse
import pandas as pd
import re
import os
import sys


def sample_bam2fq(sample_info):
    """
    For each sampleid generate fastq file(s) from bam.
    """
    # unpack arguments
    sample_data, bam_home, output_dir = sample_info
    # patttern of bam file names
    bam_pattern = '^IonXpress(RNA)?_0%s_rawlib.(basecaller.)?bam$'
    sampleid = list(sample_data['sampleid'])[0]
    sample_bams = []
    for index, row in sample_data.iterrows():
        rid = row['runid']
        bc = row['barcode'][2:]
        subdir, dirs, files = next(os.walk(os.path.join(bam_home, rid)))

        bam_regexp = re.compile(bam_pattern % bc)
        for f in files:
            if bam_regexp.match(f):
                bam_file = os.path.join(bam_home, rid, f)
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
        description="Generate fastq files from bam files.",
        epilog="Finds bam files for samples described in design file."
               "If sample has been run in muliple runs, fq files are combine to one."
               "Generated fastq files are placed in 'output_dir'.")

    parser.add_argument('-d', '--design_file',
                        type=argparse.FileType('r'),
                        help="path to design file")
    parser.add_argument('-s', '--bam_home',
                        type=str,
                        help="path to original bam files")
    parser.add_argument('-a', '--output_dir',
                        type=str,
                        help="path to target directory with analysis")
    parser.add_argument('-n', '--numthreads',
                        default=1,
                        type=int,
                        help="number of threads, max 10")

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
    colnames[:3] = [s.lower() for s in colnames[:3]]
    df.columns = colnames

    # get unique sample names
    # sample name can repeat it design file
    samples = list(set(df['sampleid']))

    # create dir for fastq files
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    numthreads = int(args.numthreads)
    if numthreads > 1:
        # max threads is 10
        n = (10 if numthreads > 10 else numthreads)
        import multiprocessing
        arguments = []
        for sampleid in samples:
            # rows corresponding to particular sample
            # more than one if sample has been run in multiple runs
            sample_data = df[df['sampleid'] == sampleid]
            arguments.append((sample_data,
                              args.bam_home,
                              args.output_dir))

        pool = multiprocessing.Pool(processes=n)
        pool.map(sample_bam2fq, arguments)
    else:
        for sampleid in samples:
            sample_data = df[df['sampleid'] == sampleid]
            sample_bam2fq((sample_data, args.bam_home, args.output_dir))

if __name__ == "__main__":
    main(sys.argv[1:])
