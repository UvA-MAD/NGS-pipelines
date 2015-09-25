import click
import pysam
import os
import csv
import itertools
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from BCBio import GFF
from itertools import groupby
from functools import reduce
from collections import Counter
import sys 

def filter_and_count(bamfile,spike_type):
    """
    Filters and counts reads in one bam file.
    Returns dictionary with the counts
    key -- spike
    value -- count of reads aligned to that spike
    """
    samfile = pysam.Samfile(bamfile, "rb")
    # get names of size spikein that reads have been align to
    spike_refs = samfile.references

    # make a dictionary to lookup min lenght of aligned reads
    min_read_size = {}
    if spike_type=="ERCC": 
       for spike in spike_refs:
           min_read_size[spike] = 32
    else: 
       for spike in spike_refs:
          size = int(spike.split("_")[2])
          if size <= 40:
              min_size = round(0.8 * size)
          else:
              min_size = 32
          min_read_size[spike] = min_size

    # iterate over references and count reads that have min_read_size
    spike_counts = defaultdict(int)
    for spike in spike_refs:
        alignments = samfile.fetch(spike)
        for aln in alignments:
            if aln.rlen >= min_read_size[spike]:
                spike_counts[spike] += 1
    return(spike_counts)


# create a click group to use subcommands
@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam-dir',
              type=click.Path(exists=True, resolve_path=True),
              help='directory with bam files')
@click.option('--count-dir',
              type=click.Path(),
              help='directory to save count table in')
def count_pirna(bam_dir, count_dir):
    """
    Makes count table for piRNAs
    """
    # get the list of files that match the basname
    bams = [f for f in os.listdir(bam_dir) if f.endswith("_sorted.bam")]
    samples = [b.split("_pirna_aln_sorted.bam")[0] for b in bams]

    ref_lists = [pysam.Samfile(os.path.join(bam_dir, bam),
                               'rb').references
                 for bam in bams]
    refs = list(set(itertools.chain(*ref_lists)))

    count_dict = {}
    for bam, sample in zip(bams, samples):
        samfile = pysam.Samfile(os.path.join(bam_dir, bam), 'rb')
        alns = samfile.fetch()
        counts = [0] * len(refs)
        for aln in alns:
            counts[aln.rname] += 1

    df = pd.DataFrame.from_dict(count_dict)
    df.index = refs
    csv_path = os.path.join(count_dir, "CountTable_pirna.txt")
    df.to_csv(csv_path, sep="\t")


# subcommand used for filtering and counting spikes
@cli.command()
@click.option('--basename', help='basename of the bam files with alignments')
@click.option('--spike_type', help='Type of spikes. Default: size and norm spikes ')
@click.option('--bam-dir',
              type=click.Path(exists=True, resolve_path=True),
              help='directory with bam files')
@click.option('--count-dir',
              type=click.Path(),
              help='directory to save count table in')
def count_spikes(basename, spike_type, bam_dir, count_dir):
    """
    Applies filter rules and counts reads witch adhere to ther rules.
    """
    # get the list of files that match the basname
    bams = [f for f in os.listdir(bam_dir) if basename in f]
    bams = [f for f in bams if f.endswith(".bam")]
    # extract sample names from bam file names
    samples = ["".join(b.split(basename)).strip() for b in bams]

    count_dicts = []
    for b in bams:
        bamfile = os.path.join(bam_dir, b)
        counts = filter_and_count(bamfile,spike_type)
        count_dicts.append(counts)

    # put count dicts into a pandas datatable
    count_table = pd.DataFrame(count_dicts, index=samples)
    # that gives table with spikes as columns and we need samples coulmns
    # and spikes rows, therefore we transpoze the table
    count_table = count_table.T

    # check if output dir is available, if not create one
    if not os.path.exists(count_dir):
        os.makedirs(count_dir)

    # fill NaN's with zeros (zero counts)
    count_table = count_table.fillna(value=0)
    if spike_type == "ERCC": 
        count_file = os.path.join(count_dir, "CountTable_spike.txt")
        with open(count_file, 'w') as fh:
            count_table.to_csv(fh, sep="\t")
    else:
            
    # split into normalisation and size spike count tables
        is_size_spike = count_table.index.map(lambda x: x.startswith('SSPK_'))
        is_norm_spike = count_table.index.map(lambda x: x.startswith('NSPK_'))
        size_spike_table = count_table[is_size_spike]
        norm_spike_table = count_table[is_norm_spike]

    # save count table as tab delimited file
        size_count_file = os.path.join(count_dir, "CountTable_size_spike.txt")
        with open(size_count_file, 'w') as fh:
            size_spike_table.to_csv(fh, sep="\t")

        norm_count_file = os.path.join(count_dir, "CountTable_norm_spike.txt")
        with open(norm_count_file, 'w') as fh:
            norm_spike_table.to_csv(fh, sep="\t")


# subcommand used for counting number of reads
@cli.command()
@click.option('--fq-dir',
              type=click.Path(exists=True, resolve_path=True),
              help='directory with fastq files')
@click.option('--count-dir',
              type=click.Path(),
              help='directory to save count table in')
def count_fq_reads(fq_dir, count_dir):
    """
    Count reads in fastq files and save csv file with counts per sample.
    """
    fastqs = [f for f in os.listdir(fq_dir) if f.endswith(".fastq")]
    with open(os.path.join(count_dir, "total_reads.csv"), 'w') as fout:
        csv_writer = csv.writer(fout, delimiter="\t")
        for f in fastqs:
            fq = SeqIO.parse(os.path.join(fq_dir, f), "fastq")
            sample_name = f.split(".fastq")[0]
            count = sum([1 for i in fq])
            csv_writer.writerow([sample_name, count])


@cli.command()
@click.option('--input',
              help="fastq file to trim",
              type=click.File('rt'))
@click.option('--output',
              help='filtered reads file name',
              type=click.File('wt'))
@click.option('--min-len',
              type=int,
              help='min lenght of the read')
@click.option('--report-file-name',
              help="File name to save the results",
              type=str) 
@click.option('--sample-id',
              help="Sample ID to put in the report",
              type=str)              
@click.option('--report-stats',
              type=bool)                           
def filter_short_reads(input, output, min_len, report_file_name="report.txt", report_stats=False, sample_id=''):
    """
    Removes reads shorter than min-len from fastq file
    Reports the number and percentage of reads that pass length filter 
    """
    #reads = [r for r in SeqIO.parse(input, 'fastq') if len(r) >= min_len]
    allReads = list(SeqIO.parse(input, 'fastq'))
    reads =  [r for r in allReads if len(r) >= min_len]
    # if reportStats is set TRUE, calculate and print filter related statistics
    if report_stats:
        allReadsCount = len(allReads)
        filteredReadsCount = len(reads)
        percentAccepted = filteredReadsCount/allReadsCount
        reportArray = ['Filtering on Sample: ' + sample_id,
                       'Total Number of Reads: ' + str(allReadsCount),
                       '# of Reads Longer Than ' + str(min_len) + ': ' + str(filteredReadsCount),
                       'Percent of Reads Longer Than ' + str(min_len) + ': ' + str(percentAccepted*100) + '%'
                       ]
        reportString = '\n'.join(reportArray)
        build_report('Filter Short Reads of '+sample_id, report_file_name, reportString)
    SeqIO.write(reads, output, format='fastq')
   
@cli.command()
@click.option('--report_header',
              help="Header describing the lines to follow",
              type=str)
@click.option('--report_file_name',
              help="File name to save the results",
              type=str)       
# call the build report function from the command line
def call_build_report(report_header, report_file_name):
    """
    Calls the build_report function 
    """              
    build_report(report_header, report_file_name)


# A function for appending text to a report file          
def build_report(report_header, report_file_name, string='empty'):
    """
    Appends a string to a report file meant to keep track of the results
    """    
    # if file does not exist, create it 
    if os.path.isfile(report_file_name)==False:
        file = open(report_file_name, 'w')
        file.write('In this file a report on pipeline run is stored \n')
        file.close()
    # if no string to write is provided to the function, take one from the stdin
    if string=="empty":
        string=get_input(sys.stdin)
        
    file = open(report_file_name, 'a')
    textArray = [ '\n ################## &&&&&&&&&&&&&&&&& ##################',
                 '>>>> '+report_header+'\n',
                 string+'\n'
    ]    
    file.write('\n'.join(textArray))
    file.close()

@cli.command()
def test_click():
    build_report("this is the header", "example.txt")

@cli.command()
@click.option('--dat',
              help="mirna annotation in embl format",
              type=click.File('rt'))
@click.option('--org',
              help="organism (ie. dre, dme, cel)",
              type=str)
@click.option('--gff',
              help="output gff file",
              type=click.File('wt'))
def embl2gff(dat, org, gff):
    """
    Parse embl file and estract mature miRNA location information.
    """
    # extract records
    dat_parser = SeqIO.parse(dat, "embl")
    # extract organism specific miRNAs
    org_mirnas = [mirna for mirna in dat_parser if mirna.name.startswith(org)]
    for mirna in org_mirnas:
        mirna.id = mirna.name
    GFF.write(org_mirnas, gff)


@cli.command()
@click.option('--dat',
              help="mirna annotation in embl format",
              type=click.File('rt'))
@click.option('--bamfile',
              help="sorted and indexed alignments in bam format",
              type=click.STRING)
@click.option('--out',
              help="output counts file",
              type=click.File('wt'))
def count_miRNAs(dat, bamfile, out):
    """
    Count mirRNAs from alignments to pre-miRNAs.

    Read alignment manipulation is done with pysam.
    Documentation of pysam can be found:
    http://pysam.readthedocs.org/en/latest/api.html
    """

    def make_group_aln_func(hairpin):
        """
        Create a function that will group alignments.

        Alignments are asigned into two groups:
        'mature': those overlaping with mature miRNA
        'pre': those mapped to hairpin outside mature annotations
        """
        mature_mirnas = hairpin.features
        # build the mature ranges list
        # check that feature is of "miRNA" type 
        # features of type "modified base" mess the script
        mature_ranges = [(int(mirna.location.start), int(mirna.location.end))
                         for mirna in mature_mirnas if mirna.type=='miRNA' ]
                         
            
        def group_aln(aln):
            # First check if it falls inside the the mature ranges
            # So loop over the mature ranges (a hairpin can have 1 or 2 matures)
            # If read is inside either of the matures, 
            # "return" will break the function execution (and "for" loop)
            # and it will never reach the premature quality check block
            for mr in mature_ranges:
                is_in = aln.pos >= mr[0]-3 and aln.pos < mr[1]
                extend_3 = (aln.pos + aln.inferred_length) < (mr[1] + 3)
                if (is_in and extend_3):
                    # mature miRNA
                    if is_quality_mature(aln):
                        return('mature')
                    else:
                        return('low_qual')
                        
            # If it is not inside neither of the mature ranges, it must be pre-miRNA
            # Check the quality and return accordingly            
            if is_quality_pre(aln):
                return('pre')
            else:
                return('low_qual')              
                
        return(group_aln)

    def is_quality_mature(aln):
        """
        Discard reads aligning to mature miRNAs
        according to specific rules.
        """
        # allow 3' softcliping but not 5'
        # sofcliping is a tupple (4,...) in a cigar tupple list
        # check if it is as a first typple of the alignment cigar
        cig = aln.cigar
        if cig[0][0] == 4:
            return(False)

        # drop alignments with more than 3 nucleotides softclipped on 3'
        if (cig[-1][0] == 4 and cig[-1][1] > 3):
            return(False)
            
        # allow for 10%  of mismatches + indels
        # NM tag defines the edit distance, that is the number of substitutions, 
        # insertions and deletions summed up   
        n_mismatch = [tag[1] for tag in aln.tags if tag[0] == 'NM'][0]
        if (n_mismatch > 0.1*aln.rlen):
            return(False)
        return(True)

    def is_quality_pre(aln):
        """
        Discard reads aligned to hairpin
        according to specific rules.
        """
        cig = aln.cigar
        # allow for 10%  of mismatches + indels
        # NM tag defines the edit distance, that is the number of substitutions, 
        # insertions and deletions summed up 
        n_mismatch = [tag[1] for tag in aln.tags if tag[0] == 'NM'][0]
        n_softclip = sum([c[1] for c in cig if c[0] == 4])
        if (n_mismatch + n_softclip > 0.1*aln.rlen):
            return(False)
        return(True)

    def make_count_alns_func(hairpin):
        def count_alns(countlist, group_aln):
            mature_mirnas = [feature for feature in hairpin.features]
            group, g_aln = group_aln
            for aln in g_aln:
                start = aln.pos
                end = aln.pos+1 + aln.inferred_length
                length = aln.inferred_length
                if group == 'pre':
                    # +1 to start so that we fix the coordinates 
                    # (sam standards count form 1, while python from 0)
                    seqid = "%s_%s_%d_%d_%d" % (hairpin.name,
                                                "pre", start+1, end,
                                                length)
                    countlist += [seqid]
                elif group == 'mature':
                    # check to which mature miRNA the read aligns
                    # +3 to the starting position to compensate for the flexibility we allow 
                    # when we decide if it is mature or not
                    # otherwise it will not fit into the coordinates of the reference
                    # Further, added condition type=="miRNA" to ensure only this type is put in the list 
                    mirna = [mirna for mirna in mature_mirnas if (start+3 >= int(mirna.location.start) and start < int(mirna.location.end) and mirna.type=='miRNA')][0]
                    name = mirna.qualifiers['product'][0]
                    # +1 to start so that we fix the coordinates 
                    # (sam standards count form 1, while python from 0)
                    seqid = "%s_%s_%d_%d_%d" % (hairpin.name,
                                                name, start+1, end,
                                                length)
                    countlist += [seqid]
                    # for the case of annotation issues
                elif group == 'annotationIssues':
                    seqid = "%s_%s_%d_%d_%d" % (hairpin.name,
                                                "annotationProblem!!!", start+1, end,
                                                length) 
                    countlist += [seqid]
                else:
                    continue
            return(countlist)
        return(count_alns)

    # read in miRNA annotations from embl formated dat file
    annots = SeqIO.parse(dat, 'embl')
    alignments = pysam.Samfile(bamfile, 'rb')
    references = alignments.references

    # fetch annotation for reference sequences
    # prsent in the bam file
    # for some reason name is valid seqid fo the sequence
    hairpin_annots = [a for a in annots if a.name in references]

    # loop over all SeqRecords in annotations
    counts_tupples = []
    for hairpin in hairpin_annots:
        # mature miRNAs are the features of hairpin SeqRecords
        # create grouping function which knows miRNA annotations
        group_aln = make_group_aln_func(hairpin)
        # get reads aligning to this hairpin
        hairpin_alns = alignments.fetch(hairpin.name)
        # group alignments into those mapped to mature miRNA
        # and those mapping autside (pre-miRNA)
        grouped_alns = groupby(hairpin_alns, group_aln)
        # count
        count_alns = make_count_alns_func(hairpin)
        counts_list = reduce(count_alns, grouped_alns, [])
        # cont same occurences of seqid
        counts_dict = Counter(counts_list)
        # store them as list of tupples, because it's difficult
        # to join dictionaries
        counts_tupples += [i for i in counts_dict.items()]

    # create dir for counts if needed
    counts_dir = os.path.dirname(out.name)
    if not os.path.exists(counts_dir):
        os.makedirs(counts_dir)

    # export counts as csv file
    csv_writer = csv.writer(out)
    for row in counts_tupples:
        csv_writer.writerow(row)


@cli.command()
@click.option('--dir',
              help="count tables directory",
              type=click.STRING)
@click.option('--suffix',
              help="common suffix of count files, used to extract smapleids",
              type=click.STRING)
@click.option('--out',
              help="output counts file",
              type=click.File('wt'))
def merge_count_tables(dir, suffix, out):
    """
    Merge count tables on seqid's into single table.
    """
    # list files in dir
    count_files = [f for f in os.listdir(dir) if f.endswith(".csv")]
    samples = [f.split(suffix)[0] for f in count_files]
    dfs = [pd.read_csv(os.path.join(dir, f), index_col=0) for f in count_files]
    merged_dfs = reduce(lambda df1, df2: pd.merge(df1, df2,
                                                  left_index=True,
                                                  right_index=True,
                                                  how="outer"), dfs)
    merged_dfs = merged_dfs.fillna(0)
    merged_dfs.columns = samples
    merged_dfs.to_csv(out, sep="\t")

# With this function, I get the stdin
def get_input(stdin):
    reportString = ""
    for line in iter(stdin.readline, ''):
        reportString+=line
    stdin.close()
    return(reportString)

if __name__ == '__main__':
    cli()
