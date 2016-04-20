"""
RNA-Seq mapping related functions:

    - run_tophat
    - run_bam2fastx
    - sort_and_index
    - filter_mapped_reads
    - map_reads

"""
import os
import sys
from util import run_command, num_lines

def run_tophat(output_dir, genome, log_handle, r1, r2="", gff=None,
               num_threads=1, read_mismatches=2, max_multihits=20,
               sort_by_name=False, extra_args=""):
    """
    Uses Tophat to map reads with the specified settings.

    Several useful Tophat parameters have been included as keyword arguments,
    along with their default settings for Tophat 2.

    The documentation for these and other useful parameters for this script are
    included below for convenience.

    -N/--read-mismatches
        Final read alignments having more than these many mismatches are
        discarded. The default is 2.
        
    --no-novel-juncs
        Only look for reads across junctions indicated in the supplied GFF or 
        junctions file. (ignored without -G/-j)
        
    -G/--GTF <GTF/GFF3 file>	
        Supply TopHat with a set of gene model annotations and/or known
        transcripts, as a GTF 2.2 or GFF3 formatted file. If this option is
        provided, TopHat will first extract the transcript sequences and use
        Bowtie to align reads to this virtual transcriptome first. Only the
        reads that do not fully map to the transcriptome will then be mapped on
        the genome. 

    -x/--transcriptome-max-hits
        Maximum number of mappings allowed for a read, when aligned to the
        transcriptome (any reads found with more then this number of mappings
        will be discarded).

    -p/--num-threads <int>
        Use this many threads to align reads. The default is 1.

    -g/--max-multihits <int>
        Instructs TopHat to allow up to this many alignments to the reference
        for a given read, and choose the alignments based on their alignment
        scores if there are more than this number. The default is 20 for read
        mapping. Unless you use --report-secondary-alignments, TopHat will
        report the alignments with the best alignment score. If there are more
        alignments with the same score than this number, TopHat will randomly
        report only this many alignments. In case of using
        --report-secondary-alignments, TopHat will try to report alignments up
        to this option value, and TopHat may randomly output some of the
        alignments with the same score to meet this number.
    """
    # create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # FASTA filename without extension
    genome_basename = os.path.splitext(genome)[0]

    # @NOTE 2014/02/06 -- Tophat 2.0.10 fails for some reads when attempting
    # to use more than one thread. For now, just perform the mapping
    # using a single-thread to be safe...

    # build command
    if gff is not None:
        cmd = 'tophat --no-novel-juncs -G %s ' % gff
    else:
        cmd = 'tophat '

    cmd += "--num-threads %d --max-multihits %d --read-mismatches %d %s -o %s %s %s %s" % (
           num_threads, max_multihits, read_mismatches, extra_args,
           output_dir, genome_basename, r1, r2)

    # run tophat
    ret = run_command(cmd, log_handle)

    # check to see if tophat succeeded and stop execution otherwise
    if ret != 0:
        log_handle.error("# Error running tophat (%s)!" % genome)
        print("# Error running tophat (%s)!" % genome)
        sys.exit()

    # sort and index bam output using samtools
    log_handle.info("# Sorting and indexing Tophat output")
    sort_and_index(os.path.join(output_dir, 'accepted_hits'), log_handle,
                   sort_by_name)
    sort_and_index(os.path.join(output_dir, 'unmapped'), log_handle, sort_by_name)

    return 0

def sort_and_index(base_output, log_handle, sort_by_name=False):
    """Sorts and indexes .bam files using samtools.
   
    Arguments
    ---------
    base_output: str
        Output filepath without the file extension portion.
    log_handle: logging.Handle
        Handler to use for logging.
    """
    # sort bam
    sort_cmd = 'samtools sort %s %s -o %s' % (
        "-n" if sort_by_name else "",
        base_output + ".bam", base_output + "_sorted.bam")
    run_command(sort_cmd, log_handle)

    # delete unsorted verion
    #os.remove(base_output + ".bam")

    # index bam
    index_cmd = 'samtools index %s' % (base_output + '_sorted.bam')

    # 2014/04/01
    # samtools index sometimes gets stuck during execution, even after the
    # indexing has finished; since it isn't necessary for downstream processes
    # the indexing will be done asynchronously for now.
    run_command(index_cmd, log_handle, wait=False)
    log_handle.info("# Done sorting and indexing")

def run_bam2fastx(bam, fastq, log_handle):
    """
    Uses the tophat bam2fastx tool to convert a bam file to fastq

    Note that bam2fastx will detect the .gz extension at the end of the output
    filename and generate gzip-compressed output.

    Arguments
    ---------
    bam: str
        Input bam file
    fastq: str
        Output fastq file
    log_handle: logging.Handle
        Handler to use for logging.
    """
    # if output is expected to be xz-compressed, first generate uncompressed
    # fastq output
    fastq = fastq.replace('.xz', '')

    #bam2fastx parameters:
    #    -q fastq
    #    -A all reads
    #    -P pair-end
    #    -o output filepath
    cmd = "bam2fastx -q -A -P -o %s %s" % (fastq, bam)

    return run_command(cmd, log_handle)

def filter_mapped_reads(r1, r2, genome, tophat_dir, output_fastq, log_handle,
                        gff=None, read_mismatches=2, num_threads_tophat=1):
    """
    Maps reads using tophat and discards any that align to the genome.

    Arguments
    ---------
    r1: str
        Filepath to R1 reads to be filtered.
    r2: str
        Filepath to R2 reads to be filtered.
    genome: str
        Filepath to genome FASTA to map reads against
    tophat_dir: str
        Directory to save tophat output to
    output_fastq: str
        Filepath to save filtered FASTQ files to
    log_handle: logging.Handle
        Handler to use for logging.
    gff: string
        Filepath to GFF to use during mapping (optional)
    read_mismatches: int
        Number of mismatches to allow when scanning with Tophat.
    num_threads_tophat: int
        Number of threads to use when running Tophat (default: 1)
    """
    # bam / fastq filepaths
    bam_input = os.path.join(tophat_dir, 'unmapped_sorted.bam')

    # check to see if tophat output already exists
    if not os.path.exists(bam_input):
        # map reads to genome using tophat
        log_handle.info("# Mapping against %s" % os.path.basename(genome))

        tophat_args = '--no-mixed --no-discordant'

        ret = run_tophat(tophat_dir, genome, log_handle, r1, r2,
                         read_mismatches=read_mismatches, gff=gff,
                         num_threads=num_threads_tophat, 
                         sort_by_name=True,
                         extra_args=tophat_args)

        # number of reads before filtering
        num_reads_total = num_lines(r1) / 4
        num_unmapped = num_lines(bam_input) / 2

        # delete uneeded accepted_hits files
        #os.remove(os.path.join(tophat_dir, "accepted_hits_sorted.bam"))

        log_handle.info(
            "# Ignoring %d reads which mapped to specified genome (%d remaining)" %
            (num_reads_total - num_unmapped, num_unmapped)
        )
    else:
        log_handle.info("# Skipping %s: output already exists" % tophat_dir)

    # convert remaining unmapped reads from bam to fastq
    r1_fastq = output_fastq.replace('removed.fastq', 'removed.1.fastq')
    r2_fastq = output_fastq.replace('removed.fastq', 'removed.2.fastq')

    if not os.path.exists(r1_fastq):
        # run bam2fastx to convert unmapped bam file to fastq
        ret = run_bam2fastx(bam_input, output_fastq, log_handle)

        if ret != 0:
            log_handle.error("# Error running bam2fastx (%s)!" % genome)
            print("# Error running bam2fastx (%s)!" % genome)
            sys.exit()

        # in the case of xz-compressed files, uncompressed .fastq files
        # are outputted from the above command so we will now compress those
        if output_fastq.endswith('.xz'):
            uncompressed_r1 = r1_fastq.replace('.xz', '')
            uncompressed_r2 = r2_fastq.replace('.xz', '')
            ret = run_command('xz %s' % uncompressed_r1, log_handle)
            ret = run_command('xz %s' % uncompressed_r2, log_handle)

def map_reads(feature_name, build_dir, sample_id, read_num, num_threads_tophat, log_handle):
    """
    Maps the filtered reads back to the genome. If the trimmed version of
    the read successfully maps this is indicative of a possible trans-spliced
    or polyadenylated read.

    Arguments
    ---------
    feature_name: str
        Name of UTR feature.
    build_dir: str
        Base diretory to save mapped reads to.
    sample_id: str
        Name of sample currently processing.
    read_num: str
        Read number currently being processed.
    num_threads_tophat: int
        Number of threads to use when running Tophat.
    log_handle: logging.Handle
        Handler to use for logging.
    """
    output_dir = '%s/%s/tophat/%s_filtered_trimmed' % (
        build_dir, sample_id, read_num
    )
    log_handle.info("# Mapping filtered reads back to genome")

    # input read base directory
    basedir = '%s/%s/fastq' % (build_dir, sample_id)

    # R1 input filepath (including matched sequence)
    if read_num == '1':
        r2_suffix = '2'

        r1_filepath = (
            '%s/unfiltered/%s_%s_1_%s_trimmed.fastq.gz' %
            (basedir, sample_id, read_num, feature_name)
        )

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('1_%s_untrimmed' % feature_name, '2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2 input filepath
    else:
        r1_suffix = '1'
        r2_filepath = (
            '%s/unfiltered/%s_%s_2_%s_trimmed.fastq.gz' %
            (basedir, sample_id, read_num, feature_name)
        )
        r1_filepath = r2_filepath.replace('2_%s_untrimmed' % feature_name, '1')

    # Map reads using Tophat
    tophat_args = '--transcriptome-max-hits 1 --no-mixed --no-discordant'

    ret = run_tophat(output_dir, args.target_genome, log_handle, r1_filepath,
                     r2_filepath, gff=args.target_gff, max_multihits=1,
                     num_threads=num_threads_tophat,
                     extra_args=tophat_args)

    log_handle.info("# Finished mapping hits to genome")

