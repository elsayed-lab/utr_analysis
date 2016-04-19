"""
Helper function:

    - run_command
    - gzip_str
    - get_next_file_name
    - setup_loggers
    - num_lines

"""
import logging
import subprocess
import glob
import gzip
import re
import sys
import os

def run_command(cmd, log_handle, wait=True):
    """Runs a command and logs the output to a specified log handle"""
    log_handle.info("# " + cmd)

    # asynchronous
    if wait is False:
        process = subprocess.Popen(cmd.split())
        return 0

    # synchronous
    process = subprocess.Popen(cmd.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stdout:
        log_handle.info(stdout)
    if stderr:
        log_handle.error(stderr)

    return process.returncode


def gzip_str(filepath, strbuffer, chunk_size=65536):
    """Takes a StringIO buffer and writes the output to a gzip-compressed
    file."""
    # go to beginning of string buffer
    strbuffer.seek(0)

    # output path
    if filepath.endswith('.gz'):
        outfile = filepath
    else:
        outfile = filepath + '.gz'

    # write contents to a gzip-compressed file
    fp = gzip.open(outfile, 'wb')

    # to avoid overflow errors, we will read from the stream in chunks
    contents = strbuffer.read(chunk_size)

    while contents != '':
        fp.write(contents)
        contents = strbuffer.read(chunk_size)

    fp.close()

def get_next_file_name(base_name):
    """Returns a filepath for the next highest file number, e.g.
       file.5.log."""
    if not os.path.exists(base_name):
        return base_name
    else:
        file_nums = [int(x.split('.').pop())
                        for x in glob.glob("%s.*" % base_name)]
        next_file_num = max([0] + file_nums) + 1
        return "%s.%d" % (base_name, next_file_num)

def setup_loggers(root_dir, build_dirs, sample_ids):
    """Sets up master and task-specific loggers"""
    # setup master logger
    log_format = '%(asctime)s %(message)s'
    date_format = '%Y-%m-%d %I:%M:%S %p'
    formatter = logging.Formatter(log_format, datefmt=date_format)

    # determine log name to use
    master_log = get_next_file_name(os.path.join(root_dir, 'build.log'))

    logging.basicConfig(filename=master_log, level=logging.INFO,
                        format=log_format, datefmt=date_format)

    # log to console as well
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)

    # version information
    from ruffus import __version__ as ruffus_version
    from Bio import __version__ as biopython_version

    logging.info("# Starting UTR Analysis")
    logging.info("# Python %s" % sys.version)
    logging.info("# Ruffus %s" % ruffus_version)
    logging.info("# Biopython %s" % biopython_version)
    logging.info("# Command:\n%s" % " ".join(sys.argv))

    # create dictionary of log handlers for sample-specific info
    loggers = {}

    # add previously processed samples
    input_globstr = (
        '%s/*/tophat/*_sl_reads/accepted_hits_sorted.bam' % build_dirs['sl']
    )
    for filepath in glob.glob(input_globstr):
        # get sample id
        sample_ids.append(re.match('.*/(.*)/.*', filepath).groups()[0])

    # setup sample-specific loggers
    for sample_id in sample_ids:
        loggers[sample_id] = {}

        for analysis in ['sl', 'rsl', 'polya', 'polyt']:
            loggers[sample_id][analysis] = {}

            for read_num in ['1', '2']:
                if analysis == 'sl':
                    bdir = build_dirs['sl']
                elif analysis == 'rsl':
                    bdir = build_dirs['rsl']
                elif analysis == 'polya':
                    bdir = build_dirs['polya']
                else:
                    bdir = build_dirs['polyt']

                sample_log_name = get_next_file_name(
                    os.path.join(bdir, sample_id, 'log', '%s_%s_%s.log' % (
                        sample_id, analysis, read_num
                    ))
                )
                loggers[sample_id][analysis][read_num] = logging.getLogger(
                    sample_id + analysis + read_num
                )
                handler = logging.FileHandler(sample_log_name)
                handler.setFormatter(formatter)
                loggers[sample_id][analysis][read_num].addHandler(handler)

    return loggers

def num_lines(filepath):
    """Returns the number of lines in a specified file"""
    if filepath.endswith('.gz'):
        fp = gzip.open(filepath, 'rb')
    elif filepath.endswith('.xz'):
        import backports.lzma as lzma
        fp = lzma.open(filepath, 'rb')
    else:
        fp = open(filepath)

    # count number of lines
    for i, line in enumerate(fp, 1):
        pass

    fp.close()
    return i

    return loggers

