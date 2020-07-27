#!/usr/bin/env python
# coding: utf-8

import sys, os
from glob import glob
import subprocess
import csv
import re

import pandas as pd


def get_fpath(dpath, globpattern):
    flist = glob(os.path.join(dpath, globpattern))
    if len(flist) != 1:
        raise Exception(f'get_fpath: more than one file matches {globpattern}: {flist}')
    return flist[0]


def read_log(fpath):
    fastq1 = None
    fastq2 = None
    snps = None
    reads1 = None
    reads2 = None
    unmappedS = None
    unmapped0 = None
    unmapped1 = None
    unmapped2 = None
    unmappedreadsS = None
    unmappedreads0 = None
    unmappedreads1 = None
    unmappedreads2 = None

    with open(fpath) as fp:
        for cnt, line in enumerate(fp):
            m = re.search(r'/bin/snippy.*--pe1\s*(\S+).*--pe2\s*(\S+)')
            if m is not None:
                fastq1 = m.group(1)
                fastq2 = m.group(2)
                reads1 = count_fastq_reads(fastq1)
                reads2 = count_fastq_reads(fastq2)

            # Converted 2 SNPs to TAB format.
            m = re.search(r'Converted (\d+) SNPs to TAB format.')
            if m is not None:
                snps = m.group(1)

            # ### samtools fastq -f 12 -v 20 --threads 3 -c 5 -N -s snp_1A3_all.unmapped_SE.fq.gz -0 snp_1A3_all.unmapped_R0.fq.gz -1 snp_1A3_all.unmapped_R1.fq.gz -2 snp_1A3_all.unmapped_R2.fq.gz snp_1A3_all.bam
            m = re.search(r'samtools fastq.*--s\s*(\S+).*--0\s*(\S+).*--1\s*(\S+).*--2\s*(\S+)')
            if m is not None:
                unmappedS = m.group(1)
                unmapped0 = m.group(2)
                unmapped1 = m.group(3)
                unmapped2 = m.group(2)
                unmappedreadsS = count_fastq_reads(unmappedS)
                unmappedreads0 = count_fastq_reads(unmapped0)
                unmappedreads1 = count_fastq_reads(unmapped1)
                unmappedreads2 = count_fastq_reads(unmapped2)


            # READ 829529 WRITTEN 821523
            # EXCLUDED 23426 EXAMINED 806103
            # PAIRED 801192 SINGLE 4911
            # DULPICATE PAIR 6018 DUPLICATE SINGLE 1988
            # DUPLICATE TOTAL 8006
    return {
        'fastq1' : fastq1,
        'fastq2' : fastq2,
        'snps' : snps,
        'reads_fastq1' : reads1,
        'reads_fastq2': reads2,
        'unmappedS' : unmappedS,
        'unmapped0' : unmapped0,
        'unmapped1' : unmapped1,
        'unmapped2' : unmapped2,
        'unmappedreadsS' : unmappedreadsS,
        'unmappedreads0' : unmappedreads0,
        'unmappedreads1': unmappedreads1,
        'unmappedreads2': unmappedreads2,
    }


def count_fastq_reads(fpath):
    catcmd = 'cat'
    if fpath.endswith('.gz'):
        catcmd = 'zcat'
    cmd = f'echo $({catcmd} {fpath} | wc -l) / 4 | bc'
    result = subprocess.run(cmd, shell=True,
                   check=True, capture_output=True, text=True, stderr=subprocess.STDOUT)
    print(cmd)
    print(result.stdout)
    return result.stdout

def read_bam(fpath):
    total_pass = None
    total_fail = None
    mapped_pass = None
    mapped_fail = None
    result = subprocess.run(['samtools', 'flagstat', fpath],
                   check=True, capture_output=True, text=True, stderr=subprocess.STDOUT)
    print(f'samtools flagstat {fpath}')
    print(result.stdout)
    # 821523 + 0 in total (QC-passed reads + QC-failed reads)
    m = re.search(r'(\d+) + (\d+) in total', result.stdout)
    if m is not None:
        total_pass = m.group(1)
        total_fail = m.group(2)
    # 798115 + 0 mapped (97.15% : N/A)
    m = re.search(r'(\d+) + (\d+) mapped', result.stdout)
    if m is not None:
        mapped_pass = m.group(1)
        mapped_fail = m.group(2)
    return {
        'total' : total_pass,
        'mapped' : mapped_pass,
        'total_fail' : total_fail,
        'mapped_fail': mapped_fail,
    }


    # 821523 + 0 in total (QC-passed reads + QC-failed reads)
    # 18 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 798115 + 0 mapped (97.15% : N/A)
    # 818655 + 0 paired in sequencing
    # 409272 + 0 read1
    # 409383 + 0 read2
    # 789480 + 0 properly paired (96.44% : N/A)
    # 795174 + 0 with itself and mate mapped
    # 163 + 0 singletons (0.02% : N/A)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models.')
    parser.add_argument('--snpdir', help='snippy dir', required=True)
    parser.add_argument('--out', help='output file', required=True)

    args = parser.parse_args()
    if not os.path.isdir(args.snpdir):
        raise (Exception(f'Directory not found: {args.snpdir}'))

    dpath = args.snpdir
    # need the following files:
    #report_fpath = get_fpath(dpath, '*.report.txt')
    log_fpath = get_fpath(dpath, '*.log')
    bam_fpath = get_fpath(dpath, '*.bam')

    result = read_log(log_fpath)
    result.update(read_bam(bam_fpath))
    result['snpdir'] = args.snpdir

    with open(args.out, 'w', newline='') as csvfile:
        fieldnames = result.keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(result)


# reference                      snp_1A3_all.log
# ref.fa                         snp_1A3_all.raw.vcf
# ref.fa.fai                     snp_1A3_all.report.txt
# snp_1A3_all.aligned.fa         snp_1A3_all.subs.vcf
# snp_1A3_all.bam                snp_1A3_all.tab
# snp_1A3_all.bam.bai            snp_1A3_all.txt
# snp_1A3_all.bed                snp_1A3_all.unmapped_R0.fq.gz
# snp_1A3_all.consensus.fa       snp_1A3_all.unmapped_R1.fq.gz
# snp_1A3_all.consensus.subs.fa  snp_1A3_all.unmapped_R2.fq.gz
# snp_1A3_all.csv                snp_1A3_all.unmapped_SE.fq.gz
# snp_1A3_all.filt.vcf           snp_1A3_all.vcf
# snp_1A3_all.gff                snp_1A3_all.vcf.gz
# snp_1A3_all.html               snp_1A3_all.vcf.gz.csi
