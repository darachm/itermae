#!/usr/bin/env python3

# Tiny little script just to turn a fastq into a SAM file

import sys
import gzip
import argparse
from Bio import SeqIO

def format_sam_record(record_id, sequence, qualities, tags,
        flag='0', reference_name='*', 
        mapping_position='0', mapping_quality='255', cigar_string='*',
        reference_name_of_mate='=', position_of_mate='0', template_length='0'
    ):
    return "\t".join([
            record_id,
            flag,
            reference_name,
            mapping_position,
            mapping_quality,
            cigar_string,
            reference_name_of_mate,
            position_of_mate,
            template_length,
            sequence,
            qualities,
            tags
        ])


# Name and description of this program
parser = argparse.ArgumentParser("") 
parser.add_argument("--input",default="STDIN")
parser.add_argument("-z","--gzipped",action="store_true")
args = parser.parse_args()


if args.input == "STDIN":
    if args.gzipped:
        with gzip.open(sys.stdin,"rt") as input_file_gz:
            input_seqs = SeqIO.parse(input_file_gz,"fastq")
    else:
        input_seqs = SeqIO.parse(sys.stdin,"fastq")
else:
    if args.gzipped:
        with gzip.open(args.input,"rt") as input_file_gz:
            input_seqs = SeqIO.parse(input_file_gz,"fastq")
    else:
        input_seqs = SeqIO.parse(args.input_file,"fastq")

for i in input_seqs:
    print(
        format_sam_record(
            i.id, str(i.seq),
            ''.join( [chr(i+33) for i in i.letter_annotations['phred_quality']] ),
            ''
        )
    )


