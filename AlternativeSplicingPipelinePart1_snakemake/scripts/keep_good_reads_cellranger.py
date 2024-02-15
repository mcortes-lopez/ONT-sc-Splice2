import numpy as np
import pysam
import pandas as pd
import os
import argparse
import subprocess

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group('required arguments')
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    parser.add_argument("-c", "--cellranger", help="raw counts from cellranger", required=True)
    parser.add_argument("-th", "--threshold", help="threshold on the number of UMIs per gene", required=True)
    parser.add_argument("-thi", "--threshold_intergenic", help="threshold on the mapping quality to be kept", required=True) 
    args = par.parse_args()
    return args

def set_umi_tag(be, bs, seq):
    if bs < be:
        umi = seq[be:be + 12]
    else:
        umi = seq[be-12:be]
    return(umi)

args = parser_arguments()
input_file = args.input
output_file = args.output
cellranger_file = args.cellranger
thresh = args.threshold
thresh_intergenic = args.threshold_intergenic

good_bc_file = pysam.AlignmentFile(input_file, 'rb')
good_reads_file = pysam.AlignmentFile(output_file, 'wb', template=good_bc_file)

illumina = pd.read_csv(cellranger_file, index_col=[0])
print(illumina.head())

all_reads = 0
bad_reads = 0
reads_saved = 0

for read in good_bc_file.fetch():
    # bc = read.get_tag('BC')
    be = read.get_tag('BE')
    bs = read.get_tag('BB')
    seq = read.get_forward_sequence()
    all_reads += 1
    try :
        gene = read.get_tag('GE')
        res = illumina.loc[gene,:].sum()
    except KeyError as err:
        if float(read.mapping_quality) >= float(thresh_intergenic):
            umi = set_umi_tag(be,bs,seq)
            read.set_tag('U8', umi) 
            good_reads_file.write(read)
            reads_saved +=1
        else:
            bad_reads += 1
        continue
    if res >= int(thresh) and float(read.mapping_quality) > 0:
            umi = set_umi_tag(be, bs, seq)
            read.set_tag('U8', umi)
            good_reads_file.write(read)
            reads_saved += 1

print(f'We had {bad_reads} bad reads and we saved {reads_saved} reads over {all_reads}.')
good_reads_file.close()
good_bc_file.close()

