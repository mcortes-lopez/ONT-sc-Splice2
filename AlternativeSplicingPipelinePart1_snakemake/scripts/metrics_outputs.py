import numpy as np
import subprocess
import os
import argparse
import pysam
import pandas as pd 

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group('required arguments')
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    args = par.parse_args()
    return args

args = parser_arguments()
input_file = args.input
output_file = args.output

samfile = pysam.AlignmentFile(input_file, 'rb')
comp_df = pd.DataFrame(columns=['tag', 'gene', 'fragment_length', 'mapq', 'tso_found'])

for i, read in enumerate(samfile.fetch()):
    seq = read.get_forward_sequence()
    mapq = read.mapping_quality
    # get type of tag
    tag_list = [tag[0] for tag in read.tags]
    if 'BT' in tag_list:
        tag = 'BT'
    elif 'BF' in tag_list:
        tag = 'BF'
    else:
        tag = 'BR'
    if 'TN' in tag_list:
        tso = 'no'
    else:
        tso = 'yes'
    try:
        gene = read.get_tag('XF')
        comp_df.loc[i, :] = [tag, gene, len(seq), mapq, tso]
    except KeyError as err:
        comp_df.loc[i, :] = [tag, 'none', len(seq), mapq, tso]

samfile.close()

comp_df.to_csv(output_file)
