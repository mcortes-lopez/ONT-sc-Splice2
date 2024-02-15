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
    parser.add_argument("-th", "--thresh", help="threshold for filtering", required=True)
    args = par.parse_args()
    return args

def comp_dist(seq_bc, bc):
    return(np.sum([seq_bc[i] == bc[i] for i in range(min(len(bc), len(seq_bc)))]))

args = parser_arguments()
input_file = args.input
output_file = args.output
thresh = args.thresh

samfile = pysam.AlignmentFile(input_file, 'rb')
good_reads = pysam.AlignmentFile(output_file, 'wb', template=samfile)

comp_df = pd.DataFrame(columns=['tag', 'edit_distance_best_match', 'edit_distance_sec_match', 'gene', 'fragment length', 'mapping quality'])

for i, read in enumerate(samfile.fetch()):
    seq = read.get_forward_sequence()
    bc = read.get_tag('BC')
    bc_s = read.get_tag('BB')
    bc_e = read.get_tag('BE')
    quality = float(read.mapping_quality)
    ### BE and BB are 1-based
    if bc_s < bc_e:
        seq_bc = seq[bc_s - 1 : bc_e]
    else:
        seq_bc = seq[bc_e - 1 : bc_s][::-1]
    comp_read = read.get_tag('B1')
    try: 
        sec_read = read.get_tag('B2')
    except KeyError as err:
        sec_read = 'none'
    # comp_read = comp_dist(seq_bc, bc)
    # comp_read = comp(seq_bc, bc)
    # read.set_tag('MI', str(comp_read))
    # comp_read = comp(seq_bc, bc)
    if comp_read <= int(thresh):
        good_reads.write(read)
    # get type of tag
    tag_list = [tag[0] for tag in read.tags]
    if 'BT' in tag_list:
        tag = 'BT'
    elif 'BF' in tag_list:
        tag = 'BF'
    else:
        tag = 'BR'
    try:
        cat = read.get_tag('XF')
        comp_df.loc[i, :] = [tag, comp_read, sec_read, cat, len(seq), quality]
    except KeyError as err:
        comp_df.loc[i, :] = [tag, comp_read, sec_read, 'none', len(seq), quality]
    

good_reads.close()
samfile.close()

comp_df.to_csv(input_file[:-4] + '_stats.txt')
