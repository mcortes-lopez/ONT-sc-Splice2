import os
import argparse

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group('required arguments')
    parser.add_argument("-i", "--input", help="directory/path to input sub fastq files", required=True)
    parser.add_argument("-o", "--output", help="outputdir", required=True)
    args = par.parse_args()
    return args


args = parser_arguments()
input_dir = args.input
print(input_dir)
output_dir = args.output

list_files = os.listdir(input_dir)
for f in list_files:
    if os.path.getsize(input_dir + '/' +  f) >0:
        os.rename(input_dir + '/' + f, output_dir + '/' + f)
        print('File successfuly moved')
