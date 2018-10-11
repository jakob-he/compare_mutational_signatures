#!/usr/bin/python3

import pandas as pd
import argparse
import numpy as np
import datetime
import sys

def parse_arguments():
    '''function that parses command line arguments'''
    parser = argparse.ArgumentParser(description = 'Convert a txt file containing mutations to a VCF file')
    parser.add_argument('--mutations', type = str, help = 'path to the mutation txt file', required = True)
    parser.add_argument('--output', type = str, help = 'path where the final VCF is saved', required = True)
    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()

    #read data into dataframe
    df = pd.read_csv(args.mutations, sep = '\t', header = None, comment = '#')
    #extract columns needed for mutation position format
    df.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    df = df[['FILTER','#CHROM','POS','REF','ALT']]
    #add 'chr'
    df['#CHROM'] = [f"chr{chromsome}" for chromsome in df['#CHROM']]
    with open(f'{args.output}final_subs_mutation_pos.txt','w') as output:
        df.to_csv(output,sep = '\t', index = False, header = False)


if __name__ == '__main__':
    main()
