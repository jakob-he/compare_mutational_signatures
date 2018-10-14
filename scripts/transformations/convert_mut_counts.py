#!/usr/bin/python3

import pandas as pd
import argparse
import numpy as np
import datetime
import os
import re


def parse_arguments():
    '''function that parses command line arguments'''
    parser = argparse.ArgumentParser(description = 'Convert a txt file containing mutations to a VCF file')
    parser.add_argument('--countmatrix', type = str, help = 'path to the mutation count file')
    parser.add_argument('--matrixdir', type = str, help = 'directory of mutation count data')
    args = parser.parse_args()
    return args



def transform_mut_counts(matrixfile):
    data = pd.read_csv(matrixfile, sep='\t').transpose()
    data.columns = data.iloc[0]
    data.drop(data.index[0],inplace=True)
    new_format = []
    for mutation in data.columns:
        motif = re.search(r"\[(\w\>\w)\]", mutation)
        substitution = motif.group(1)
        new_format.append(f'{substitution}:{mutation[motif.start()-1]}{substitution[0]}{mutation[motif.end()]}')
    data.columns = new_format
    return data

def save_data(path,data):
    with open(path,'w') as output:
        data.to_csv(output,index = True)

def main():
    args = parse_arguments()

    if args.matrixdir:
        files = os.listdir(args.matrixdir)
        for countmatrix in files:
            if countmatrix.endswith(".txt"):
                transformed_data = transform_mut_counts(os.path.join(args.matrixdir,countmatrix))
                save_data(os.path.basename(countmatrix),transformed_data)

    else:
        transformed_data = transform_mut_counts(args.countmatrix)
        save_data(os.path.basename(args.countmatrix),transformed_data)


if __name__ == '__main__':
    main()
