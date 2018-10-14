#!/usr/bin/python3

import pandas as pd
import argparse
import numpy as np
import datetime
import os

def parse_arguments():
    '''function that parses command line arguments'''
    parser = argparse.ArgumentParser(description = 'Convert a txt file containing mutations to a VCF file')
    parser.add_argument('--mutations', type = str, help = 'path to the mutation txt file')
    parser.add_argument('--mutationdir', type = str, help = 'directory of mutation data')
    args = parser.parse_args()
    return args

def transform_mut_pos(filepath: str):
    #read data into dataframe
    df = pd.read_csv(filepath, sep = '\t', header = None)

    df.columns = ['INFO','#CHROM','POS','REF','ALT']
    qual = np.repeat(".",df.shape[0])
    df['QUAL'] = qual
    df['ID'] = qual
    df['FILTER'] = qual

    #sort VCF by chromsome
    df['#CHROM'] = [chromosome.replace("chr","") for chromosome in df['#CHROM']]
    df['INFO'] = ["SID="+sampleid for sampleid in df['INFO']]
    df.replace({'X':23,'Y':24, 'MT':25}, inplace = True)
    df['#CHROM'] = pd.to_numeric(df['#CHROM'])
    df.sort_values(by = ['#CHROM'], inplace = True)
    df['#CHROM'].replace({23:'X',24:'Y', 25:'MT'}, inplace = True)
    df['#CHROM'] = ["chr"+str(chromosome) for chromosome in df['#CHROM']]

    #rearrange dataframe
    df = df[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']]

    #add VCF header
    now = datetime.datetime.now()
    with open(f'{os.path.basename(filepath).split(".")[0]}.vcf','w') as output:
        output.write(f'##fileformat=VCFv4.1\n##fileDate={now.strftime("%Y%m%d")}\n##source=TxttoVCFconversion\n##reference=hg19\n##INFO=<ID=SID,Number=1,Type=String,Description="Sample ID\n')
        df.to_csv(output,sep = '\t', index = False)


def main():
    args = parse_arguments()

    if args.mutationdir:
        files = os.listdir(args.mutationdir)
        for mutationdata in files:
            if mutationdata.endswith(".txt"):
                transform_mut_pos(os.path.join(args.mutationdir,mutationdata))

    else:
        transform_mut_pos(args.mutations)


if __name__ == '__main__':
    main()
