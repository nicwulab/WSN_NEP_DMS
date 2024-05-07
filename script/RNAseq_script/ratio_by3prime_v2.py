import argparse
import pysam
from glob import glob
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
from functools import partial
import multiprocessing
def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence
def count_sequences_with_pattern(bam_file, cRNA_pattern,mRNA_pattern):

    cRNA_count = {}
    mRNA_count = {}
    # Open BAM file for reading
    # Open BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    idx=0
    for chrom in cRNA_pattern.keys():
    # Iterate over aligned reads
        for read in bamfile.fetch(chrom):
            if read.is_unmapped or read.is_secondary:
                continue
            if read.is_read1:
                seq = read.query_sequence
                name = read.query_name
                pattern = cRNA_pattern[chrom]
                mrna_pattern = mRNA_pattern[chrom]
                if pattern in seq:
                    cRNA_count.setdefault(chrom, set()).add(name)
                if mrna_pattern in seq:
                    mRNA_count.setdefault(chrom, set()).add(name)
        
            if read.is_read2:
                seq = read.query_sequence
                name = read.query_name
                pattern = reverse_complement(cRNA_pattern[chrom])
                mrna_pattern = reverse_complement(mRNA_pattern[chrom])
                if pattern in seq:
                    cRNA_count.setdefault(chrom, set()).add(name)
                if mrna_pattern in seq:
                    mRNA_count.setdefault(chrom, set()).add(name)
    bamfile.close()
    # Create a new dictionary with keys as the lengths of the sets
    mRNA_new = {inner_key: len(mRNA_count[inner_key]) for inner_key in mRNA_count.keys()}
    cRNA_new = {inner_key: len(cRNA_count[inner_key]) for inner_key in cRNA_count.keys()}
    
    name = bam_file.split('/')[-1].split('-')[0][:-9]
    # Convert each dictionary to a DataFrame
    temp_mRNA_df = pd.DataFrame([mRNA_new], index=[name])
    temp_mRNA_df.to_csv(f'data/tmp/{name}_mRNA.csv')
    temp_cRNA_df = pd.DataFrame([cRNA_new], index=[name])
    temp_cRNA_df.to_csv(f'data/tmp/{name}_cRNA.csv')
    
def main():
    parser = argparse.ArgumentParser(description="count FASTQ sequences containing a specified pattern.")
    parser.add_argument("-i","--input_file", default='result/*/*Aligned.sortedByCoord.out.bam', type=str, 
                        help="Input FASTQ file path")
    parser.add_argument("-o","--output_file", default='result/NEP_RNA_ratio.tsv', type=str, 
                        help="Input FASTQ file path")
    args = parser.parse_args()
    nprocs=64
    cRNA_pattern_dict={'LC333182.1':'GTCGAATAGTTTAAAAACGA','LC333183.1':'CTTGTCCTTCATGAAAAAATGC','LC333184.1':'CCATACTGTCCAAAAAAGTA','LC333185.1':'CAGAAATATAAGGAAAAACAC','LC333186.1':'GTACGACAATTAAAGAAAAATAC','LC333187.1':'AAGTAGTTTGTTCAAAAAACT','L25818.1':'ATAGAGCTGGAGTAAAAAACTA','CY034136.1':'GCTTATTTAATAATAAAAAACAC'}
    # mRNA_pattern_dict={'LC333182.1':'GTCGAATAGTTTAAAAAAAA','LC333183.1':'CTTGTCCTTCATGAAAAAAAAA','LC333184.1':'CCATACTGTCCAAAAAAAAA','LC333185.1':'CAGAAATATAAGGAAAAAAAA','LC333186.1':'GTACGACAATTAAAGAAAAAAAA','LC333187.1':'AAGTAGTTTGTTCAAAAAAAA','L25818.1':'ATAGAGCTGGAGTAAAAAAAAA','CY034136.1':'GCTTATTTAATAATAAAAAAAAA'}
    # cRNA_pattern_dict={'LC333182.1':'AAAAACGA','LC333183.1':'AAAAAATGC','LC333184.1':'AAAAAAGTA','LC333185.1':'AAAAACAC','LC333186.1':'AAAAATAC','LC333187.1':'AAAAAACT','L25818.1':'AAAAAACTA','CY034136.1':'AAAAAACAC'}
    mRNA_pattern_dict={'LC333182.1':'AAAAAAA','LC333183.1':'AAAAAAAA','LC333184.1':'AAAAAAAA','LC333185.1':'AAAAAAA','LC333186.1':'AAAAAAA','LC333187.1':'AAAAAAAA','L25818.1':'AAAAAAAA','CY034136.1':'AAAAAAAA'}
    files =glob(args.input_file)
    
    
    if nprocs > 4:
        pool = multiprocessing.Pool(processes=nprocs)
        pool.map(partial(count_sequences_with_pattern,cRNA_pattern=cRNA_pattern_dict,mRNA_pattern=mRNA_pattern_dict),
                 files)
    else:
        for file in files:
            count_sequences_with_pattern(file,cRNA_pattern_dict,mRNA_pattern_dict)
    
    # Create an empty DataFrame
    cRNA_df = pd.DataFrame()
    mRNA_df = pd.DataFrame()
    for f in tqdm(files, desc="Processing", unit="file"):
        # mRNA, cRNA = count_sequences_with_pattern(f,cRNA_pattern_dict,mRNA_pattern_dict)
        name = f.split('/')[-1].split('-')[0][:-9]
        # Convert each dictionary to a DataFrame
        temp_mRNA_df = pd.read_csv(f'data/tmp/{name}_mRNA.csv')
        temp_cRNA_df = pd.read_csv(f'data/tmp/{name}_cRNA.csv')
        # Append the temporary DataFrame to the main DataFrame
        mRNA_df = mRNA_df.append(temp_mRNA_df, ignore_index=False)
        cRNA_df = cRNA_df.append(temp_cRNA_df, ignore_index=False)
    mRNA_df=mRNA_df.reset_index()
    cRNA_df=cRNA_df.reset_index()
    mRNA_df.fillna(0, inplace=True)
    cRNA_df.fillna(0, inplace=True)
    mRNA_df.to_csv('result/mRNA_count_by3prime_new.csv')
    cRNA_df.to_csv('result/cRNA_count_by3prime_new.csv')

if __name__ == "__main__":
    main()