import pysam
import pandas as pd
from glob import glob
from tqdm import tqdm

def bin_splice_length(splice_length):
    # Adjust the bin size based on your requirements
    bin_size = 20
    return bin_size * (splice_length // bin_size)

def bin_splice_position(splice_position):
    # Adjust the bin size based on your requirements
    bin_size = 500
    return bin_size * (splice_position // bin_size)

def quantify_splice_with_pysam(bam_file):
    # Open BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # Dictionary to store splice ranges and positions for each chromosome
    splice_ranges = {}
    splice_positions = {}

    # Iterate over aligned reads
    for read in samfile.fetch():
        if read.is_unmapped or read.is_secondary:
            continue
        if read.is_forward and read.is_read2 and read.is_paired:
            # Extract CIGAR string
            cigar_tuples = read.cigartuples

            # Process CIGAR to identify introns
            intron_start = None
            for op, length in cigar_tuples:
                if op == 3:  # 'N' in CIGAR indicates intron
                    intron_length = length
                    reference_length = samfile.lengths[samfile.get_tid(read.reference_name)]
                    # using precentage to represent
                    splice_length = (intron_length/reference_length)*100

                    splice_position = read.reference_start + sum([l for o, l in cigar_tuples if o not in [2, 3]])

                    # Bin the splice length and position
                    binned_splice_length = bin_splice_length(splice_length)
                    binned_splice_position = bin_splice_position(splice_position)

                    # Get chromosome information
                    chromosome = read.reference_name

                    # Create dictionaries for each chromosome if not exists
                    if chromosome not in splice_ranges:
                        splice_ranges[chromosome] = {}
                    if chromosome not in splice_positions:
                        splice_positions[chromosome] = {}

                    # Update the dictionaries
                    if binned_splice_length not in splice_ranges[chromosome]:
                        splice_ranges[chromosome][binned_splice_length] = set()
                    else:
                        splice_ranges[chromosome][binned_splice_length].add(read.query_name) 

                    if binned_splice_position not in splice_positions[chromosome]:
                        splice_positions[chromosome][binned_splice_position] = set()
                    else:
                        splice_positions[chromosome][binned_splice_position].add(read.query_name)
        if read.is_reverse and read.is_read1 and read.is_paired:
            # Extract CIGAR string
            cigar_tuples = read.cigartuples

            # Process CIGAR to identify introns
            intron_start = None
            for op, length in cigar_tuples:
                if op == 3:  # 'N' in CIGAR indicates intron
                    intron_length = length
                    reference_length = samfile.lengths[samfile.get_tid(read.reference_name)]
                    # using precentage to represent
                    splice_length = (intron_length/reference_length)*100

                    splice_position = read.reference_start + sum([l for o, l in cigar_tuples if o not in [2, 3]])

                    # Bin the splice length and position
                    binned_splice_length = bin_splice_length(splice_length)
                    binned_splice_position = bin_splice_position(splice_position)

                    # Get chromosome information
                    chromosome = read.reference_name

                    # Create dictionaries for each chromosome if not exists
                    if chromosome not in splice_ranges:
                        splice_ranges[chromosome] = {}
                    if chromosome not in splice_positions:
                        splice_positions[chromosome] = {}

                    # Update the dictionaries
                    if binned_splice_length not in splice_ranges[chromosome]:
                        splice_ranges[chromosome][binned_splice_length] = set()
                    else:
                        splice_ranges[chromosome][binned_splice_length].add(read.query_name) 

                    if binned_splice_position not in splice_positions[chromosome]:
                        splice_positions[chromosome][binned_splice_position] = set()
                    else:
                        splice_positions[chromosome][binned_splice_position].add(read.query_name) 

    # Close BAM file
    samfile.close()
    # Create a new dictionary with keys as the lengths of the sets
    splice_ranges_new = {outer_key: {inner_key: len(splice_ranges[outer_key][inner_key]) for inner_key in splice_ranges[outer_key]} for outer_key in splice_ranges}
    splice_positions_new = {outer_key: {inner_key: len(splice_positions[outer_key][inner_key]) for inner_key in splice_positions[outer_key]} for outer_key in splice_positions}

    return splice_ranges_new, splice_positions_new



if __name__ == "__main__":
    
    result = 'result/'
    bam_file_pattern = "result/*/*Aligned.sortedByCoord.out.bam" 
    
    files = glob(bam_file_pattern)
    for f in tqdm(files, desc="Processing", unit="file"):
        splice_ranges, splice_positions = quantify_splice_with_pysam(f)
        name = f.split('/')[-1].split('Aligned')[0]
        # Convert dictionaries to dataframes
        df_splice_ranges = pd.DataFrame.from_dict(splice_ranges, orient='index')
        df_splice_positions = pd.DataFrame.from_dict(splice_positions, orient='index')

        # Sort columns by ranges
        df_splice_ranges = df_splice_ranges.reindex(sorted(df_splice_ranges.columns), axis=1)
        df_splice_positions = df_splice_positions.reindex(sorted(df_splice_positions.columns), axis=1)

        # Save dataframes as TSV files
        df_splice_ranges.to_csv(f'{result}{name}/{name}_deletion_length.tsv', sep='\t')
        df_splice_positions.to_csv(f'{result}{name}/{name}_deletion_positions.tsv', sep='\t')
