import os
import pandas as pd
from glob import glob
# Directory where the count files are located
count_dir = '/home/yiquan2/NEP/result/*/*_counts.txt'
# Glob all count files in the result directory
count_files = glob(count_dir)

# Initialize an empty DataFrame to store the counts
merged_counts = pd.DataFrame()

# Loop through each count file and merge it into the DataFrame
for count_file in count_files:
    # Extract the sample name from the file name
    sample_name = count_file.split('/')[-1].split('-')[0][:-9]
    # Read the count file into a DataFrame
    counts_df = pd.read_csv(count_file, sep='\t', header=None, names=['gene', sample_name])
    # Merge the counts into the main DataFrame
    if merged_counts.empty:
        merged_counts = counts_df
    else:
        merged_counts = pd.merge(merged_counts, counts_df, how='outer', on='gene')
# Sort columns by name
merged_counts = merged_counts.reindex(sorted(merged_counts.columns), axis=1)


# Save the merged counts to a TSV file
merged_counts.to_csv('result/NEP_count_table.tsv', sep='\t', index=False)
