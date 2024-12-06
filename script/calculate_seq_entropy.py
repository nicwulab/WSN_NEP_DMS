import matplotlib.pyplot as plt
from Bio import AlignIO
from collections import Counter
import math
import numpy as np


def calculate_entropy(column):
    total = len(column)
    counts = Counter(column)
    entropy = 0.0
    for residue, count in counts.items():
        if residue != 'X':  # ignore gaps
            frequency = count / total
            entropy -= frequency * math.log2(frequency)
    return entropy


# Function to filter out duplicate sequences from the MSA
def filter_unique_sequences(msa):
    unique_seqs = []
    seen_seqs = set()

    for record in msa:
        seq_str = str(record.seq)
        if seq_str not in seen_seqs:
            unique_seqs.append(record)
            seen_seqs.add(seq_str)

    return unique_seqs


# Function to calculate entropy for the MSA and specific positions range (e.g., 2 to 113)
def calculate_entropy_for_msa(msa_file, start_pos, end_pos):
    alignment = AlignIO.read(msa_file, "fasta")  # Read the MSA file
    unique_alignment = filter_unique_sequences(alignment)  # Remove duplicate sequences
    num_columns = alignment.get_alignment_length()

    # Make sure the positions are within the alignment length
    if end_pos > num_columns:
        end_pos = num_columns

    entropy_scores = []

    # Adjust indexing for the range [start_pos, end_pos], Python indices are 0-based
    for i in range(start_pos - 1, end_pos):
        column = [str(record.seq[i]) for record in unique_alignment]
        entropy = calculate_entropy(column)
        entropy_scores.append(entropy)

    return np.array(entropy_scores)


# Bar plot for entropy across the selected region
def plot_entropy_bar(entropy_scores, start_pos, end_pos, output_file):
    plt.figure(figsize=(12, 1))

    # Adjust x-axis for the selected region and create a bar plot
    positions = np.arange(start_pos, end_pos + 1)
    plt.bar(positions, entropy_scores, color='black', width=0.8)


    # Set border (spine) linewidth to 1.5 and remove ticks
    ax = plt.gca()  # Get current axis
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    ax.tick_params(left=False, bottom=False)  # Remove ticks
    # Optional: Add grid and layout adjustments
    plt.tight_layout()
    plt.xlim(1, 114)
    # Save the plot as a PNG file
    plt.savefig(output_file, dpi=1200, format='png')
    plt.show()

# Line plot for entropy across the selected region
def plot_entropy_line(entropy_scores, start_pos, end_pos, output_file):
    plt.figure(figsize=(12, 1))

    # Adjust x-axis for the selected region
    plt.plot(np.arange(start_pos, end_pos + 1), entropy_scores, color='black', linewidth=2)

    # Set labels and title

    # Set border (spine) linewidth to 1.5 and remove ticks
    ax = plt.gca()  # Get current axis
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    ax.tick_params(left=False, bottom=False)  # Remove ticks
    # Optional: Add grid and layout adjustments
    plt.tight_layout()
    plt.xlim(1, 114)
    # Save the plot as a PNG file
    plt.savefig(output_file, dpi=1200, format='png')
    plt.show()
# Heatmap plotting function with x-axis starting at 1 and saving the plot
def plot_entropy_heatmap(entropy_scores, output_file):
    plt.figure(figsize=(12, 2))  # Adjust width and height as needed
    entropy_matrix = np.expand_dims(entropy_scores, axis=0)  # Create a matrix for the heatmap
    plt.imshow(entropy_matrix, aspect='auto', cmap='viridis', interpolation='nearest')

    # Set x-axis labels to start from 1
    plt.xticks(ticks=np.arange(len(entropy_scores)), labels=np.arange(1, len(entropy_scores) + 1))

    plt.colorbar(label='Shannon Entropy')
    plt.title("Entropy Across MSA Positions")
    plt.xlabel("Position")
    plt.yticks([])  # Remove y-axis ticks since it’s a 1D heatmap
    plt.tight_layout()

    # Save the plot as a PNG file with 300 ppi (dots per inch)
    plt.savefig(output_file, dpi=1200, format='png')
    plt.show()
# Usage
msa_file = "../Data/NEP_95.fasta"  # Replace with your MSA file
start_pos = 2  # Start position
end_pos = 113  # End position
# Calculate entropy only for the positions 2 to 113
entropy_scores = calculate_entropy_for_msa(msa_file, start_pos, end_pos)

# Save heatmap to a PNG file with 300 ppi
output_file = "../graph/NEP_entropy_barplot.png"
# plot_entropy_heatmap(entropy_scores, output_file)
plot_entropy_bar(entropy_scores, start_pos, end_pos, output_file)
