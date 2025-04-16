import argparse
import os
import numpy as np
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

def compute_gap_fractions(msa_array):
    col_gap_counts = np.sum(msa_array == "-", axis=0)
    fractions = col_gap_counts / msa_array.shape[0]
    return [np.mean(fractions >= threshold) for threshold in [0.1, 0.25, 0.5, 0.75, 1.0]]

def compute_column_pw_scores(msa_array):
    matrix = substitution_matrices.load("BLOSUM62")
    col_pw_scores = np.zeros(msa_array.shape[1])
    
    for col_idx, col in enumerate(msa_array.T):
        unique, counts = np.unique(col[col != "-"], return_counts=True)
        if np.sum(counts) < 2:
            continue;
        norm_counts = counts / np.sum(counts)

        # Initialize
        curr_col_pw_score = 0

        # Compute weighted pairwise scores
        for i, aa1 in enumerate(unique):
            for j, aa2 in enumerate(unique):
                score = matrix.get((aa1, aa2), matrix.get((aa2, aa1), 0))
                weight = norm_counts[i] * norm_counts[j]  # Product of norm_counts
                curr_col_pw_score += score * weight

        col_pw_scores[col_idx] = curr_col_pw_score

    #print(f"col_pw_scores: {col_pw_scores}")
    return np.percentile(col_pw_scores, [25, 50, 75, 90, 100])


def main(input_fasta, output_dir):
    if not os.path.exists(output_dir):
        print(f"Error: Output directory '{output_dir}' does not exist.")
        exit(1)
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    msa_array = np.array([list(str(seq.seq)) for seq in sequences])
    
    unaligned_file = os.path.join(output_dir, "unaligned.fasta")
    with open(unaligned_file, "w") as f:
        for seq in sequences:
            f.write(f">{seq.id}\n{str(seq.seq).replace('-', '')}\n")
    
    num_seqs = len(sequences)
    msa_length = msa_array.shape[1]
    gap_fractions = compute_gap_fractions(msa_array)
    quartiles = compute_column_pw_scores(msa_array)
    
    stats_file = os.path.join(output_dir, "msa_statistics.txt")
    with open(stats_file, "w") as f:
        f.write(f"MSA statistics for {os.path.basename(input_fasta)}:\n")
        f.write(f"Num. seqs: {num_seqs}\n")
        f.write(f"MSA length: {msa_length}\n")
        for perc, frac in zip([10, 25, 50, 75, 100], gap_fractions):
            f.write(f"Frac. columns with at least {perc}% gaps: {frac:.2f}\n")
        f.write(f"Columns' sum of pairwise scores by BLOSUM62 weigthed by AA freqs.\n")
        f.write(f"Column sum weighted pw scores at 25%: {quartiles[0]:.2f}\n")
        f.write(f"Column sum weighted pw scores at 50%: {quartiles[1]:.2f}\n")
        f.write(f"Column sum weighted pw scores at 75%: {quartiles[2]:.2f}\n")
        f.write(f"Column sum weighted pw scores at 90%: {quartiles[3]:.2f}\n")
        f.write(f"Column sum weighted pw scores at 100% (highest): {quartiles[4]:.2f}\n")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute MSA statistics and generate unaligned sequences.")
    parser.add_argument("input_fasta", help="Input FASTA file of an amino-acid MSA")
    parser.add_argument("output_dir", help="Output directory for the results")
    args = parser.parse_args()
    main(args.input_fasta, args.output_dir)
