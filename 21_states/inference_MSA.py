import numpy as np
import argparse

# Define standard amino acid and nucleotide alphabets
AA_ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"  # 21 states (gap + 20 AA)
DNA_ALPHABET = "-ACGT"  # 5 states (gap + 4 bases)

def convert_msa_to_fasta(input_path, output_path, alphabet="AA"):
    """
    Converts a NumPy .npy MSA file to FASTA format.
    
    Parameters:
        input_path (str): Path to the input .npy file.
        output_path (str): Path to save the output .fasta file.
        alphabet (str): "AA" for amino acids, "DNA" for nucleotides.
    """
    # Select the alphabet
    if alphabet == "AA":
        ALPHABET = AA_ALPHABET
    elif alphabet == "DNA":
        ALPHABET = DNA_ALPHABET
    else:
        raise ValueError("Invalid alphabet choice. Use 'AA' for amino acids or 'DNA' for nucleotides.")
    
    # Load the MSA matrix
    msa = np.load(input_path)

    # Check validity of indices
    max_index = max(np.unique(msa))
    if max_index >= len(ALPHABET):
        raise ValueError(f"MSA contains unexpected values. Max index found: {max_index}, but alphabet size is {len(ALPHABET)}.")

    # Write to FASTA file
    with open(output_path, "w") as f:
        for i, sequence in enumerate(msa):
            seq_str = "".join(ALPHABET[idx] for idx in sequence)  # Convert indices to characters
            f.write(f">seq_{i}\n{seq_str}\n")  # FASTA format

    print(f"FASTA file saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a NumPy .npy MSA file to FASTA format.")
    parser.add_argument("input_file", type=str, help="Path to the input .npy file")
    parser.add_argument("output_file", type=str, help="Path to save the output .fasta file")
    parser.add_argument("--alphabet", type=str, choices=["AA", "DNA"], default="AA", help="Choose 'AA' for amino acids or 'DNA' for nucleotides (default: AA)")

    args = parser.parse_args()
    convert_msa_to_fasta(args.input_file, args.output_file, args.alphabet)

