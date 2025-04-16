import os
import sys

# Ensure the script can find cython_code and field_bmDCA no matter where it is run from
# Get the absolute path of the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the script's directory to sys.path
sys.path.insert(0, script_dir)

import argparse
import numpy as np
import cython_code.generation_sequence as ge
import field_bmDCA.import_msa as im
from Bio import Phylo
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Process MSA, tree, and parameter files and generate outputs.")
    parser.add_argument("--msa", required=True, help="Path to input MSA FASTA file")
    parser.add_argument("--tree", required=True, help="Path to input tree file in Newick format")
    parser.add_argument("--bmDCA", required=True, help="Path to bmDCA parameter file")
    parser.add_argument("--output_dir", required=True, help="Directory where output files will be saved")
    parser.add_argument("--eq_flips", type=int, default=10000, help="Number of flips to perform on root seq towards equilibrium")
    
    args = parser.parse_args()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    print("Starting script...", flush=True)

    # Load the MSA FASTA file from which the bmDCA parameters were inferred:
    with open(args.msa, "r") as msa_fasta:
        msa_records = list(SeqIO.parse(msa_fasta, "fasta"))
    
    # Create a mapping from original sequence names to their order
    msa_name_to_index = {record.id: idx for idx, record in enumerate(msa_records)}
    mapping_file = os.path.join(output_dir, "bitbol_gen_name_to_index.tsv")
    with open(mapping_file, "w") as f:
        for original_name, index in msa_name_to_index.items():
            f.write(f"{original_name}\t{index}\n")
    print(f"Saved MSA name-to-index mapping to {mapping_file}", flush=True)

    # Load and process the tree file
    print(f"Reading tree file: {args.tree}", flush=True)
    tree = Phylo.read(args.tree, "newick")

    # Rename the tips on the input tree as numericals 0, 1, 2, ... such that they match the order of the sequences in the input MSA
    # for example, if the input MSA is:
    # >homo
    # MG-
    # >fish
    # MKK
    # >mouse
    # MGK
    #
    # and the tree is: ((homo:0.3,mouse:0.3):0.4,fish:0.5);
    #
    # then the new tree labels should be:
    # homo - 0, fish - 1, mouse - 2
    # resulting in the tree:
    # ((0:0.3,2:0.3):0.4,1:0.5);
    print("Renaming all nodes as numerical (input tree)...", flush=True)
    for clade in tree.get_terminals():
        if clade.name in msa_name_to_index:
            clade.name = str(msa_name_to_index[clade.name])
        else:
            raise ValueError(f"Tree tip '{clade.name}' not found in MSA!")
    
    tree_output_file = os.path.join(output_dir, "bitbol_gen_input_tree_numerical_names.newick")
    Phylo.write(tree, tree_output_file, "newick")
    print(f"Saved numerical input tree to {tree_output_file}", flush=True)
    
    # Root the tree
    print("Rooting tree at midpoint...", flush=True)
    tree.root_at_midpoint()
    
    # Save rooted tree
    tree_output_file = os.path.join(output_dir, "bitbol_gen_rooted_tree_numerical_names.newick")
    Phylo.write(tree, tree_output_file, "newick")
    print(f"Saved numerical rooted tree to {tree_output_file}", flush=True)

    # Load bmDCA parameters
    print(f"Loading field and coupling parameters from {args.bmDCA}", flush=True)
    Field, Coupling = im.import_msa_bmDCA_dynamic(args.bmDCA)
    
    # Initialize sequence generator
    print("Initializing sequence generator...", flush=True)
    msa_gen = ge.Creation_MSA_Generation(Field, Coupling)

    # Generate simulated MSA
    print("Generating a simulated MSA along the input tree with equilibrium flips...", flush=True)
    simulated_msa = msa_gen.msa_tree_phylo(tree.clade, args.eq_flips)
    print("Finished simulating an MSA along the input tree.", flush=True)

    # Save simulated MSA
    simulated_msa_file = os.path.join(output_dir, "bitbol_gen_sim_eq_msa_along_tree.npy")
    np.save(simulated_msa_file, simulated_msa)
    print(f"Saved simulated MSA to {simulated_msa_file}", flush=True)
    
    print("Script completed successfully!", flush=True)

if __name__ == "__main__":
    main()
