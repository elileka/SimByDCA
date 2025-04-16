# SimByDCA

**_Simulating synthetic amino-acid sequences along a tree while respecting structural constraints._**

This code is adapted from the Bitbol's group simulator (forked from [here](https://github.com/Bitbol-Lab/Phylogeny-Partners/tree/v2.0). Please refer to their papers by [Gerardos et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010147) and by [Lupo et al.](https://www.nature.com/articles/s41467-022-34032-y) for a full description of the simulator. Please cite their work if you use this adaptation.

Briefly, given a multiple sequence alignment of biological interest (**bio-MSA**), from which Potts' couplings (J) and fields (h) parameters were inferred using [bmDCA](https://github.com/ranganathanlab/bmDCA.git) and for which a tree was inferred (**bio-T**), the simulator will generate amino-acid (AA) sequences along bio-T by accepting/rejecting proposed mutations according to the Metropolis criterion defined by h and J.

The main differences between this fork and the original:
* Small bug fix in cython_code/**generation_sequence.pyx** based on this [fix](https://github.com/Bitbol-Lab/Phylogeny-ESM2/blob/4d75497116427948de2bb1d7722483e3b95f3781/MSAGenerator/MSAGenerator.py#L52-L55).
* Modified field_bmDCA/**import_msa.py** to determine the MSA length dynamically.
* Generalized the manager-script **generation_sequence.py** so it takes in arguments, converts tree labels to integers and can be run from outside the directory.
* Added auxiliary scripts to convert from npy to FASTA and unalign the simulated MSA.
* Cleaned up: removal of anything that is not used for AA simulations.

### Installation

Requirements: cython, biopython

After cloning this repo: 
```
cd SimByDCA/21_states/cython_code
cythonize -i generation_sequence.pyx
```

### Required input for running the simulator

* An MSA of real (biological) amino-acid sequences (bio-MSA) in FASTA format
* A phylogenetic tree bio-T which was inferred from bio-MSA in Newick format. The simulation will occur along bio-T
* A file `parameters_final.txt` with the results of bmDCA inference on bio-MSA
* A Path to write the output
* The number of flips to generate a root sequence that resembles the bio-MSA: (in our experience: `--eq_flips 100000` works fine)

To infer parameters with bmDCA, we refer you to [the bmDCA repository](https://github.com/ranganathanlab/bmDCA.git) and provide the following tips:
* Our step-by-step installation commands of armadillo and bmDCA can be found in the file `armadillo_and_bmDCA_installation_commands.txt`
* Let _L_ be the number of columns in bio-MSA, then the number of sequences in bio-MSA should be in the order of _L<sup>2</sup>_ for effective inference
* Inference (even when multithreaded) may take several days (5-7 in our experience)
* Once finished, convert the result files as follows:
  `arma2ascii -p parameters_h_final.bin -P parameters_J_final.bin`. This will generate a big file `parameters_final.txt` of the following structure:
```
J 0 1 0 0 0.513048
J 0 1 0 1 0.0502829
J 0 1 0 2 -0.0384053
...
h 0 0 21.1422
...
```
### Running the simulator

```
python SimByDCA/21_states/generation_sequence.py --msa path/to/bio-MSA.fasta --tree path/to/bio-T.newick --bmDCA path/to/bmDCA-on-bio-MSA/parameters_final.txt --output_dir path/to/results/dir --eq_flips 100000
```
### The simulator's output
For technical reasons, the tip labels of bio-T need to be converted to integers for the simulation process. The simulator saves the mapping from the original labels to the numerical labels as well as the bio-T with numerical labels. Additionally, bio-T will be rooted before simulation begins so its rooted version will also be saved.

In all, the simulator will save 4 files in the results directory:
* **bitbol_gen_name_to_index.tsv** - map of original label to number
* **bitbol_gen_input_tree_numerical_names.newick** - bio-T with numerical labels
* **bitbol_gen_rooted_tree_numerical_names.newick** - rooted version of bio-T (the actual input to the simulator)
* **bitbol_gen_sim_eq_msa_along_tree.npy** - the simulated MSA

### Processing the output

To obtain a FASTA file from the npy output file, run:
```
python SimByDCA/21_states/inference_MSA.py bitbol_gen_sim_eq_msa_along_tree.npy bitbol_gen_sim_eq_msa_along_tree.fasta
```

To unalign the MSA and summarize it, run:
```
python SimByDCA/21_states/summarize_and_unalign_sim_msa.py bitbol_gen_sim_eq_msa_along_tree.fasta path/to/results/dir
```
This will create a file `unaligned.fasta` with the unaligned simulated sequences and a file `msa_statistics.txt`
