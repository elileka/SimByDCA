# SimByDCA

**_Simulating synthetic amino-acid sequences along a tree while respecting structural constraints._**

This code is adapted from the Bitbol's group simulator (forked from [here](https://github.com/Bitbol-Lab/Phylogeny-Partners/tree/v2.0). Please refer to their papers by [Gerardos et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010147) and by [Lupo et al.](https://www.nature.com/articles/s41467-022-34032-y) for a full description of the simulator. Please cite their work if you use this adaptation.

Briefly, given a multiple sequence alignment of biological interest (**bio-MSA**), from which Potts' couplings (J) and fields (h) parameters were inferred using bmDCA and for which a tree was inferred (**bio-T**), the simulator will generate amino-acid (AA) sequences along bio-T by accepting/rejecting proposed mutations according to the Metropolis criterion defined by h and J.

Major differences of this fork to the original:
* Small bug fix in cython_code/**generation_sequence.pyx** based on this [fix](https://github.com/Bitbol-Lab/Phylogeny-ESM2/blob/4d75497116427948de2bb1d7722483e3b95f3781/MSAGenerator/MSAGenerator.py#L52-L55).
* Modified field_bmDCA/**import_msa.py** to determine the MSA length dynamically.
* Generalized the manager-script **generation_sequence.py** so it takes in arguments, converts tree labels to integers and can be run from outside the directory.
* Added auxiliary scripts to convert from npy to FASTA and unalign the simulated MSA.
* Cleaned up: removal of anything that is not used for AA simulations.

### Installation

Requirements: cython, biopython

After cloning: 
```
cd SimByDCA/21_states/cython_code
cythonize -i generation_sequence.pyx
```

### Input requirements

The simulator requires:
* MSA of real (biological) amino-acid sequences (bio-MSA) in FASTA format
* a phylogenetic tree bio-T which was inferred from bio-MSA in Newick format. The simulation will occur along bio-T
* a txt file `parameters_final.txt` with the results of bmDCA inference on bio-MSA
* path to write the output
* number of flips to generate a root sequence that resembles the bio-MSA: (in our experience: `--eq_flips 100000` works fine)

To infer paramters with bmDCA, we refer you to [this repository](https://github.com/ranganathanlab/bmDCA.git) and provide the following tips:
* The number of sequences in bio-MSA should be in the order of square the number of columns in it.
* Inference may take several days (5-7 in our experience)
* Once finished, convert the result files as follows:
  `arma2ascii -p parameters_h_final.bin -P parameters_J_final.bin`. This will generate a big file `parameters_final.txt` of the following structure:
```
# J 0 1 0 0 0.513048
# J 0 1 0 1 0.0502829
# J 0 1 0 2 -0.0384053
...
# h ...
...
```
