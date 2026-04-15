# Hominid-Phylo-Forest
My final project for BINF6251

If I uploaded the wrong link to the Canvas assignment, follow [this link](https://github.com/Justin-Blu3berry/Hominid-Phylo-Tree/pull/2) to the peer review pull request

## Project Overview
Still working on this as of 04/15/2026 at 4:22 pm. EST, I unfortunately have class right now so I can't quite finish the report. 

## Installation & Setup

## Quick-Start
Pip-install all packages listed in `requirements.txt`
### For use with toy data
Create a directory named `data` in the root project directory  
Generate data files by running in the root project directory `python scripts/toy_data_generator.py --outdir data/`  

This file simulates the evolution of five genes inherited by six species from a common ancestor. The five genes were generated with various levels of conservation, as proxied by the number of mutations allowed per generation. The following are generated:  

- `data/expected_outputs.txt`: a plaintext listing of the toy data's gene names and the Newick representation of the tree structure that was used during the simulated evolution to determine inheritance. All of the trees here use the same topology, but the limb lengths on each gene's tree is different, as these were calculated by aligning each child's sequence to its parent during the simulated evolution process. While interesting, these limb lengths are fundamentally based on alignments that we would not have the data to do in real data. As a result, the topology is more important than the lengths here. 
- `data/toy_gene_<i>.fna`: a fasta file containing the evolved sequence for toy_gene_\<i\> in all of the extant species listed in `toy_config.json`, as well as their ancestors' sequences that were used during the simulated sequence evolution. These ancestral sequences would not be in realistic data, but were useful for debugging and are ignored by `main.py`'s analysis.  
- `data/expected_gene_<i>.png`: a png file depicting the topology for all of the Newick trees that were printed in `data/expected_outputs.txt`. As mentioned before, these trees will all share the same topology that was used to determine inheritance during the simulated evolution process that generated the sequence, and their only differences are in the limb lengths.  

After generating the toy data, run the analysis from the root of the project directory: `python scripts/main.py --config toy_config.json --indir data/ --outdir outputs/`

Note that --indir and --outdir are optional and will check the working directory for directories named `data` and `outputs` if nothing is provided, and it will attempt to make these directories if they aren't found. However, the functionality to run `toy_data_generator` and `data_download` after creating the data directory from scratch has yet to be implemented, so it is advised that you point the script to a data directory with the sequences of interest or that you put the sequences of interest into `data/`.  
This will produce the follwing in the `outputs/` directory:  

- `<gene_name>_outputs.text`: a plaintext file containing the Newick string representation of the tree constructed using the species' sequences for this gene and the distance based on that gene. One of these will exist for each gene specified in `toy_config.json`.  
- `<gene_name>_tree.png`: a png file showing the tree generated based on the family's sequences for the given gene.  
- `mean_outputs.text`: a plaintext file containing a the Newick string representation of a tree and the distance based the average of each position on all genes' distance matrices (i.e. the mean matrix's value at i,j is equal to the mean of all the distance matrices' values at position i,j). The resulting distance matrix and tree essentially account for all of the genes by coming up with what the "average" distances between two given species is when looking at all of the genes of interest. Currently, it requires all of the distance matrices to have the same shape, necessitating that each species have sequences for every gene in the analysis.  

At present, a message prints to the terminal when any species are missing sequences to flag the user to what species are missing sequences so they can be removed from the analysis or so their sequences can be obtained. In a future version, the script will dynamically determine what sequences each species is missing (if any) and exclude the corresponding matrix when trying to determine the species in question's mean distances to each other node (e.g. *Gorilla beringei* is missing a sequence for ACE2. The ACE2 matrix, lacking a column for *Gorilla beringei*, is then not included when determining the mean distance between *Gorilla beringei* and anything else). 

### For Use with Hominid Sequences
Create a directory named `data` in the root of the project directory

Download data using `python scripts/data_download.py --config config.json --outdir data/` 

Note: you can edit the gene names and species being searched for my editing `config.json` directly, but results aren't guaranteed, as the tree based on the mean of each interspecies distance will not generate if any species lack sequences for any of the genes, and the NCBI search terms were weirdly sensitive during testing. 