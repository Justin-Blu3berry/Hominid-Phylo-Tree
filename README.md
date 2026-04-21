# Hominid-Phylo-Forest
My final project for BINF6251

If I uploaded the wrong link to the Canvas assignment, follow [this link](https://github.com/Justin-Blu3berry/Hominid-Phylo-Tree/pull/2) to the peer review pull request

Last updated: 04/21/2026 at 6:00p. EST
Current report progress: FINAL_REFLECTION.md (one section left) and TESTING.md have yet to be completed

## Project Overview
When creating a phylogenetic tree using a distance-based treebuilding algorithm, the tree's topology is dependent on the distances between all of the species in the analysis; however, it is not practical to use each species's entire genome when calculating these distances. As a result, the distances between the species tend to be proxied by calculating the distances between each species's sequence for a given gene. 

The question that I sought to answer in this project was to what degree the topology of a phylogenetic tree is affected by the choice of which gene is used to calculate interspecies distances. 

The algorithm of choice for this project was the Neighbor Joining algorithm. 

My script takes sequence data as inputs in the form of fasta files, where each gene of interest gets its own file and has an entry for each species of interest. A configuration json file is used by the main driver script to point to what fasta files should be read for analysis. The main driver script will generate a series of plaintext files in the outputs directory with names following the format `<gene_name>_outputs.txt`, one for each gene, containing the pretty-printed distance matrix derived from that gene's sequences and the Newick string representation of the phylogenetic tree constructed using those distances. It will also generate an image of the aforementioned tree. Finally, it also generates a "mean" distance matrix, which contains a mean distance for each pair of species in the analysis, as calculated using that pair's distance on each gene's distance matrix. The resulting distance matrix and a diagram of the resulting tree are saved to a text file as `mean_outputs.txt` and `mean_tree.png`. 

## Installation & Setup
Step-by-step instructions to set up the environment, including:
Python version (and other relevant tools, if any).
How to install dependencies (e.g., via pip, conda).
Any system-level requirements (if applicable).

This project was written and tested using Python 3.10.12, and the package manager of choice was pip 22.0.2. I recommend pip-installing all packages listed with a version number in `requirements.txt` into your virtual environment. Anything without a listed version is a part of the standard Python library and does not need to be pip-installed. 

From the repo, download the `scripts/` directory and its contents, `config.json`, and `toy_config.json`. 

I designed this project's scripts to be run from the command line in the root of the project's directory, which ensures that `data/` and `outputs/` directories will be created outside of the `scripts/` directory. 

## Quick-Start
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

- `<gene_name>_outputs.txt`: a plaintext file containing the Newick string representation of the tree constructed using the species' sequences for this gene and the distance based on that gene. One of these will exist for each gene specified in `toy_config.json`.  
- `<gene_name>_tree.png`: a png file showing the tree generated based on the family's sequences for the given gene.  
- `mean_outputs.txt`: a plaintext file containing a the Newick string representation of a tree and the distance based the average of each position on all genes' distance matrices (i.e. the mean matrix's value at i,j is equal to the mean of all the distance matrices' values at position i,j). The resulting distance matrix and tree essentially account for all of the genes by coming up with what the "average" distances between two given species is when looking at all of the genes of interest. Currently, it requires all of the distance matrices to have the same shape, necessitating that each species have sequences for every gene in the analysis.  

At present, a message prints to the terminal when any species are missing sequences to flag the user to what species are missing sequences so they can be removed from the analysis or so their sequences can be obtained. In a future version, the script will dynamically determine what sequences each species is missing (if any) and exclude the corresponding matrix when trying to determine the species in question's mean distances to each other node (e.g. *Gorilla beringei* is missing a sequence for ACE2. The ACE2 matrix, lacking a column for *Gorilla beringei*, is then not included when determining the mean distance between *Gorilla beringei* and anything else). 

### For Use with Hominid Sequences
Create a directory named `data` in the root of the project directory

Download data using `python scripts/data_download.py --config config.json --outdir data/` 

Now, you are ready to run the analysis from the root of the project directory: `python scripts/main.py --config toy_config.json --indir data/ --outdir outputs/`. The same output files will be generated here as with the toy data in the previous section, with the main difference being the number of genes in the analysis and which species are used. 

## Usage and Options

The scripts that a user would want to run are:
- `scripts/toy_data_generator.py` to generate the testing data that can be used to demonstrate treebuilding functionality. 
- `scripts/data_download.py` to download fasta sequences that are used to build phylogenetic trees. 
- `scripts/main.py` to build the trees for the species of interest using their respective sequences for the genes of interest. 

These scritps are run from the command line while in the root of the project directory (the parent to the `scripts/` directory) to ensure that the `data/` and `outputs/` directories are not created inside of the `scripts/` directory (though if you ran all of these while your working directory was set to `scripts/`, nothing would break). Usage of these scripts is demonstrated in the previous section. 

All three of these scripts require a configuration file to run from the command line. Config files are expected to be json files with the fields: "API_key", "list_gene_names", "list_species", and "email". 

`toy_config.json` only requires "list_gene_names" and "list_species" to be filled out, as it's intended use is with the script `toy_data_generator.py` and with `main.py`. Since this script makes a hard-coded tree structure to determine which sequences are inherited by which species during the simulated sequence evolution, it is advised that the species names in this config file are never changed. Additionally, since `toy_data_generator.py` generates a predetermined number of genes with a specific naming scheme, I advise not changing the gene names in `toy_config.json` either. 

`config.json` is provided in the repo as a template for the user, and all of its fields are required when used to download sequences via `data_download.py`. The user's requests for sequence data from NCBI's nucleotide database will be rate-limited if they don't have an API key, and I haven't tested if that impedes funcitonality, so I recommend making an account on NCBI's website and generating an API key. Your email is also required by the Entrez package whenever making a request. When used with `main.py`, only "list_gene_names" and "list_species" are required, as the other two fields only matter for the NCBI queries in the download script. There needs to be at least one gene in "list_gene_names" so an output can be generated, and there need to be at least four species in "list_species" so that the Neighbor Joining algorithm can resolve the tree. 

Note: you can edit the gene names and species being searched for my editing `config.json` directly, but results aren't guaranteed, as the tree based on the mean of each interspecies distance will not generate if any species lack sequences for any of the genes, and the NCBI search terms were weirdly sensitive during testing. 

!!! IMPORTANT !!!!
Any time you add or remove either species or genes in your config file, run `data_download.py` using that config file before then running `main.py`. 

This ensures that the genes and species in the config file actually have data for `main.py` to use. 

## Limitations and Assumptions
Most of the assmuptions that I made during this project related to the files and the data that the user is running the scripts with. 
#### Assumptions made about inputs
- The user remembered to download their data files before running `main.py`.
- Each input file contains a header for every species listed in the config file
- An input file exists for each gene in the config file
- The user has not changed the key names in their config file
In the interest of time, I had to cut a lot of the file validation functions that would allow these to be handled. 
#### Assumptions made during data download
- A sequence of length longer than the cutoff (50,000 bp by default) is a contig, scaffold, or a whole assembly, and anything shorter is a gene sequence. This may be more adequately addressed by reading the header for the fasta sequence to check for terms like "full-genome", "reference", "complete", "shotgun", etc. that indicate a sequence is more than just a gene, though it comes with the downside of downloading each candiate sequence before checking if it's too long to use, which risks time and memory inefficiency. 
- All species in the config file have sequences (observed or predicted) for each of the genes of interest
- The names for the genes of interest are correct (the NCBI nucleotide database is VERY particular about gene names)
- The first result on the nucleotide database that fulfills the length requirement (<= the cuttoff length) is the right sequence for the gene + species combo that we're looking for.  
This ignores the potential of being from an experimental assembly, rather than that species's widly-accepted reference assembly. This also ignores the potential of an mRNA transcript being the first result (databases like nucleotide tend to use T instead of U in RNA sequences), which means that the potential effects of alternative splicing on the subsequent alignments are ignored. 

#### Assumptions made by the "mean" matrix
- All of the genes' distance matrices have the same shape and the same list of species on them
- A species in the config file is either missing sequences for NONE of the genes, or it's missing sequences for ALL of the genes (no tolerance for in-between)
In the interest of time, I've had to defer making `get_mean_matrix()` actually check which species are on each tree to correctly identify each species's row and column on each distance matrix. This forces it to assume that each tree's distance matrix has the same list of species on it. 
#### Assumptions made in tree-building
- The distance matrices are additive, this is an unavoidable assumption inherent to the algorithm I chose
- Distance calculations are calculated assuming the lowest possible number of mutations required for one sequence to become another. This is another somewhat unavoidable assumption that we make because mutations that occurred and were overwritten by future mutations have no way of being detected. 
- There are more than three species in the analysis. This assumption has to be made because at 3 species, the Q-matrix cannot resolve which species is the outgroup (see `IMPLEMENTATION.md` for details), and at 2 or fewer species, nothing on the tree can be resolved. 