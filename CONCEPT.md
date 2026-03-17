# Pseudocode

Wow


## Data download

Data will be downloaded by a python script that requests the sequences from NCBI, using a config file containing the user's API key (if they have one) and a list of the sequence IDs to download. 
```
Read `config.json` to grab the API key and list of sequence IDs, sorted by gene and species

Navigate to the data folder

For each gene:
  make an API request to NCBI to download the sequence data
  Write all sequences for the same orthologous gene to one FASTA file- this file now has one sequence for each speices' version of the gene
    Name these files according to <gene_symbol>.fasta

```

## File parsing

Sequence data needs to be read from the fasta files in the `/data` directory

### Identifying if a species has sequence data in the file associated with a given gene

This is necessary because each file is dedicated to a different gene, and we want to make sure that all the species that we're analyzing have sequence data. 

```
Expect to receive:
Path to the directory to be searched in, as a path
list of gene symbols
list of species

initialize an empty dict to track file statuses
go to the path
go thru every gene symbol and check if <gene_symbol>.fasta exists
  Add entry to status dict for this gene symbol, with the value of "missing file" if it can't be found
  if it exists, open the file
    initialize a dict mapping species names to booleans, with the default value of false, this will indicate if they're in the file
    read the file line by line, if any of the species names can be found on the line, set its value in the dict to true
at the end of reading the file, if any of the genes still have a value of false, they're missing from the file, so add an entry to the status dict for this file with "missing species: " listing all the missing species
```


### Parsing the sequence data from the fasta files
```
Expect to receive:
file path as a string, where we already know that the file path exists
list of species: list of strings, these are species we're looking to analyze, ignore all others

```

## Sequence Comparison and Tree-Building

### 


## Driver

```
Parse command line arguments to get the path to the config file, the data directory, and the outputs directory

Check that all of these paths exist
  if the config file doesn't exist or wasn't passed, exit and inform the user that they need to get a config file
  if the data directory doesn't exist, make one in the working directory
  if the outputs directory doesn't exist, make one in the working directory

Open the config file
  Grab the list of gene symbols to use for analysis
  Grab the list of species for later as well
  Grab the API key

Check if any files in the data directory for the gene symbols called out in the config file (expect files to be named <gene_symbol.fasta>)
Also attempt to grep species names from the fasta files to make sure every species is covered
  If they all exist and no species are missing, move on
  If any of them don't exist
    Make a list of the gene symbols that lack a fasta file
    run functions from the download script to generate the missing files or write in the sequences from species missing from the files that did exist
      Inform the user what files needed to be downloaded

Now make a tree on a gene-by-gene basis:
For every gene symbol in the list (that came from the config file)
  Open that gene's file and make a dictionary mapping species name to their respective sequence for this gene (helper function that takes the file path)
  Make a distance matrix using the genes in the dict (helper function)
  Write the distance matrix to a file (pretty print as a pandas dataframe) in the outputs/<gene_symbol>/ directory
  Use the distance matrix to make a Newick tree via Neighbor Joining (helper function for tree building), save this Newick representation to a plaintext file under the outputs/<gene_symbol>/ directory
  Use Biopython's Phylo package to visualize the Newick tree, save the image with the name <gene_symbol>.<image_format> under the outputs/<gene_symbol>/ directory

Now we need to compare distance matrices between genes

Make a copy of the distance matrix
Also initialize an empty dictionary for storing each comparisons' sample statistics
For every position (i,j) on the distance matrix (they all have the same shape and the same set of coordinates is always the same pair of species):
  Check which species are being compared at this position and check if the comparison between these two species has already been recorded in the dictonary (homo-sap to pan-tro will be the same as pan-tro to homo-sap, so don't do all these laborious steps more than you have to)
    If the comparison at this position has already been made, skip it
  For each gene:
      Get the distance value for this coordinate and add it to a list
  After getting the values from all the genes, take the mean, standard deviation, and standard error of the distances in this list
  Write the mean to this position (i,j) and to its opposite (j,i) on the mean matrix
  In the sample stats dictionary, record the mean, standard deviation, and standard error, suc that the dict is formatted like this: {<species1>-<species2> : { mean: <val>, std_dev: <val>, std_err: <val> } }

Convert the sample stats dictionary into a pandas dataframe and write it to a file in the output/ directory
  TODO: do something with this information PLEASE

Now make a tree based on the average distance matrix, save the newick representation as a plaintext file in the outputs/ directory
Visualize the tree using Biopython's Phylo package, save the image in the outputs/ directory

Now feed all of the newick trees into a helper function that compares the groups these trees formed, write output into a plaintext report in the outputs/ directory
```
