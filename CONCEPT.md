
# Project Summary

This project is meant to answer the research question of how different a family's phylogenetic tree will be when generated using different genes as the basis for distance calculations. This will be answered using the Neighbor-Joining algorithm to build phylogenetic trees for the family hominidae, using fasta sequences for individual orthologous genes retrieved from NCBI's datasets command line integration. 


# Inputs, Outputs, and Assumptions

## Expected Inputs
- FASTA sequence data for each orthologous gene, where each gene has its own FASTA file and each sequence in a given file corresponds to a different species' ortholog of the same gene

## Expected Outputs
- A series of plaintext files containing distance matrices, where each distance matrix is based on a one of the genes in the input data, and one of the distance matrices is based on the average of all the other distance matrices' values for each position
- A tab-separated file containing the sample statistics for each species-species comparison, including the mean, standard deviation, and standard error of the mean for the distance between each pair of species
- A series of image files depicting visualizations of the phylogenetic trees constructed based on each gene, including the tree generated from the "average" distance matrix
- A plaintext file listing the Newick string representation of each tree, including the tree generated from the "average" distance matrix
- A plaintext report comparing the groupings that all of the trees made

## Key Assumptions
- Distances between two species' sequences are additive (or nearly additive)
- An ortholog is under equal selective pressure in one species to the orthologs in other species (i.e. the molecular clock for a given gene runs at the same speed in all species being looked at)
- The limb length for an internal node cannot be resolved without at least one leaf on the distance matrix that doesn't descend from it
- Comomonalities between two sequences are always due to the sequences having inherited those commonalities from a common ancestor, rather than independently mutating in the same exact manner.
- The minimum number of mutations occurred between an ancestor and its descendents (i.e. a nucleotide that is the same between two species is assumed to have not mutated to a different nucleotide and then mutated back to the original in the time between their shared ancestor and now). 

# Pseudocode

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
Helper function to make a list of what species are missing from a fasta file
expect to receive:
path to the fasta file to inspect, we're assuming that it exists
list of species that should have sequences in this file

initialize a dict mapping species names to booleans, with the default value of false, this will indicate if they're in the file
  read the file line by line, if any of the species names can be found on the line, set its value in the dict to true
  at the end of reading the file, if any of the genes still have a value of false, they're missing from the file
  so return a tuple of all the species that had a value of false
If there are no genes with a value of false, return an empty tuple
```
```
Helper function to check for missing data
Expect to receive:
Path to the directory to be searched in, as a path
list of gene symbols
list of species
we're trying to see if a given gene has any fasta data in the given directory and, if the gene already has a fasta file, if it's missing any of the species we're looking at in our analysis

initialize an empty dict to track file statuses
go to the path
go thru every gene symbol and check if <directory_to_search>/<gene_symbol>.fasta exists
  if it doesn't exist, Add entry to status dict for this gene symbol, with the value of "missing file"
  if it exists, open the file
    then use a helper function to get a list of the species from our list that aren't in this file
    if the helper function returned a tuple with any species listed, add an entry to the status dict for this file with "missing species: " listing all the species with a value of false
if the status dict has any entries in it, print to stdout all the entries from the status dict
```


### Parsing the sequence data from the fasta files
```
Expect to receive:
file path as a string, where we already know that the file path exists
list of species: list of strings, these are species we're looking to analyze, ignore all others

KEY ASSUMPTION: the data file exists and has already been validated for containing the species of interest

initialize a dictionary mapping species names to their respective sequence for the gene this file is dedicated to
initialize a current_species and a current_sequence

open the file
  read the file line by line
    if the line is a header:
      check if this is the first line in the file (there will not already be a sequence ready to be packed up and added to the dict)
        If this isn't the first line in the file, make an entry in the dict mapping the current_species to the current_sequence
        reassign the current_sequence to an empty string
      extract the species name from the sequence ID (this deserves a helper function), and update the current_species
    if the line isn't a header, just append the line (stripped of whitespace) to the current_sequence

at the end, just return the dict mapping species names to their respective sequence
```

## Sequence Comparison and Tree-Building

### Building a distance matrix
```
Expect to receive:
list of species names
dictionary mapping a species's name to its sequence

initialize an empty N x N matrix, where N is number of species in the list, and each value is initialized as 0

Iterate through every species and its index (i) in the list
  Iterate through every species that comes after the one the current loop is on, along with its associated index (j)
    Use a common distance metric to determine the distance between the current two species' sequences (considering Jaccard or Hamming, or adapting Smith-Waterman alignment scores as a distance metic)
    Fill in both positions on the matrix that correspond to this pair of species (`(i,j)` and `(j,i)`) with the calculated distance metric

At the end, just return the matrix
```

### Node Object
```
Needs to have attributes:
- Name
- Limb length
- The name of its parent

Init function
Takes name and limb length as mandatory inputs, with parent name as an optional input
Assigns name and limb length to corresponding attributes
Assigns parent name if provided, otherwise uses default of "default parent"

__repr__ function
we can make the string representation of this node be condusive to the tree printing its Newick string representation
Spit out a string that concatenates the following: name + colon symbol + limb length
```

### Tree Object
```
Needs the following attributes:
- Nodes: Dict mapping the name of each node on the tree to the node object they represent
- Top_layer: List of the nodes that are in the top layer (they don't have parents that are also nodes)
- Edges: Some way of mapping parent to list of children

Methods:

initialization:
  Mandatory inputs: Nothing
  Optional inputs: node_list: List of node objects to put on the tree
                   top_layer: List of node objects that are in the top layer
                   node_mapping: Dict mapping nodes to a list of nodes that descend from them
  
  All of these need to check that they received the right types as input by the way and use the default values if they don't
  
  Assign the nodes attribute to the the provided node_list parameter, if none is provided just use an empty list
  Assign the top_layer attribute to the provided top_layer param if nothing is provided just use an empty list as well
    If it's empty, it runs the helper func to identify parentless nodes
  
  Assign the edges attribute to a defaultdict with a default value of empty list
  Attempt to unpack the values from the node_mapping dict into the defaultdict


Add leaf:
  Takes as input the node object to be added
  Assumes that the node is being added without any connections to anything else

  Add the given node to the self.nodes list and to self.top_layer list

Identify parentless nodes (aka the nodes that should be in the top layer):
  Doesn't take any params, just uses self.edges list as iterable

  Make a shallow copy of the self.nodes list

  iterate through the values of the edges dict (each value is a list of nodes that are the children of something else)
  remove any values that appear in any of these lists

  return the list of whatever remains after removing everything that was found to be the child of something else


Generate a name for a new internal node based on its children:
  Takes as input: List of child nodes

  Generate a string that concatenates "(", the comma-separated string representations of all the children (using the str.join func), and ")"
    NOTE: the string representation of a child will always be "<child.name>: <child.limb_length>"
    NOTE: since the name for an internal node is just the Newick string representaiton of the subtree that descends from that internal node, this function
          doesn't have to work any differently if a given child is internal or terminal

  
Add a parent for a set of children- This is our "add_edge": 
  Takes as input: the list of the nodes that will be the children to this new internal node
                  distance from this internal node to an outgroup (this gets calculated in the neighbor_joining func since that scope has access to the distance matrix)

  Generate a name for this new internal node using the helper function
  Create a node object for this new internal node, passing the newly-minted name as its name and its distance to an outgroup as its limb length
    NOTE: the neighbor-joining func can update this as well later when this node is grouped with something else, and if this internal node is one of nodes in the 2x2 matrix at the end,
          then its distance to the other node on the 2x2 matrix will be its limb length anyway (so it doesn't matter that it won't be updated)
  Create an entry for this new internal node in the edges dict and map it to its list of child nodes
  Update each child's entry for the name of its parent

  Attempt to remove all the children from the self.top_layer list
  Add this new internal node to the self.nodes list

  return the name for this internal node (so that the neighbor joining func can access it to update the labels for the distance matrix)
  

__repr__ function (this is how we generate the newick string of the entire tree):
  Takes no params, just uses the self.top_layer

  Just calls the function that generates a name for a new internal node, except it passes all of the nodes in the top layer as the children
```

### Tree-Building via Neighbor-Joining

#### Make a Q-matrix
```
Helper function to calculate the Q-matrix
  Takes the distance matrix as input

  Initialize an empty matrix with the same shape as the distance matrix
  Iterates through half the positions on the distance matrix that hasn't been updated yet
  (using something like for i in range(0, nrow) with for j in range(i, ncol) to avoid hitting BOTH i,j and j,i)
    Use +infinity as the value for diagonals (where i=j)
    For non-diagonal posiitons (i != j)
      Calculate Q for i,j using the formula: (n-2)*d(i,j) - sum(row_i) - sum(row_j)
      Fill in positions i,j and j,i with this Q value

  After it's done filling out the matrix, return it
```

#### Calculate Limb Lengths
```
Helper function to calculate limb lengths for two nodes on the distance matrix
  Takes as input:
    The distance matrix
    The row numbers for the nodes in question (i and j)

  get the sum across the row for both of the nodes at issue (call these sumrow_i, sumrow_j)
  for i: limb_length = ( distance(i, j) + ( (sumrow_i - sumrow_j) / (n-2) ) ) /2
    could also use this for j, but it's faster to just do the shortcut below
  for j: limb_length = distance(i,j) - limb_length_i

  return the limb lengths for i and j
```

#### Calculate Distances to an Internal Node
```
Helper function to calculate the distance from the parent to two given nodes to all others
  Takes as input:
    The distance matrix
    The row numbers for the two nodes in question, i and j

  Initialize an empty numpy array of floats for the distances
  Iterate over all of the row numbers that aren't the two given (0 to numrows - 1, excluding i and j)
    apply the formula: dist_to_curr = ( distance(i, current) + distance(j, current) - distance(i, j) ) / 2
    add the distance to the current outgroup to the numpy array

  return the numpy array of distances
    NOTE: the indices in this array will line up with the distance matrix and the label list AFTER i and j are removed from the matrix
```

#### Neighbor-Joining
```
Function to handle graph creation
  Takes distance matrix and list of labels as input

  Initialize an empty tree

  Iterate until the matrix is 2x2
    Make a Q-matrix using what the distance matrix currrently looks like
    Find the position on the Q-matrix with the lowest Q (if multiple are tied and their coordinates aren't just flipped versions of each other, i.e. more than 2 are tied, just pick one)
      this gives us coordinates i,j
    Calculate the limb lengths for i and j
    Calculate the distances to the remaining nodes for the parent to i and j
    
    Create node objects for i and j using their labels from the label list + their newly-calculated limb lengths
    Add the nodes for i and j to the graph object
    Create a parent node on the graph for nodes i and j, assign the name that this func will return to a variable
    Remove i and j from the distance matrix (+ labels list)
      To remove from the matrix, use numpy.delete and pass the distance matrix, a list containing i and j, and do this for axes 0 and 1
    Add the name for the parent to i and j to the labels list
    Add a row + column at the end of the matrix for the parent to i and j
      Reshape the numpy array of distances to 1 row + length(array) columns and append to the matrix along axis 0 (or use vstack)
      Append a 0 to the end of the numpy array to represent this node's distance to itself
      Reshape the numpy array to length(array) rows + 1 column and append to the matrix along axis 1 (or use hstack)

  Once we've broken out of the loop, we know that there are only two nodes left on the matrix
    ASSUMPTION: at least one of the nodes left on the matrix is an internal node (this will have a node object associated with it, and its limb length will have
                already been set to its distance to the only other node on the matrix), but one of the nodes on the matrix may still be a leaf. 
                This remaining leaf will NOT yet have a node object associated with it on the graph
    if any of the remaining species on the distance matrix doesn't yet have a node object for it, create a node object for it, with its limb length
    equal to its distance to the only other remaining species on the matrix
      This is checked by just seeing if it's in the list of names in graph.nodes
    Then add the node for this leaf to the graph

    Now just return the newick string representation of the tree
```

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
  In the sample stats dictionary, record the mean, standard deviation, and standard error, such that the dict is formatted like this: {<species1>-<species2> : { mean: <val>, std_dev: <val>, std_err: <val> } }

Convert the sample stats dictionary into a pandas dataframe and write it to a file in the output/ directory
  TODO: do something with this information PLEASE

Now make a tree based on the average distance matrix, save the newick representation as a plaintext file in the outputs/ directory
Visualize the tree using Biopython's Phylo package, save the image in the outputs/ directory

Now feed all of the newick trees into a helper function that compares the groups these trees formed, write output into a plaintext report in the outputs/ directory
```

# Complexity and Bottlenecks

I'm avoiding saving too many matrices to memory by overwriting the same distance matrix for every gene I use and just writing each one to its own file before I overwrite it. Writing the newick string representations of trees to files after I'm done with them. Only working with one gene/file at a time and writing its outputs to files before moving on to the next in this manner should ideally help prevent memory from being occupied by information associated with trees that are already done being made. I'm primarily concerned about the runtime of the algorithm, as Neighbor-Joining has a time complexity of $O(n^3)$, which 

I have concerns about:
- the matrices stored in memory being too large (particularly for alignment scoring to calculate distance between two large sequences), which may be a justification for using a different distance metric
- the time complexity for distance calculations potentially scaling poorly for the long sequences from NCBI
- The runtime being poor because Neighbor-Joining generates each tree too slowly and I'm generating a bulk amount of trees (upgma may shine here but makes worse trees)

# Validation and Testing Plan
I plan to test the helper functions (especially the ones for downloading data and parsing files) using unit tests via the pytest module. This allows me to use fixtures to test the helper functions against expected outputs. 

I plan to test my implementation as a whole by creating short sequences by hand for a handful of genes from an "ancestor" species and then mutating copies of the ancestral sequences one generation at a time to create the "modern" species' sequences. The key during this simulated evolution would be that the relationships between the sequences are tracked as they're generated, so the phylogenetic tree for these simulated genes is known before the tree-building algorithm even looks at them. This gives the toy data expected outputs to compare my results to as I test my code. The way I choose to mutate the sequences as I work with them can be manipulated to force certain edge cases, creating different genes to test different edge cases during tree building. For example, I could generate sequences such that two leaves (A & C) both have short limb lengths, have parents that are, themselves, not that distant from one another, and where one of the two leaves (A) has a sibling (B) with a very long limb length. In this case, the distance between the two cousins is small in magnitude, but the long limb length of the sibling should not lead to it being classified as the outgroup due to the use of the Q-matrix during Neighbor Joining. I would expect, in this example, that the sibling relationship between A and B is correctly identified and C is correctly identified as the outgroup, with B having an extremely long limb length and A C, the parent to A & B, and the parent to C all having extremely short limb lengths. 
