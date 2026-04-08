# Project Snapshot

The goal of this project is to determine the impact of gene choice on the output of a tree-building algorithm by building a series of trees for the family hominidae using a variety of orthologous sequences as the basis for inter-species distance calculations. The tree-building algorithm of choice is the Neighbor-Joining algorithm using a normalized Smith-Waterman similarity as the basis for distance calculations. Neighbor-Joining was preferred over UPGMA for this application because, despite UPGMA's advantagous speed, Neighbor-Joining generally creates more robust trees thanks to its use of a transformation of distance, Q. The project has a working demonstration of the Neighbor-Joining algorithm, which has been successfully implemented, and this demo can be seen by running `scripts/tree_building_utils.py` from the terminal in the project directory. 

# What Is Implemented - NOT DONE YET

The Neighbor-Joining function has successfully been implemented, as well as the classes for node and tree objects, the ability to extract a species's name from a sequence's fasta header, and the ability to read in the sequence data from the fasta files as a dictionary. The primary functionality that has yet to be implemented is reading the config file, downloading the sequence data from NCBI, and calculating the sample statistics for each species-species distance. So far, I have not needed to make any major deviations from the pseudocode I laid out in `CONCEPT.md`. 

# Prototype Demo Description

Currently, a small demo of the Neighbor-Joining function is implemented in the main driver for `scripts/tree_building_utils.py`, which can be accessed by running that module from the command line while in the main project directory (i.e. the parent to the `scripts/` directory). A more robust demonstration is planned that will generate a larger dataset using `scripts/toy_data_generator.py` when it is finished, and that data will be written to the `data/` directory as fasta files. In the meantime, this minimal demo does not require any input files, instead using a hard-coded example with three genes across a family of four species. It prints the resulting trees for each of the three genes to the terminal, as well as information about how the analysis proceeded, should the reader wish to look at the distance and Q matrices that each gene's sequences gave rise to. 

The expected output for all three genes in this demo was a grouping of species_A with species_D, and a grouping of species_B with species_C. A portion of the output from this demo is shown below:  
```
...
Gene_3:
species_A:CCCCAGC
species_B:CCGGCCC
species_C:CCGCC
species_D:CTTCACC
['species_A', 'species_B', 'species_C', 'species_D']
[[0.         0.57142857 0.4        0.42857143]
 [0.57142857 0.         0.4        0.57142857]
 [0.4        0.4        0.         0.4       ]
 [0.42857143 0.57142857 0.4        0.        ]]
Gene_3 Q-matrix:
[[        inf -1.8        -1.8        -1.94285714]
 [-1.8                inf -1.94285714 -1.8       ]
 [-1.8        -1.94285714         inf -1.8       ]
 [-1.94285714 -1.8        -1.8                inf]]
Tree 1:
((species_A:0.1429,species_D:0.0):0.4286,(species_B:0.1429,species_C:0.0):0.4286);
Tree 2:
((species_B:0.1071,species_C:0.1786):0.2143,(species_A:0.1071,species_D:0.1786):0.1071);
Tree 3:
((species_B:0.2857,species_C:0.1143):0.2857,(species_A:0.2143,(species_B:0.2857,species_C:0.1143):0.0714):0.2143,species_D:0.2143);
```
# Data Documentation

The sequences for the three genes in the existing minimal demo were created by hand, starting with an assumed ancestral sequence. Species A and D were given the same mutation as one-another, and a different mutation was given to both B and C. This simulated the splitting of A and D's common ancestor with that of B and C. Each of the species were then given substitutions independent of the other species, simulating further independent evolution of B and C from their common ancestor, and the same for A and D. The mutations given during this second divergence were selected to test the algorithm's resilience against different types of mutation. Gene_one was just a normal test case using substitutions, gene_two was made to test if inversions made it difficult to resolve species B and C's relationship, and gene_three was made to test the effect of a deletion event that partially deletes one of the features that makes B and C distinct from A and D. 

This acts as a small, manually-generated example of the process that the `scripts/toy_data_generator.py` script will implement. By creating an intended tree up-front and then simulating the nucleotide evolution of the sequences from the ancestor, the test data will have known relationships to act as the "ground truth" to verify the accuracy of the trees. During this simulated evolution, we can directly compare the sequence of a node to its parent to find it's "true" limb length. I am eager to see how well the limb lengths calculated during the Neighbor-Joining function will approximate the "true" limb lengths claculated during the evolutionary process. Unfortunately, these "true" limb lengths cannot be obtained for real species without sequence data from the species they evolved from, which isn't realistic when looking at families or any larger taxa that span even greater evolutionary time. So generating toy data like in this way allows me to ensure that the algorithm's process for estimating limb lengths is as strong as possible, which may instill some confidence in the figures generated for real species. I am concerned, though, that these comparisons may be flawed, especially since the sequences for the real species will not be conducive to this method of verification. 

# Initial Observations

Making an extra-small toy example work with only 3 species is disturbingly difficult because the Q matrix seems to end up being filled with all the same values. Perhaps I should have expected this based on the way that Q is calculated by looking at the distances between the two species and all outgroups, which likely falls apart when there is only one outgroup. This blunder while I was generating my miniaturized toy example may prove surprisingly helpful, though, as it has taught me that the smallest tree that the algorithm can resolve is one with four nodes. This leads me to think that I should force the user to supply four species at minimum, rather than only three. 

Aside from that, the algorithm is behaving as expected on the miniaturized demo and during testing in the python console where I imported all of my functions into the session and generated similar toy data on the spot. The runtime is unsurprisingly fast on such small sequences and such few species, and I will certainly be taking advantage of my script for generating toy data to verify that the algorithm runs well with longer seqeuences and a greater number of species in consideration, once the toy data generator is up and running that is. 

# Reflection on Changes and Challenges

My implemention has not made any major deviations from my pseudocode in `CONCEPT.md`. The only deviation I can find was that I added functionality to the tree object's `add_node()` method that lets the user pass a list of the new node's children and, if it was provided, the new node is given an entry in `Tree.edges` that is set to the given list of children. My pseudocode had placed this exact functionality in `make_parent()`, and I had chosen to move it when I was writing `add_node()` to make the method more amenable to adding internal nodes. Truthfully, it doesn't make a large difference. 

The most noteworthy challenge that I had encountered while scripting the Tree-building functionality was that of mutable default values. While I had my modules imported into the python terminal and was toying around with small sequences, I had discovered that if I ran the same line of code to generate a tree off the same data, it put an entire copy of the tree inside of itself, and another copy would appear inside the tree every time it was run. I had discovered the hard way that by initializing an empty tree at the start of `neighbor_joining()` before appending values and adding dict entries to it, I was editing the instances for the default values in my Tree class's constructor. 

Another oversight that I discovered while making my toy data for the demo in `tree_building_utils.py` was that the Q-matrix for a tree with three species on it will be a 3 by 3 matrix where all the non-diagonal positions have the same value. This was very interesting consequence of the Q equation. When there are only three species on the matrix (i, j, and k), the equation for Q starts as:  
$Q_ij = (nrow - 2) * dist(i, j) - sum(dist_matrix[row_i]) - sum(dist_matrix[row_j])$  
Which, with 3 species on the matrix, becomes:  
$Q_ij = 1 * dist(i,j) - (dist(i,j) + dist(i,k)) - (dist(i,j) + dist(j,k))$  
$Q_ij = dist(i,j) + -2*dist(i,j) + -dist(i,k) + -dist(j,k)$  
$Q_ij = -dist(i,j) + -dist(i,k) + -dist(j,k)$  
Now for the Q for j and k, it works out to:  
$Q_jk = 1 * dist(j,k) - (dist(j,j) + dist(i,j) + dist(j,k)) - (dist(k,k) + dist(k,i) + dist(j,k))$  
$Q_jk = dist(j,k) + -2*dist(j,k) + -dist(i,j) + -dist(i,k)$  
$Q_jk = -dist(j,k) + -dist(i,j) + -dist(i,k)$  
Seem familiar? The same goes for i and k:  
$Q_ik = 1 * dist(i,k) - (dist(i,i) + dist(i,j) + dist(i,k)) - (dist(k,k) + dist(i,k) + dist(j,k))$  
$Q_ik = dist(i,k) + -2*dist(i,k) + -dist(i,j) + -dist(j,k)$  
$Q_ik = -dist(i,k) + -dist(i,j) + -dist(j,k)$  
We end up with $$Q_ij = Q_jk = Q_ijk$$ no matter what the distances between any of i, j, and k are. This was really fun to discover first-hand, as I had been quite confused why my toy data would continue creating a tree that groups A and B as siblings, even as I continued to add mutations shared by B and C to widen the distance between them and A. I was worried that my algorithm was seriously flawed until I had it print out the Q matrix for each of my trees and discovered that, even as I widened the distance from the siblings B and C to the outgroup A, the Q matrix continued to maintain identical values on all the non-diagonal positions. 

To keep the project manageable, I am sticking to a simple cast of the mean, standard deviation, and standard error of the mean as the summary statistics of choice when it comes to comparing a given position's distances across all the distance matrices. I have already spent so much time on the tree building algorithm itself, and between the time I anticipate needing to debug the queries for NCBI data and file validation, I don't believe I will reasonably be able to compute anything more advanced for each of the distances or do statistical tests to see if different categories of genes (i.e. genes on one chromosome vs another, or antibodies vs tumor suppressors) generally give rise to more consistent distance metrics. 

# Next Steps

The next feature to be implemented will be the generation of sample data using `toy_data_generator.py`, which should allow for more comprehensive testing of the algorithm. This shouldn't take too long to get operational, and will be a small detour on the way to the more imprtant changes. 

Namely, I need to implement `main.py`'s ability to read the user's configuration json file to ascertain what species and genes are being analyzed. The helper functions for both validating the existing sequence data and for downloading sequences from NCBI both need to be written. I will also have to move the functions `_get_species_from_header` and `read_fasta` from the tree-building module, `tree_building_utils.py` to a separate module which will also be home to the other file parsing functions. 

The next suite of features to add all involve calculating sample statistics on the distances between every pair of species across the series of distance matrices that were generated. I didn't write much pseudocode for this because I wanted to leave open the option to do something more comprehensive, but for the time being, my plan is to just get the mean, standard deviation, and the standard error of the mean for each distance. I then wanted to make a separate distance matrix that contains the mean for each position across all the matrices, and use that matrix to construct the "mean" tree. 

I would also like to write a helper function that identifies the different ways in which each tree groups the species, but doing that might prove challenging. A fast way of accomplishing this might be to just remove the numbers from each gene's newick string and make a dict mapping each unique newick string to the genes that cause the species to be categorized in that way. This isn't very resilient, though, against two trees that make the same categorizations but list literally any of the nodes in a different order. Perhaps I could, alternatively, make functions that tell the user if two tree objects are exactly the same (same placements AND limb lengths) or have the same clusters (they make the same placements, but limb lengths are different). As annoying as it would be to store all the tree objects in memory to report on all the trees' clusters, that might be an easier solution to comparing the clusters in each tree without needing to do anything fancy with the Newick strings. These comparison functions will also make testing and validation easier, so I think I've talked myself into it. 

# Generative AI Disclosure (If Used)
I have not used generative AI on this project in any way. 
