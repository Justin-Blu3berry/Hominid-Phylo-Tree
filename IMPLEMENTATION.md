# Project Snapshot

The goal of this project is to determine the impact of gene choice on the output of a tree-building algorithm by building a series of trees for the family hominidae using a variety of orthologous sequences as the basis for inter-species distance calculations. The tree-building algorithm of choice is the Neighbor-Joining algorithm using a normalized Smith-Waterman similarity as the basis for distance calculations. Neighbor-Joining was preferred over UPGMA for this application because, despite UPGMA's advantagous speed, Neighbor-Joining generally creates more robust trees thanks to its use of a transformation of distance, Q. The project has a working demonstration of the Neighbor-Joining algorithm, which has been successfully implemented, and this demo can be seen by running `scripts/tree_building_utils.py` from the terminal in the project directory. 

# What Is Implemented - NOT DONE YET
Describe what parts of your pseudocode are fully implemented, partially implemented, or still missing.
If there are deviations from your pseudocode (e.g., additional steps, different data structures), explain why.

The Neighbor-Joining function has successfully been implemented, as well as the classes for node and tree objects, the ability to extract a species's name from a sequence's fasta header, and the ability to read in the sequence data from the fasta files as a dictionary. 

TODO: check for deviations from pseudocode

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

# Reflection on Changes and Challenges -NOT DONE YET
Discuss:
Where your implementation diverged from your original plan in Part 2.
Key challenges faced (e.g., data issues, debugging difficulties, design changes).
Decisions you made to keep the project manageable (e.g., simplifying the model, reducing input size).
# Next Steps - NOT DONE YET
Outline what you still need to do for the final part:
Additional features or refactoring.
More robust testing or validation.
Documentation and Quick Start preparation.
# Generative AI Disclosure (If Used)
I have not used generative AI on this project in any way. 
