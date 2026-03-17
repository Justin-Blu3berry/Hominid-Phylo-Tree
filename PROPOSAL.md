# Hominid Phylogenetic Forest

#### Justin Wildman, for BINF 6251 at Northeastern University in Spring 2026

## Research Question: 
Is there a gene shared by members of the family hominidae that gives rise to a different phylogenetic tree than is commonly accepted? 

Conventionally, a tree building algorithm isn't using the entire genome of the species in question as the basis for its distance matrix, which necessitates a strategic choice in what sequence is used. This choice is deeply affected by the diversity of the organisms in the group. A strong candidate needs to be sufficiently conserved across the species being evaluated to establish orthology but not so highly conserved that it lacks the variability required to resolve differences between the species being placed (Bohle & Gabaldon, 2012). The 16S rRNA gene, for instance, is highly conserved across kingdoms and is useful for cases like identifying bacterial species from environmental samples, where it is unknown if the cells in question are bacteria or archaea. However, such highly conserved sequences generally don't offer the kind of fidelity needed to resolve placements within families or within genera (Janda & Abbott, 2007). This leads us to the question of how wrong the choice of gene can get. If there are genes that give rise to an incorrect phylogenetic tree, it could be worthwhile to investigate why they give rise to the wrong tree. Does the gene fail to resolve differences between the species in question, giving rise to an ambiguous distance matrix? Could the function of the funciton of the gene provide an explanation for why two distantly related species have highly similar sequences for it? 

## Algorithm and Algorithm Class: 
This project will be using the Neighbor Joining algorithm, which is a member of the neighbor-joining family of phylogenetic tree building algorithms. It was chosen for its relative simplicity to implement, and while UPGMA is somewhat easier to implement and runs faster, Neighbor Joining generally makes better trees. 

The distance matrix required to build a tree with Neighbor Joining could be constructed using a Smith-Waterman alignment, which gives me the freedom to use a more robust set of scoring rules than a static mismatch penalty that equally penalizes all substitutions. As much as I'd like to use BLOSUM, the substitution matrix that assigns scores to substitutions based on how likely the resulting amino acid mutation would be, using BLOSUM would require me to identify the ORFS in the sequences, which is a potentially time-consuming feature to add that might not be worth it given the timeline of this project. 
As such, I intend to take advantage of the `textdistance` package's `smith_waterman.distance` function (or the similar `smith_waterman.normalized_distance` function). This function calculates the distance by setting up the alignment scoring matrix for the Smith-Waterman algorithm, grabbing the maximum value from the scoring matrix to use as the similarity score, and then subtracting the similarity score from the length. The normalized distance and similarity functions just divide these distance and similarity scores by the length, respectively, as the two normalized scores add to 1. 


## Data plan: 
This project requires DNA sequences for a number of genes shared amongst memebers of the family hominidae. Whole genomes aren't necessary, as I'm only looking to calculate distances based on individual genes. Instead, I intend to access the sequences for each of the genes individually, which are easily accessible on NCBI's genome browser as FASTA files. The phylogenetic tree I will be referencing will be the tree Rocatti and colleagues (2019). This tree is only a cladogram, so it won't help evaluate the branch lengths on the trees each gene gives rise to; however, it provides placements for multiple species of *Pongo*, multiple species of *Gorilla*, and multiple subspecies of *Pan troglodytes*. The tree by Pontzer (2012) in Nature was the only source I could find while doing research that had a hominid phylogenetic tree with labeled branch lengths, and it lacks the fidelity within each genus that the tree by Rocatti and colleagues (2019) has. Because I plan to use sequences from these additional species that aren't on the tree from Pontzer's Nature article (2012), there's no way to know from this phylogeny what the branch lengths for each species of *Pongo* or *Gorilla* should be or whether the stated lengths are for any one species, an average of the two species in each genus, or some alternative. Unfortunately, this means that I can only verify the trees I generate on the basis of how they group the species in the family hominidae. 

### Brief "prototype data" plan: 
To test and debug the algorithm, I will generate short sequences to use as toy data. 
These synthesized sequences only need to represent a few species (maybe 3 to 6), represent a few genes (maybe around 3), and have short sequences (maybe ~20 bp). I could then generate these sequences by coming up with a single ancestor sequence for each gene and then splitting it off and mutating it to generate the sequences for all the extant species. By doing this, I will know what the phylogenetic relationships between the species will be for the toy data before making a tree, which lets me come up with an expected output to compare to. This simulates the evolution of a gene's sequence that the real hominid genes would have underwent as family hominidae evolved. 

## Success Criteria: NOT DONE
The algorithm will report the distance matrix it constructed for each orthologous gene and a visualization for the tree based on that distance matrix. A tree whose groupings of species match the phylogeny from Rocatti and colleagues (2019), ignoring the other species on the cladogram, will be considered to have grouped the family correctly. Unfortunately, as previously mentioned, I was not able to find a source that allows me to verify the distance matrices or branch lengths against a ground truth. Therefore, I will have to settle for comparing these values across all of the trees that I generate. While this doesn't allow for the trees to be directly fact-checked, the distribution of each node's branch length could be used to make a statistical estimate of their "true" value. 

## Pitfall scan: 
Gene selection will be difficult, as I'll need a way to identify genes to include in the dataset that are at risk of generating poor trees, but the whole point of the project is to identify genes the make poor trees. Another method might be to hand-pick multiple genes at a time that fall under a common category such as "highly conserved housekeepers that are vital to basic cellular functions," including ATP-synthase, 16S rRNA, and DNA-polymerase III. Another such group might include genes responsible for coordinating skeletal development, which would be a group of genes expected to show great variability within the family due to the diversity of skeletal robustness seen across the family. A simpler approach might be to randomly select a large enough sample of genes that the number of genes that give rise to poor trees and the number that give rise to good trees are likely to be observable (greater than 0). 

There are additional concerns about genes not being represented across all the species in the family (e.g. a gene is selected, and all species except *Gorilla gorilla* have sequence data for it). Ideally the refseq data on NCBI is robust enough that this shouldn't affect too many genes, so this may never even manifest as a problem. 

A potential issue may be memory blowup due to holding a dozen large alignment matrices because the sequences are long. Luckily, only one distance calculation between sequences is occuring at any given time, so this could be resolved by assigning each alignment matrix to the same variable every time an alignment needs to be made. When a new alignment is started, the variable is then reassigned to the new alignment matrix, and the old one is elligible to be reclaimed by Python's native garbage collection, and it will no longer be taking up memory. 

The runtime might also be a concern due to needing to align N sequences M times and building M graphs, where N is the number of species represented and M is the number of genes being looked at. 
I also need a way to quantify or evaluate how similar trees are. For comparing the groupings that each tree makes, it might be necessary to count how many species are grouped together that "shouldn't" be, though this isn't a very robust metric since it double-counts each misplaced species. Another approach might be a sort of Hamming distance analog, where the differences in two trees' placements are quantified by counting how many species would need to move for one tree to become the other. As for branch lengths, there isn't a great source to use as a ground truth for the phylogeny's branch lengths, so one way to describe how different a given branch length on two different trees might be via fold-change. 

## Planned repository structure: 
```
├── python scripts
│   ├── .venv
│   ├── main.py
│   ├── unit_tests
│   │   ├── tree_utils_test.py
│   │   └── main_test.py
│   └── tree_building_utils.py
├── outputs
│   ├── <gene_name>
│   │   ├── <graph image>
│   │   ├── <plaintext file containing distance matrix for this gene>
│   │   └── <plaintext file for the Newick representations of this gene's graph structure>
│   ├── ...
│   └── <this is done for all genes>
├── data (this will not be on the repo because this many fasta files ought not be stored on github)
│   ├── <species name>
│   │   ├── gene_1.fasta
│   │   ├── ...
│   │   └── <this is done for all genes>
│   ├── ...
│   └── <this is done for all species>
├── PROPOSAL.md
├── CONCEPT.md
├── IMPLEMENTATION.md
├── FINAL_REFLECTION.md
├── README.md
├── requirements.txt
└── data.txt (explaining where to find download the data, this might eventually be replaced by a python script that downloads the data for the user)
```

## Generative AI disclosure
Generative AI was not used in the writing of this proposal nor for the formulation of this project idea. 

## References: 
Ari, E., Ittzés, P., Podani, J., Thi, Q. C., & Jakó, E. (2012). Comparison of Boolean analysis and standard phylogenetic methods using artificially evolved and natural mt-tRNA sequences from great apes. Molecular phylogenetics and evolution, 63(1), 193–202. https://doi.org/10.1016/j.ympev.2012.01.010  
Bohle, H. M., & Gabaldón, T. (2012). Selection of marker genes using whole-genome DNA polymorphism analysis. Evolutionary bioinformatics online, 8, 161–169. https://doi.org/10.4137/EBO.S8989  
Janda, J. M., & Abbott, S. L. (2007). 16S rRNA gene sequencing for bacterial identification in the diagnostic laboratory: pluses, perils, and pitfalls. Journal of clinical microbiology, 45(9), 2761–2764. https://doi.org/10.1128/JCM.01228-07  
Pontzer, H. (2012). Overview of hominin evolution. Nature.com. https://www.nature.com/scitable/knowledge/library/overview-of-hominin-evolution-89010983/  
Rocatti, G., & Perez, S. I. (2019). The Evolutionary Radiation of Hominids: a Phylogenetic Comparative Study. Scientific reports, 9(1), 15267. https://doi.org/10.1038/s41598-019-51685-w  
