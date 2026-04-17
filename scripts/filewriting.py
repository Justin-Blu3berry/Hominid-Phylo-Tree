# Module for printing and saving outputs to files

# import modules
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO
from Bio import Phylo

# import from my own modules
from tree_objects import Tree


def write_fasta(seq_dict: dict[str, str], fasta_path: Path) -> None:
    """
    Function to write orthologous sequences from various species to a fasta file
    Used in toy_data_generator.py
    ASSUMES: all sequences in the dict are for the same gene, as seen across all the species in the dict

    Parameters
    ----------
    seq_dict : dict[str, str]
        dict mapping species names to their respective sequences for this shared gene
    fasta_path : Path
        path to the fasta file being written
    """
    # write the header and sequence for this species
    text_to_write = (f"> {species}\n{sequence}\n" for species, sequence in seq_dict.items())
    text_to_write = "".join(text_to_write)
    fasta_path.write_text(text_to_write)


def write_outputs(out_path: Path, tree: Tree, seq_labels: list[str], dist_matrix: np.ndarray) -> None:
    """
    Function that writes the tree object and the distance matrix to a text file

    Parameters
    ----------
    gene_name: str
        Name of the gene that was used to generate tree and dist_matrix
    tree : Tree
        Tree being recorded
    seq_labels: list[str]
        List of names for the species on the distance matrix, in the order they appear on the matrix
    dist_matrix : np.ndarray
        Additive distance matrix between species based on their sequences for a given gene
    """
    # figure out indent for header line
    name_lengths = [len(species_name) for species_name in seq_labels]
    longest_name = max(name_lengths)
    left_space = " " * longest_name

    # make a copy of the list so we can format its elements to show along the top of the matrix
    labels = [species_name for species_name in seq_labels]
    # add spaces to make species names align with distances, rounded to nearest 3 decimals
    for i in range(len(labels)):
        if len(labels[i]) < 5:
            # add space to account for the leading 0 and . in each printed distance
            labels[i] += " " * (5- len(labels[i]))
    # join the labels with tabs
    labels = f"\t".join(labels)
    labels = f"{left_space}\t{labels}"

    # now make each line after that
    rows = []
    # iterate through the rows of the distance matrix
    for i, species_name in enumerate(seq_labels):
        # add the indent to the left side where the species labels will go
        indent = " " * (longest_name - len(species_name))
        # format the strings for the distances on the matrix
        cells = [f"{distance:.3f}" for distance in dist_matrix[i]]
        cells = f"\t".join(cells)
        # combine to create the string for this row on the file
        rows.append(f"{species_name}{indent}\t{cells}")
    # now join all the rows of the distance matrix together
    rows = f"\n".join(rows)

    # write the output to the file
    out_path.write_text(f"{str(tree)}\n{labels}\n{rows}")


def plot_tree(tree: Tree, image_path: Path, gene_name: str) -> None:
    """
    Function to generate a plot of the given tree and save it to the given path

    Parameters
    ----------
    tree : Tree
        Phylogenetic tree being plotted
    image_path : Path
        the file path where the image will be saved
    gene_name: str
        The name of the gene used as the basis of the tree
    """
    # convert the string for the text structure into something file-like so phylo.draw will use it
    handle = StringIO(str(tree))

    # read the tree so Phylo can properly display it 
    phylo_tree = Phylo.read(handle, "newick")  # type: ignore

    # start the plot so we can add details to it after throwing the tree on
    fig, ax = plt.subplots(1,1)

    # stick the tree onto the plt plot
    Phylo.draw(phylo_tree, axes=ax, do_show=False,  # type: ignore
               branch_labels=lambda clade: f"{clade.branch_length:.3f}" if clade.branch_length is not None else clade.branch_length)
    
    # title the tree diagram
    plt.title(f"Phylogenetic tree based on {gene_name}")

    plt.savefig(image_path)
    # plt.show()
