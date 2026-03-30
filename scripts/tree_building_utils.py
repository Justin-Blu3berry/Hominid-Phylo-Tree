"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import from modules
import numpy as np
import textdistance as td
from pathlib import Path
from collections import defaultdict

# import from my modules
from tree_objects import Node, Tree


def _get_species_from_header(header: str, verbose: bool = False) -> str:
    """
    Function to parse a fasta sequeence header to pull out just the species name

    Expects header format of:
    >NC_<refseq ID>:<genomic location> <genus> <species> <subspecies> isolate <isolate ID> chromosome <number>, 
    <assembly name>, <space-delimited description of sequencing type>
    This is space-delimited, where the refseq-ID and location are the first element, followed by separate elements
    for the genus, species, subspecies, the word "isolate", the isolate ID, and more after that

    Parameters
    ----------
    header : str
        fasta header to be parsed
    verbose: bool, default is False
        indicates if debug messages should be displayed

    Returns
    -------
    str
        species name, as it appears in the fasta header
    """
    # remove whitespace characters from the ends, snap to lowercase
    header_list = header.strip().lower().split()

    if verbose:
        print(header_list)

    # identify the position for the word "isolate" to accomodate entries that list the subspecies and entries that don't
    if "isolate" in header_list:
        # set the upper bound to the location of the word "isolate" (it gets cut off because slices are exclusive on the right)
        upper = header_list.index("isolate")
    else:
        # if it doesn't appear, just set a default behavior
        upper = 3

    # expect the species name to be the slice of the header list from index 1 to the index where "isolate" occurs
    # (for species like Gorilla gorilla gorilla, using a slice from just 1:3 won't work, since it cuts off the subspecies)
    full_name = header_list[1:upper]

    return "_".join(full_name)


def read_fasta(infile: Path) -> dict[str, str]:
    """
    Function to read the fasta for a given gene and make a dict mapping species names to
    their respective sequences for this gene

    IMPORTANT: this function assumes that the fasta file only contains ONE sequence per species. 
               If one species has multiple sequences in this file, all but the last are overwritten

    This function assumes that the file path already exists and has been validated for 
    containing the species of interest

    Parameters
    ----------
    infile : Path
        pathlib Path object pointing to the fasta file to be read

    Returns
    -------
    dict[str, str]
        dictionary mapping species names to their respective sequence for a given gene
    """
    # initialize sequences dict
    seq_dict = {}
    curr_species = ""
    curr_seq = ""

    # open the file
    with infile.open(mode="r", encoding="utf-8") as fasta:
        
        # iterate over the lines
        for line in fasta:
            
            # check if the line is a header
            if ">" in line:
                # attempt to pack up the current species and current sequence into the sequences dict
                if curr_seq and curr_species:
                    seq_dict[curr_species] = curr_seq
                
                # update the species using the current header
                curr_species = _get_species_from_header(line)
                # reset the current sequence
                curr_seq = ""

            else:
                # this means we're not looking at a header, append the current line to the current sequence
                curr_seq += line.strip()

    return seq_dict


def make_dist_matrix(species_names: list[str], seq_dict: dict[str, str]) -> np.ndarray:
    """
    Function to generate an N x N distance matrix, where N is the number of species and the distance
    is calculated based on the normalized Smith-Waterman similarity between each pair of sequences

    Parameters
    ----------
    species_names : list[str]
        list of species whose sequences are being compared
        The ordering of names in this list will match the orders of the rows and columns of the dist matrix
    seq_dict : dict[str, str]
        dictionary mapping species names to their respective sequences

    Returns
    -------
    np.ndarray
        a 2d numpy matrix with zeros along the diagonals
    """
    # initialize an empty N x N matrix with valeus of 0 at every position
    num_species = len(species_names)
    dist_matrix = np.zeros(shape=(num_species, num_species))

    # iterate through every species and its index in the species list
    for i, curr_one in enumerate(species_names):
        # iterate through every species taht comes after the current one
        # note that because we use a list slice here, slice_idx is less than 
        # the actual column number, and the difference equals i
        for slice_idx, curr_two in enumerate(species_names[i:]):
            
            # get the row number
            j = slice_idx + i

            # skip distance calculation if the comparison is equal to itself
            if i == j:
                continue

            # calculate the distance between the two species' sequences
            # TODO: decide later if I want to use a different distance metric
            curr_dist = td.smith_waterman.normalized_distance(seq_dict[curr_one], seq_dict[curr_two])
            # this produces distances between 0 and 1 by normalizing for the length of the aligned string

            # now update the distance matrix at position i,j and at j,i
            dist_matrix[i][j] = curr_dist
            dist_matrix[j][i] = curr_dist
    
    return dist_matrix


def make_q_matrix(dist_matrix: np.ndarray) -> np.ndarray:
    """
    Function to calculate the Q-matrix for neighbor-joining
    Q is a transformation of the distance between two nodes that more or less indicates
    how much more related two species are to each other than they are to the rest of the tree

    Parameters
    ----------
    dist_matrix : np.ndarray
        2d numpy matrix indicating evolutionary distance between species
        Has dimensions N by N, where N is the number of species being compared
    Returns
    -------
    np.ndarray
        A 2d numpy matrix of identical shape to the distance matrix showing the Q-value for each pair of species
        on the distance matrix
    """
    # initialize an N x N matrix with np.inf as the default value (this way when we argmin the q-matrix we
    # don't get a comparison between a species and itself)
    q_matrix = np.full(dist_matrix.shape, np.inf)
    nrow = q_matrix.shape[0]
    ncol = q_matrix.shape[1]

    if nrow != ncol:
        print("Number of rows on distance matrix isn't equal to number of columns!!!\n"
              f"Distance matrix at issue:\n{dist_matrix}")

    # iterate through the non-redundant comparisons on the distance matrix
    for i in range(0, nrow):
        for j in range(i, ncol):
            # skip the diagonals (comparisons to self)
            if i == j:
                continue

            # calcluate q for the two species being compared at i,j on the dist_matrix
            q_ij = (nrow - 2) * dist_matrix[i, j] - sum(dist_matrix[i,:]) - sum(dist_matrix[j,:])
            # now update the q-matrix with this value at both i,j and j,i
            q_matrix[i,j] = q_ij
            q_matrix[j,i] = q_ij
    
    return q_matrix
