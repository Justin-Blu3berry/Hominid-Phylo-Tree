"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import from modules
import numpy as np
import textdistance as td
from pathlib import Path

# import from my modules
from tree_objects import Node, Tree



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


def calc_limb_lengths(dist_matrix: np.ndarray, i: int, j: int) -> tuple[int, int]:
    """
    Function to calculate the limb lengths for two nodes on the distance matrix
    given the row/column number for each of them

    Parameters
    ----------
    dist_matrix : np.ndarray
        The matrix of distances between each pair of species, has shape (N, N) where
        N is the number of species
    i: int
        the row/column number for one of the species of interest
    j: int
        the row/column number for the other species of interest
    
    Returns
    -------
    Tuple of ints
        the limb length for the species at row i, followed by the limb lenght for 
        the species at row j
    """

    # get the sum across rows for species i and j
    sumrow_i = np.sum(dist_matrix[i,:])
    sumrow_j = np.sum(dist_matrix[j,:])

    num_species = dist_matrix.shape[0]

    # calculate the limb lengths for i and j
    length_i = ( dist_matrix[i,j] + ( (sumrow_i - sumrow_j) / (num_species - 2) ) ) / 2
    length_j = dist_matrix[i,j] - length_i

    return length_i, length_j


def calc_internal_dist(dist_matrix: np.ndarray, i: int, j: int) -> np.ndarray:
    """
    Function to calculate the distance from a parent of two given nodes to
    all others on the distance matrix

    Parameters
    ----------
    dist_matrix : np.ndarray
        the matrix displaying the distances between all combinations of two species
        all distances are floats between 0 and 1
    i : int
        the row number for one of the child nodes
    j : int
        the column number for the other child node

    Returns
    -------
    np.ndarray
        an array of the distances between this internal node to all other nodes on
        the distance matrix that don't descend from this one
        NOTE: the indices in this array will line up with the distance matrix and 
        the label list AFTER i and j are removed from the matrix)
    """

    # initialize the output
    num_species = dist_matrix.shape[0]
    distances = np.zeros(shape=num_species - 2, dtype=float)

    # save the distance from i to j for later
    dist_ij = dist_matrix[i,j]

    # iterate over the rows, excluding the rows for the two given species (i and j)
    rows_to_hit = [idx for idx in range(num_species) if idx not in (i,j)]
    for new_idx, curr_row in enumerate(rows_to_hit):

        # calc distance from internal node to this species in the dist matrix
        dist_to_curr = ( dist_matrix[i, curr_row] + dist_matrix[j, curr_row] - dist_ij ) / 2.0
        
        # store the distance for later
        distances[new_idx] = dist_to_curr
    
    return distances


def _is_symmetric(dist_matrix: np.ndarray) -> bool:
    """
    Function to check if the distances on the distance matrix are symmetric
    i.e. that matrix[i,j] is equal to matrix[j,i]

    Parameters
    ----------
    dist_matrix : np.ndarray
        2d numpy matrix of distances between species

    Returns
    -------
    bool
        boolean indicating True if the matrix is symmetric and False if not
    """
    # iterate through the positions on the matrix
    for i in range(dist_matrix.shape[0]):
        for j in range(dist_matrix.shape[1]):

            # check if the value at i,j is equal to the value at j,i
            if dist_matrix[i,j] != dist_matrix[j,i]:
                # if not, print out an error message and quit
                print(f"Distance matrix is not symmetric! \n"
                      f"Expected matrix[{i},{j}] to be equal to matrix[{j},{i}].\n"
                      f"matrix[{i},{j}] = {dist_matrix[i,j]}, while matrix[{j},{i}] = {dist_matrix[j,i]}.")
                return False
    # if we made it this far, we never encountered any asymmetries
    return True


def neighbor_joining(dist_matrix: np.ndarray, seq_labels: list[str], verbose: bool = False) -> Tree | None:
    """
    Function to create the tree structure for a group of sequences based
    on their distance matrix and the list of labels for the rows/columns
    of the distance matrix

    Parameters
    ----------
    dist_matrix : np.ndarray
        numpy matrix of distances between every pair of species
        has shape N by N, where N is the number of species in seq_labels
        Distances in this matrix should be symmetric dist_matrix[i,j] should equal dist_matrx[j,i]
    seq_labels : list[str]
        list of names for the species on each row/column of the sequence matrix

    Returns
    -------
    Tree
        the phylogenetc tree structure for the sequences in the distance matrix
    """

    # validate that the matrix has the right dimensions
    num_species = len(seq_labels)
    if dist_matrix.shape != (num_species, num_species):
        # this also catches if number of rows is not equal to number of columns
        # print for debugging
        print(f"Distance matrix does not match species labels! Expected shape: "
              f"{num_species} by {num_species}, and got: {dist_matrix.shape}.")
        # quit the func if this core assumption is volated
        return None
    
    # check if the distances are symmetric
    if not _is_symmetric(dist_matrix):
        return None
    
    # now that we're sure the matrix is valid, we can start
    phylo_tree = Tree()
    working_matrix = dist_matrix.copy()
    working_labels = seq_labels.copy()
    if verbose:
        print(f"Starting, empty tree:\n{phylo_tree}")
        print(f"Starting, working_matrix: \n{working_matrix}")
        print(f"Starting, working_labels: {working_labels}")

    # iterate until the distance matrix has a shape of 2x2
    while working_matrix.shape != (2,2):

        # make a q-matrix
        curr_q_matrix = make_q_matrix(working_matrix)
        if verbose:
            print(f"Q_matrix:\n{curr_q_matrix}")

        # grab the coordinates for the FIRST minimum to appear on the q-matrix, reading across rows
        i, j = np.unravel_index(np.argmin(curr_q_matrix), shape=curr_q_matrix.shape)
        if verbose:
            print(f"Identified minimum on q-matrix at [{i}, {j}]. These are {working_labels[i]}, {working_labels[j]}")

        # get the limb lengths for the species on rows i and j
        # NOTE: i and j are converted here from np.int64 to int, which can change their values
        # in cases of extreme overflow
        limblength_i, limblength_j = calc_limb_lengths(working_matrix, int(i), int(j))
        if verbose:
            print(f"Limb lengths are i:{limblength_i}, j:{limblength_j}")

        # get the distances from the parent to i and j to all other nodes on the dist matrix
        parental_distances = calc_internal_dist(working_matrix, int(i), int(j))

        # only create node objects for i and j if these species are not already on the graph structure
        if working_labels[i] in phylo_tree.nodes:
            node_i = phylo_tree.nodes[working_labels[i]]
        else:
            node_i = Node(working_labels[i], limblength_i)
            phylo_tree.add_node(node_i)
        
        if working_labels[j] in phylo_tree.nodes:
            node_j = phylo_tree.nodes[working_labels[j]]
        else:
            node_j = Node(working_labels[j], limblength_j)
            phylo_tree.add_node(node_j)

        if verbose:
            print(f"Current tree before adding parent: {phylo_tree}")

        # create a parent node on the graph for nodes i and j
        parent_name = phylo_tree.make_parent([node_i, node_j], min(parental_distances))
        if verbose:
            print(f"Current tree after adding parent: {phylo_tree}")

        # remove species i and j from the distance matrix AND the labels list
        # iterate through the row numbers for the two species from greatest to least so
        # the indices for rows to be removed don't need to be recalculated
        for idx in sorted([i,j], reverse=True):
            # delete column and then delete row
            working_matrix = np.delete(working_matrix, idx, axis=1)
            working_matrix = np.delete(working_matrix, idx, axis=0)
            del working_labels[idx]
        if verbose:
            print(f"After removing rows i and j, working matrix:\n{working_matrix}\n"
                f"With labels: {working_labels}")
        
        # add the name for the parent to the labels list
        working_labels.append(parent_name)

        # add a row at the end of the matrix for this row
        new_row = np.reshape(parental_distances, shape=(1,len(parental_distances)))
        if verbose:
            print(f"New row: {new_row}")
            print(f"Working matrix: \n{working_matrix}")
        working_matrix = np.vstack((working_matrix, new_row))
        
        # add a 0 in the last position to represent this node's distance to itself
        parental_distances = np.append(parental_distances, 0)

        # add a new column for this node as well
        new_col = np.reshape(parental_distances, shape=(len(parental_distances), 1))
        working_matrix = np.hstack((working_matrix, new_col))
    
    if verbose:
        print(f"Working matrix at end:\n{working_matrix}, with labels:{working_labels}")
    # now there's only two nodes left on the matrix, and >=1 is internal
    for idx, species_name in enumerate(working_labels):
        if species_name not in phylo_tree.nodes:
            # this works because the matrix is guaranteed to be 2x2 now, so the only options are idx=0 and idx=1
            other_idx = 1 - idx
            # create a node object for this species, using its distance to the only other species 
            # as its limb length
            new_node = Node(species_name, limb_length=working_matrix[idx,other_idx])
            phylo_tree.add_node(new_node)

    # now every node is on the graph, and we can't make any more parents because there is no outgroup for
    # those parents to calculate their limb length
    return phylo_tree


def get_mean_dist_matrix(matrix_list: list[np.ndarray]) -> np.ndarray:
    """
    Function to make a distance matrix where each position contains the mean 
    of that position's corresponding values in all the distance matrices provided
    E.g. if there are 3 matrices, position (i,j) on the mean matrix will be mean of the values on
    position (i,j) in the three matrices provided

    Parameters
    ----------
    matrix_list : list[np.ndarray]
        list of additive distance matrices between species
        NOTE: must all have same shape and be two-dimensional

    Returns
    -------
    np.ndarray
        2d additive distance matrix with same shape as each of the distance matrices provided
    """
    # get the shapes of all the matrices
    numrows, numcols = matrix_list[0].shape

    # initialize the output
    mean_matrix = np.zeros(shape=(numrows, numcols))

    # iterate through the positions in all the distance matrices
    for i in range(numrows):
        for j in range(numcols):
            # get the mean for this position on all matrices
            corresponding_values = [dist_mat[i,j] for dist_mat in matrix_list]
            curr_mean = np.mean(corresponding_values)
            # write this to the mean matrix
            mean_matrix[i,j] = curr_mean
    
    return mean_matrix


if __name__ == "__main__":
    # little demonstration with quick toy data
    # intended tree has A and D as siblings under one common ancestor with B and C as siblings under a different common ancestor
    print("Expected output: ((species_B,species_C),(species_A,species_D));")
    species = ["species_A", "species_B", "species_C", "species_D"]
    
    # gene whose sequences show A and D's common ancestor maintaining ancestral state, B and C's ancestor having a mutation, and
    # A+D and B+C each developing another mutation independent of each other (assume ancestral state is ATAAAAA)
    gene_one = {"species_A": "ATTGAAA",
                "species_B": "AAAACGA",
                "species_C": "ATAACGA",
                "species_D": "ATAGAAA"}
    # gene whose purpose is to see what happens if a mutation in B+C's ancestor is inverted in either of B or C
    # assume ancestral state is GAGGGGG
    gene_two = {"species_A": "GAGGTGG",
                "species_B": "GACGGGG",
                "species_C": "GCAGGGG",
                "species_D": "GAAATGG"}
    # gene designed to see if a deletion affects B and C being resolved as siblings
    # assume ancestral state is CCCCCCC
    gene_three = {"species_A": "CCCCAGC",
                  "species_B": "CCGGCCC",
                  "species_C": "CCGCC",
                  "species_D": "CTTCACC"}
    
    genes = [gene_one, gene_two, gene_three]

    trees = []

    for i, curr_seq_dict in enumerate(genes):
        print(f"Gene_{i+1}:")
        for k, v in curr_seq_dict.items():
            print(f"{k}:{v}")

        # make a distance matrix for this gene
        curr_dist_matrix = make_dist_matrix(species, curr_seq_dict)
        print(f"{species}\n{curr_dist_matrix}\nGene_{i+1} Q-matrix:\n{make_q_matrix(curr_dist_matrix)}")

        # make the tree
        curr_tree = neighbor_joining(curr_dist_matrix, species, False)
        trees.append(curr_tree)
        
    for i in range(len(genes)):
        print(f"Tree {i+1}:\n{trees[i]}")
