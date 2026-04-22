# Module that simulates sequence evolution to generate "realistic" toy sequences at varying levels of conservation

# imports
import argparse
import numpy as np
from textdistance import jaccard
from pathlib import Path

from tree_objects import Node, Tree
from filewriting import write_fasta, plot_tree


def get_cli_args() -> argparse.Namespace:
    """
    Function to get the arguments from the command line

    Returns
    -------
    argparse.Namespace
        instance of argparse arguments
    """
    # create parser
    parser = argparse.ArgumentParser(description="Provide a file path to write toy data to.")
    
    # make an argument to get the path to the data directory
    parser.add_argument("--outdir", "--data", "-o",
                        dest="outputs_dir",
                        type=str,
                        help="Path to the directory where FASTA files will be written..",
                        required=True)
    
    # make an argument to get the path to the data directory
    parser.add_argument("--nseq", "-n",
                        dest="num_seqs",
                        type=int,
                        help="Number of genes to generate sequences for.",
                        required=False)
    
    return parser.parse_args()


def make_fake_tree() -> Tree:
    """
    Due to time contraints, this function  generates a hard-coded tree instead of letting you parameterize
    the number of generations this should run for (+ number of extant species)

    Returns
    -------
    Tree
        A tree object for the phylogeny of the species in the toy data
    """

    # make a pseudo-graph structure mapping species names to their children's names
    name_graph = {"ancestor": ("ABCDE", "great_grandparent_F"),
                    "ABCDE": ("ABC", "parent_DE"),
                    "great_grandparent_F": ("grandparent_F",),
                    "ABC": ("parent_A", "BC"),
                    "parent_DE": ("DE",),
                    "grandparent_F": ("parent_F",),
                    "parent_A": ("A",),
                    "BC": ("B", "C"),
                    "DE": ("D", "E"),
                    "parent_F": ("F",)}
    
    # make a list of the names for all the nodes that will go onto the graph
    species_names = ("ancestor", "ABCDE", "great_grandparent_F", "ABC", "parent_DE", "grandparent_F","parent_A", "BC", "DE", "parent_F",
                     "A", "B", "C", "D", "E", "F")
    
    # make node objects for each of these species, with a placeholder limb length
    nodes_dict = {species_name: Node(species_name, 0) for species_name in species_names}

    node_mapping = {}

    # iterate through all of the species that have children
    for parent_name, child_names in name_graph.items():
        # access this parent's node
        curr_node = nodes_dict[parent_name]

        # get the node objects for this species's children
        children = [nodes_dict[child_name] for child_name in child_names]

        # map this node obejct to its children
        node_mapping[curr_node] = children
    
    # top layer
    top_layer = [nodes_dict["ancestor"]]

    return Tree(nodes_dict, top_layer, node_mapping)


def generate_ancestral_seqs(rng: np.random.Generator,
                            desired_len: int = 15,
                            num_genes: int = 5,
                            ) -> list[str]:
    """
    Function to generate num_genes random DNA sequences of desired_len

    Parameters
    ----------
    rng : np.random.Generator
        numpy.random Generator object, provided for reproducibility if a default seed is used
    desired_len : int, optional
        desired length of DNA sequences, by default 15
    num_genes : int, optional
        desired number of DNA sequences, by default 5

    Returns
    -------
    list[str]
        list of length num_genes, where each item is a DNA sequence of length desired_len
    """

    # generate num_genes random gene sequences
    ancestral_seqs = []
    for i in range(num_genes):
        random_seq = "".join(rng.choice(("A", "C", "G", "T"), desired_len))
        ancestral_seqs.append(str(random_seq))
    
    return ancestral_seqs


def mutated(sequence: str,
            rng: np.random.Generator,
            mutation_range: tuple[int, int] = (0,3),
            size_range: tuple[int, int] = (1, 3),
            mut_types: tuple[str, str, str, str, str, str] | None = None) -> str:
    """
    Function to generate a randomly-mutated version of the given sequence

    Parameters
    ----------
    sequence : str
        DNA sequence being mutated
    rng : np.random.Generator
        numpy.random Generator object, provides RNG seed and allows for reproducibility
    mutation_range: tuple[int, int], optional
        indicates the minimum (at index 0) and maximum (at index 1) number of mutations to perform
        Inclusive on both sides
        Default is 0 to 3
    size_range: tuple[int, int], optional
        indicates the minimum (at index 0) and maximum (at index 1) number of nucleotides that each mutation can affect
        Inclusive on both sides
        Default is 1 to 3
    mut_types: tuple[str], optional
        indicates the types of mutations that are allowed
        By default, it includes: "insert", "delete", "duplicate", "invert", "substitute", "translocate"

    Returns
    -------
    str
        mutated DNA sequence
    """
    # handle default value for mut_types (too long to put in definition line)
    mut_types = mut_types if mut_types else ("insert", "delete", "duplicate", "invert", "substitute", "translocate")

    # select the number of mutations to perform
    num_mutations = rng.integers(mutation_range[0], mutation_range[1], endpoint=True)

    # make a working copy of the sequence
    working_seq = sequence

    nucleotides = ("A", "C", "G", "T")

    # repeat the following for however many mutations we're choosing to do
    for i in range(num_mutations):

        # randomly select a mutation
        curr_mutation = rng.choice(mut_types)
        
        mut_size = rng.integers(size_range[0], size_range[1], endpoint=True)
        mut_location = rng.integers(0, len(working_seq) - mut_size)

        # do a switch case to make the mutation happen (different mutations need different params)
        if curr_mutation == "insert":
            inserted_seq = rng.choice(nucleotides)
            working_seq = working_seq[:mut_location] + inserted_seq + working_seq[mut_location:]

        elif curr_mutation == "delete":
            # make a deletion of the selected size at the selected location
            working_seq = working_seq[:mut_location] + working_seq[mut_location + 1:]

        elif curr_mutation == "duplicate":
            # use a list slice to duplicate the substring
            duplicated = working_seq[mut_location: mut_location + mut_size] * 2
            working_seq = working_seq[:mut_location] + duplicated + working_seq[mut_location + mut_size:]

        elif curr_mutation == "invert":
            # use a list slice to invert a section of the sequence
            inverted = working_seq[mut_location + mut_size: mut_location: -1]
            working_seq = working_seq[:mut_location] + inverted + working_seq[mut_location + mut_size:]

        elif curr_mutation == "substitute":
            # substitute with a random nucleotide
            new_nuc = str(rng.choice(nucleotides, size=1)[0])
            working_seq = working_seq[:mut_location] + new_nuc + working_seq[mut_location + 1:]

        elif curr_mutation == "translocate":
            # get the location the substring will be moved to
            new_location = rng.integers(0, len(sequence))
            # grab the part that needs to be moved
            mobile_element = working_seq[mut_location: mut_location + mut_size]
            # remove the mobile element from the working seq
            working_seq =  working_seq[:mut_location] + working_seq[mut_location + mut_size:]
            # insert the mobile element back in at its new location
            working_seq = working_seq[:new_location] + mobile_element + working_seq[new_location:]
    
    return working_seq


def fix_internal_node_names(tree: Tree,
                            layers: tuple[tuple[str, ...], ...]) -> None:
    """
    Function to change the names of the internal nodes in tree to be the Newick representations
    of the subtree that descends from them (this makes the shortcut for newick string generation
    continue to hold up without having to use recursion)

    Parameters
    ----------
    tree : Tree
        tree object whose internal nodes are having their names changed
    layers : tuple[tuple[str, ...]]
        2d array of node names in order from top layer (tuple in position 0) to leaves (tuple in position -1)
        A node whose name is in layer i will find its parent on layer i-1 and its children on layer i+1
    """
    # iterate through the layers from bottom to top (ignoring leaves)
    for layer in layers[len(layers)-2:: -1]:
        for species_name in layer:
            # grab the children for this species
            curr_node = tree.nodes[species_name]
            children = tree.edges[curr_node]

            # update this node object's name based on its children
            new_name = tree._generate_parent_name(children)
            tree.change_node_name(curr_node, new_name)


def evolve_sequences(tree: Tree, 
                     ancestral_seq: str,
                     rng: np.random.Generator,
                     mutation_range: tuple[int, int] = (0, 3)) -> tuple[str, dict[str, str]]:
    """
    Function to generate a tree structure with a hard-coded number of generations worth of nodes
    with sequences for each leaf
    Unfortunately, due to time constraints, this function has to use a hard-coded graph

    Parameters
    ----------
    tree: Tree
        The tree object containing the species that are having sequences generated
    ancestral_seqs : str
        DNA sequence for the species ancestral to the whole rest of the tree
    rng : np.random.Generator
        numpy.random Generator object, provided for reproducibility if a default seed is used
    mutation_range: tuple[int, int], optional
        minimum and maximum number of mutations to do on each child seq (default is 0 to 3)

    Returns
    -------
    dict[str: float]
        Limb lengths calculated for each node on the tree
    dict[str,str]
        Dictionary mapping the names of the leaves to their generated sequences
    """
    
    # initialize empty sequences
    sequences = {species: "" for species in tree.nodes}
    sequences["ancestor"] = ancestral_seq
    limb_lengths = {species: 0.0 for species in tree.nodes}
    
    # make tuples to represent the species names on each layer of the graph
    layers = (("ancestor",), 
              ("ABCDE", "great_grandparent_F"), 
              ("ABC", "parent_DE", "grandparent_F"),
              ("parent_A", "BC", "DE", "parent_F"),
              ("A", "B", "C", "D", "E", "F"))
    
    # iterate though the species that will have (go layer by layer because dicts are technically unordered)
    for layer in layers[:-1]:
        for curr_species in layer:
            # access this species's node so its children can be accessed
            parent_node = tree.nodes[curr_species]
            parent_seq = sequences[curr_species]
            for child in tree.edges[parent_node]:
                # make one mutated sequence for each of this species's children
                child_seq = mutated(parent_seq, rng, mutation_range)
                
                sequences[child.name] = child_seq

                # calculate the "true" limb length for the child based on the Smith Waterman alignment with its parent
                curr_limb_length = jaccard.normalized_distance(child_seq, parent_seq)
                # overwrite the limb length recorded for this node on the tree
                child_node = tree.nodes[child.name]
                child_node.limb_length = curr_limb_length

    # fix the internal node names so that str(tree) is actually the newick string instead of "(ancestor:0);"
    fix_internal_node_names(tree, layers)

    return str(tree), sequences


if __name__ == "__main__":
    # get CLI options
    cli_args = get_cli_args()

    data_path = Path(cli_args.outputs_dir)

    num_seqs = cli_args.num_seqs if cli_args.num_seqs else 5

    rng = np.random.default_rng(42)

    # generate ancestral sequences
    ancestral_seqs = generate_ancestral_seqs(rng, desired_len=30, num_genes=num_seqs)

    # initialize
    newick_strings = []

    # evolve sequences
    for i, sequence in enumerate(ancestral_seqs):
        # make/reset the fake tree
        tree = make_fake_tree()

        # evolve the sequences on the tree (make each gene progressively less conserved)
        newick_string, sequences = evolve_sequences(tree, sequence, rng, mutation_range=(i, i+1))

        # write the sequences to a fasta file
        fasta_path = data_path / f"toy_gene_{i}.fna"
        write_fasta(sequences, fasta_path)

        newick_strings.append(newick_string)

        # save the image for this expected tree
        image_path = data_path / f"expected_gene_{i}.png"
        plot_tree(tree, image_path, f"expected for toy_gene_{i}")
    
    # write the expected newick string outputs
    trees_path = data_path / "expected_outputs.txt"
    text_to_write = [f"toy_gene{i}:\n{newick_strings[i]}\n" for i in range(len(newick_strings))]
    text_to_write = "".join(text_to_write)
    trees_path.write_text(text_to_write)

    print(f"Finished writing toy data. Check {data_path} for outputs. "
           "Gene 0 is expected to be the most conserved, and gene 4 is expected to be the least conserved")
