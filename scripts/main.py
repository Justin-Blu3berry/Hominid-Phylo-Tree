"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import modules
import argparse
import json
import numpy as np
from pathlib import Path

# import from my own modules
from file_parsing import str_to_dir, read_fasta
from tree_objects import Tree
from tree_building_utils import make_dist_matrix, neighbor_joining, get_mean_dist_matrix

def get_cli_args() -> argparse.Namespace:
    """
    Function to get the arguments from the command line

    Returns
    -------
    argparse.Namespace
        instance of argparse arguments
    """
    # create parser
    parser = argparse.ArgumentParser(description="Provide a JSON config file to generate phylogenetic trees.")

    # make an argument to take the config file
    parser.add_argument("--config", "-c",
                        dest="config_file",
                        type=str,
                        help="Path to the JSON config file including NCBI API Key, gene names, and species names",
                        required=True)
    
    # make an argument to get the path to the data directory
    parser.add_argument("--indir", "--data", "-d",
                        dest="data_dir",
                        type=str,
                        help="Path to the directory containing FASTA sequences to be analyzed.",
                        required=False)
    
    # make an argument to get the path to the outputs directory
    parser.add_argument("--outputs", "--outdir", "-o",
                        dest="outputs_dir",
                        type=str,
                        help="Path to directory where outputs will be saved.",
                        required=False)
    
    return parser.parse_args()


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

    # make a header line showing species names along the top of the matrix
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
        cells = [f"{distance:.3f}" for distance in dist_matrix[i,]]
        cells = f"\t".join(cells)
        # combine to create the string for this row on the file
        rows.append(f"{species_name}{indent}\t{cells}")
    # now join all the rows of the distance matrix together
    rows = f"\n".join(rows)

    # write the output to the file
    out_path.write_text(f"{str(tree)}\n{labels}\n{rows}")


if __name__ == "__main__":
    
    # get CLI options
    cli_args = get_cli_args()

    config_path = Path(cli_args.config_file)
    
    # make the output directory if it doesn't exist or wasn't provided
    output_dir = str_to_dir(path=cli_args.outputs_dir, default="outputs")

    # make the data directory if it doesn't exist or wasn't provided
    data_path = str_to_dir(path=cli_args.data_dir, default="data")

    # load in the config info
    with config_path.open() as config_file:
        config = json.load(config_file)
    
    api_key = config["API_key"]
    genes = [gene.lower() for gene in config["list_seq_ids"]]
    list_species = [species.lower() for species in config["list_species"]]

    matrices = []

    for gene in genes:
        print(f"Progress, {gene}: Opening fasta file for {gene}...")
        # read the sequences for all species' copies of this gene
        curr_fasta = data_path / f"{gene}.fna"
        seq_dict = read_fasta(curr_fasta)

        print(f"Progress, {gene}: building distance matrix...")

        # make a distance matrix for this gene
        dist_matrix = make_dist_matrix(list_species, seq_dict)
        matrices.append(dist_matrix)

        print(f"Progress, {gene}: building phylogenetic tree...")

        # make a tree
        curr_tree = neighbor_joining(dist_matrix, list_species)

        print(f"Progress, {gene}: writing Newick string and distance matrix to file...")

        # make a path for the output file
        output_path = output_dir / f"{gene}_outputs.text"
        write_outputs(output_path, curr_tree, list_species, dist_matrix)  # type: ignore

        print(f"Finished {gene}!")
    
    print("Progress: creating mean distance matrix...")
    # get the mean distance matrix
    mean_matrix = get_mean_dist_matrix(matrices)

    print("Progress: creating tree based on mean distance matrix...")
    # make a tree
    curr_tree = neighbor_joining(mean_matrix, list_species)

    print("Progress: writing matrix and resulting tree to file...")
    # make a path for the output file
    output_path = output_dir / f"mean_outputs.text"
    write_outputs(output_path, curr_tree, list_species, mean_matrix)  # type: ignore

    print("Finished mean matrix and tree!")

        


    # TODO: use biopython Phylo package to visualize the tree
    # TODO: save image to outputs
    # TODO: repeat for all genes
    # TODO later: get summary stats for every combo of species on dist matrix
    # TODO later: make dist matrix + tree + image of mean dist matrix
    # TODO later: compare groupings of species
    # TODO later: generate report of comparisons
