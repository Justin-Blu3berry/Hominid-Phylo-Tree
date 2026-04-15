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
from filewriting import write_outputs, plot_tree

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
    genes = [gene.strip().lower() for gene in config["list_gene_names"]]
    list_species = [species.strip().lower() for species in config["list_species"]]
    list_species = ["_".join(full_name.split()) for full_name in list_species]

    matrices = []

    # initialize a dict tracking which species have sequences for which genes (important for later)
    seq_labels = {}
    excluded_species = {}

    for gene in genes:
        print(f"Progress, {gene}: Opening fasta file for {gene}...")
        # read the sequences for all species' copies of this gene
        curr_fasta = data_path / f"{gene}.fna"
        seq_dict = read_fasta(curr_fasta)
        print({key: len(value) for key, value in seq_dict.items()})

        print(f"Progress, {gene}: building distance matrix...")

        # if any species failed to pull up a sequence for this gene during data download, we want to safely ignore
        # the species in question
        seq_labels[gene] = [species for species in list_species if species in seq_dict]
        
        # tell the user if a species lacks a sequence for a given gene
        curr_excluded = set(list_species).difference(seq_labels[gene])
        if curr_excluded:
            print(f"Following species were excluded from distance matrix for {gene} due to not having a sequence: {curr_excluded}")
            excluded_species[gene] = curr_excluded

        # make a distance matrix for this gene
        dist_matrix = make_dist_matrix(seq_labels[gene], seq_dict)
        matrices.append(dist_matrix)

        print(f"Progress, {gene}: building phylogenetic tree...")

        # make a tree
        curr_tree = neighbor_joining(dist_matrix, seq_labels[gene])

        print(f"Progress, {gene}: writing Newick string and distance matrix to file...")

        # make a path for the output file
        output_path = output_dir / f"{gene}_outputs.text"
        write_outputs(output_path, curr_tree, seq_labels, dist_matrix)  # type: ignore

        print(f"Progress, {gene}: saving image of generated tree...")

        # visualize the tree
        image_path = output_dir / f"{gene}_tree.png"
        plot_tree(curr_tree, image_path, gene)  # type: ignore

        print(f"Finished {gene}!")
    
    # tell the user there are unrepresented species that they requested in their config
    if excluded_species:
        print(f"The following genes were missing for the following species. \n{excluded_species}\n"
              f"Rerun data_download.py with {config_path}, pointing to {data_path}, or remove the named species from {config_path}"
              "to ensure accuracy of the mean tree.")

    print("Progress: creating mean distance matrix...")
    # get the mean distance matrix
    mean_matrix = get_mean_dist_matrix(matrices)

    print("Progress: creating tree based on mean distance matrix...")
    # make a tree
    curr_tree = neighbor_joining(mean_matrix, seq_labels[genes[0]])

    print("Progress: writing matrix and resulting tree to file...")
    # make a path for the output file
    output_path = output_dir / "mean_outputs.text"
    write_outputs(output_path, curr_tree, seq_labels[genes[0]], mean_matrix)  # type: ignore

    print("Progress: saving image of tree generated with mean distance matrix...")

    # visualize the tree
    image_path = output_dir / "mean_tree.png"
    plot_tree(curr_tree, image_path, "mean distances")  # type: ignore

    print("Finished mean matrix and tree!")

        
    # TODO later: get summary stats for every combo of species on dist matrix
    # TODO later: compare groupings of species
    # TODO later: generate report of comparisons
