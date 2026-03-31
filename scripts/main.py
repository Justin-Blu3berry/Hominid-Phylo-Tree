"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import modules
import argparse
from pathlib import Path

# import from my own modules
from file_validation import str_to_dir

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
    output_path = str_to_dir(path=cli_args.outputs_dir, default="outputs")

    # make the data directory if it doesn't exist or wasn't provided
    data_path = str_to_dir(path=cli_args.data_dir, default="data")
    
    print(config_path)
    print(data_path)
    print(output_path)

    # TODO: read the config file (leave addressing assumption that the file exists and is filled out directly for later)
    # TODO: make ordered list of species and ordered list of genes from config
    # TODO LATER: validate that all the genes' files exist and have all the species we care about
    # TODO: iterate over gene names in the config file
    # TODO: open current gene's fasta file
    # TODO: make dict of species_name: sequence for current gene
    # TODO: make distance matrix
    # TODO: make tree from distance matrix (maybe just save the str(tree) bc we don't need the full object)
    # TODO: write distance matrix & tree string to file
    # TOOD: use biopython Phylo package to visualize the tree
    # TODO: save image to outputs
    # TODO: repeat for all genes
    # TODO later: get summary stats for every combo of species on dist matrix
    # TODO later: make dist matrix + tree + image of mean dist matrix
    # TODO later: compare groupings of species
    # TODO later: generate report of comparisons
