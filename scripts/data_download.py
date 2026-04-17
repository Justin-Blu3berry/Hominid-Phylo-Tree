# this is the dreaded "query ncbi refseq data, format it perfectly, and don't make any mistakes" script

# imports
import argparse
import json
from pathlib import Path
from Bio import Entrez

from filewriting import write_fasta


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
    
    # make an argument to take the config file
    parser.add_argument("--config", "-c",
                        dest="config_file",
                        type=str,
                        help="Path to the JSON config file including NCBI API Key, gene names, and species names",
                        required=True)
    
    # make an argument to get the path to the data directory
    parser.add_argument("--outdir", "--data", "-o",
                        dest="data_dir",
                        type=str,
                        help="Path to the directory where FASTA files will be written..",
                        required=True)
    
    return parser.parse_args()


def search_sequence(species_name: str, gene_name: str, max_attempts: int = 20, max_length: int = 50000, verbose: bool = False) -> str:
    """
    Function to query NCBI's gene database for a given species's sequence
    for the named gene

    Parameters
    ----------
    species_name : str
        species whose sequence is being searched for 
    gene_name : str
        name of the gene whose sequence we're finding
    max_attempts: int, optional
        Maximum number of search results to search through before giving up on finding the sequence
        Will not do more attempts than there are search results
        Default is 20
    max_length: int, optional
        filters out a search result for having a sequence that is too long
        This avoids a sequence returning an entire chromosome or more than just the gene
        Default is 50,000 bp
    verbose: bool, optional
        Indicates if debugging messages are printed
        Default is False

    Returns
    -------
    str
        nucleotide sequence for the gene
    """
    if verbose:
        print(f"Starting search for {gene_name} in {species_name}")
    # make search string using nucleotide's advanced search syntax
    search_string = f"{gene_name}[GENE] AND {species_name}[Organism]"
    # generate the handle for the search
    handle = Entrez.esearch(db="nucleotide", term=search_string, retmax=max_attempts)
    
    # convert handle to a dict
    ids = Entrez.read(handle)["IdList"]  # type: ignore

    # make sure we don't do more attempts than we have search results
    if len(ids) < max_attempts:
        max_attempts = len(ids)
        if verbose:
            print(f"Max attempts adjusted to {max_attempts} to match number of search results")

    # attempt the following until we run out of attempts or find a reasonable seq
    for i in range(max_attempts):
    
        # grab the ID for the first search result
        curr_id = ids[i]  # type: ignore

        if verbose:
            print(f"Candidate ID for {gene_name} in {species_name}: {curr_id}")

        # make another search to validate length without having to do string processing
        newhandle = Entrez.efetch(db="nucleotide",id=curr_id, rettype="gb", retmode="text")
        genbank_info = newhandle.readline().strip().split()
        if verbose:
            print(f"genbank info for this gene: {genbank_info}")
        
        # grab the sequence's length from the genbank info
        length = int(genbank_info[2])
        
        if verbose:
            print(f"Length of sequence in current search result: {length}, max length is {max_length}\n"
                  f"Length > max_length: {length > max_length}")

        # if the sequence is safe, we can export this sequence
        if length <= max_length:
            # search up the fasta sequence
            newhandle = Entrez.efetch(db="nucleotide",id=curr_id, rettype="FASTA", retmode="text")
            
            # convert the result to a string (format is handle -newline- sequence -newline-, which is how we want to write it into a file)
            outstring = newhandle.read()

            return outstring
    
    # report that we failed to find 
    print(f"Could not find sequence for gene {gene_name} in species {species_name} in {max_attempts} search results")
    return f">null: {species_name} {gene_name} NOT FOUND\n\n"


if __name__ == "__main__":

    # read the command line args
    cli_args = get_cli_args()
    config_path = Path(cli_args.config_file)
    data_path = Path(cli_args.data_dir)

    # read config file
    with config_path.open() as config_file:
        config = json.load(config_file)
    Entrez.api_key = config["API_key"]
    genes = [gene.strip().lower() for gene in config["list_gene_names"]]
    list_species = [species.strip().lower() for species in config["list_species"]]
    Entrez.email = config["email"].strip().lower()

    # iterate over the genes
    for gene_name in genes:
        # make a path for the output
        fasta_path = data_path / f"{gene_name}.fna"

        with fasta_path.open(mode="w", encoding="utf-8") as outfile:
        
            # iterate through the species of interest
            for species_name in list_species:

                # query NCBI to find a DNA sequence in fasta format for the named gene and species
                curr_entry = search_sequence(species_name, gene_name)
                # the formatting we get from NCBI automatically separates header and sequences for us
                outfile.write(curr_entry)
            