# imports
from pathlib import Path


def _attempt_mkdir(pathname: str) -> None:
    """
    Funtion to make a relative path from the present working directory to the named path, prints to stdout if
    the path already exists

    Parameters
    ----------
    pathname : str
        the path to be made (if it doesn't already exist)
    """
    # create a new path object relative to the present working directory
    relative = Path(pathname)
    # attempt to make the directory
    try:
        Path.mkdir(relative)
    except FileExistsError:
        print(f"Tried to make directory: {relative}, but failed because it already exists")


def str_to_dir(path: str, default: str) -> Path:
    """
    Function that takes a relative path to a directory as a string and attempts to make the named directory
    If the directory's name is empty (none or empty string), a path is created using the default name instead

    Parameters
    ----------
    path : str
        _description_
    default : str
        _description_

    Returns
    -------
    Path object pointing to the directory of interest
        Either points to ./<path> or ./<default> if path is None or an empty string
        These files get pointed to regardless of if the path already existed or was just created
    """
    if path is None or path == "":
        # if no path was provided, attempt to create one using the default name
        _attempt_mkdir(default)
    
        # now we know the "data" directory exists, reassign the data path
        return Path(default)
    
    else:
        # if the data path was provided, just attempt to make the directory as named 
        _attempt_mkdir(path)

        # convert it to a path object
        return Path(path)


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

    # grab the sequence for the last species
    seq_dict[curr_species] = curr_seq

    return seq_dict
