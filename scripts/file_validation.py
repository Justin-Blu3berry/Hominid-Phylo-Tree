"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import from modules
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
    if path is None or path is "":
        # if no path was provided, attempt to create one using the default name
        _attempt_mkdir(default)
    
        # now we know the "data" directory exists, reassign the data path
        return Path(default)
    
    else:
        # if the data path was provided, just attempt to make the directory as named 
        _attempt_mkdir(path)

        # convert it to a path object
        return Path(path)
