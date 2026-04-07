# I'm feeling too lazy to make the toy data by hand

# imports
import random
from tree_objects import Node, Tree

def make_fake_tree(num_generations: int, 
                   leaf_names: list[str], 
                   ancestral_seq: str|None = None, 
                   default_length: int = 15) -> tuple[str, dict[str, str]]:
    """
    Function to generate a tree structure with num_generations generations worth of nodes
    with sequences for each leaf

    Parameters
    ----------
    num_generations : int
        number of generations through which to preceed with the tree
    leaf_names: list[str]
        list of names for the leaves the tree needs to end up with
    ancestral_seq : str | None, optional, 
        DNA sequence for the species ancestral to the whole rest of the tree
        by default None
        If nothing is passed, we generate a random DNA sequence
    default_length: int
        Default length for a sequence to be generated if no ancestral_seq is provided

    Returns
    -------
    str
        Newick String representation for the tree generated
    dict[str,str]
        Dictionary mapping the names of the leaves to their generated sequences
    """
    # make the ancestral sequence
    if ancestral_seq is None:
        # generate a random sequence, look at this beautiful syntax </3
        ancestral_seq = "".join(random.choices(("A", "C", "G", "T"), k=default_length))

    # now make a node object for the ancestor
    ancestor = Node("ancestor", limb_length=0)

    # 