"""
Justin Wildman, BINF6251 Final Project: Hominin Phylogenetic Forest
"""

# import from modules
from collections import defaultdict


class Node:
    """
    Object representing nodes on a phylogenetic tree. Works for both terminal (leaves) and internal nodes on the tree.
    A name for the node and the distance to its parent are mandatory to create an object of this class

    Attributes:
        name: str, required; this node's name
              If it's internal, this node's name will be the Newick string representation of the subtree that this node
              is ancestral to
        limb_length: float, required; this node's distance to its parent
        parent_name: str, optional; the name of this node's parent
                     This is only updated if this node's parent can be resolved and made into a node object
    """
    def __init__(self, name: str, limb_length:float, parent_name: str = "unnamed_parent"):
        self.name = name
        self.limb_length = limb_length
        self.parent_name = parent_name
    
    def __repr__(self) -> str:
        """
        Makes a string representation of this ndoe in a way that's as conducive to making a Newick
        string as possible

        Returns
        -------
        str
            Newick-adjacent string representing this node
        """
        return f"{self.name}:{round(self.limb_length,4)}"
    

class Tree:
    """
    Object representing a phylogenetic tree as a graph structure

    Attributes:
        nodes: dict[str, Node]
            Dict mapping node names to their Node objects
        top_layer: list[Node]
            List of node objects in the tree that don't have parents
        edges: dict[Node, Node]
            Dict mapping Node objects to the Nodes they are parental to
    Methods:
        add_node(self, new_node: Node, children: list[Node] = []) -> None
            Function to add a new node to the graph without connections to anything else
        get_parentless_nodes(self) -> list[Node]
            Function to identify nodes on the graph that don't have any parents 
            (and therefore belong on self.top_layer)
        make_parent(self, children: list[Node], dist_to_outgroup: float, verbose: bool = True) -> str
            Function to make a new internal node that is parental to all of the Nodes in children
        __repr__(self) -> str
            Function to create the Newick String representation of the entire graph
    """
    def __init__(self, node_dict: dict[str, Node] | None = None,
                 top_layer: list[Node] | None = None,
                 node_mapping: dict[Node, list[Node]] | None = None):
        
        # set attributes
        self.nodes = node_dict if node_dict else {}
        self.top_layer = top_layer if top_layer else []
        self.edges = defaultdict(list)

        # check if the top layer is empty but self.nodes isn't
        if self.top_layer is [] and self.nodes:
            self.top_layer = self.get_parentless_nodes()
        
        # attempt unpack its pairings into the defaultdict (does nothing if node_mapping is empty)
        if node_mapping:
            for node, children in node_mapping.items():
                self.edges[node] = children
        
    
    def add_node(self, new_node: Node, children: list[Node] = []) -> None:
        """
        Function to add a new node to the graph without connections to anything else

        Parameters
        ----------
        new_node : Node
            Node object being added to the graph
        children: list[Node], optional
            list of Node objects that this new node is the parent to
        """
        # add the new node to the nodes list and to the top layer
        self.nodes[new_node.name] = new_node
        self.top_layer.append(new_node)
        
        # add this node to the edges dictionary if any children are specified
        if children:
            self.edges[new_node] = children

    
    def get_parentless_nodes(self) -> list[Node]:
        """
        Function to identify nodes on the graph that don't have any parents (and therefore belong on self.top_layer)
        """

        # make a shallow copy of the self.nodes list
        parentless_nodes = [node for node in self.nodes.values()]

        # iterate through the lists of edges in self.edges
        for curr_edges in self.edges.values():
            # iterate through the children in the current list of edges
            for child_node in curr_edges:
                # check if the current node from the edge list is in the list of parentless nodes
                if child_node in parentless_nodes:
                    # remove this node from the list of parentless nodes because this one clearly has a parent
                    parentless_nodes.remove(child_node)
        
        # return list of whatever nodes remain after removing every node that's a child to something else
        return parentless_nodes
    

    def make_parent(self, children: list[Node], dist_to_outgroup: float, verbose: bool = True) -> str:
        """
        Function to make a new internal node that is parental to all of the Nodes in children

        Parameters
        ----------
        children : list[Node]
            list of Nodes that this new internal node is ancestral to
        dist_to_outgroup : float
            the distance from this internal node to a leaf it doesn't descend from is used as a proxy for its limb length.
            This gets updated later if another node is made that is parental to this.

        Returns
        -------
        str
            the name for this parental node so the Neighbor Joining function can access this new Node object
        """

        # generate a name for this new internal node
        parent_name = _generate_parent_name(children)
        
        # create a node object for this internal node
        new_node = Node(parent_name, dist_to_outgroup)
        self.add_node(new_node, children)

        # iterate through the child nodes
        for child_node in children:
            # attempt to remove this child from the top layer
            try:
                self.top_layer.remove(child_node)
            except ValueError:
                if verbose:
                    print(f"Attempted to remove {child_node} from tree's top layer, but it isn't already in the top layer.\n"
                        f"Top layer at time of attempted removal: {self.top_layer}.\n"
                        f"Nodes dict at time of attempted removal: {self.nodes}.\n"
                        f"Edges at time of attempted removal: {self.edges}.")

            # also update the children's parental names too while we're here
            child_node.parent_name = new_node.name
        
        # even though we've already made the node object and that's what's important, we return the name
        # so the neighbor joining function can access this internal node later on
        return parent_name
    

    def __repr__(self) -> str:
        """
        Function to create the Newick String representation of the entire graph

        Returns
        -------
        str
            the Newick String representation of the tree
        """
        # get the name for a hypothetical parent to all the nodes in the top layer
        # cap it off with semicolon to indicate a completed tree
        return _generate_parent_name(self.top_layer) + ";"



def _generate_parent_name(children: list[Node]) -> str:
    """
    Function to generate a name for an internal node based on what its children are
    This node's name is going to be the Newick string representation of the subtree that descends from this node

    NOTE: because internal nodes are resolved from bottom to top AND because all internal nodes are getting named
            like this, there is actually no need to recurse through the tree to make the newick string representation
            for either the main tree or any subtree inside of it

    Parameters
    ----------
    children : list[Node]
        list of Nodes that this one is directly ancestral to

    Returns
    -------
    str
        the Newick String representation of the subtree that descends from this tree
    """
    # initialize the name for this parent
    parent_name = "("

    # concatenate on the string representations of the children
    for child_node in children:
        # this concatenates "{child_name}:{limb_length}" thanks to how we defined __repr__ for Node objects
        if child_node == children[-1]:
            parent_name += str(child_node)
        else:
            parent_name += str(child_node) + ","

    # cap off at the end with a closing parenthesis
    return parent_name + ")"
