"""
DataStructures.py
Module of custom data structures for statistics computation.
"""

def partition_of_multi_indices(multiset,N_of_var=None):
    """
    Converts a partition of multisets into a partition of multi index
    See K.D. Smith's paper "A Tutorial on Multivariate k-statistics ..." for definitions
    
    Parameters
    ----------
    multiset : list of tuples
    N_of_var : int
        The number of variables to be used in the multi index representation
        Ex: 
            In a multiset representation 
                <x^2> is (0,0)  , <xy> is (0,1)   and <y^2> is (1,1)
            Their equivalent in a multi-indices representation are
                <x^2> is (2,0)  , <xy> is (1,1)   and <y^2> is (0,2)
            If N_of_var = 3 then the multi-indices representation are
                <x^2> is (2,0,0), <xy> is (1,1,0) and <y^2> is (0,2,0)
    """
    if N_of_var is None :
        m = 0
        N_of_var = max([max(tpl) for tpl in multiset]) + 1
    mlt_idx = []
    for tpl in multiset :
        idx_tpl=()
        for i in range(N_of_var):
            idx_tpl += (tpl.count(i),)
        mlt_idx += [idx_tpl,]
    return mlt_idx
    
def multi_idxs_to_multiset(multi_idxs):
    """
    Converts a multi index into a multisets
    See K.D. Smith's paper "A Tutorial on Multivariate k-statistics ..." for definitions
    
    Parameters
    ----------
    multi_idxs : tuple
    
    Ex: 
        In the 1D case the multi indices for <x^2> is (2,)  and the multiset is       (0,0)
        In the 2D case the multi indices for <x^2> is (2,0) and the multiset is still (0,0)
                       the multi indices for <xy>  is (1,1) and the multiset is still (0,1)
    """
    mltset = ()
    for i,idx in enumerate(multi_idxs) :
        for j in range(idx):
            mltset += (i,)
    return mltset

class IntPartitionTree:
    """
    Data structure for storing coefficients associated with integer partitions.
    The data structure looks up an integer partition in the following manner:
        1. Sort the parts in ascending order.
        2. Look up the child node of the root node associated with the first part.
        3. Recursively apply step 2 until the node associated with the full partition is reached.

    Attributes
    ----------
    min_branch : int
        Minimum of the parts associated with each child.
        Since parts are increasing along paths in the tree, we can assume a minimum value.
    value : object
        Value associated with the int partition corresponding to the root node.
    assume_sorted : bool
        Whether or not to assume the int partition arguments come pre-sorted.
    n_children : int
        Number of children of the root node.
    children : list of IntPartitionTree or None
        Children of the root node. None represents an empty child.

    Methods
    -------
    get_coef(int_partition)
        Look up the coefficient associated with an integer partition.
        Return None if the value hasn't been defined yet.
    set_coef(int_partition, value)
        Set the value of the coefficient associated with the integer partition.
    n_nodes()
        Count the number of nodes in the tree.
    depth()
        Compute the depth of the tree.

    Notes
    -----
    This tree does not have a significant memory footprint, even when evaluating
    high-order k-statistics. When evaluating a k-statistic of order i, the
    number of values stored in this tree is the partition number p(i), while the
    number of power sum products that must be computed is the Bell number B(i).
    B(i) scales much, much faster than p(i). For example, computing a 20th
    order k-statistic would save p(20) = 627 entries in this tree, while also
    requiring B(k) ~ 10^14 computations of power sum products. In short: the
    k-statistic computation becomes intractable long before this tree grows to
    a concerning size.
    """

    def __init__(self, min_branch=1, value=None, assume_sorted=False):
        """
        Initialize either a new IntPartitionTree or a child in an IntPartitionTree.

        Parameters
        ----------
        min_branch : int, optional
            Minimum possible value of parts in partitions starting at this node.
            Used for creating child nodes.
            To create an IntPartitionTree, leave at the default value of 1.
        value : object
            Value associated with partitions terminating at this node.
        assume_sorted : bool
            Whether or not to assume the int partition arguments come pre-sorted.
        """
        self.min_branch = min_branch
        self.value = value
        self.assume_sorted = assume_sorted
        self.n_children = 0
        self.children = []

    def get_coef(self, int_partition):
        """
        Look up the coefficient associated with an integer partition.

        Parameters
        ----------
        int_partition : sequence of ints
            Integer partition to look up.
        """
        if len(int_partition) == 0:
            return self.value
        sorted_int_partition = int_partition if self.assume_sorted else sorted(int_partition)
        child_idx = sorted_int_partition[0] - self.min_branch
        if self.n_children <= child_idx:
            return None
        if self.children[child_idx] is None:
            return None
        return self.children[child_idx].get_coef(sorted_int_partition[1:])

    def set_coef(self, int_partition, value):
        """
        Set the coefficient associated with a given integer partition.

        Parameters
        ----------
        int_partition : sequence of ints
            Integer partition to look up.
        value : object
            Value to associate with the integer partition.
        """
        if len(int_partition) == 0:
            self.value = value
            return
        sorted_int_partition = int_partition if self.assume_sorted else sorted(int_partition)
        child_idx = sorted_int_partition[0] - self.min_branch
        if self.n_children <= child_idx:
            self.children += [None] * (child_idx - self.n_children + 1)
            self.n_children = child_idx + 1
        if self.children[child_idx] is None:
            self.children[child_idx] = IntPartitionTree(sorted_int_partition[0], assume_sorted=True)
        self.children[child_idx].set_coef(sorted_int_partition[1:], value)

    def n_nodes(self):
        """
        Count the number of nodes in the tree.

        Returns
        -------
        n : int
            Number of nodes in the tree.
            Note that None placeholders for child nodes are not included in this count.
        """
        n_ancestors = 0
        for child in self.children:
            if child is not None:
                n_ancestors += child.n_nodes()
        return n_ancestors + 1

    def depth(self):
        """
        Calculate the depth of the tree.

        Returns
        -------
        d : int
            Depth of the tree, i.e., length of the longest path from the root node.
            Note that None placeholders for child nodes are not included in this count.
        """
        max_ancestor_depth = 0
        for child in self.children:
            if child is not None:
                max_ancestor_depth = max(max_ancestor_depth, child.depth())
        return max_ancestor_depth + 1
