"""
Combinatorics.py
Module of combinatorics-related generators and methods.
"""

from math import factorial
from more_itertools import distinct_permutations
from numpy import r_,zeros,full

from itertools import combinations,product,permutations,groupby
from collections import defaultdict
from math import perm,comb,prod

def integer_partitions(n, min_value=1, max_value=None):
    """
    Generate all partitions of n where each part is at least min_value.
    
    Parameters:
    - n: the integer to partition.
    - min_value: the minimum allowed value for any part (default is 1).
    - max_value: the maximum allowed value for the next part (default is n).
    
    The partitions are yielded in non-increasing order.
    """ 
    if max_value is None:
        max_value = n
    if n == 0:
        yield []
        return
    # Only consider values i that are between min_value and min(n, max_value)
    for i in range(min(n, max_value), min_value - 1, -1):
        for partition in integer_partitions(n - i, min_value, i):
            yield [i] + partition

def simplex_iter(s, max_vals):
    """
    Generator over all tuples of integers (i1, i2, ..., id) with the following properties:
        * i1 + i2 + ... + id = s
        * i1, i2, ..., id >= 1
        * ij <= max_vals[j]

    Parameters
    ----------
    s : int
        Sum to which the indices are constrained.
    max_vals : sequence of positive ints
        Maximum value that each mode of the index can take on.

    Yields
    ------
    idx : tuple of ints
        Index tuple satisfying the three properties above.
    """
    if s <= 0:
        return
    elif len(max_vals) == 0:
        return
    elif len(max_vals) == 1:
        if max_vals[0] >= s:
            yield (s,)
        else:
            return
    else:
        for i in range(1, 1 + min(max_vals[0], s)):
            for indices in simplex_iter(s - i, max_vals[1:]):
                yield (i,) + indices

def _set_partitions_slow(set):
    """
    Generator over all partitions of the given set.

    Parameters
    ----------
    set : sequence of objects
        Set to partition. Note that repeated elements are treated as distinct.

    Yields
    ------
    pi : list of tuples of elements from set
        Partition of the set.
    """
    if len(set) == 0:
        return
    elif len(set) == 1:
        yield [(set[0],)]
    else:
        head = set[0]
        for tail_parts in _set_partitions_slow(set[1:]):
            for i, part in enumerate(tail_parts):
                new_part = (head,) + part
                yield tail_parts[:i] + [new_part] + tail_parts[i + 1:]
            yield [(head,)] + tail_parts
                    
def disjoint_product(A, B):
    result = []
    for (count_A, *partition_A), (count_B, *partition_B) in product(A, B):
        new_count = count_A * count_B  # Multiply the leading counts
        new_sublists = tuple(partition_A + partition_B)  # Concatenate the sublists
        result.append((new_count, *new_sublists))
    return result

def limited_combinations(limits, k, start=0, current_counts=None):
    """
    A generator that yields valid combinations of length k using indexed elements,
    where each index has a limited number of repetitions.
    
    The output is a tuple where each value represents how many times the corresponding element is taken.
    """
    if current_counts is None:
        current_counts = [0] * len(limits)

    if sum(current_counts) == k:
        yield tuple(current_counts)
        return

    for i in range(start, len(limits)):
        if current_counts[i] < limits[i]:
            new_counts = current_counts[:]
            new_counts[i] += 1
            yield from limited_combinations(limits, k, i, new_counts)
            
def multi_index_to_combination(elements, multi_index):
    """
    Converts a multi-index into the actual combination.
    
    elements: List of items to choose from.
    multi_index: Tuple representing the multiplicity of each element.
    
    Returns a list containing the actual combination.
    """
    combination = []
    for count, element in zip(multi_index, elements):
        combination.extend([element] * count)
    return combination
            
def limited_permutations(limits, k, current_perm=None, remaining_counts=None):
    """
    A generator that yields valid permutations of length k using indexed elements,
    where each index has a limited number of repetitions.

    The output is a tuple representing a permutation, where each value corresponds
    to an index being chosen.
    """
    if current_perm is None:
        current_perm = []
    if remaining_counts is None:
        remaining_counts = list(limits)

    if len(current_perm) == k:
        yield tuple(current_perm)
        return

    for i in range(len(limits)):
        if remaining_counts[i] > 0:
            remaining_counts[i] -= 1
            yield from limited_permutations(limits, k, current_perm + [i], remaining_counts)
            remaining_counts[i] += 1  # Backtrack
            
def index_to_permutation(elements,perm_index):
  return [elements[i] for i in perm_index] 

def remove_duplicates_and_count(set):
    """
    Count the repeating elements of a set or multiset and converts it to (number, (element))
    """
    count_dict = defaultdict(int)
    for sublist in set:
        key = tuple(sorted(sublist))
        count_dict[key] += 1
            # Construct the result where count is the first element
    return [tuple([count] + list(sublist)) for sublist, count in count_dict.items()]
      
def reduce_partitions(partition):
    """
    Counts the identical block and reduces the representation.
    """
    counts = []
    red_part = []
    for key, group in groupby(partition):
        counts   += ( len(list(group)) ,)
        red_part += ( tuple(key) ,)
    return counts,red_part
    
def sum_duplicates(input_list):
    """
    Reduces the (number, *elements ) representation by summing over identical elements
    """
    result_dict = defaultdict(int)
    for num, *tuples in input_list:
        
        sorted_tuples = tuple(sorted(tuples))
        result_dict[sorted_tuples] += num
    result_list = [(num, *key) for key, num in result_dict.items()]
    
    return result_list
    
def count_occurrences(data: tuple, max_value: int) -> tuple:
    """Counts occurrences of each value from 0 to max_value-1 in the given tuple."""
    counts = [0] * max_value  # Initialize a list of zeros
    for num in data:
        if 0 <= num < max_value:  # Ensure values are within range
            counts[num] += 1
    return tuple(counts)
    
def generate_fused_combinations(partition_A, partition_B):
    """
    Generates all valid unions of elements from partition_A and partition_B, using indices and permutations for partition_B.
    """
    min_len = min(len(partition_A), len(partition_B))
    fused_results = []
    
    red_count_A,partition_red_A = reduce_partitions(partition_A)
    red_count_B,partition_red_B = reduce_partitions(partition_B)
    
    # Try fusing 1 to min_len elements together using indices
    for i in range(1, min_len + 1):
        for idx_comb_A in limited_combinations(red_count_A, i): 
            comb_A      = multi_index_to_combination( partition_red_A,idx_comb_A )
            rem_idx_A = [i-j for i,j in zip(red_count_A,idx_comb_A)]
            remaining_A = multi_index_to_combination( partition_red_A, rem_idx_A )
            multiplicity_A = prod([comb(n,k) for n,k in zip(red_count_A,idx_comb_A)])
            for idx_perm_B in limited_permutations(red_count_B, i):
                perm_B = index_to_permutation(partition_red_B,idx_perm_B)
                idx_comb_B = count_occurrences(idx_perm_B,len(partition_red_B))
                rem_idx_B = [i-j for i,j in zip(red_count_B,idx_comb_B)]
                remaining_B = multi_index_to_combination( partition_red_B, rem_idx_B )
                multiplicity_B = prod([perm(n,k) for n,k in zip(red_count_B,idx_comb_B)])
                # Fuse the elements corresponding to the indices
                fused = [tuple(sorted(a+b)) for a, b in zip(comb_A, perm_B)]    
                # Construct new sublist arrangement
                new_sublists = tuple(fused + remaining_A + remaining_B)
                fused_results.append((multiplicity_A*multiplicity_B,*new_sublists))
    return fused_results
    
def conjoint_product(A, B):
    result = []
    for (count_A, *partition_A), (count_B, *partition_B) in product(A, B):
        new_count = count_A * count_B  # Multiply the leading counts
        # Generate all fused combinations
        fused_combinations = generate_fused_combinations(partition_A, partition_B)
        # Store each valid combination in the result list
        for sublists in fused_combinations:
            multiplicity = sublists[0]
            result.append((new_count*multiplicity, *sublists[1:]))
    return result
    
def partitions_composition(A, B):
    """
    Composes/Multiplies the set partitions A and B together by first taking the union of the outputs of disjoint_product(A, B) and conjoint_product(A, B).
    """
    return disjoint_product(A, B) + conjoint_product(A, B)
    
def set_partitions(set):
    if len(set) == 0:
        return []
    
    set = sorted(set)  # Ensure order
    
    # Group identical elements together
    groups = [list(group) for _, group in groupby(set)]
    
    # Initialize result with partitions of the first group
    current_multiset = remove_duplicates_and_count(_set_partitions_slow(groups[0]))

    for i in range(1, len(groups)):
        next_multiset = remove_duplicates_and_count(_set_partitions_slow(groups[i]))
        current_multiset = partitions_composition(current_multiset, next_multiset)

        current_multiset = sum_duplicates(current_multiset)
    return iter(current_multiset)
        
def mu_partitions(set):
    if len(set) == 0:
        return
    else :
        gen = set_partitions(set)  # Assuming set_partitions is a function that returns a generator
        while True:
            nxt = next(gen, None)  # Get the next element or None if generator is exhausted
            if nxt is None:
                break  # End the loop if generator is exhausted
            
            # Check if one of the elements has length 1
            if all(len(idx_tpl) != 1 for idx_tpl in nxt[1:]):
                yield nxt  # Yield the valid tuple
            
def ff(n, i):
    """
    Returns the falling factorial (n)_i.

    Parameters
    ----------
    n : int
        Argument to the falling factorial.
    i : int
        Number of terms to include.

    Returns
    -------
    ff : int
        Falling factorial (n)_i.

    Notes
    -----
    The falling factorial is computes the product n(n-1)...(n-1+1).
    For example, (4)_2 = 4(3) = 12.
    """
    if i <= 0:
        return 1
    elif i == 1:
        return n
    else:
        return n * ff(n-1, i-1)
