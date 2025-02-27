"""
Combinatorics.py
Module of combinatorics-related generators and methods.
"""

from math import factorial

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


def set_partitions(set):
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
        for tail_parts in set_partitions(set[1:]):
            for i, part in enumerate(tail_parts):
                new_part = (head,) + part
                yield tail_parts[:i] + [new_part] + tail_parts[i + 1:]
            yield [(head,)] + tail_parts
            
def mu_partitions(set):
    """
    Returns partitions for centered moments 
    A.k.a. removes all set partition with a block of size 1.
    
    It works by :
    1. Generating an iterator if the integer partitions.
    2. Getting distinct permutations of each integer partition.
    3. Producing growth strings based on the permutations.
    4. Converting the growth strings to the corresponding elements of the set.

    Yields:
        - Each element of the partitions for centered moments .
    """
    if not set :
        return 
    for int_partition in integer_partitions(len(set), 2):
        for block_shape in distinct_permutations(int_partition):
           for gs in growth_string_from_blocks_shape(block_shape, set):
               yield growth_string_to_partition(gs,set)
            
def retricted_combinations(iterable, r):
    """
    Very similar to itertools.combinations.
    Produces the combinations placing the next digit in a restricted growth string
    Always returns the first element first and then combines the rest
    retricted_combinations('ABCD', 2)   → AB AC AD
    retricted_combinations(range(4), 3) → 012 013 023
    Use in a restricted growth string
    
    Say the growths string's state is currently
        0010x10xx 
        and say that next digits to be place are two 2s.
        There are 3 spots available to be placed in.
        The first x must always be occupied by a 2 for the restricted growth string 
        to respect its growth inequality. Therefore the options are
        00102102x
        0010210x2
        This generator is given the indices for the possible position of the two 2s 
        will return the right combinations of possible placements.
        retricted_combinations([5,8,9],2) -> (5,8) and (5,9).
    """
    if r > len(iterable):
        return
    
    first = iterable[0]
    
    r-=1
    pool = tuple(iterable[1:])
    n = len(pool)
    
    indices = list(range(r))

    yield (first,) + tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield (first,) +tuple(pool[i] for i in indices)

class growth_string_from_blocks_shape:
    """
    Generates all valid growth string corresponding to a certain blocks shape (aka a tuple of blocks sizes).
    """
    def __init__(self, blocks_shape, set):
        self.GS_vals = r_[:len(blocks_shape)]       # Ex: (0,1,2) 
        self.blocks_shape = blocks_shape            # Ex: (2,1,1)
        self.set = set                              # Ex: (A,B,C,D)
        self.max_vals = r_[:len(set)]               # Ex: (0,1,2,3)
        self.max_val = len(set)-1
        
        self.GS    = zeros( (len(set),), int )      # Ex: (0,0,0,0)
        
        self.available = full(  (len(set),), True ) # Ex: (T,T,T,T)
        
        # Initializing a list for each number to be placed in the growth string
        self.it = []
        self.GS_idx = []
        
        self._first_call_= True
        
    def _init_from_(self,start=0) :
        for val,size in zip(self.GS_vals[start:-1],self.blocks_shape[start:-1]) :
            w    = self.available  
            comb = retricted_combinations( self.max_vals[w], size)
            nxt  = next(comb)
            self.GS_idx += (nxt,) 
            self.GS[ [*nxt] ] = val
            self.available[ [*nxt] ] = False
            self.it += [comb]
        
        # The last indices can always be deduced.
        self.GS[self.available] = self.GS_vals[-1]

    def __iter__(self):
        return self
        
    def __next__(self):
        if self._first_call_ :
            self._init_from_(0)
            self._first_call_ = False
            return self.GS
        else :
            for i,(it,idx) in enumerate(zip(self.it[::-1],self.GS_idx[::-1])) :
                #pdb.set_trace()
                i += 1
                self.available[[*idx]] = True
                del self.GS_idx[-1]
                try :                 
                    nxt = next(it)
                    self.GS_idx += (nxt,) 
                    val = self.GS_vals[-i-1]
                    self.GS[ [*nxt] ] = val
                    self.available[ [*nxt] ] = False
                    
                    start = len(self.GS_vals)-i
                    self._init_from_(start)
                    
                    return self.GS
                except StopIteration : 
                    del self.it[-1] 
                    continue
            raise StopIteration

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


def binom(n, k):
    """
    Computes the binomial coefficient n choose k.
    I.e., computs the number of ways to choose k-element subsets from a collection with n elements.

    Parameters
    ----------
    n : int
        Number of elements in the collection.
    k : int
        Size of subsets of the collection to choose.

    Returns
    -------
    b : int
        Binomial coefficient n choose k.

    Notes
    -----
    The binomial coefficient is given by n! / k! (n-k)!
    Thus, the computation is simplified using falling factorials:
        (n, k) = n! / k! (n-k)! = (n)_k / (k)_(k-1)
    """
    return ff(n, k) // ff(k, k-1)
