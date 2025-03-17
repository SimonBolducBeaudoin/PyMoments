def growth_string_to_partition(growth_string,set):
    ret = []
    for i in range(max(growth_string)+1):
        ret += [[]]
    for idx,elem in zip(growth_string,set) :
        ret[idx] += [elem]
    return ret

def restricted_combinations(iterable, r):
    """
    Similar to itertools.combinations, but always returns the first element 
    first and then combines the rest.
    retricted_combinations('ABCD', 2)   → AB AC AD
    retricted_combinations(range(4), 3) → 012 013 023
    
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
    r -= 1
    pool = tuple(iterable[1:])
    
    for combo in combinations(pool, r):
        yield (first,) + combo

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
            comb = restricted_combinations( self.max_vals[w], size)
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
            if len(self.set)==0 :
                return []
            self._init_from_(0)
            self._first_call_ = False
            return self.GS
        else :
            for i,(it,idx) in enumerate(zip(self.it[::-1],self.GS_idx[::-1])) :
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