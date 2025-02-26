#!/bin/env/python
#! -*- coding: utf-8 -*-

from Combinatorics import set_partitions, mu_partitions
from DataStructures import partition_of_multi_indices,multi_idxs_to_multiset
from math import factorial

class _conversion_base:
    def __new__(cls, *args, as_iterator=False, **kwargs):
        """Decides whether to return an iterator or a list."""
        instance = super().__new__(cls)  # Create instance normally
        if as_iterator:
            return instance  # Return the iterator
        else :
            instance.__init__(*args, **kwargs)  # Explicitly call __init__
            return list(instance)  # Collect all elements into a list
    def __init__(self,multiset,LaTex=False,multi_index_in=False,multi_index_out=False,**kwargs):
        self.multiset_in = multiset if not(multi_index_in) else multi_idxs_to_multiset(multiset)
        self._LaTex=LaTex
        self._multi_index=multi_index_out
    def __iter__(self):
        return self
    def _next_term_frmt(self):
        if self._multi_index :
            nxt = self.next_term()
            coef = nxt[0]
            mltset = partition_of_multi_indices(nxt[1:])
            return [coef,mltset]
        else :
            return self.next_term()
    def next_term_latex(self,var="x"):
        nxt = self._next_term_frmt()
        return _conversion_base.term_to_latex(nxt,var=var)
    def next_term(self):
        raise StopIteration
    def __next__(self):
        if self._LaTex :
            return self.next_term_latex()
        else :
            return self._next_term_frmt()
    @staticmethod
    def term_to_latex(list_of_prod,var="x"):
        coef = list_of_prod[0]
        if coef == 1 :
            s = ""
        elif coef == -1 :
            s = "-"
        else :
            s = f"{coef}"
        s += " "
        for tup in list_of_prod[1:] :
            subscripts = ','.join([f'{i}' for i in tup])
            s += f"{var}_{{{subscripts}}}"
            s += " "
        return s

class _to_moment(_conversion_base):
    def __init__(self, *args,**kwargs):
        super().__init__(*args,**kwargs)
        self.gen = set_partitions(self.multiset_in)
        
class _to_centered_moment(_conversion_base):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.gen = mu_partitions(self.multiset_in)
                          
class cumulants_to_moment(_to_moment):
    """
    This class implements Equation (1.3) from K.D. Smith:

    .. math::
        m_{[i_1, i_2, \dots, i_k]} = \sum_{\pi \in \Pi_k} \prod_{B \in \pi} \kappa_{[i_j | j \in B]}

    By default, calling this class returns a full list of terms representing the moment.
    If `as_iterator=True` is passed, it behaves as an iterator and yields terms one by one.

    Parameters
    ----------
    multiset : list of tuples
        A multiset (a set that allows repetitions) representing the moment
        to be computed from cumulants. Each tuple corresponds to a random
        variable's index. Examples:

        <X_0^2>=<X_0 X_0> = is (0,0)
        <X_1^2>=<X_1 X_1> = is (1,1)
                <X_0 X_1> = is (0,1)
                <X_0 X_1 X_0> = is (0,1,0) or (0,0,1) is equivalent
    as_iterator : bool, optional
        If True, returns an iterator instead of a full list. Default is False.

    Returns
    -------
    list or iterator
        If `as_iterator=False` (default), returns a list where each element represents a term
        of the expansion. Each term is a list containing a coefficient (always 1) followed by tuples
        representing cumulant terms.
        If `as_iterator=True`, returns an iterator that yields terms one by one.

    Examples
    --------
    The second cumulant is the covariance:

    .. math::
        <X_i X_j> = \kappa_{i,j} + \kappa_{i} \kappa_{j}

    Example usage:

    >>> i, j = 0, 1
    >>> terms = cumulants_to_moment((0,1))
    >>> print(terms)
    [[1, (0,1)], [1, (0,), (1,)]]

    >>> it = cumulants_to_moment((0,1), as_iterator=True)
    >>> next(it)
    [1, (0,1)]
    >>> next(it)
    [1, (0,), (1,)]
    """
    def next_term(self):
        return [1] + next(self.gen)
    def next_term_latex(self):
        return super().next_term_latex(var="\kappa")
               
class cumulants_to_mu(_to_centered_moment):
    """
    Same as cumulants_to_moment but for centered moments
    """
    def next_term(self):    
        return [1] + next(self.gen)
    def next_term_latex(self):
        return super().next_term_latex(var="\kappa")
                
class moments_to_cumulant(_to_moment):
    """
    This class implements Equation (1.6) from K.D. Smith:

    .. math::
        \kappa_{[i_1, i_2, \dots, i_k]} = \sum_{\pi \in \Pi_k} (-1)^{|\pi| - 1} (|\pi| - 1)! \prod_{B \in \pi} m_{[i_j \mid j \in B]}

    By default, calling this class returns a full list of terms representing the cumulant.
    If `as_iterator=True` is passed, it behaves as an iterator and yields terms one by one.

    Parameters
    ----------
    multiset : list of tuples
        A multiset (a set that allows repetitions) representing the cumulant
        to be computed from moments. Each tuple corresponds to a random
        variable's index. Examples:

        <k_{2,0}> = is (0,0)
        <k_{0,2}> = is (1,1)
        <k_{1,1}> = is (0,1)
        <k_{2,1}> = is (0,1,0) or (0,0,1) is equivalent
    as_iterator : bool, optional
        If True, returns an iterator instead of a full list. Default is False.

    Returns
    -------
    list or iterator
        If `as_iterator=False` (default), returns a list where each element represents a term
        of the expansion. Each term is a list containing a coefficient followed by tuples
        representing moment terms.
        If `as_iterator=True`, returns an iterator that yields terms one by one.

    Examples
    --------
    The second cumulant is the covariance:

    .. math::
        \kappa_{i,j} = <X_i X_j> - <X_i><X_j>

    Example usage:

    >>> i, j = 0, 1
    >>> terms = moments_to_cumulant((0,1))
    >>> print(terms)
    [[1, (0,1)], [-1, (0,), (1,)]]

    >>> it = moments_to_cumulant((0,1), as_iterator=True)
    >>> next(it)
    [1, (0,1)]
    >>> next(it)
    [-1, (0,), (1,)]
    """
    
    def next_term(self):
        nxt = next(self.gen)
        pi = len(nxt)
        coef = (-1)**(pi-1)*factorial(pi-1)
        return [coef] + nxt
    def next_term_latex(self):
        return super().next_term_latex(var="m")
        
class mus_to_cumulant(_to_centered_moment):
    """
    Same as moments_to_cumulant but for centered moments
    """
    def next_term(self):
        nxt = next(self.gen)
        pi = len(nxt)
        coef = (-1)**(pi-1)*factorial(pi-1)
        return [coef] + nxt
    def next_term_latex(self):
        return super().next_term_latex(var="\mu")