#!/bin/env/python
#! -*- coding: utf-8 -*-

import unittest

try : # Absolute import
    # Works with calling 
    # python -m unittest discover -s PyMoments/tests -p "test_*.py
    # From PyMoments's parent directorie
    from PyMoments.Conversions import cumulants_to_moment,cumulants_to_mu,moments_to_cumulant,mus_to_cumulant
    
except ModuleNotFoundError : # Relative import
    # Works with 
    # pytest PyMoments/tests
    from ..Conversions import cumulants_to_moment,cumulants_to_mu,moments_to_cumulant,mus_to_cumulant


class TestCumulantsToMoment(unittest.TestCase):
    
    def test_default_behavior(self):
        # Test that the default behavior returns a list of terms
        
        # Example input: (0, 1)
        terms = cumulants_to_moment((0, 1))
        expected_terms = [
            [1, (0,), (1,)],  # Second term: corresponds to <X_0> <X_1>
            [1, (0, 1)]       # First term: corresponds to <X_0 X_1>
        ]
        
        # Check if the output matches the expected terms
        self.assertEqual(terms, expected_terms)
        
    def test_iterator_behavior(self):
        # Test iterator behavior with as_iterator=True
        
        # Example input: (0, 1)
        it = cumulants_to_moment((0, 1), as_iterator=True)
        
        # Expected output: the iterator should yield terms one by one
        first_term = next(it)
        second_term = next(it)
        
        self.assertEqual(first_term , [1, (0,), (1,)])
        self.assertEqual(second_term, [1, (0, 1)]    )
        
        # Ensure that the iterator raises StopIteration after yielding all terms
        with self.assertRaises(StopIteration):
            next(it)

    def test_empty_input(self):
        # Test behavior with empty input, should handle gracefully
        
        # Empty input (should return an empty list of terms)
        terms = cumulants_to_moment(())
        self.assertEqual(terms, [])
        
        # Test iterator with empty input (should not yield any terms)
        it = cumulants_to_moment((), as_iterator=True)
        with self.assertRaises(StopIteration):
            next(it)

    def test_multiple_terms(self):
        # Test input with more complex structure
        
        # Example input: (0, 1, 2)
        terms = cumulants_to_moment((0, 1, 2))
        expected_terms = [[1, (0,), (1,), (2,)],
            [1, (0, 1), (2,)],
            [1, (0, 2), (1,)],
            [1, (0,), (1, 2)],
            [1, (0, 1, 2)]]
        
        self.assertEqual(terms, expected_terms)
        
    def test_all_1d_cumulants_to_moment_on_wiki(self):
        """
        This one uses multi_index notation
        """
        terms = cumulants_to_moment((1,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(1,)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((2,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(2,)]], [1, [(1,), (1,)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((3,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(3,)]], [3, [(1,), (2,)]], [1, [(1,), (1,), (1,)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((4,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(4,)]],
            [4, [(1,), (3,)]],
            [3, [(2,), (2,)]],
            [6, [(1,), (1,), (2,)]],
            [1, [(1,), (1,), (1,), (1,)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((5,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(5,)]],
            [5, [(1,), (4,)]],
            [10, [(2,), (3,)]],
            [10, [(1,), (1,), (3,)]],
            [15, [(1,), (2,), (2,)]],
            [10, [(1,), (1,), (1,), (2,)]],
            [1, [(1,), (1,), (1,), (1,), (1,)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((6,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(6,)]],
            [6, [(1,), (5,)]],
            [15, [(2,), (4,)]],
            [15, [(1,), (1,), (4,)]],
            [10, [(3,), (3,)]],
            [60, [(1,), (2,), (3,)]],
            [20, [(1,), (1,), (1,), (3,)]],
            [15, [(2,), (2,), (2,)]],
            [45, [(1,), (1,), (2,), (2,)]],
            [15, [(1,), (1,), (1,), (1,), (2,)]],
            [1, [(1,), (1,), (1,), (1,), (1,), (1,)]]]
        self.assertEqual(terms, expected_terms)
        
        
    def test_bivariate_cumulants_to_moment_in_ref1(self):
        """
        ref1 Bi-variate k-statistics and Cumulants of Their Joint Sampling Distribution
        M.B. Cook 
        Biometrika, Jun. 1951, Vol. 38, No. 1/2 (Jun., 1951), pp. 179-195
        """
        terms = cumulants_to_moment((4,0,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(4, 0)]],
            [4, [(1, 0), (3, 0)]],
            [3, [(2, 0), (2, 0)]],
            [6, [(1, 0), (1, 0), (2, 0)]],
            [1, [(1, 0), (1, 0), (1, 0), (1, 0)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((2,2,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(2, 0), (0, 2)]],
            [1, [(2, 0), (0, 1), (0, 1)]],
            [1, [(1, 0), (1, 0), (0, 2)]],
            [1, [(1, 0), (1, 0), (0, 1), (0, 1)]],
            [1, [(2, 2)]],
            [2, [(2, 1), (0, 1)]],
            [2, [(1, 0), (1, 2)]],
            [4, [(1, 0), (1, 1), (0, 1)]],
            [2, [(1, 1), (1, 1)]]]
        self.assertEqual(terms, expected_terms)
        
        terms = cumulants_to_moment((3,3,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(3, 0), (0, 3)]],
            [3, [(3, 0), (0, 1), (0, 2)]],
            [1, [(3, 0), (0, 1), (0, 1), (0, 1)]],
            [3, [(1, 0), (2, 0), (0, 3)]],
            [9, [(1, 0), (2, 0), (0, 1), (0, 2)]],
            [3, [(1, 0), (2, 0), (0, 1), (0, 1), (0, 1)]],
            [1, [(1, 0), (1, 0), (1, 0), (0, 3)]],
            [3, [(1, 0), (1, 0), (1, 0), (0, 1), (0, 2)]],
            [1, [(1, 0), (1, 0), (1, 0), (0, 1), (0, 1), (0, 1)]],
            [1, [(3, 3)]],
            [3, [(3, 1), (0, 2)]],
            [3, [(3, 2), (0, 1)]],
            [3, [(3, 1), (0, 1), (0, 1)]],
            [3, [(2, 0), (1, 3)]],
            [3, [(1, 0), (2, 3)]],
            [9, [(2, 0), (1, 1), (0, 2)]],
            [9, [(2, 0), (1, 2), (0, 1)]],
            [9, [(1, 0), (2, 1), (0, 2)]],
            [9, [(1, 0), (2, 2), (0, 1)]],
            [9, [(2, 2), (1, 1)]],
            [9, [(2, 1), (1, 2)]],
            [9, [(2, 0), (1, 1), (0, 1), (0, 1)]],
            [9, [(1, 0), (2, 1), (0, 1), (0, 1)]],
            [18, [(2, 1), (1, 1), (0, 1)]],
            [3, [(1, 0), (1, 0), (1, 3)]],
            [9, [(1, 0), (1, 0), (1, 1), (0, 2)]],
            [9, [(1, 0), (1, 0), (1, 2), (0, 1)]],
            [18, [(1, 0), (1, 1), (1, 2)]],
            [9, [(1, 0), (1, 0), (1, 1), (0, 1), (0, 1)]],
            [18, [(1, 0), (1, 1), (1, 1), (0, 1)]],
            [6, [(1, 1), (1, 1), (1, 1)]]]
        self.assertEqual(terms, expected_terms)
        
    def test_bivariate_cumulants_to_mu_in_ref1(self):
        """
        ref1 Bi-variate k-statistics and Cumulants of Their Joint Sampling Distribution
        M.B. Cook 
        Biometrika, Jun. 1951, Vol. 38, No. 1/2 (Jun., 1951), pp. 179-195
        """    
        terms = cumulants_to_mu((5,1,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(5, 1)]],
            [5, [(4, 0), (1, 1)]],
            [10, [(3, 0), (2, 1)]],
            [10, [(2, 0), (3, 1)]],
            [15, [(2, 0), (2, 0), (1, 1)]]]
        self.assertEqual(terms, expected_terms)
        
    def test_bivariate_mus_to_cumulant_in_ref1(self):
        """
        ref1 Bi-variate k-statistics and Cumulants of Their Joint Sampling Distribution
        M.B. Cook 
        Biometrika, Jun. 1951, Vol. 38, No. 1/2 (Jun., 1951), pp. 179-195
        """    
        terms = mus_to_cumulant((5,1,),multi_index_in=True,multi_index_out=True)
        expected_terms = [[1, [(5, 1)]],
            [-5, [(4, 0), (1, 1)]],
            [-10, [(3, 0), (2, 1)]],
            [-10, [(2, 0), (3, 1)]],
            [30, [(2, 0), (2, 0), (1, 1)]]]
        self.assertEqual(terms, expected_terms)
        
if __name__ == "__main__":
    unittest.main()
