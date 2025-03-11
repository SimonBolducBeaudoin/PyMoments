#!/bin/env/python
#! -*- coding: utf-8 -*-

import unittest

try : # Absolute import
    # Works with calling 
    # python -m unittest discover -s PyMoments/tests -p "test_*.py
    # From PyMoments's parent directorie
    from PyMoments.Combinatorics import simplex_iter,set_partitions,integer_partitions,mu_partitions,ff,binom
    from PyMoments.Combinatorics import fuse_sublists,list_disjoint_product,generate_fused_combinations,list_conjoint_product
except ModuleNotFoundError : # Relative import
    # Works with 
    # pytest PyMoments/tests
    from ..Combinatorics import simplex_iter,set_partitions,integer_partitions,mu_partitions,ff,binom
    from ..Combinatorics import fuse_sublists,list_disjoint_product,generate_fused_combinations,list_conjoint_product

class TestCombinatorics(unittest.TestCase):

    def test_simplex_iter(self):

        # Test empty case
        empty_case_1 = list(simplex_iter(1, [1, 1]))
        empty_case_2 = list(simplex_iter(0, [1]))
        self.assertEqual(len(empty_case_1), 0)
        self.assertEqual(len(empty_case_2), 0)

        # Test one-element cases
        singleton_case_1 = list(simplex_iter(1, [1]))
        singleton_case_2 = list(simplex_iter(5, [1, 2, 3, 4, 5]))
        self.assertListEqual(singleton_case_1, [(1,)])
        self.assertListEqual(singleton_case_2, [(1, 1, 1, 1, 1)])

        # Test some simple cases
        simple_case_1 = set(simplex_iter(4, [1, 2, 3]))
        simple_case_2 = set(simplex_iter(5, [1, 2, 3]))
        simple_case_1_true = {(1, 1, 2), (1, 2, 1)}
        simple_case_2_true = {(1, 1, 3), (1, 2, 2)}
        self.assertSetEqual(simple_case_1, simple_case_1_true)
        self.assertSetEqual(simple_case_2, simple_case_2_true)

    def test_set_partitions(self):

        # Test partitions of an empty set
        s0_parts_est = list(set_partitions([]))
        self.assertEqual(len(s0_parts_est), 0)

        # Test partitions of a 1-element set
        s1 = ['dirigible']
        s1_parts_true = {
            frozenset([('dirigible',)])}
        s1_parts_est = set(map(frozenset, set_partitions(s1)))
        self.assertSetEqual(s1_parts_est, s1_parts_true)

        # Test partitions of a 2-element set
        s2 = [1, -1]
        s2_parts_true = {
            frozenset([(1,), (-1,)]),
            frozenset([(1, -1)])}
        s2_parts_est = set(map(frozenset, set_partitions(s2)))
        self.assertSetEqual(s2_parts_est, s2_parts_true)

        # Test partitions of a 3-element set
        s3 = ['apple', 'banana', 1.4]
        s3_parts_true = {
            frozenset([('apple',), ('banana',), (1.4,)]),
            frozenset([('apple', 'banana'), (1.4,)]),
            frozenset([('apple', 1.4), ('banana',)]),
            frozenset([('apple',), ('banana', 1.4)]),
            frozenset([('apple', 'banana', 1.4)])}
        s3_parts_est = set(map(frozenset, set_partitions(s3)))
        self.assertSetEqual(s3_parts_est, s3_parts_true)

        # Test partitions of a 4-element set
        s4 = [1, 2, 3, 4]
        s4_parts_true = {
            frozenset([(1,), (2,), (3,), (4,)]),
            frozenset([(1, 2), (3,), (4,)]),
            frozenset([(1, 3), (2,), (4,)]),
            frozenset([(1, 4), (2,), (3,)]),
            frozenset([(2, 3), (1,), (4,)]),
            frozenset([(2, 4), (1,), (3,)]),
            frozenset([(3, 4), (1,), (2,)]),
            frozenset([(1, 2), (3, 4)]),
            frozenset([(1, 3), (2, 4)]),
            frozenset([(1, 4), (2, 3)]),
            frozenset([(1, 2, 3), (4,)]),
            frozenset([(1, 2, 4), (3,)]),
            frozenset([(1, 3, 4), (2,)]),
            frozenset([(2, 3, 4), (1,)]),
            frozenset([(1, 2, 3, 4)])}
        s4_parts_est = set(map(frozenset, set_partitions(s4)))
        self.assertSetEqual(s4_parts_est, s4_parts_true)

        # Validate sizes of larger sets
        bell_numbers = [52, 203, 877, 4140, 21147, 115975]
        for n in range(5, 11):
            count = 0
            for _ in set_partitions(list(range(n))):
                count += 1
            self.assertEqual(count, bell_numbers[n-5])

    def test_ff(self):

        self.assertEqual(ff(0, 0), 1)

        for n in range(1, 11):
            self.assertEqual(ff(n, 0), 1)
            ff_true = n
            for i in range(1, n + 1):
                self.assertEqual(ff(n, i), ff_true)
                ff_true *= (n - i)

    def test_binom(self):

        # First 5 rows of Pascal's triangle
        self.assertEqual(binom(0, 0), 1)
        self.assertEqual(binom(1, 0), 1)
        self.assertEqual(binom(1, 1), 1)
        self.assertEqual(binom(2, 0), 1)
        self.assertEqual(binom(2, 1), 2)
        self.assertEqual(binom(2, 2), 1)
        self.assertEqual(binom(3, 0), 1)
        self.assertEqual(binom(3, 1), 3)
        self.assertEqual(binom(3, 2), 3)
        self.assertEqual(binom(3, 3), 1)
        self.assertEqual(binom(4, 0), 1)
        self.assertEqual(binom(4, 1), 4)
        self.assertEqual(binom(4, 2), 6)
        self.assertEqual(binom(4, 3), 4)
        self.assertEqual(binom(4, 4), 1)

        # Selected larger values
        self.assertEqual(binom(5, 2), 10)
        self.assertEqual(binom(7, 2), 21)
        self.assertEqual(binom(8, 2), 28)
        self.assertEqual(binom(8, 4), 70)
        self.assertEqual(binom(20, 14), 38760)

class TestIntegerPartitions(unittest.TestCase):
    def test_zero(self):
        """The only partition of 0 should be the empty partition."""
        partitions = list(integer_partitions(0))
        self.assertEqual(partitions, [[]])

    def test_one(self):
        """Partitioning 1 should produce only [[1]]."""
        partitions = list(integer_partitions(1))
        self.assertEqual(partitions, [[1]])

    def test_five_default(self):
        """Test partitions of 5 with default min_value=1."""
        partitions = list(integer_partitions(5))
        expected = [
            [5],
            [4, 1],
            [3, 2],
            [3, 1, 1],
            [2, 2, 1],
            [2, 1, 1, 1],
            [1, 1, 1, 1, 1]
        ]
        # Order is not important so we compare sorted lists.
        self.assertEqual(sorted(partitions), sorted(expected))

    def test_eleven_default(self):
        """Test partitions of 11 with default min_value=1."""
        partitions = list(integer_partitions(11))
        expected_count = 56  # There are 56 partitions of 11
        self.assertEqual(len(partitions), expected_count)

    def test_five_min_value_2(self):
        """Test partitions of 5 with parts not smaller than 2."""
        partitions = list(integer_partitions(5, min_value=2))
        expected = [
            [5],
            [3, 2]
        ]
        self.assertEqual(sorted(partitions), sorted(expected))
    
    def test_five_max_value_3(self):
        """Test partitions of 5 with parts not larger than 3."""
        partitions = list(integer_partitions(5, max_value=3))
        expected = [
            [3, 2],
            [3, 1, 1],
            [2, 2, 1],
            [2, 1, 1, 1],
            [1, 1, 1, 1, 1]
        ]
        self.assertEqual(sorted(partitions), sorted(expected))
        
    def test_eleven_max_value_5(self):
        """Test partitions of 11 with parts not larger than 5."""
        partitions = list(integer_partitions(11, max_value=5,min_value=2))
        expected = [
            [5,4,2],[5,3,3],[5,2,2,2],[4,4,3],[4,3,2,2],[3,3,3,2],[3,2,2,2,2]
        ]
        self.assertEqual(sorted(partitions), sorted(expected))

    def test_five_min_2_max_3(self):
        """Test partitions of 5 with parts between 2 and 3 (inclusive)."""
        partitions = list(integer_partitions(5, min_value=2, max_value=3))
        expected = [
            [3, 2]
        ]
        self.assertEqual(sorted(partitions), sorted(expected))
            
    def test_non_increasing_order(self):
        """Ensure that each partition is in non-increasing order."""
        partitions = list(integer_partitions(5))
        for partition in partitions:
            with self.subTest(partition=partition):
                self.assertTrue(all(partition[i] >= partition[i+1] for i in range(len(partition)-1)),
                                f"{partition} is not in non-increasing order")

class TestMuPartitions(unittest.TestCase):
    def test_empty_set(self):
        """mu_partitions should return nothing for an empty set."""
        partitions = list(mu_partitions([]))
        self.assertEqual(partitions, [])

    def test_single_element(self):
        """mu_partitions should return nothing for a single-element set."""
        partitions = list(mu_partitions(["A"]))
        self.assertEqual(partitions, [])

    def test_two_elements(self):
        """mu_partitions should return a single partition: both elements together."""
        partitions = list(mu_partitions(["A", "B"]))
        expected = [[["A", "B"]]]  # Example expected output (assuming binary representation)
        self.assertEqual(partitions, expected)

    def test_three_elements(self):
        """mu_partitions should remove partitions that contain singletons."""
        partitions = list(mu_partitions(["A", "B", "C"]))
        expected=[[['A', 'B', 'C']]]
        self.assertEqual(partitions, expected)

    def test_four_elements(self):
        """mu_partitions should generate correct partitions for four elements."""
        partitions = list(mu_partitions(["A", "B", "C", "D"]))
        # Check that none of the partitions contain singleton blocks
        expected = [[['A', 'B', 'C', 'D']],
        [['A', 'B'], ['C', 'D']],
        [['A', 'C'], ['B', 'D']],
        [['A', 'D'], ['B', 'C']]]
        self.assertEqual(partitions, expected)

class TestFuseSublists(unittest.TestCase):
    def test_basic_fusion(self):
        self.assertEqual(fuse_sublists(('a', 'b'), ('c', 'd')), ('a', 'b', 'c', 'd'))

    def test_fusion_with_duplicates(self):
        self.assertEqual(fuse_sublists(('a', 'b', 'a'), ('b', 'a', 'c')), ('a', 'a', 'a', 'b', 'b', 'c'))

    def test_fusion_with_empty_sublist(self):
        self.assertEqual(fuse_sublists((), ('x', 'y', 'z')), ('x', 'y', 'z'))
        self.assertEqual(fuse_sublists(('x', 'y', 'z'), ()), ('x', 'y', 'z'))

class TestGenerateFusedCombinations(unittest.TestCase):
    def test_basic_fusion(self):
        sublists_A = [('a',)]
        sublists_B = [('b',)]
        expected_output = [(('a', 'b'),)]
        self.assertEqual(generate_fused_combinations(sublists_A, sublists_B), expected_output)

    def test_different_lengths(self):
        sublists_A = [('a',), ('b',)]
        sublists_B = [('x',)]
        expected_output = [(('a', 'x'), ('b',)), (('b', 'x'), ('a',))]
        self.assertEqual(generate_fused_combinations(sublists_A, sublists_B), expected_output)

    def test_permutation_check(self):
        sublists_A = [('a',), ('b',)]
        sublists_B = [('x',), ('y',)]
        expected_output = [(('a', 'x'), ('b',), ('y',)),
            (('a', 'y'), ('b',), ('x',)),
            (('b', 'x'), ('a',), ('y',)),
            (('b', 'y'), ('a',), ('x',)),
            (('a', 'x'), ('b', 'y')),
            (('a', 'y'), ('b', 'x'))]
        result = generate_fused_combinations(sublists_A, sublists_B)
        self.assertCountEqual(result, expected_output)  # Order doesn't matter

    def test_empty_sublists(self):
        sublists_A = []
        sublists_B = [('x',)]
        expected_output = []
        self.assertEqual(generate_fused_combinations(sublists_A, sublists_B), expected_output)

        sublists_A = [('a',)]
        sublists_B = []
        expected_output = []
        self.assertEqual(generate_fused_combinations(sublists_A, sublists_B), expected_output)

    def test_duplicate_elements(self):
        sublists_A = [('a', 'a'), ('b',)]
        sublists_B = [('x', 'x'), ('y',)]
        expected_output = [(('a', 'a', 'x', 'x'), ('b',), ('y',)),
             (('a', 'a', 'y'), ('b',), ('x', 'x')),
             (('b', 'x', 'x'), ('a', 'a'), ('y',)),
             (('b', 'y'), ('a', 'a'), ('x', 'x')),
             (('a', 'a', 'x', 'x'), ('b', 'y')),
             (('a', 'a', 'y'), ('b', 'x', 'x'))]
        result = generate_fused_combinations(sublists_A, sublists_B)
        self.assertCountEqual(result, expected_output)

class TestListConjointProduct(unittest.TestCase):
    def test_basic_case(self):
        A = [(1, ('a',))]
        B = [(2, ('b',))]
        expected_output = [(2, ('a', 'b'))]  # Count is multiplied, sublists are fused
        self.assertEqual(list_conjoint_product(A, B), expected_output)

    def test_different_lengths(self):
        A = [(1, ('a',)), (3, ('c',))]
        B = [(2, ('b',))]
        expected_output = [
            (2, ('a', 'b')),
            (6, ('b', 'c'))
        ]
        self.assertEqual(list_conjoint_product(A, B), expected_output)

    def test_multiplication_of_counts(self):
        A = [(2, ('x',)), (3, ('y',))]
        B = [(4, ('z',))]
        expected_output = [
            (8, ('x', 'z')),
            (12, ('y', 'z'))
        ]
        self.assertEqual(list_conjoint_product(A, B), expected_output)

    def test_empty_input(self):
        A = []
        B = [(2, ('b',))]
        self.assertEqual(list_conjoint_product(A, B), [])

        A = [(1, ('a',))]
        B = []
        self.assertEqual(list_conjoint_product(A, B), [])

        A = []
        B = []
        self.assertEqual(list_conjoint_product(A, B), [])

    def test_onesublist_fusion(self):
        A = [(1, ('a', 'b'))]
        B = [(2, ('x', 'y'))]
        expected_output = [ (2, ('a', 'b', 'x', 'y')) ]
        result = list_conjoint_product(A, B)
        self.assertCountEqual(result, expected_output)  
        
    def test_complex_fusion(self):
        A = [(1, ('a', 'b')),(2,('c',),('d',))]
        B = [(3, ('e', 'f')),(4,('g',),('h',))]
        expected_output = [(3, ('a', 'b', 'e', 'f')),
             (4, ('a', 'b', 'g'), ('h',)),
             (4, ('a', 'b', 'h'), ('g',)),
             (6, ('c', 'e', 'f'), ('d',)),
             (6, ('d', 'e', 'f'), ('c',)),
             (8, ('c', 'g'), ('d',), ('h',)),
             (8, ('c', 'h'), ('d',), ('g',)),
             (8, ('d', 'g'), ('c',), ('h',)),
             (8, ('d', 'h'), ('c',), ('g',)),
             (8, ('c', 'g'), ('d', 'h')),
             (8, ('c', 'h'), ('d', 'g'))]
        result = list_conjoint_product(A, B)
        self.assertCountEqual(result, expected_output)  
        
if __name__ == '__main__':
    unittest.main()