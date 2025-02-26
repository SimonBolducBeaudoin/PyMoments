#!/bin/env/python
#! -*- coding: utf-8 -*-


import unittest

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from Combinatorics import integer_partitions  # Now an absolute import


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

if __name__ == '__main__':
    unittest.main()
