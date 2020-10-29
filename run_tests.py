#!/usr/bin/env python

import unittest
import trf_utils

class TestTrfUtils(unittest.TestCase):

    def test_filter_overlaps(self):
        recs = {'contig1': [{'start':100,'end':140,'aln_score':10}, # same start diff end
                            {'start':100,'end':150,'aln_score':10},
                            {'start':155,'end':170,'aln_score':10}, # test sorting
                            {'start':120,'end':140,'aln_score':10}, # entirely contained
                            {'start':140,'end':160,'aln_score':40},
                            {'start':140,'end':160,'aln_score':70},
                            {'start':156,'end':158,'aln_score':10}, # entirely contained
                            {'start':180,'end':200,'aln_score':10}]}
        rec2 = {'contig1': [{'start':100,'end':150,'aln_score':10},
                            {'start':140,'end':160,'aln_score':70},
                            {'start':155,'end':170,'aln_score':10},
                            {'start':180,'end':200,'aln_score':10}]}
        self.assertEqual(
            trf_utils.filter_overlaps(recs),
            rec2)

if __name__ == "__main__":
    unittest.main()
