from unittest import TestCase
import seqtools


class TestHam_dist(TestCase):

    def test_correct_distance(self):
        assert seqtools.ham_dist('abra', 'acba') == 2
