from unittest import TestCase
import hgvs

__author__ = 'Alex H Wagner'


class TestVariant(TestCase):
    def setUp(self):
        self.p_sub = hgvs.Variant('FGFR3:p.R248C')

    def test_p_sub_to_c(self):
        assert self.p_sub.c.hgvs == 'ENST00000352904:c.742C>T'

    def test_p_sub_to_g(self):
        assert self.p_sub.g.hgvs
