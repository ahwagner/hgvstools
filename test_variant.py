from unittest import TestCase
import hgvs

__author__ = 'Alex H Wagner'


class TestVariant(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p_sub = hgvs.Variant('FGFR3:p.R248C', reference_assembly=37)
        cls.p_sub_neg = hgvs.Variant('ALK:p.F1174I', reference_assembly=37)

    def test_p_sub_to_c(self):
        assert self.p_sub.c.hgvs == 'ENST00000352904:c.742C>T'

    def test_p_sub_rep(self):
        assert self.p_sub.hgvs == ('4:g.1803564C>T',
                                   'ENST00000352904:c.742C>T',
                                   'ENSP00000231803:p.R248C')

    def test_p_sub_to_g(self):
        assert self.p_sub.g.hgvs == '4:g.1803564C>T'

    def test_sub_ucsc(self):
        assert self.p_sub.g.ucsc == 'chr4:1803564-1803564'

    def test_sub_ensembl(self):
        assert self.p_sub_neg.g.ensembl == '2:29443697-29443697'

    def test_negative_strand_rep(self):
        assert self.p_sub_neg.hgvs == ('2:g.29443697A>T',
                                       'ENST00000389048:c.3520T>A',
                                       'ENSP00000373700:p.F1174I')
