import re
import requests
import string
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.SeqUtils import seq1
from Bio.Seq import reverse_complement

__author__ = 'Alex H Wagner'


class HgvsBase:

    def __init__(self, ref_seq_id, start, stop, ref, alt, edit_type, predicted=False):
        self.ref_seq_id = ref_seq_id
        self.start = str(start)
        self.stop = str(stop)
        self.ref = str(ref)
        self.alt = str(alt)
        self.edit_type = edit_type
        self.predicted = predicted
        self._validate()

    ref_regex = re.compile(r'[ACTG]+$')
    alt_regex = ref_regex
    start_regex = re.compile(r'[\d+-]+$')
    stop_regex = start_regex

    def _validate(self):
        if self.edit_type not in ['substitution']:
            raise ValueError('Edit type {} not currently supported.'.format(self.edit_type))
        if not self.ref_seq_id.isalnum():
            raise ValueError('Expected only alphanumerics in ref_seq_id ({}).'.format(self.ref_seq_id))
        if not self.ref_regex.match(self.ref):
            raise ValueError('Unexpected ref sequence format.'.format(self.ref))
        if not self.start_regex.match(self.start):
            raise ValueError('Unexpected start position format.'.format(self.start))
        if not self.stop_regex.match(self.stop):
            raise ValueError('Unexpected stop position format.'.format(self.stop))
        if not self.alt_regex.match(self.alt):
            raise ValueError('Unexpected alt sequence format.'.format(self.alt))

    def __repr__(self):
        return self.hgvs

    @property
    def info(self):
        d = {'id': self.ref_seq_id,
             'start': self.start,
             'stop': self.stop,
             'ref': self.ref,
             'alt': self.alt,
             'edit_type': self.edit_type,
             'predicted': self.predicted}
        return d


class P(HgvsBase):

    ref_regex = re.compile('[{}]+$'.format(aa1))
    three_letter_regex = re.compile('({})+$'.format('|'.join(['({})'.format(x) for x in aa3])), re.I)
    alt_regex = ref_regex

    def __init__(self, ref_seq_id, start, stop, ref, alt, edit_type, predicted=False):

        if self.three_letter_regex.match(ref) and self.three_letter_regex.match(alt):
            ref = seq1(ref)
            alt = seq1(alt)
        super().__init__(ref_seq_id, start, stop, ref, alt, edit_type, predicted)


    @property
    def hgvs(self):
        if self.edit_type == 'substitution':
            return "{}:p.{}{}{}".format(self.ref_seq_id, self.ref, self.start, self.alt)

    ref_regex = re.compile(r'')

class C(HgvsBase):

    prefix = 'c'

    def __init__(self, ref_seq_id, start, stop, strand, ref, alt, edit_type, predicted=False):
        self.strand = str(strand)
        super().__init__(ref_seq_id, start, stop, ref, alt, edit_type, predicted=False)

    def _validate(self):
        super()._validate()
        if self.strand not in ('1', '-1'):
            raise ValueError("Expected strand info ({}) to be 1 or -1.".format(self.strand))
            # TODO: Add support for '+' / '-' as well


    @property
    def hgvs(self):
        if self.edit_type == 'substitution':
            return "{}:{}.{}{}>{}".format(self.ref_seq_id, self.prefix, self.start, self.ref, self.alt)

    @property
    def info(self):
        d = super().info
        d['strand'] = self.strand
        return d


class G(C):

    prefix = 'g'

    def __init__(self, chromosome, start, stop, ref, alt, edit_type):
        self.ref_seq_id = chromosome
        self.start = str(start)
        self.stop = str(stop)
        self.ref = str(ref)
        self.alt = str(alt)
        self.edit_type = edit_type
        self._validate()

    def _validate(self):
        HgvsBase._validate(self)
        chromosome_pattern = re.compile(r'(?:CHR)?([MTXY\d]+)')
        if chromosome_pattern.match(self.ref_seq_id):
            self.ref_seq_id = chromosome_pattern.match(self.ref_seq_id).group(1)
        else:
            raise ValueError("Expected chromosome ({}) to match regex /(chr)?[mtxy\d]+/i.".format(self.ref_seq_id))

    @property
    def info(self):
        d = {'start': self.start,
             'stop': self.stop,
             'ref': self.ref,
             'alt': self.alt,
             'chromosome': self.ref_seq_id,
             'edit_type': self.edit_type}
        return d

    @property
    def ucsc(self):
        return str('{}:{}-{}'.format('chr' + self.ref_seq_id, self.start, self.stop))

    @property
    def ensembl(self):
        return str('{}:{}-{}'.format(self.ref_seq_id, self.start, self.stop))

    # The chromosome property allows users to intuitively grab the expected chromosome string from G objects.
    # self.ref_seq_id is maintained to seamlessly inherit the hgvs property from the C class.

    @property
    def chromosome(self):
        return self.ref_seq_id

    @chromosome.setter
    def chromosome(self, value):
        self.ref_seq_id = value


class Variant:

    _primary_regex = re.compile(r'(?P<id>\S+):(?P<prefix>[pgc])\.(?P<edits>.*)')
    _p_sub_regex = re.compile(r'(?P<res1>[a-zA-Z]{1,3})(?P<pos>\d+)(?P<res2>[a-zA-Z*=]{1,3})')
    _p_ins_regex = re.compile(r'abracadabra')  # TODO: Change this to something sensible.
    _session = requests.Session()

    def __repr__(self):
        return str(self.hgvs)

    @property
    def ucsc(self):
        return self.g.ucsc

    @property
    def ensembl(self):
        return self.g.ensembl

    @property
    def hgvs(self):
        return self.g.hgvs, self.c.hgvs, self.p.hgvs

    @property
    def info(self):
        return self.g.info, self.c.info, self.p.info

    def __init__(self, hgvs=None, reference_assembly='current', species='human'):
        self.local = False  # TODO: Add support for local annotation files
        self.species = species
        self.p = None
        self.c = None
        self.g = None
        if hgvs is None:
            raise ValueError("Variants without hgvs nomenclature are currently not supported.")
        if reference_assembly == 'current':
            self._subdomain = ''
        elif isinstance(reference_assembly, int):
            self._subdomain = 'grch' + str(reference_assembly) + '.'
        elif reference_assembly.isnumeric():
            self._subdomain = 'grch' + reference_assembly + '.'
        elif reference_assembly.lower().startswith('grch'):
            self._subdomain = reference_assembly.lower() + '.'
        else:
            raise ValueError("reference_build should be of the format XX or grchXX (e.g. 'grch37')")
        m = Variant._primary_regex.match(hgvs)
        if not m:
            raise ValueError("hgvs string ({}) not recognized.".format(hgvs))
        self.edit = m.group('edits')
        if '[' in self.edit:
            raise ValueError("hgvs alleles are currently not supported.")
        prefix = m.group('prefix')
        if prefix == 'p':
            if Variant._p_sub_regex.match(self.edit):
                self.edit_type = 'substitution'
                p_edit_match = Variant._p_sub_regex.match(self.edit)
                start = p_edit_match.group('pos')
                stop = start
                ref = p_edit_match.group('res1')
                alt = p_edit_match.group('res2')
            elif Variant._p_ins_regex.match(self.edit):
                self.edit_type = 'insertion'
                p_edit_match = Variant._p_ins_regex.match(self.edit)
                # TODO: Finish filling out start, stop, ref, alt
            else:
                raise ValueError("This type of protein change ({}) is not currently supported.")
            self.p = P(m.group('id'),
                       start,
                       stop,
                       ref,
                       alt,
                       self.edit_type)
            self._p_fill()
        else:
            raise ValueError("Prefix type '{}' not supported.".format(prefix))

    def _vep_hgvs_rest(self, hgvs):
        url = 'http://{}rest.ensembl.org/vep/{}/hgvs/{}?content-type=application/json'\
            .format(self._subdomain, self.species, hgvs)
        resp = self._session.get(url)
        resp.raise_for_status()
        return resp.json()[0]

    def _p_fill(self):
        try:
            t = self._select_best_vep_transcript(self.p.hgvs)

            self.c = C(t['transcript_id'],
                       t['cds_start'],
                       t['cds_end'],
                       t['strand'],
                       t['codons'].split('/')[0].translate(str.maketrans('', '', string.ascii_lowercase)),
                       t['codons'].split('/')[1].translate(str.maketrans('', '', string.ascii_lowercase)),
                       self.edit_type)
        except ValueError as e:
            raise e  # TODO: Add in call to Ensembl overlap if no transcripts are returned by VEP query

        self._c_to_g()
        self.p.ref_seq_id = self._get_transcript_info(self.c.ref_seq_id)['Translation']['id']

    def _c_to_g(self):
        t = self._map_cds_to_genome()

        if self.c.strand == '-1':
            ref = reverse_complement(self.c.ref)
            alt = reverse_complement(self.c.alt)
        else:
            ref = self.c.ref
            alt = self.c.alt
        self.g = G(t['seq_region_name'],
                   t['start'],
                   t['end'],
                   ref,
                   alt,
                   self.edit_type)

    def _c_to_p(self):
        t_info = self._get_transcript_info(self.c.ref_seq_id)
        p_id = t_info['Translation']['id']
        self.p = None  # TODO: Go from Ensembl ID and coordinates to a protein change. Use BioPython and http://goo.gl/40sMGm

    def _get_transcript_info(self, transcript_id):
        """Return a dictionary conforming to the responses here: http://rest.ensembl.org/documentation/info/lookup"""
        if not self.local:
            try:
                url = 'http://{}rest.ensembl.org/lookup/id/{}?content-type=application/json;expand=1'\
                        .format(self._subdomain, transcript_id)
                resp = self._session.get(url)
                resp.raise_for_status()
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('Ensembl ID Lookup query failed: {}'.format(msg))
            return resp.json()
        else:
            raise ValueError('Local annotations not yet supported.')

    def _select_best_vep_transcript(self, hgvs):
        if not self.local:
            try:
                r = self._vep_hgvs_rest(hgvs)
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('VEP REST query failed: {}'.format(msg))
            transcripts = r['transcript_consequences']
        else:
            raise ValueError('Local annotations not yet supported.')
        result = list()
        backup = list()
        for i, transcript in enumerate(transcripts):
            try:
                if (transcript['biotype'] != 'protein_coding' or
                    str(transcript['protein_start']) != self.p.start or
                    str(transcript['protein_end']) != self.p.stop):
                    continue
                backup.append(i)
                if 'polyphen_score' not in transcript.keys():
                    continue
            except KeyError:
                continue
            result.append(i)
        if len(result) > 1:
            result = sorted(result, key=lambda x: transcripts[x]['transcript_id'], reverse=True)
            result = sorted(result, key=lambda x: transcripts[x]['polyphen_score'])
            best = result[-1]
        elif len(result) == 1:
            best = result[0]
        elif len(result) == 0 and len(backup) > 0:
            result = sorted(backup, key=lambda x: transcripts[x]['transcript_id'], reverse=True)
            best = result[-1]
        else:
            raise ValueError("No matching transcripts in VEP response!")
        return transcripts[best]

    def _map_cds_to_genome(self):
        if not self.local:
            try:
                url = 'http://{}rest.ensembl.org/map/cds/{}/{}..{}?content-type=application/json'\
                    .format(self._subdomain, self.c.ref_seq_id, self.c.start, self.c.stop)
                resp = self._session.get(url)
                resp.raise_for_status()
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('Ensembl CDS Map query failed: {}'.format(msg))
            result = resp.json()['mappings']
            if len(result) != 1:
                raise ValueError('Expected exactly one mapping for transcript {}, received {}.'
                                 .format(self.c.ref_seq_id, len(result)))
        else:
            raise ValueError('Local annotations not yet supported.')
        return result[0]


if __name__ == '__main__':
    v = Variant('ALK:p.F1174I', reference_assembly=37)
    print(v)
    print(v.g.ensembl)
