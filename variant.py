import re
import requests
import string
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.SeqUtils import seq1
from Bio.Seq import reverse_complement, Seq
from Bio.Alphabet import generic_protein

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
        if self.edit_type not in ['substitution', 'insertion']:
            raise ValueError('Edit type {} not currently supported.'.format(self.edit_type))
        if not self.ref_seq_id.isalnum():
            raise ValueError('Expected only alphanumerics in ref_seq_id ({}).'.format(self.ref_seq_id))
        if not self.ref_regex.fullmatch(self.ref):
            raise ValueError('Unexpected ref sequence format.'.format(self.ref))
        if not self.start_regex.fullmatch(self.start):
            raise ValueError('Unexpected start position format.'.format(self.start))
        if not self.stop_regex.fullmatch(self.stop):
            raise ValueError('Unexpected stop position format.'.format(self.stop))
        if not self.alt_regex.fullmatch(self.alt):
            raise ValueError('Unexpected alt sequence format.'.format(self.alt))

    def __repr__(self):
        return self.hgvs

    @property
    def start_pos(self):
        return self.start

    @property
    def stop_pos(self):
        return self.stop

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

    ref_regex = re.compile('([{}]+)|-'.format(aa1))
    three_letter_regex = re.compile('({})+$'.format('|'.join(['({})'.format(x) for x in aa3])), re.I)
    alt_regex = ref_regex
    start_regex = re.compile('[a-zA-Z]{1,3}\d+')
    stop_regex = re.compile('[a-zA-Z]{1,3}\d*')

    def __init__(self, ref_seq_id, start, stop, ref, alt, edit_type, predicted=False):

        if self.three_letter_regex.match(ref) and self.three_letter_regex.match(alt):
            ref = seq1(ref)
            alt = seq1(alt)
        super().__init__(ref_seq_id, start, stop, ref, alt, edit_type, predicted)

    @property
    def start_pos(self):
        return re.sub('[^0-9]', '', self.start)

    @property
    def stop_pos(self):
        return re.sub('[^0-9]', '', self.stop)

    @property
    def hgvs(self):
        if self.edit_type == 'substitution':
            return "{}:p.{}{}".format(self.ref_seq_id, self.start, self.alt)
        elif self.edit_type == 'insertion':
            return "{}:p.{}_{}ins{}".format(self.ref_seq_id, self.start, self.stop, self.alt)

class C(HgvsBase):

    prefix = 'c'

    def __init__(self, ref_seq_id, start, stop, strand, ref, alt, edit_type, predicted=False):
        self.strand = str(strand)
        if self.strand == '+':
            self.strand = '1'
        elif self.strand == '-':
            self.strand = '-1'
        super().__init__(ref_seq_id, start, stop, ref, alt, edit_type, predicted=False)

    def _validate(self):
        super()._validate()
        if self.strand not in ('1', '-1'):
            raise ValueError("Expected strand info ({}) to be 1 or -1.".format(self.strand))

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
    _p_ins_regex = \
        re.compile(r'(?P<res1>[a-zA-Z]{1,3}\d+)_(?P<res2>[a-zA-Z]{1,3}\d+)ins(?P<ins>\w+)')
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

    def __init__(self, hgvs=None, reference_assembly='current', species='human', **kwargs):
        self.local = False  # TODO: Add support for local annotation files
        self.species = species
        self.p = kwargs.get('P', None)
        self.c = kwargs.get('C', None)
        self.g = kwargs.get('G', None)
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
        if hgvs is None:
            raise ValueError("Variants without hgvs nomenclature are currently not supported.")
            # TODO: add checks for self.c, .g, .p, and fill in missing.
        else:
            self._parse_hgvs(hgvs)

    def _parse_hgvs(self, hgvs):
        m = self._primary_regex.match(hgvs)
        if not m:
            raise ValueError("hgvs string ({}) not recognized.".format(hgvs))
        self.edit = m.group('edits')
        if '[' in self.edit:
            raise ValueError("hgvs alleles are currently not supported.")
        prefix = m.group('prefix')
        if prefix == 'p':
            if Variant._p_sub_regex.fullmatch(self.edit):
                self.edit_type = 'substitution'
                p_edit_match = self._p_sub_regex.match(self.edit)
                ref = p_edit_match.group('res1')
                alt = p_edit_match.group('res2')
                start = ref + p_edit_match.group('pos')
                stop = start
            elif Variant._p_ins_regex.fullmatch(self.edit):
                self.edit_type = 'insertion'
                p_edit_match = Variant._p_ins_regex.match(self.edit)
                start = p_edit_match.group('res1')
                stop = p_edit_match.group('res2')
                ref = '-'
                alt = p_edit_match.group('ins')
            else:
                raise ValueError("This type of protein edit ({}) is not currently supported.".format(self.edit))
            self.p = P(m.group('id'),
                       start,
                       stop,
                       ref,
                       alt,
                       self.edit_type)
            self._p_fill()
        else:
            raise ValueError("Parsing of prefix type '{}' currently not supported.".format(prefix))

    def _vep_hgvs_rest(self, hgvs):
        url = 'http://{}rest.ensembl.org/vep/{}/hgvs/{}?content-type=application/json'\
            .format(self._subdomain, self.species, hgvs)
        resp = self._session.get(url)
        resp.raise_for_status()
        return resp.json()[0]

    def _p_fill(self):
        t = self._infer_best_transcript(self.p.hgvs)

        self.c = C(t['transcript_id'],
                   t['cds_start'],
                   t['cds_end'],
                   t['strand'],
                   t['ref'],
                   t['alt'],
                   self.edit_type)

        self._c_to_g()
        self.p.ref_seq_id = self._lookup_rest_by_id(self.c.ref_seq_id)['Translation']['id']

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
        t_info = self._lookup_rest_by_id(self.c.ref_seq_id)
        p_id = t_info['Translation']['id']
        self.p = None  # TODO: Go from Ensembl ID and coordinates to a protein change. Use BioPython and http://goo.gl/40sMGm

    def _lookup_rest_by_id(self, lookup_id, expand=True):
        """Return a dictionary conforming to the responses here: http://rest.ensembl.org/documentation/info/lookup"""
        if not self.local:
            if expand:
                expand = ';expand=1'
            else:
                expand = ''
            try:
                url = 'http://{}rest.ensembl.org/lookup/id/{}?content-type=application/json{}'\
                        .format(self._subdomain, lookup_id, expand)
                resp = self._session.get(url)
                resp.raise_for_status()
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('Ensembl ID Lookup query failed: {}'.format(msg))
            return resp.json()
        else:
            raise ValueError('Local annotations not yet supported.')

    def _lookup_rest_by_symbol(self, symbol):
        if not self.local:
            try:
                url = 'http://{}rest.ensembl.org/lookup/symbol/{}/{}?content-type=application/json;expand=1'\
                        .format(self._subdomain, self.species, symbol)
                resp = self._session.get(url)
                resp.raise_for_status()
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('Ensembl Symbol Lookup query failed: {}'.format(msg))
            return resp.json()
        else:
            raise ValueError('Local annotations not yet supported.')

    def _lookup_transcripts_rest(self, hgvs):
        lookup_id = hgvs.split(':')[0]
        try:
            r = self._lookup_rest_by_id(lookup_id)
        except ValueError as e:
            if str(e).count('Expand option only available for Genes and Transcripts'):
                r = self._lookup_rest_by_id(lookup_id, expand=False)
            else:
                r = self._lookup_rest_by_symbol(lookup_id)
        if r['object_type'] == 'Gene':
            return r['Transcript']
        elif r['object_type'] == 'Transcript':
            return [r]
        elif r['object_type'] == 'Translation':
            lookup_id = r['Parent']
            r = self._lookup_rest_by_id(lookup_id)
            return [r]
        else:
            raise ValueError('Unexpected response type ({}) from lookup (expected gene, transcript, or translation).'
                             .format(r['object_type']))

    def _infer_best_transcript(self, hgvs):
        """
        Infer the best transcript that corresponds to the mutation specified by the HGVS string.
        Return a dictionary including keys transcript_id, cds_start, cds_end, strand, ref, and alt.
        """
        if not self.local:
            try:
                r = self._vep_hgvs_rest(hgvs)
                transcripts = r['transcript_consequences']
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                try:
                    r = self._lookup_transcripts_rest(hgvs)
                    transcripts = self._select_p_compatible_transcripts(r)
                except requests.HTTPError as e:
                    msg2 = e.response.content.decode()
                    raise ValueError('VEP REST query and lookup failed:\n'
                                     '1){}\n2){}'.format(msg, msg2))
        else:
            raise ValueError('Local annotations not yet supported.')
        result = list()
        backup = list()
        for i, transcript in enumerate(transcripts):
            try:
                if (transcript['biotype'] != 'protein_coding' or
                    str(transcript['protein_start']) != self.p.start_pos or
                    str(transcript['protein_end']) != self.p.stop_pos):
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
        t = transcripts[best]
        if 'codons' in t:
            (t['ref'], t['alt']) = t['codons'].translate(str.maketrans('', '', string.ascii_lowercase)).split('/')
        else:
            # infer ref and alt
            raise ValueError('Not implemented!')
        return t

    def _select_p_compatible_transcripts(self, transcripts):
        out = list()
        for transcript in transcripts:
            if transcript['biotype'] != 'protein_coding':
                continue
            protein_id = transcript['Translation']['id']
            protein_seq = self._get_sequence(protein_id)

    def _get_sequence(self, seq_id):
        if not self.local:
            try:
                url = 'http://{}rest.ensembl.org/sequence/id/{}?content-type=application/json'\
                        .format(self._subdomain, seq_id)
                resp = self._session.get(url)
                resp.raise_for_status()
            except requests.HTTPError as e:
                msg = e.response.content.decode()
                raise ValueError('Ensembl sequence query failed: {}'.format(msg))
            return Seq(resp.json()['seq'], generic_protein)
        else:
            raise ValueError('Local annotations not yet supported.')

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
    v = Variant('ERBB2:p.P780_Y781insGSP', reference_assembly=37)
    print(v)
    print(v.g.ensembl)
