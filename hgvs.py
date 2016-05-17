import re
import requests
from Bio.Data import CodonTable
from Bio.SeqUtils import seq1, seq3

__author__ = 'Alex H Wagner'


class HgvsBase:

    def __init__(self, ref_seq_id, start, stop, ref, alt, mut_type, predicted=False):
        self.ref_seq_id = ref_seq_id
        self.start = str(start)
        self.stop = str(stop)
        self.ref = str(ref)
        self.alt = str(alt)
        self.mut_type = mut_type
        self.predicted = predicted
        if not self._validate():
            raise ValueError

    def _validate(self):
        if self.mut_type == 'substitution':
            return all((self.ref_seq_id.isalnum(),
                        self.ref.isalpha(),
                        self.start.isnumeric(),
                        self.stop.isnumeric(),
                        self.alt.isalpha()))
        else:
            return False

    def __repr__(self):
        return self.hgvs

    @property
    def hgvs(self):
        raise NotImplementedError


class P(HgvsBase):

    @property
    def hgvs(self):
        if self.mut_type == 'substitution':
            return "{}:p.{}{}{}".format(self.ref_seq_id, self.ref, self.start, self.alt)

class C(HgvsBase):

    prefix = 'c'

    @property
    def hgvs(self):
        if self.mut_type == 'substitution':
            return "{}:{}.{}{}>{}".format(self.ref_seq_id, self.prefix, self.start, self.ref, self.alt)


class G(C):

    prefix = 'g'

    def __init__(self, ref_seq_id, chromosome, start, stop, strand, ref, alt, mut_type):
        self.ref_seq_id = ref_seq_id
        self.start = str(start)
        self.stop = str(stop)
        self.ref = str(ref)
        self.alt = str(alt)
        self.chromosome = chromosome
        self.strand = str(strand)
        self.mut_type = mut_type
        self._validate()


class Variant:

    regex = re.compile(r'(?P<id>\S+):(?P<prefix>[pgc])\.(?P<edits>.*)')
    p_sub_regex = re.compile(r'(?P<res1>[a-zA-Z]{1,3})(?P<pos>\d+)(?P<res2>[a-zA-Z*=]{1,3})')
    p_ins_regex = re.compile(r'abracadabra')  # TODO: Change this to something sensible.

    def __init__(self, hgvs=None, reference_assembly='current', species='human'):
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
        m = Variant.regex.match(hgvs)
        if not m:
            raise ValueError("hgvs string ({}) not recognized.".format(hgvs))
        self.edit = m.group('edits')
        if '[' in self.edit:
            raise ValueError("hgvs alleles are currently not supported.")
        prefix = m.group('prefix')
        if prefix == 'p':
            if Variant.p_sub_regex.match(self.edit):
                self.edit_type = 'substitution'
                p_edit_match = Variant.p_sub_regex.match(self.edit)
                start = p_edit_match.group('pos')
                stop = start
                ref = p_edit_match.group('res1')
                alt = p_edit_match.group('res2')
            elif Variant.p_ins_regex.match(self.edit):
                self.edit_type = 'insertion'
                p_edit_match = Variant.p_ins_regex.match(self.edit)
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
        resp = requests.get(url)
        resp.raise_for_status()
        return resp.json()[0]

    def _p_fill(self):
        try:
            r = self._vep_hgvs_rest(self.p.hgvs)
        except requests.HTTPError as e:
            msg = e.response.content.decode()
            raise ValueError('VEP REST lookup failed: {}'.format(msg))
        candidate_transcripts = r['transcript_consequences']
        try:
            transcript = self._select_best_vep_transcript(candidate_transcripts)

            self.c = C(transcript['transcript_id'],
                       transcript['cds_start'],
                       transcript['cds_end'],
                       r['allele_string'].split('/')[0],
                       r['allele_string'].split('/')[1],
                       self.edit_type)

            self.g = G(transcript['transcript_id'],
                       r['start'],
                       transcript['cds_end'])
        except ValueError as e:
            raise e  # TODO: Add in call to Ensembl overlap if no transcripts are returned by VEP query

    def _select_best_vep_transcript(self, transcripts):
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


if __name__ == '__main__':
    v = Variant('FGFR3:p.R248C', reference_assembly=37)
    print(v.c)
