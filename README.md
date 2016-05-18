# HGVS Tools <img src=https://travis-ci.org/ahwagner/hgvs_tools.svg?branch=master />

**HGVS Tools** enables rapid forward and reverse translation between protein, cDNA, and genomic HGVS variant
descriptions. In addition, it provides reasonable rules for inferring transcripts from ambiguous identifiers, automating the process of concordant, unambiguous transcript selection.

## Usage

### Construct a variant
You may construct a variant object directly from an HGVS string, using Chromosome names, HGNC Gene Symbols, or Ensembl
Transcript or Protein IDs.
```python
>>> from hgvs_tools import Variant
>>> v1 = Variant('FGFR3:p.R248C')
>>> v2 = Variant('9:g.22125504G>C')
>>> v3 = Variant('ENST00000003084:c.1431_1433delTTC')
```

### View variant HGVS strings
```python
>>> v1.g.hgvs
ENST00000352904:g.1803564C>T
>>> v1.c.hgvs
ENST00000352904:c.742C>T
>>> v1.p.hgvs
ENSP00000231803:p.R248C
```

### View variant genomic coordinates
```python
>>> v1.g.ucsc
chr4:1803564-1803564
>>> v1.g.ensembl
4:1803564-1803564
```

### Get a dictionary with info about a protein, cdna, or genomic description
```python
>>> v1.g.info
{'edit_type': 'substitution', 'strand': '1', 'ref': 'C', 'start': '1803564', 'chromosome': '4', 'alt': 'T',
 'id': 'ENST00000352904', 'stop': '1803564'}
>>> v1.c.info
{'alt': 'T', 'stop': '742', 'edit_type': 'substitution', 'id': 'ENST00000352904', 'ref': 'C', 'predicted': False, 
 'start': '742'}
>>> v1.p.info
{'predicted': False, 'start': '248', 'edit_type': 'substitution', 'alt': 'C', 'stop': '248', 'ref': 'R',
 'id': 'ENSP00000231803'}
```

## Example
```python
>>> from hgvs_tools import Variant
>>> v = Variant('FGFR3:p.R248C', reference_assembly=37)
>>> print(v)
('ENST00000352904:g.1803564C>T', 'ENST00000352904:c.742C>T', 'ENSP00000231803:p.R248C')
```
