# HGVS Tools

The HGVS Tools enable rapid forward and reverse translation between protein, cDNA, and genomic HGVS variant
descriptions. HGVS Tools provides reasonable rules for inferring transcripts from ambiguous identifiers, simplifying the
process of 

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

## Example
```python
>>> from hgvs_tools import Variant
>>> v = Variant('FGFR3:p.R248C', reference_assembly=37)
>>> print(v)
('ENST00000352904:g.1803564C>T', 'ENST00000352904:c.742C>T', 'ENSP00000231803:p.R248C')
```
