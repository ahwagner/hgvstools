# HGVS Tools

The HGVS Tools enable rapid forward and reverse translation between protein, cDNA, and genomic HGVS variant
descriptions.

## Example Usage

```python
>from hgvs_tools import Variant
>v = Variant('FGFR3:p.R248C', reference_assembly=37)
>print(v)

('ENST00000352904:g.1803564C>T', 'ENST00000352904:c.742C>T', 'ENSP00000231803:p.R248C')
```
