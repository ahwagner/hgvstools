from collections import defaultdict
from Bio.Data import CodonTable

reverseCodonTable = defaultdict(list)
for k, v in CodonTable.unambiguous_dna_by_name['Standard'].forward_table.items():
    reverseCodonTable[v].append(k)

# Human codon usage frequencies from: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
codonUsageTable = None