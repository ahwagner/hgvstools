from collections import defaultdict
from Bio.Data import CodonTable

reverseCodonTable = defaultdict(list)
for k, v in CodonTable.unambiguous_dna_by_name['Standard'].forward_table.items():
    reverseCodonTable[v].append(k)