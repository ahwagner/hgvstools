import numpy as np


def ham_dist(seq1, seq2):
    """Returns the hamming distance between equal length sequences seq1 and seq2."""
    if len(seq1) != len(seq2):
        raise ValueError('Expected seq1 and seq2 to be same length!')
    dist = 0
    for char1, char2 in zip(seq1, seq2):
        if char1 != char2:
            dist += 1
    return dist


if __name__ == '__main__':
    ham_dist('aaa', 'aba')