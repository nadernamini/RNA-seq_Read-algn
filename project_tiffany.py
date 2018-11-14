""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
from shared import *

ALPHABET = [TERMINATOR] + BASES
RADIX = 100


def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
    
    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """
    sorted_suffixes = sort_suffixes(s, [(s[i:i+RADIX], i) for i in range(len(s))])
    sa = [pair[1] for pair in sorted_suffixes]
    return sa


def sort_suffixes(s, suffixes, ranks=None):
    if ranks is None:
        ranks = {}
        sorted_suffixes = sorted(suffixes, key=lambda x: x[0])
    else:
        sorted_suffixes = sorted(suffixes, key=lambda x: ranks[x])
    conflicts = {}
    bucket = []
    rank = 0
    curr_su = None
    for i in range(len(sorted_suffixes)):
        su, idx = sorted_suffixes[i]
        if su != curr_su:
            if len(bucket) > 1:
                conflicts[rank-1] = bucket
            bucket = []
            rank = i + 1
            curr_su = su
        ranks[(su, idx)] = rank
        bucket.append((su, idx))
    if len(bucket) > 1:
        conflicts[rank - 1] = bucket
    resolve_conflicts(s, sorted_suffixes, conflicts, ranks)
    return sorted_suffixes


def resolve_conflicts(s, sorted_suffixes, conflicts, ranks):
    for i, c in conflicts.items():
        resolved = [(s[idx-RADIX:idx], idx-RADIX) for su, idx in
                    sort_suffixes(s, [(s[idx+RADIX:idx+RADIX+RADIX], idx+RADIX) for su, idx in c], ranks)]
        for k in range(len(resolved)):
            r = resolved[k]
            ranks[r] += k
        sorted_suffixes[i:i+len(resolved)] = resolved


def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    L = ''
    for i in sa:
        L += s[i - 1]
    return L


def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted
    """
    return ''.join(sorted(L))


def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    """
    M = {c: -1 for c in ALPHABET}
    prev_c = None
    for i in range(len(F)):
        c = F[i]
        if c != prev_c:
            M[c] = i
            prev_c = c
    return M


def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    occ = {c: [0 for _ in range(len(L))] for c in ALPHABET}
    for i in range(len(L)):
        for char in ALPHABET:
            if char == L[i]:
                occ[char][i] = occ[char][i - 1] + 1
            else:
                occ[char][i] = occ[char][i - 1]
    return occ


def set_ep(c, M, occ):
    indices = sorted(M.values())
    try:
        return indices[indices.index(M[c]) + 1] - 1
    except IndexError:
        return len(occ[c]) - 1


def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    """
    length = len(p)
    sp = M[p[length - 1]]
    ep = set_ep(p[length - 1], M, occ)
    for i in range(length - 2, -1, -1):
        curr_range = (sp, ep + 1)
        sp = M[p[i]] + occ[p[i]][sp - 1]
        ep = M[p[i]] + occ[p[i]][ep] - 1
        if sp > ep:
            return curr_range, length - i - 1
    if sp > ep:
        return curr_range, length - i - 1
    return (sp, ep + 1), length


MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000


class Aligner:
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        s = genome_sequence + TERMINATOR
        self.sa = get_suffix_array(s)
        L = get_bwt(s, self.sa)
        F = get_F(L)
        self.M = get_M(F)
        self.occ = get_occ(L)
        self.genes = known_genes


    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        pass
