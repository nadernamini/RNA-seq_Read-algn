"""
    RNA Alignment Assignment

    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.

    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys  # DO NOT EDIT THIS
from shared import *
import numpy as np
from time import time

ALPHABET = [TERMINATOR] + BASES


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
    el, k = len(s), 1000
    return radix_sorting(s, [(i * 100, (i + 1) * 100) for i in range(k // 100)])


def lex_order(bs):
    sa = []
    for lst in bs:
        sa.extend(lst)
    return sa


def radix_sorting(ss, kk):
    src = ss
    buckets = [list(range(len(src)))]
    buckets_sort = []

    for s_i, e_i in kk:
        # print(buckets)
        len_check, t = True, -time()
        for idx in range(len(buckets)):
            new_buckets = {}
            for iddx in range(len(buckets[idx])):
                indx = buckets[idx][iddx]
                ky = src[indx + s_i: indx + e_i]
                if ky in new_buckets:
                    new_buckets[ky].append(buckets[idx][iddx])
                else:
                    new_buckets[ky] = [buckets[idx][iddx]]

            for k in sorted(new_buckets):
                if len_check and len(new_buckets[k]) > 1:
                    len_check = False
                buckets_sort.append(new_buckets[k])

        buckets, buckets_sort = buckets_sort, []
        # print(e_i, time() + t, len(buckets), len_check)
        if len_check:
            return lex_order(buckets)


def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    L, l = [], len(s)
    for i in range(l):
        L.append(s[(sa[i] - 1) % l])
    return ''.join(L)


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
    d, last_ch = {c: -1 for c in ALPHABET}, -1
    for i in range(len(F)):
        ch = F[i]
        if ch != last_ch:
            d[ch] = i
            last_ch = ch
    return d


def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    d = {c: [0 for _ in range(len(L))] for c in ALPHABET}
    for let in ALPHABET:
        for i in range(len(L)):
            if let == L[i]:
                d[let][i] = d[let][i - 1] + 1
            else:
                d[let][i] = d[let][i - 1]
    return d


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
    rng, length = None, 0

    sp, ep = 0, -1

    for i in range(len(p) - 1, -1, -1):
        sp, ep = M[p[i]] + occ[p[i]][sp - 1] if i != (len(p) - 1) else M[p[i]], M[p[i]] + occ[p[i]][ep] - 1
        if sp > ep:
            break
        else:
            rng, length = (sp, ep + 1), len(p) - i

    return rng, length


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
        self.sa = get_suffix_array(genome_sequence + TERMINATOR)
        L = get_bwt(genome_sequence, self.sa)
        self.occ, self.M = get_occ(L), get_M(get_F(L))
        self.known_genes = known_genes
        self.known_transcriptome = {}
        for gene in known_genes:
            self.known_transcriptome[gene] = {}
            for iso in known_genes[gene].isoforms:
                isoform, exons, total_len = '', [], 0
                for ex in iso.exons:
                    isoform += genome_sequence[ex.start: ex.end]
                    exons.append((total_len, ex.start, ex.end - ex.start))
                    total_len += ex.end - ex.start
                self.known_transcriptome[gene][iso] = [isoform, exons]

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces.
        Each piece consists of a start index in the read, a start index in the genome, and a length
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the Project.pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        pass

    def align_to_transcriptome(self, read_sequence):
        match = (-1, float("inf"))
        for gene in self.known_transcriptome:
            for iso in self.known_transcriptome[gene]:
                alignment, mismatches = Aligner.align_to_isoform(read_sequence, self.known_transcriptome[gene][iso][0],
                                                                 self.known_transcriptome[gene][iso][1])
                if alignment != -1 and mismatches <= match[1]:
                    match = (alignment, mismatches)
        return match

    @staticmethod
    def align_to_isoform(read, isoform, exons):
        match = (-1, float("inf"))
        for i in range(len(isoform) - len(read) + 1):
            j, mismatches = 0, 0
            while j < len(read):
                if isoform[i + j] != read[j]:
                    mismatches += 1
                if mismatches > MAX_NUM_MISMATCHES:
                    break
                j += 1
            if j == len(read) and mismatches < match[1]:
                match = (i, mismatches)
        return match if match[0] == -1 else (Aligner.find_alignment(len(read), match[0], exons), match[1])

    @staticmethod
    def find_alignment(read_len, align_start, exons):
        # maximum number of un-gapped alignments
        k = 3

        def find_start_location(lo, hi):
            mid = (lo + hi) // 2
            if exons[mid][0] <= align_start < exons[mid][0] + exons[mid][2]:
                return mid
            elif exons[mid][0] + exons[mid][2] <= align_start:
                return find_start_location(mid + 1, hi)
            else:
                return find_start_location(lo, mid - 1)

        idx = find_start_location(0, len(exons))
        align = [(0, exons[idx][1] + (align_start - exons[idx][0]),
                  read_len if read_len + (align_start - exons[idx][0]) <= exons[idx][2]
                  else exons[idx][2] - (align_start - exons[idx][0]))]
        read_len -= exons[idx][2] - (align_start - exons[idx][0])
        while read_len > 0:
            idx += 1
            align.append((align[-1][0] + align[-1][2], exons[idx][1],
                          read_len if read_len <= exons[idx][2]
                          else exons[idx][2]))
            read_len -= exons[idx][2]
        return align if len(align) <= k else -1
