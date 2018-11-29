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
        s = genome_sequence[::-1] + TERMINATOR
        self.genome_length = len(genome_sequence)
        self.sa = get_suffix_array(s)
        L = get_bwt(s, self.sa)
        self.occ, self.M = get_occ(L), get_M(get_F(L))
        self.known_genes = known_genes
        self.known_transcriptome = {}
        for gene in known_genes:
            self.known_transcriptome[gene] = {}
            for iso in gene.isoforms:
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
        alignment, mismatches = self.align_to_transcriptome(read_sequence)
        if alignment != -1:
            return alignment
        else:
            aligned, alignment = self.align_to_genome(read_sequence)
            if aligned:
                return alignment
            else:
                return []

    def align_to_transcriptome(self, read_sequence):
        """
        Performs an alignment to all the known isoforms (i.e. the annotated transcriptome) and returns the one with the
        smallest number of mismatches.

        :param read_sequence: RNA read sequence
        :return: a tuple (alignment, mismatches)
            alignment: alignment to genome if one exists else -1
            mismatches: number of mismatches in the alignment if one exists otherwise float("inf")
        """
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
        """
        Performs alignment to the isoform element-by-element as long as it does not exceed MAX_NUM_MISMATCHES.

        :param read: RNA read sequence
        :param isoform: isoform
        :param exons: exons in the isoform
        :return: a tuple (alignment, mismatches)
            alignment: alignment to genome if one exists else -1
            mismatches: number of mismatches in the alignment if one exists otherwise float("inf")
        """
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
        return match if match[0] == -1 else (Aligner.find_transcriptome_alignment(len(read), match[0], exons), match[1])

    @staticmethod
    def find_transcriptome_alignment(read_len, align_start, exons):
        """
        Finds the alignment to the genome sequence of length read_len starting at align_start in the isoform created by
        the exons.

        :param read_len: length of the RNA read
        :param align_start: start of the alignment in the isoform created by the exons
        :param exons: the exons
        :return: alignment in the format mentioned in Aligner.align if less than or equal to number of alignment parts,
                 else -1
        """
        # maximum number of un-gapped alignments
        k = 3

        def find_start_location(lo, hi):
            """
            Recursive binary search function that returns the index of the exon in which align_start lies.

            :param lo: lower index of the search
            :param hi: higher index of the search
            :return: the exon in which the start of the alignment lies
            """
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

    def exact_match(self, string):
        """
        Returns whether or not an exact match of the string is found in the reversed genome sequence and if found, the
        range in the suffix array indicating where they start in the reversed genome.

        :param string: string to check
        :return: a tuple (match, range)
            match: a boolean indicating whether an exact match of the string exists in the reversed genome sequence
            range: None if not match, otherwise the range returned by exact_suffix_matches
        """
        _range, length = exact_suffix_matches(string, self.M, self.occ)
        if length == len(string):
            return True, _range
        else:
            return False, None

    @staticmethod
    def replace_base_at_index(read_sequence, idx, base):
        """
        Returns a copy of read_sequence with the idx-th element replaced by base.

        :param read_sequence: reversed RNA read sequence
        :param idx: the index of the replacement
        :param base: the base character to replace the idx-element of read_sequence
        :return: read_sequence with the idx-th element replaced with base
        """
        return read_sequence[:idx] + base + read_sequence[idx + 1:]

    def check_seeds(self, _range, length, high_bound=None):
        """
        Returns a valid set of seeds, i.e. starting locations in the original genome sequence.

        :param _range: range of the seed locations in the suffix array, self.sa
        :param length: length of the match for the genome suffixes starting at the index locations in the range in the
                       suffix array
        :param high_bound: bound which every index location in _range + length has to be less than
        :return: list of valid seeds
        """
        valid_seeds = []
        for i in range(_range[0], _range[1]):
            location = self.sa[i] + length
            if not high_bound or location < high_bound and MIN_INTRON_SIZE <= high_bound - location <= MAX_INTRON_SIZE:
                valid_seeds.append(self.genome_length - location)
        return valid_seeds

    def find_max_so_far(self, max_so_far, len_read, mismatches):
        """
        A recursive function that finds (if it exists) a **single** ungapped alignment of the reversed RNA read to the
        reversed genome with the mismatches constraint. It does so by considering all the permutations of the base at
        the index of maximum alignment (second element) in max_so_far and recursively calling the function on the new
        max_so_far if there exists a match in the genome (by calling exact_match(string)). The recursive call stack
        terminates if the number of mismatches (fourth element) of max_so_far exceeds number of allowed mismatches
        (no alignment) or the length (second element) of max_so_far becomes equal to len_read (alignment!).

        :param max_so_far: the current best alignment of the reversed read sequence and the reversed genome
        :param len_read: length of the original RNA read
        :param mismatches: maximum number of mismatches allowed
        :return: a tuple (aligned, mx_so_far)
            aligned: a boolean indicating whether or not a **single** ungapped alignment was found
            mx_so_far: best max_so_far found when the terminating conditions are satisfied as documented above
        """
        if max_so_far[3] < mismatches and max_so_far[1] < len_read:
            misses = max_so_far[3]
            changes = {}
            matched, new_range = self.exact_match(max_so_far[2][:max_so_far[1] + 1])
            if matched:
                changes[''] = (new_range, max_so_far[1] + 1, max_so_far[2], misses)
            possible_changes = [b for b in BASES if b != max_so_far[2][max_so_far[1]]]
            for change in possible_changes:
                new_length = max_so_far[1] + 1
                new_read = self.replace_base_at_index(max_so_far[2], max_so_far[1], change)
                matched, new_range = self.exact_match(new_read[:new_length])
                if matched:
                    changes[change] = (new_range, new_length, new_read, misses + 1)
            for change in changes:
                test = self.find_max_so_far(changes[change], len_read, mismatches)
                if test:
                    aligned, maxi = test
                    if aligned:
                        return aligned, maxi
            return False, max_so_far
        elif max_so_far[1] == len_read and max_so_far[3] <= mismatches:
            return True, max_so_far
        else:
            return False, max_so_far

    @staticmethod
    def replace_bases_at_indices(read_sequence, _range, bases):
        """
        Returns a copy of read_sequence with the range elements replaced by bases.

        :param read_sequence: reversed RNA read sequence
        :param _range: the range of the replacement
        :param bases: the base characters to replace the idx-element of read_sequence
        :return: read_sequence with the range elements replaced with bases
        """
        return read_sequence[:_range[0]] + bases + read_sequence[_range[1]:]

    def fast_find_max_so_far(self, max_so_far, len_read, mismatches):
        if max_so_far[3] < mismatches and max_so_far[1] < len_read:
            misses = max_so_far[3]
            changes = {}
            possible_changes = [b for b in BASES if b != max_so_far[2][len_read - max_so_far[1] - 1]]
            for change in possible_changes:
                new_read = self.replace_bases_at_indices(max_so_far[2], (len_read - max_so_far[1] - len(change),
                                                                         len_read - max_so_far[1]), change)
                new_range, new_length = exact_suffix_matches(new_read, self.M, self.occ)
                changes[change] = (new_range, new_length, new_read, misses + 1)
            for change in changes:
                test = self.fast_find_max_so_far(changes[change], len_read, mismatches)

                if test:
                    aligned, maxi = test
                    if aligned:
                        return aligned, maxi
            return False, ()
        elif max_so_far[1] == len_read and max_so_far[3] <= mismatches:
            return True, max_so_far
        else:
            return False, max_so_far

    def seed_finder(self, read_sequence, high_bound=None, mismatches=MAX_NUM_MISMATCHES, k=1):
        """
        Note: this function only looks for a single ungapped alignment of the READ_SEQUENCE in the genome sequence
        provided in __init__.

        Returns a tuple as documented below indicating whether an alignment has been found and if found, return a list
        of **single** tuple containing the start index of the alignment in the *actual (unreversed)* read_sequence,
        the start index of the alignment in the genome sequence, and the length of the alignment (which should be equal
        to len(read_sequence)).

        :param read_sequence: reversed RNA read sequence
        :param high_bound: the bound that all alignments need to be less than
        :param mismatches: maximum number of mismatched allowed
        :param k: maximum number of parts of alignment allowed
        :return: a tuple (aligned, seeds)
            aligned: a boolean indicating whether or not a **single** ungapped alignment was found
            seeds: [] if not aligned, otherwise a list of a single 3-element tuple containing the start of the alignment
            in *actual (unreversed)* read_sequence and the genome sequence, and the length of the alignment (which is
            equal to the length of the read_sequence)
        """
        len_read, seeds_len, seeds = len(read_sequence), 0, []
        if k == 1:
            max_so_far, updated = (None, 0, read_sequence, 0), True

            test, max_so_far = self.find_max_so_far(max_so_far, len_read, mismatches)
            if max_so_far:
                if max_so_far[1] == len_read:
                    seeds = self.check_seeds(max_so_far[0], max_so_far[1], high_bound)
                    if seeds:
                        return True, (seeds, max_so_far[1])
                    else:
                        return False, max_so_far[2]
                else:
                    return False, max_so_far[2]
            return False, ()

    def fast_seed_finder(self, read_sequence, high_bound=None, mismatches=MAX_NUM_MISMATCHES, k=1):
        len_read, seeds_len, seeds = len(read_sequence), 0, []
        if k == 1:
            _range, length = exact_suffix_matches(read_sequence, self.M, self.occ)
            max_so_far, updated = (_range, length, read_sequence, 0), True

            test, max_so_far = self.fast_find_max_so_far(max_so_far, len_read, mismatches)
            if max_so_far:
                if max_so_far[1] == len_read:
                    seeds = self.check_seeds(max_so_far[0], max_so_far[1], high_bound)
                    if seeds:
                        return True, (seeds, max_so_far[1])
                    else:
                        return False, max_so_far[2]
                else:
                    return False, max_so_far[2]
            return False, ()

    def align_to_genome(self, read_sequence):
        """
        Note: this function only looks for a single ungapped alignment of the READ_SEQUENCE in the genome sequence
        provided in __init__.

        Returns a tuple as documented below indicating whether an alignment has been found and if found, return a list
        of **single** tuple containing the start index of the alignment in the READ_SEQUENCE, the start index of the
        alignment in the genome sequence, and the length of the alignment (which should be equal to len(read_sequence)).

        :param read_sequence: RNA read sequence
        :return: a tuple (aligned, seeds)
            aligned: a boolean indicating whether or not a **single** ungapped alignment was found
            seeds: [] if not aligned, otherwise a list of a single 3-element tuple containing the start of the alignment
            in read_sequence and the genome sequence, and the length of the alignment (which is equal to the length of
            the read_sequence)
        """
        reverse_read_sequence = read_sequence[::-1]
        aligned, seeds = self.fast_seed_finder(reverse_read_sequence)
        if aligned:
            return True, [(0, seeds[0][0], 50)]
        else:
            return False, []
