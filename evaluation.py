import itertools

from shared import *


def index_isoform_locations(known_isoforms, unknown_isoforms):
    """
    known_isoforms: a python set of all known isoforms (listed in genes.tab)
    unknown_isoforms: a python set of all unknown isoforms (listed in genes.tab)
    """
    genome_isoform_offsets = {}
    for isoform in itertools.chain(known_isoforms, unknown_isoforms):
        isoform_offset = 0
        for exon in isoform.exons:
            for i in range(exon.start, exon.end):
                genome_isoform_offsets.setdefault(i, []).append((isoform, isoform_offset))
                isoform_offset += 1
    return genome_isoform_offsets


def evaluate_alignment(genome_sequence, read_sequence, alignment, unknown_isoforms, genome_isoform_offsets):
    """
    genome_sequence: sequence of the genome, doesn't matter if '$'-terminated or not
    read_sequence: sequence of the read being evaluated
    alignment: alignment returned by Aligner.align
    unknown_isoforms: a python set of all unknown isoforms (listed in genes.tab)
    genome_isoform_offsets: output from calling index_isoform_locations(known_isoforms, unknown_isoforms)
    """
    isoform_offset_to_count = {}
    read_index = 0
    r_end_prev = -1
    for r_start, g_start, length in alignment:
        assert r_start >= r_end_prev, 'alignment segment %s, %s, %s overlaps with previous alignment segment' % (
        r_start, g_start, length)
        for i, g_i in enumerate(range(g_start, g_start + length)):
            if genome_sequence[g_i] == read_sequence[r_start + i]:
                for isoform, isoform_offset in genome_isoform_offsets.get(g_i, []):
                    key = (isoform, isoform_offset - read_index)
                    isoform_offset_to_count[key] = isoform_offset_to_count.get(key, 0) + 1
            read_index += 1
        r_end_prev = r_start + length
    if len(isoform_offset_to_count) == 0:
        return CASE_UNALIGNED, 0
    num_matches, isoform = max([(count, isoform) for (isoform, _), count in isoform_offset_to_count.items()])

    num_mismatches = len(read_sequence) - num_matches
    if num_mismatches <= MAX_NUM_MISMATCHES:
        if isoform in unknown_isoforms:
            return CASE_HIDDEN_GENE, num_mismatches
        else:
            return CASE_GENE, num_mismatches
    return CASE_UNALIGNED, 0
