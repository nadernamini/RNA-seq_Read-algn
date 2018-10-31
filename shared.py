# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #

BASES = ['A', 'C', 'G', 'T']
TERMINATOR = '$'
MAX_NUM_MISMATCHES = 6

# Where a read is aligned to
CASE_GENE = 'gene'
CASE_HIDDEN_GENE = 'hidden_gene'
CASE_UNALIGNED = 'unaligned'


# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #

class IdElement(object):
    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id


# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #

class Gene(IdElement):
    def __init__(self, gene_id, isoforms):
        """
        gene_id: a string representing the unique id of this gene
        isoforms: a list of Isoform objects representing this gene's isoforms
        """
        assert type(gene_id) == str and type(isoforms) == list
        for isoform in isoforms:
            assert isinstance(isoform, Isoform)
        self.id = gene_id
        self.isoforms = isoforms

    def __repr__(self):
        return 'gene\t%s\t%s' % (self.id, ';'.join(isoform.id for isoform in self.isoforms))


# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #

class Isoform(IdElement):
    def __init__(self, isoform_id, exons):
        """
        isoform_id: a string representing unique id of this isoform
        exons: a list of Exon objects representing this isoform's exons 
        """
        assert type(isoform_id) == str and type(exons) == list
        for exon in exons:
            assert isinstance(exon, Exon)
        self.id = isoform_id
        self.exons = exons

    def __repr__(self):
        return 'isoform\t%s\t%s' % (self.id, ';'.join(exon.id for exon in self.exons))


# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #

class Exon(IdElement):
    def __init__(self, exon_id, start, end):
        """
        exon_id: a string representing unique id of this exon
        start: index of start of exon inclusive
        end: index of end of exon exclusive
        """
        assert type(exon_id) == str and type(start) == int and type(end) == int
        self.id = exon_id
        self.start = start
        self.end = end

    def __repr__(self):
        return 'exon\t%s\t%s\t%s' % (self.id, self.start, self.end)

# DO NOT MODIFY THIS FILE. WE WILL USE OUR OWN COPY DURING GRADING #
