# Genomic Tools

__Author:__ David Piñeyro
__Date:__ 2020-11-05
__Version:__ 0.0.1
__License:__ GPL-3

This package is a collection of variables, classes and functions meant
to process genomic data.

As now, it contains three handy classes:

- GeneAnnotated: this class is menat to represent a gene, generating
an object with all the infromation as well as relevant sequences that
can be obtained from a given annotation (GTF format). It also have some
relevant functionality.

- FragmentedSeq: a class to represent fragmented sequences belonging to the
same unit. To be more specific, it's meant to represent exons and introns
from a given gene.

- CodonCounter: making use of a GeneAnnotated object, this class is meant to
produce codon usage counting and some other codon related functionalities.


Contact: dapineyro.dev@gmail.com
