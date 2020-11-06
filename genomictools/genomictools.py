#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genomic Tools: handy classes to manage genomic data.
"""

#===============================================================================
#
# Author: David Piñeyro
# Contact info: dapineyro.dev@gmail.com
# Date: 2020-10-06
# Version: 0.0.1
#
# License: GPL-3.
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#===============================================================================

#===============================================================================
# IMPORTS
import argparse
import os
import pyranges as pr
import numpy as np
import pandas as pd
import multiprocessing
from Bio import SeqIO


#===============================================================================
# GLOBAL VARIABLES
SCRIPT_DESCRIPTION = """
This is a collection of useful classes to handle genomic sequences:
    FragmentedSeq: a class to represent fragmented sequences, e.g. exons of a
        gene.
    GeneAnnotated: a class to represent annotated genes.
    CodonCounter: a class to perform codon counts.

Note: the sample code that is executed with this script demonstrate some
    GeneAnnotated functionality for the PRF-1 project (Lida). It was developed
    for a previous version of GeneAnnotated class and is not tested.

Version: 0.0.1
Author: David Piñeyro
Date: 2020-10-06
License: GPL-3
"""
# This global var is not used but I conserved it just in case it will.
SELENOPROTEINS_HUMAN = [
    'DIO1',
    'DIO2',
    'DIO3',
    'GPX1',
    'GPX2',
    'GPX3',
    'GPX4',
    'GPX6',
    'SELENOF',
    'SELENOH',
    'SELENOI',
    'SELENOK',
    'SELENOM',
    'SELENON',
    'SELENOO',
    'SELENOP',
    'MSRB1',
    'SELENOS',
    'SELENOT',
    'SELENOV',
    'SELENOW',
    'SEPHS2',
    'TXNRD1',
    'TXNRD2',
    'TXNRD3'
]


#===============================================================================
# FUNCTIONS

#===============================================================================
# CLASSES

class FragmentedSeq():
    """
    Container class for a fragmented seq. It's meant to contain introns, exons,
    CDS and other sequences not linearly sited in the genome. By the moment it
    only accepts sequences with an exon_id.

    Attributes:
        chromosome [string]: chromosome id. E.g: 'chr14', 'chrM'.
        start [int]: start position, reversed if strand == '-'.
        end [int]: end position, reversed if strand == '-'.
        strand [string]: sequence strandness ('+' or '-').
        frame [int/NaN]: an int indicating the CDS frame of the fragment or
            NaN for non-coding fragments.
        exon_id [string]: enembl exon identifier.
        seq [Bio.SeqRecord.SeqRecord]: the transcribed sequence of the fragment.
            It will be reverse complemented if strand == '-'.
    Methods:
    """
    def __init__(self, annot_row, genome):
        """
        Arguments:
            annot_row [pd.Series]: a Pandas Series object, usually comming from
                a DataFrame row of an annotation DataFrame, ideally generated
                using PyRanges package by reading a GTF file.
            genome [list]: list of Bio.SeqRecord.SeqRecord objects, one per
                chromosome.
        """
        self.chromosome = annot_row.loc['Chromosome']
        self.start = annot_row.loc['Start']
        self.end = annot_row.loc['End']
        self.strand = annot_row.loc['Strand']
        self.frame = annot_row.loc['Frame']
        self.exon_id = annot_row.loc['exon_id']
        self.seq = GeneAnnotated.get_genomic_seq(self.chromosome,
                                                 self.start,
                                                 self.end,
                                                 genome)
        self.seq.seq = self.seq.seq.transcribe()
        if self.strand == '-':
            self.seq = self.seq.reverse_complement()
            if self.start < self.end:
                self.start, self.end = self.end, self.start

    @property
    def frame(self):
        return self._frame
    @frame.setter
    def frame(self, val):
        if val == '.':
            self._frame = np.nan
        else:
            self._frame = int(val)


class GeneAnnotated():
    """
    Store an annotated version of a gene, with its sequences.

    Attributes:
        gene_name [string]: common gene name. E.g: 'ADAM20'.
        annotations [DataFrame]: an annotation DataFrame, ideally generated
            using PyRanges package by reading a GTF file.
        gene_ids [numpy.ndarray]: an array of unique gene_ids for a given
            gene_name. Usually one, but there are exceptions.
        genomic [list]: a list of gene dicts generated by self._get_gene.
            Usually of length 1, but there are cases with more than one
            gene_id for gene_name.
        transcripts [list]: a list of list of transcript dicts generated by
            self._get_transcripts. Usually of length 1, but there are cases
            with more than one gene_id for gene_name.
    Methods:
        GeneAnnotated.get_genomic_seq [staticmethod]
        self._get_gene
        self._get_transcripts
    """
    def __init__(self, gene_name, gtf, genome, verbose = False):
        """
        Arguments:
            gene_name [string]: common gene name. E.g: 'ADAM20'.
            gtf [pyranges]: a pyranges object generated by using
                pyranges.read_gtf. It's tested wth Gencode GTF annotations.
            genome [list]: list of Bio.SeqRecord.SeqRecord objects, one per
                chromosome.
            verbose [bool]: whether to print warning messages during
                self._get_transcripts or not. Default = False.
        """

        self.gene_name = gene_name
        self.annotations = gtf[gtf.gene_name == gene_name].df
        # In case this gene is not found in the annotations.
        if self.annotations.empty:
            print(f'[WNG] Gene {gene_name} is not found in the annotations.')
            self.gene_ids = None
            self.genomic = None
            self.transcripts = None

        else:
            # There could be different annotated genes with the same gene_name.
            self.gene_ids = pd.unique(self.annotations.loc[:, 'gene_id'].dropna())
            self.genomic = [self._get_gene(
                self.gene_name,
                self.annotations.loc[self.annotations.loc[:, 'gene_id'] == g_id, :],
                genome) for g_id in self.gene_ids]
            self.transcripts = [self._get_transcripts(
                self.gene_name,
                self.annotations.loc[self.annotations.loc[:, 'gene_id'] == g_id, :],
                genome, verbose) for g_id in self.gene_ids]
    @staticmethod
    def get_genomic_seq(chromosome, start, end, genome):
        """
        Collect a genomic sequence.

        Arguments:
            chromosome [string]: the chromosome id.
            start [int]: start position of the gene.
            end [int]: end position of the gene.
            genome [list]: list of Bio.SeqRecord.SeqRecord objects, one per
                chromosome.
        Return:
            seq [Bio.SeqRecord.SeqRecord]: the genomic sequence of the provided
                position.
        """
        # Check start is the smallest number.
        if start > end:
            start, end = end, start
        genome_ids = [g.id for g in genome]
        seq = genome[genome_ids.index(chromosome)][start:end].upper()
        return seq

    def _get_gene(self, gene_name, annot, genome):
        """
        Creates a gene dictionary.

        Arguments:
            gene_name [string]: common gene name.
            annot [DataFrame]: annotations derived from a GTF object imported
                using PyRanges and subsetted by gene_name.
            genome [list]: list of Bio.SeqRecord.SeqRecord objects, one per
                chromosome.

        Return:
            gene [dict]: a dictionary with gene features:
                'chromosome': chromosome [string],
                'start': start [int]: it is reversed for - strand seqs,
                'end': end [int]: it is reversed for - strand seqs,
                'gene_id': gene_id [string]: ensembl gene identifier,
                'gene_name': common gene name,
                'gene_type': gene_type [string]: protein_coding, etc,
                'seq': seq [Bio.SeqRecord.SeqRecord]: reversed for - strand seq
        """
        chromosome = pd.unique(annot.loc[annot.loc[:, 'Feature'] == 'gene',
                                        'Chromosome'])
        assert (len(chromosome) == 1), \
            f'[ERROR] GTF malformed for gene: {gene_name}.'

        chromosome = chromosome[0]
        start = int(annot.loc[annot.loc[:, 'Feature'] == 'gene', 'Start'])
        end = int(annot.loc[annot.loc[:, 'Feature'] == 'gene', 'End'])
        gene_id = annot.loc[annot.loc[:, 'Feature'] == 'gene', 'gene_id'].iloc[0]
        gene_type = annot.loc[annot.loc[:, 'Feature'] == 'gene', 'gene_type'].iloc[0]
        strand = annot.loc[annot.loc[:, 'Feature'] == 'gene', 'Strand'].iloc[0]
        seq = self.get_genomic_seq(chromosome, start, end, genome)
        if strand == '-':
            seq = seq.reverse_complement()
            if start < end:
                start, end = end, start

        gene = {'chromosome':chromosome,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'gene_type': gene_type,
                'seq': seq}
        return gene

    def _get_transcripts(self, gene_name, annot, genome, verbose):
        """
        Creates a list of Transcript objects, one for each gene transcript.

        Arguments:
            gene_name [string]: common gene name.
            annot [DataFrame]: annotations derived from a GTF object imported
                using PyRanges and subsetted by gene_name.
            genome [list]: list of Bio.SeqRecord.SeqRecord objects, one per
                chromosome.
            verbose [bool]: whether to print warning messages during
                self._get_transcripts or not. Default = False.

        Return:
            transcripts [list]: a list of dictionaries, each with the following
                transcript features:
                'chromosome': chromosome [string],
                'start': start [int]: it is reversed for - strand seqs,
                'end': end [int]: it is reversed for - strand seqs,
                'strand': strand [string],
                'transcript_id': transcript_id [string],
                'transcript_name': transcript_name [string],
                'transcript_type': transcript_type [string],
                'transcript_tag': transcript_tag [string],
                'seq_unprocessed': seq [Bio.SeqRecord.SeqRecord]: reversed
                    for - strand seq,
                'seq_processed': seq [Bio.SeqRecord.SeqRecord]:
                    merge_fragments_seq(exons). May contain UTRs.
                'seq_cds': seq [Bio.SeqRecord.SeqRecord]:
                    merge_fragments_seq(cds). Shouldn't contain UTRs.
                'exons': exons [list]: list of FragmentedSeq for exons,
                'cds': cds [list]: list of FragmentedSeq for CDS,
                'start_codon': start_codon [list ] list of FragmentedSeq with
                    the start codons for the transcript. Normally just one,
                'stop_codon': stop_codon [list ] list of FragmentedSeq with
                    the stop codons for the transcript. Normally just one,
                'utrs': utrs [list]: list of FragmentedSeq for UTR
        """
        possible_start_codons = ['AUG']
        possible_stop_codons = ['UAG', 'UAA', 'UGA']
        transcript_ids = pd.unique(annot.loc[:, 'transcript_id'].dropna())
        transcripts = []
        for t in transcript_ids:
            annot_t = annot.loc[annot.loc[:, 'transcript_id'] == t, :]
            transcript_id = t
            transcript_name = annot_t.loc[:, 'transcript_name'].values[0]
            transcript_type = annot_t.loc[:, 'transcript_type'].values[0]
            transcript_tag = annot_t.loc[:, 'tag'].values[0]
            chromosome = annot_t.loc[:, 'Chromosome'].values[0]
            start = annot_t.loc[annot_t.loc[:, 'Feature'] == 'transcript',
                                'Start'].values[0]
            end = annot_t.loc[annot_t.loc[:, 'Feature'] == 'transcript',
                              'End'].values[0]
            strand = annot_t.loc[annot_t.loc[:, 'Feature'] == 'transcript',
                                'Strand'].values[0]
            seq_unprocessed = self.get_genomic_seq(chromosome, start, end, genome)
            seq_unprocessed.seq = seq_unprocessed.seq.transcribe()
            if 'exon' in annot_t['Feature'].values:
                exons = list(annot_t.loc[annot_t.loc[:, 'Feature'] == 'exon',
                                         :].apply(lambda row: FragmentedSeq(row,
                                                                            genome),
                                                  axis = 1)
                            )
            else:
                exons = []
            if 'CDS' in annot_t['Feature'].values:
                cds = list(annot_t.loc[annot_t.loc[:, 'Feature'] == 'CDS',
                                         :].apply(lambda row: FragmentedSeq(row,
                                                                            genome),
                                                  axis = 1)
                           )
            else:
                cds = []
            if 'start_codon' in annot_t['Feature'].values:
                start_codon = list(annot_t.loc[annot_t.loc[:, 'Feature'] == 'start_codon',
                                         :].apply(lambda row: FragmentedSeq(row,
                                                                            genome),
                                                  axis = 1)
                                   )
                # Check correct start codons.
                if not all(str(s.seq.seq) in possible_start_codons for s in start_codon):
                    idx_wrong = [str(s.seq.seq) in possible_start_codons for s in start_codon]
                    wrong_codon = [str(s.seq.seq) for s in start_codon if str(s.seq.seq) not in possible_start_codons]
                    if verbose:
                        print(f'[WNG] One start codon of {gene_name}',
                          f'is splitted and/or non-canonical: {wrong_codon}')
            else:
                start_codon = []
            if 'stop_codon' in annot_t['Feature'].values:
                stop_codon = list(annot_t.loc[annot_t.loc[:, 'Feature'] == 'stop_codon',
                                         :].apply(lambda row: FragmentedSeq(row,
                                                                            genome),
                                                  axis = 1)
                           )
                # Check correct stop codons.
                if not all(str(s.seq.seq) in possible_stop_codons for s in stop_codon):
                    idx_wrong = [str(s.seq.seq) in possible_stop_codons for s in stop_codon]
                    wrong_codon = [str(s.seq.seq) for s in stop_codon if str(s.seq.seq) not in possible_stop_codons]
                    if verbose:
                        print(f'[WNG] One stop codon of {gene_name}',
                          f'is splitted and/or non-canonical: {wrong_codon}')
            else:
                stop_codon = []
            if 'UTR' in annot_t['Feature'].values:
                utrs = list(annot_t.loc[annot_t.loc[:, 'Feature'] == 'UTR',
                                         :].apply(lambda row: FragmentedSeq(row,
                                                                            genome),
                                                  axis = 1)
                            )
            else:
                utrs = []
            if strand == '-':
                seq_unprocessed = seq_unprocessed.reverse_complement()
                if start < end:
                    start, end = end, start
            transcript = {'chromosome': chromosome,
                          'start': start,
                          'end': end,
                          'strand': strand,
                          'transcript_id': transcript_id,
                          'transcript_name': transcript_name,
                          'transcript_type': transcript_type,
                          'transcript_tag': transcript_tag,
                          'seq_unprocessed': seq_unprocessed,
                          'seq_processed': self.merge_fragments_seq(exons),
                          'seq_cds': self.merge_fragments_seq(cds),
                          'exons': exons,
                          'cds': cds,
                          'start_codon': start_codon,
                          'stop_codon': stop_codon,
                          'utrs': utrs}
            transcripts.append(transcript)
        return transcripts
    def merge_fragments_seq(self, fragments_list):
        """
        Merge a list of FragementedSeq objects into a single sequence. Note
        that the fragments cannot be overlapping.

        Arguments:
            fragments_list [list]: a list of FragmentedSeq objects to merge.
        Results:
            seq_merged [Bio.SeqRecord.SeqRecord]: the merged sequence.
        """
        if len(fragments_list) == 0:
            return None
        elif len(fragments_list) == 1:
            return fragments_list[0].seq
        else:
            positions = [f.start for f in fragments_list]
            order_idx = np.argsort(positions)
            if fragments_list[0].strand == '-':
                order_idx = order_idx[::-1]
            seq_merged = fragments_list[order_idx[0]].seq
            for i in order_idx[1:]:
                seq_merged += fragments_list[i].seq
            return(seq_merged)
    def get_dist_prf_stop_splicing(self, motif):
        """
        Specific function for Lida's project. She was interested in -1 PRF
        (-1 Predicted Ribosomal Frameshift, see: http://prfdb.umd.edu/). This
        method will search for a given heptamer in all the CDS for the gene.
        The heptamer should be in frame. Once found, it will calculate the
        generated incorrect STOP codon due to slippery and reports this new
        STOP codon as well as the distance to the closest downstream splicing
        junction.

        Arguments:
            motif [string]: an heptamer motif to search for. It should be
                in RNA.
        Return:
            heptamer_list_gene [list]: a list[gene_ids] of list[transcripts],
                each with the following dict:
                heptamer_positions [list]: list of ints with the start
                    positions of the heptamer. Genomic coordinates.
                new_stop_positions [list]: list of ints with the  start
                    positions of the new stop codons generated by the -1
                    frameshifting, one for each heptamer found. Genomic
                    coordinates.
                closest_splicing [list]: position of the closest splicing
                    junction for each new stop codon. Genomic coordinates.
                distance_closest_splicing [int]: nucleotide distance between
                    each new stop codon and its closest splicing junction.
        """
        possible_stop_codons = ['UAG', 'UAA', 'UGA']
        heptamer_list_gene  = []
        for genes in range(len(self.gene_ids)):
            heptamer_list_transcripts = []
            for transcript in self.transcripts[genes]:
                if len(transcript['cds']) == 0:
                    heptamer_positions = []
                else:
                    cds = str(transcript['seq_cds'].seq)
                    heptamer_positions = GeneAnnotated.get_motif_pos_in_minus_1_frame(cds, motif)
                if len(heptamer_positions) == 0:
                    new_stop_positions = []
                    closest_splicings = []
                    distance_closest_splicing = []
                else:
                    new_stop_positions = []
                    # Get splicing sites relative coordinates, sorted.
                    splicings = []
                    for cds_exon in transcript['cds']:
                        splicings.append(cds_exon.start)
                        splicings.append(cds_exon.end)
                    if transcript['strand'] == '-':
                        splicings = sorted(splicings, reverse = True)
                        splicings = [splicings[i] -  splicings[i+1] for i in range(0, len(splicings), 2)]
                    else:
                        splicings = sorted(splicings)
                        splicings = [splicings[i+1] -  splicings[i] for i in range(0, len(splicings), 2)]
                    # Recursive sum to obtain relative splicing sites.
                    splicings = [sum(splicings[:i+1]) for i in range(len(splicings))]
                    for h in heptamer_positions:
                        # New frame start exactly at the heptamer start.
                        cds_from_slippery = cds[h:]
                        all_new_stops = [t for t in range(0, len(cds_from_slippery), 3) if cds_from_slippery[t:t+3] in possible_stop_codons]
                        if len(all_new_stops) == 0:
                            new_stop_position = np.nan
                        else:
                            new_stop_position = min(all_new_stops)
                        new_stop_positions.append(new_stop_position)
                    # Note: control the case in which the new stop is in the
                    # last exon.
                    new_stop_positions = list(np.array(heptamer_positions) + np.array(new_stop_positions))
                    if len(splicings) == 1:  # Only 1 exon.
                        closest_splicings = [np.nan]
                        distance_closest_splicing = [np.nan]
                    else:
                        closest_splicings = [min([s for s in splicings if s >= p]) if p <= max(splicings) else np.nan for p in new_stop_positions]
                        distance_closest_splicing = list(np.array(closest_splicings) - np.array(new_stop_positions))
                d = {'heptamer_positions': heptamer_positions,
                     'new_stop_positions': new_stop_positions,
                     'closest_splicings': closest_splicings,
                     'distance_closest_splicing': distance_closest_splicing}
                heptamer_list_transcripts.append(d)
            heptamer_list_gene.append(heptamer_list_transcripts)
        return heptamer_list_gene

    @staticmethod
    def get_motif_pos_in_minus_1_frame(query, motif):
        """
        Get the start position of a -1 motif . Both argument sequences
        should be of the same type (DNA or RNA) and both are expected in
        uppercase.

        Arguments:
            query [string]: the query string to search for the motif. Typically
                a CDS (RNA) sequence.
            motif [string]: a sequence to search for. It should be <= than the
                query and should be of the same type (DNA or RNA).

        Return:
            pos [list]: a list of ints with the relative positions found for
                the motif.
        """
        if len(motif) > len(query):
            print('[WNG] Query sequence smaller that the searched motif.')
            return []
        else:
            # Note that the motif should be at -1 respect translation frame.
            pos = [t for t in range(2, len(query), 3) if motif == query[t:t+len(motif)]]
            return pos

class CodonCounter():
    """Codon counter object.

    This object, initializated with a list of codons (SeqRecords), is designed
    to count each appearance of its codons in a given annotated gene
    (GeneAnnotated).

    Parameters
    ----------
    codons : list
        A list of Bio.SeqRecord.SeqRecord objects with the triplet (codon)
        sequences to search for.

    Attributes
    ----------
    codons : list
        A list of Bio.SeqRecord.SeqRecord objects with the triplet (codon)
        sequences to search for.

    """
    def __init__(self, codons):
        self.codons = codons

    @property
    def codons(self):
        return self._codons
    @codons.setter
    def codons(self, codons):
        if len(codons) < 1:
            raise ValueError('[ERROR] Codons file is empty.')
        else:
            c = [c.upper() for c in codons]
            self._codons = c

    def count_codons(self, gene):
        """Count each codon appearance in each gene CDS.

        This function finds each codon in each CDS sequence. It uses
        gencode tags to select a selection of transcripts (see
        valid_tags variable).

        Parameters
        ----------
        gene : GeneAnnotated
            A GeneAnnotated object corresponding to a gene with its transcripts
            and sequences (see GeneAnnotated documentation).

        Returns
        -------
        codon_count : list
            A list of lists, one for each protein coding cds, containing:
                [gene_name, tr_name, tr_len, counts]. Being:
                gene_name : string
                    A string with the common gene name.
                tr_name : string
                    A string with the transcript name. Is the CDS name.
                tr_len : int
                    The lenght (in nucleotides) of the CDS searched.
                counts : dict
                    A dict with the codons as keys and absolute index
                    appearances of this coding, in the particular cds. Thus,
                    if 'AUG' is searched in 'AUGAAUGCCAUG' will be:
                    {'AUG': [0, 9]}
        """
        valid_tags = ['basic', 'CCDS', 'seleno', 'appris_principal_1',
                      'appris_principal', 'appris_candidate',
                      'appris_candidate_ccds']
        codon_count = []
        gene_name = gene.gene_name
        # Each annotated gene usually contains only one transcript dict, but
        # there are a few of cases with more than one.
        for i in range(len(gene.transcripts)):
            t = gene.transcripts[i]
            # Select only 'protein_coding' transcripts.
            for tr in t:
                tr_name = tr['transcript_name']
                if ((tr['transcript_type'] == 'protein_coding') and
                    (tr['transcript_tag'] in valid_tags)):
                    tr_len = len(tr['seq_cds'])
                    counts = {}
                    for c in self.codons:
                        if ((tr['chromosome'] == 'chrM' and c.name[:2] == 'mt') or
                           (tr['chromosome'] != 'chrM' and c.name[:2] != 'mt')):
                            count = {c.name: self.find_codons(c.seq, tr['seq_cds'].seq)}
                        else:
                            count = {c.name: np.nan}
                        counts.update(count)
                    codon_count.append([gene_name, tr_name, tr_len, counts])
        return codon_count

    def find_codons(self, codon, cds):
        """Find codon location in a given cds.

        This function takes a single codon and a single cds and returns a
        list with the index positions of each codon found, or an empty list
        if the codon was not found in the given cds.

        Parameters
        ----------
        codon : {Bio.SeqRecord.SeqRecord,  string}
            Three letters string with the codon to search.
        cds : {Bio.SeqRecord.SeqRecord, string}
            The CDS sequence in wich look for codon appearance. It should be
            in uppercase letters and as RNA seq (A,G,U,C letters).

        Return
        ------
        count : list
            A list of index positions of the codons found. There are the CDS
            positions of the first nucleotide of the triplet (codon). An
            empty list is returned when no codons are found.
        """
        count = [n * 3
                    for n, c in enumerate(
                        [cds[i:i+3] for i in range(0, len(cds), 3)])
                            if str(c) == str(codon)]
        return count

    @staticmethod
    def codon_counts_2_df(codon_count):
        """Convert codon_count to a pandas.DataFrame.

        This function uses the codon_count list returned by self.count_codons
        method and returns a pandas.DataFrame from it.

        Parameters
        ----------
        codon_count : list
            A list of lists, one for each protein coding cds, containing:
                [gene_name, tr_name, tr_len, counts]. Being:
                gene_name : string
                    A string with the common gene name.
                tr_name : string
                    A string with the transcript name. Is the CDS name.
                tr_len : int
                    The lenght (in nucleotides) of the CDS searched.
                counts : dict
                    A dict with the codons as keys and absolute index
                    appearances of this coding, in the particular cds. Thus,
                    if 'AUG' is searched in 'AUGAAUGCCAUG' will be:
                    {'AUG': [0, 9]}.

        Return
        ------
        count_df : pandas.DataFrame
            A DataFrame with the info from codon_count.
        """
        # Generate a DataFrame with all the codon_counts. Note that the
        # genes without codon_counts (non-coding RNAs normally) will be
        # lost, due to the inner loop of the list comprehension. Each
        # row will be a transcript.
        count_df_pre = pd.DataFrame([i for c in codon_count for i in c],
                                    columns = ['Gene',
                                               'Transcript',
                                               'Trx_len',
                                               'Profile'])
        profiles = pd.DataFrame(list(count_df_pre['Profile']))
        count_df = pd.concat([count_df_pre.iloc[:, : -1], profiles],
                             axis = 1)

        return count_df

    @staticmethod
    def counts_by_gene(count_df, collapsing_method):
        """Collapses counts per transcript to counts per gene.

        This function takes a count_df pandas.DataFrame generated by
        self.codon_counts_2_df and collapses all codon counts to a
        single row per gene.

        Parameters
        ----------
        count_df : pandas.DataFrame
            A DataFrame generated by self.codon_counts_2_df. Essentially,
            it should have the following columns:
                'Gene': the gene name by which it will be collapsed.
                'Transcript': transcript name.
                'Trx_len': length of the counted CDS.
                ...: the rest of columns, one per each codon counted.
        collapsing_method : string
            One of the two options:
                'longer': the resulting count will come from the longer
                    cds available for that gene.
                'canonical': the resulting count will come from the transcript
                    with the first id. For instance, from a list of 'AAAS-202',
                    'AAAS-201' and 'AAAS-208' transcripts for AAAS gene, the
                    selected will be AAAS-201. All transcripts should have a
                    unique annotation.

        Return
        ------
        count_df_by_gene : pandas.DataFrame
            The collapsed DataFrame, with the following columns.
                'Gene': the gene name by which it will be collapsed.
                'Trx_len': the CDS length of the corresponding count (see
                    `collapsing_method` parameter help.
                ...: the rest of columns, one per each codon counted.
        """
        count_df_gr = count_df.groupby(by = 'Gene')
        # Collapse each df using the appropiate method.
        def collapse(df, method):
            if method == 'longer':
                # Selecting the longest transcript. idxmax will select only
                # one of the transcripts, avoiding problems with genes with
                # more than one longer transcript. It returns a series,
                # so it should be converted to a DataFrame.
                gene = df.loc[df['Trx_len'].idxmax(), :].to_frame().T
            elif method == 'canonical':
                # Selecting the longest transcript. At least in Gencode
                # human annotation, it corresponds to the one with the
                # first transcript name (in alpha-numerical order),i.e.
                # from a list of 'AAAS-202', 'AAAS-201', 'AAAS-208'
                # transcripts for AAAS gene, the selected will be
                # AAAS-201. All transcripts should have a unique annotation.
                gene = df.loc[df['Transcript'] == df['Transcript'].sort_values().values[0], :]
            else:
                raise ValueError('[ERROR] collapsing_method should' +
                                 ' be one of the following: \'longer\',' +
                                 ' \'canonical\'.')
            return gene
        count_df_gr_col = [collapse(v,
                                    collapsing_method) for k,v in count_df_gr]
        # Concatenate all df in groups.
        count_df_by_gene = pd.concat(count_df_gr_col)

        return count_df_by_gene

    @staticmethod
    def simplify(count_df_by_gene):
        """Simplify counts DataFrame.

        This function takes a counts DataFrame, previously collapsed by gene
        using self.counts_by_gene and generates a simplified DataFrame with
        the required statistics to further processing.

        Parameters
        ----------
        count_df_by_gene : pandas.DataFrame
            The collapsed DataFrame, with the following columns.
                'Gene': the gene name by which it will be collapsed.
                'Trx_len': the CDS length of the corresponding count (see
                    `collapsing_method` parameter help.
                ...: the rest of columns, one per each codon counted.

        Return
        ------
        count_simplified : pandas.DataFrame
            A DataFrame with the same row number of `count_df_by_gene`,but
            containing the required statistics only. It should contain the
            following columns:
                'Gene': the gene name.
                'Trx_codon_len': number of codons of the transcript.
                ...: for each codon counted:
                    '<codon_id>_n': number of codons found in the transcript.
                    '<codon_id>_prop': `<codon_id>_n` / `Trx_codon_len`.
                    '<codon_id>_max_stretch: longer stretch (codons found in
                        a row) in the transcript.
        """
        # Vectorized version of max_stretch that can accept a list of
        # iterables, or a pd.Series of lists, as will be the case.
        v_max_stretch = np.vectorize(CodonCounter.max_stretch,
                                     otypes = [float])

        def simplificator(x):
            """Function to create the new DataFrame from x (a row of
            the original DataFrame.
            """
            gene = x['Gene']
            trx_codon_len = x['Trx_len'] / 3
            n = x[3:].str.len()
            prop = n / trx_codon_len
            max_str = pd.Series(v_max_stretch(x[3:]))
            n.index = n.index + '_n'
            prop.index = prop.index + '_prop'
            max_str.index = x[3:].index + '_max_stretch'
            gene_s = pd.Series(gene, index = ['Gene'])
            trx_codon_len_s = pd.Series(trx_codon_len, index = ['Trx_codon_len'])
            s = pd.concat([gene_s, trx_codon_len_s, n, prop, max_str])
            return s

        count_simplified = count_df_by_gene.apply(simplificator, axis = 1)

        return count_simplified


    @staticmethod
    def max_stretch(s):
        """Count max codon stretch length.

        This function counts the longest consecutive stretch of codons. This
        function was designed to be used by self.simplify, but can be used
        standalone.

        Parameters
        ----------
        s : {list, np.array, pd.Series, etc}
            An iterable with the codon start positions. Then, the
            minimum difference between two consecutive elements is 3.
            The iterable is supposed to be ascending sorted.

        Return
        ------
        str_final : int
            The length of the longest stretch.
        """
        try:
            l = len(s)
        except TypeError:
            return np.nan
        if l == 0:
            return 0
        else:
            str_final = 1
            str_temp = 1
            for i in range(len(s)-1):
                if s[i] + 3 == s[i+1]:
                    str_temp += 1
                else:
                    if str_temp > str_final:
                        str_final = str_temp
                        str_temp = 1
                    else:
                        str_temp = 1
            return str_final

    @staticmethod
    def stats(count_simplified):
        """Compute summatory statistics.

        This function generates summatory statistics from a count_simplified
        DataFrame generated by self.simplify.

        Parameters
        ----------
        count_simplified : pandas.DataFrame
            A DataFrame with the same row number of `count_df_by_gene`,but
            containing the required statistics only. It should contain the
            following columns:
                'Gene': the gene name.
                'Trx_codon_len': number of codons of the transcript.
                ...: for each codon counted:
                    '<codon_id>_n': number of codons found in the transcript.
                    '<codon_id>_prop': `<codon_id>_n` / `Trx_codon_len`.
                    '<codon_id>_max_stretch: longer stretch (codons found in
                        a row) in the transcript.

        Return
        ------
        stats_dict : dict
            A dictionary with the following key, values:
            'total_transcripts': int
                Total number of transcripts in count_simplified.
            'total_codons': int | float
                Sum of all transcript codons in count_simplified. It should
                be an int, but a float is possible when truncated
                transcripts are included.
            'max_codons': float
                Max number of codons found in a transcript.
            'mean_codons': float
                Mean codon count of the transcripts in count_simplified.
            'median_codons': float
                Median codon count of the transcripts in count_simplified.
            'total_codon': pandas.Series
                Total codon count per analyzed codon.
            'max_codons_per_transcript': pandas.Series
                Max codon count per analyzed transcript.
            'mean_codons_per_transcript': pandas.Series
                Mean codon count per analyzed transcript.
            'median_codons_per_transcript': pandas.Series
                Median codon count per analyzed transcript.
            'max_codon_prop': pandas.Series
                Max codon proportion per analyzed transcript.
            'mean_codon_prop': pandas.Series
                Mean codon proportion per analyzed transcript.
            'median_codon_prop': pandas.Series
                Meadian codon proportion per analyzed transcript.
            'max_codon_stretch': pandas.Series
                Max codon stretch length (in number of codons) per analyzed
                transcript.
            'mean_codon_stretch': pandas.Series
                Mean codon stretch length (in number of codons) per
                analyzed transcript.
            'median_codon_stretch': pandas.Series
                Median codon stretch length (in number of codons) per
                analyzed transcript.

        """
        total_transcripts = count_simplified.shape[0]
        total_codons = count_simplified['Trx_codon_len'].sum()
        max_codons = count_simplified['Trx_codon_len'].max()
        mean_codons = count_simplified['Trx_codon_len'].mean()
        median_codons = count_simplified['Trx_codon_len'].median()

        total_codon = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_n$', regex = True)
        ].sum(skipna = True, axis = 0)
        max_codons_per_transcript = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_n$', regex = True)
        ].max(skipna = True, axis = 0)
        mean_codons_per_transcript = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_n$', regex = True)
        ].mean(skipna = True, axis = 0)
        median_codons_per_transcript = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_n$', regex = True)
        ].median(skipna = True, axis = 0)


        max_codon_prop = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_prop$', regex = True)
        ].max(skipna = True, axis = 0)
        mean_codon_prop = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_prop$', regex = True)
        ].mean(skipna = True, axis = 0)
        median_codon_prop = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_prop$', regex = True)
        ].median(skipna = True, axis = 0)


        max_codon_stretch = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_max_stretch$', regex = True)
        ].max(skipna = True, axis = 0)
        mean_codon_stretch = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_max_stretch$', regex = True)
        ].mean(skipna = True, axis = 0)
        median_codon_stretch = count_simplified.loc[
            :,
            count_simplified.columns.str.contains('_max_stretch$', regex = True)
        ].median(skipna = True, axis = 0)

        stats_dict = {
            'total_transcripts': total_transcripts,
            'total_codons': total_codons,
            'max_codons': max_codons,
            'mean_codons': mean_codons,
            'median_codons': median_codons,
            'total_codon': total_codon,
            'max_codons_per_transcript': max_codons_per_transcript,
            'mean_codons_per_transcript': mean_codons_per_transcript,
            'median_codons_per_transcript': median_codons_per_transcript,
            'max_codon_prop': max_codon_prop,
            'mean_codon_prop': mean_codon_prop,
            'median_codon_prop': median_codon_prop,
            'max_codon_stretch': max_codon_stretch,
            'mean_codon_stretch': mean_codon_stretch,
            'median_codon_stretch': median_codon_stretch
        }

        return stats_dict

    def rand_test(self, stats_dict, gtf, genome, iterations,
                  collapsing_method, cpus):
        """Randomization test to assess stats_dict.

        This provides a multiprocessing enhanced randomization test
        to calculate experimental p.Values and adj.P-values for some of the
        statistics of stats_dict.

        Parameters
        ----------
        stats_dict : dict
            The returned dict from self.stats (see help(self.stats)).
        gtf : PyRanges
            The annotation DataFrame generated using pyranges.read_gtf.
            It's tested with Gencode GTF annotations.
        genome : list
            List of Bio.SeqRecord.SeqRecord objects, one per chromosome.
        iterations : int
            Number of random iterations for the randomization test.
        collapsing_method : string
            The collapsing_method to be used in self.counts_by_gene.
            It could be 'longer' or 'canonical'. See
            help(self.counts_by_gene) for a more detailed explanation.
        cpus : int
            Number of CPUs to use in parallel calculations.

        Return
        ------
        A tuple with the following elements:
        real_data_df : pandas.DataFrame
            stats_dict with the real results converted to a DataFrame and
            with all the np.NaN values converted to 0.0.
        random_data_df : pandas.DataFrame
            A DataFrame with the same columns as real_data_df with the
            result of each randomization iteration per row.
        stats_pval : pandas.DataFrame
            A DataFrame with the info of stats_dict, and the associated pValues
            obtained using randomization tests. The collumns are the items in
            stats_dict, the rows are:
                'Sample_data': real values obtained in the sample.
                'Random_data_mean': random data mean.
                'Random_data_sd': random data standard deviation.
                'Number_equal_or_smaller': number of random samples equal to or
                    smaller than sample.
                'P_val_equal_or_smaller': p_value associated to the number of
                    equal or smaller. Essentially:
                    'Number_equal_or_smaller' / iterations.
                'Number_equal_or_greater': number of random samples equal to or
                    greater than sample.
                'P_val_equal_or_greater': p_value associated to the number of
                    equal or greater. Essentially:
                    'Number_equal_or_greater' / iterations.
        """
        # Extract a list of all the annotated genes from the GTF.
        all_genes = gtf.gene_name.unique()

        ### Run in parallel or in serial depending on Python version.
        # In Python < 3.8 there is a strict limit for an object to
        # be pickeable. So, passing gtf and genome in each iteration
        # is too much. Then, for python < 3.8 run it serial. otherwise
        # run it in parallel.
        #
        # Parallel version.
        # -----------------
        # Annotate all the genes in the genome. This is a very intesive task,
        # then using the multiprocessing version.
        a_g_zip = [(g, gtf, genome) for g in all_genes]
        with multiprocessing.Pool(processes = cpus) as pool:
            all_genes_annot = pool.starmap(GeneAnnotated, a_g_zip)
        #
        # Serial version.
        # ---------------
        #all_genes_annot = [GeneAnnotated(g, gtf, genome) for g in all_genes]
        #
        ### End of this section.

        # Count all codons from all genes. Again, very intesive.
        with multiprocessing.Pool(processes = cpus) as pool:
            codon_count = pool.map(self.count_codons, all_genes_annot)
        # Tidy up results and save a DataFrame.
        count_df = self.codon_counts_2_df(codon_count)
        # Summaryze codon counts by gene.
        count_df_by_gene = self.counts_by_gene(count_df,
                                               collapsing_method)
        # Generate a new df with simplified statistics.
        count_simplified = self.simplify(count_df_by_gene)

        # Randomization loop.
        # The number of genes to select in each iteration.
        n = stats_dict['total_transcripts']
        # A list with the reduced version of count_simplified for each
        # iteration.
        count_list = [count_simplified.sample(n = n,
                                              replace = False,
                                              axis = 0
                                             ) for _ in range(iterations)]
        with multiprocessing.Pool(processes = cpus) as pool:
            stats_list = pool.map(self.stats, count_list)

        # Covert stats_list to a pandas DataFrame.
        # Create a list of dics were each element is a single value.
        stats_list_to_df = []
        for e in stats_list:
            new_dict = {}
            for k, v in e.items():
                if isinstance(v, pd.Series):
                    new_dict.update({
                        k + '_' + v.index.values[i]: v.values[i]
                            for i in range(len(v))})
                else:
                    new_dict[k] = v
            stats_list_to_df.append(new_dict)
        # Convert to df.
        random_data_df = pd.DataFrame(stats_list_to_df)
        # Convert all NaN to 0.0.
        random_data_df.fillna(0.0, inplace = True)

        # Convert stats_dict into a DataFrame to easily compare with stats_df
        stats_dict_new = {}
        for k, v in stats_dict.items():
            if isinstance(v, pd.Series):
                stats_dict_new.update({
                    k + '_' + v.index.values[i]: v.values[i]
                        for i in range(len(v))})
            else:
                stats_dict_new[k] = v
        real_data_df = pd.DataFrame([stats_dict_new])
        # Convert all nan to 0.0.
        real_data_df.fillna(0.0, inplace = True)

        # Calculate statistics based on random test.
        n_eq_smaller = pd.DataFrame(
            random_data_df.apply(lambda col:
                np.where(col.sort_values() <= float(real_data_df[col.name].values))[0].size,
                axis = 0)
        ).T
        n_eq_greater= pd.DataFrame(
            random_data_df.apply(lambda col:
                np.where(col.sort_values() >= float(real_data_df[col.name].values))[0].size,
                axis = 0)
        ).T
        p_val_eq_smaller = n_eq_smaller / iterations
        p_val_eq_greater = n_eq_greater / iterations
        random_data_mean = pd.DataFrame(random_data_df.mean(axis = 0)).T
        random_data_sd = pd.DataFrame(random_data_df.std(axis = 0)).T

        # Generate final df.
        stats_pval = pd.concat([real_data_df,
                                random_data_mean,
                                random_data_sd,
                                n_eq_smaller,
                                p_val_eq_smaller,
                                n_eq_greater,
                                p_val_eq_greater])
        stats_pval.index = ['Sample_data',
                            'Random_data_mean',
                            'Random_data_sd',
                            'Number_equal_or_smaller',
                            'P_val_equal_or_smaller',
                            'Number_equal_or_greater',
                            'P_val_equal_or_greater']

        return random_data_df, real_data_df, stats_pval

