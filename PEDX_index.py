from pysam import VariantFile
import numpy as np
import pickle
import logging
import os


class PhaseSet:

    def __init__(self, set_name, chrom, positions):
        self.set_name = set_name
        self.chrom = chrom
        self.positions = positions  # positions of phased HET not filtered variants
        self.num_unphased = 0
        self.length_bp = 0
        self.length_var = 0
        self.trio_overlap = 0
        self.matTOhap1_patTOhap2 = 0
        self.matTOhap2_patTOhap1 = 0
        self.log2ratio = 0
        self.rephased = False

    def set_length(self):
        self.positions.sort()
        self.length_bp = self.positions[-1] - self.positions[0] + 1
        self.length_var = len(self.positions)

    def set_log2ratio(self):
        self.log2ratio = np.log2(self.matTOhap1_patTOhap2 / float(self.matTOhap2_patTOhap1))


def create_index_positions(filename, create_indexed_phase_sets=False):
    """ Creates dict with chroms as keys and list of positions as values.
    If chosen also creates dict with phase set information. """
    indexed_records = dict()
    indexed_records_per_chrom = dict()
    vcf_in = VariantFile(filename)

    # Create indices
    if not create_indexed_phase_sets:
        for rec in vcf_in.fetch():
            # Create positions index
            if rec.filter.keys()[0] == 'PASS':
                indexed_records[rec.chrom + ':' + str(rec.pos)] = rec
        return indexed_records

    else:
        indexed_phase_sets = dict()
        for rec in vcf_in.fetch():
            # Create positions index
            if rec.filter.keys()[0] == 'PASS':
                indexed_records[rec.chrom + ':' + str(rec.pos)] = rec

            # Create phase set index
            if rec.filter.keys()[0] == 'PASS':
                rec_sample = rec.samples[0]
                if rec_sample['GT'][0] != rec_sample['GT'][1]:
                    phase_set_id = rec.chrom + ':' + str(rec_sample['PS'])
                    if rec_sample.phased:
                        if phase_set_id in indexed_phase_sets:
                            indexed_phase_sets[phase_set_id].positions.append(rec.pos)
                        else:
                            indexed_phase_sets[phase_set_id] = PhaseSet(rec_sample['PS'], rec.chrom, [rec.pos])
                    else:
                        if phase_set_id in indexed_phase_sets:
                            indexed_phase_sets[phase_set_id].num_unphased += 1
        return indexed_records, indexed_phase_sets
