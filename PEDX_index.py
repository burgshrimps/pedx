from pysam import VariantFile
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
        self.hd_1 = 0
        self.hd_2 = 0
        self.log2ratio = 0

    def set_length(self):
        self.positions.sort()
        self.length_bp = self.positions[-1] - self.positions[0] + 1
        self.length_var = len(self.positions)

    def set_log2ratio(self):
        self.log2ratio = np.log2(self.hd_1 / float(self.hd_2))


def get_index_root_name(filename):
    """ Helper function to get filenames of index files """
    if filename.endswith('.vcf'):
        offset = -4
    elif filename.endswith('.vcf.gz'):
        offset = -7
    index_pos_file = filename[:offset] + '_pos_idx.pickle'
    index_phase_sets_file = filename[:offset] + '_phase_set_idx.pickle'
    return index_pos_file, index_phase_sets_file


def create_index_positions(filename, phase_set=False):
    """ Creates dict with chroms as keys and list of positions as values.
    If chosen also creates dict with phase set information. """
    positions = dict()
    vcf_in = VariantFile(filename)

    # Prepare output filename
    index_pos_file, index_phase_sets_file = get_index_root_name(filename)

    # Create indices
    if not phase_set:
        for rec in vcf_in.fetch():
            # Create positions index
            if rec.chrom in positions:
                positions[rec.chrom].append(rec.pos)
            else:
                positions[rec.chrom] = [rec.pos]
                logging.info(' -> ' + rec.chrom)

    else:
        phase_sets = dict()
        for rec in vcf_in.fetch():
            # Create positions index
            if rec.chrom in positions:
                positions[rec.chrom].append(rec.pos)
            else:
                positions[rec.chrom] = [rec.pos]
                logging.info(' -> ' + rec.chrom)

            # Create phase set index
            if rec.filter.keys()[0] == 'PASS':
                rec_sample = rec.samples[0]
                if rec_sample['GT'][0] != rec_sample['GT'][1]:
                    phase_set_id = rec.chrom + ':' + str(rec_sample['PS'])
                    if rec_sample.phased:
                        if phase_set_id in phase_sets:
                            phase_sets[phase_set_id].positions.append(rec.pos)
                        else:
                            phase_sets[phase_set_id] = PhaseSet(rec_sample['PS'], rec.chrom, [rec.pos])
                    else:
                        if phase_set_id in phase_sets:
                            phase_sets[phase_set_id].num_unphased += 1
        pickle.dump(phase_sets, open(index_phase_sets_file, 'wb'))

    pickle.dump(positions, open(index_pos_file, 'wb'))
