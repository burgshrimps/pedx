import pickle

from PEDX_index import get_index_root_name

def hamming_dist(x,y):
    """ Calculates Hamming Distance between two vectors """
    hd = 0
    for i in range(len(x)):
        hd += x[i] != y[i]
    return hd


def compare_genotypes(tenx_vcf, trio_vcf, sample_idx):
    """ Compares genotypes of 10X variants to genotypes of trio variants
    and calculates the log2-ratio of two possible scenarios:
    1. tenx.hap1 = trio.hap1 AND tenx.hap2 = trio.hap2
    2. tenx.hap1 = trio.hap2 AND tenx.hap2 = trio.hap1 """
    trio_pos_idx_file, _ = get_index_root_name(trio_vcf)
    tenx_pos_idx_file, tenx_phase_sets_idx_file = get_index_root_name(tenx_vcf)
    trio_positions = pickle.load(open(trio_pos_idx_file, 'rb'))
    tenx_positions = pickle.load(open(tenx_pos_idx_file, 'rb'))
    tenx_phase_sets = pickle.load(open(tenx_phase_sets_idx_file, 'rb'))
    for ps in tenx_phase_sets:
        print(ps)
