import pickle
from pysam import VariantFile
import os
import logging


def hamming_dist(x,y):
    """ Calculates Hamming Distance between two vectors """
    hd = 0
    for i in range(len(x)):
        hd += x[i] != y[i]
    return hd


def compare_genotypes(trio_records, tenx_records, tenx_phase_sets, sample_idx):
    """ Compares genotypes of 10X variants to genotypes of trio variants
    and calculates the log2-ratio of two possible scenarios:
    1. tenx.hap1 = trio.hap1 AND tenx.hap2 = trio.hap2
    2. tenx.hap1 = trio.hap2 AND tenx.hap2 = trio.hap1 """
    for ps_id in tenx_phase_sets:
        chrom = ps_id.split(':')[0]
        alleles_trio_hap1 = []
        alleles_trio_hap2 = []
        alleles_tenx_hap1 = []
        alleles_tenx_hap2 = []
        for pos in tenx_phase_sets[ps_id].positions:
            rec_id = chrom + ':' + str(pos)
            if rec_id in trio_records:
                trio_rec_sample = trio_records[rec_id].samples[sample_idx]
                if trio_rec_sample.phased:
                    tenx_rec_sample = tenx_records[rec_id].samples[0]
                    alleles_trio_hap1.append(trio_rec_sample.alleles[0])
                    alleles_trio_hap2.append(trio_rec_sample.alleles[1])
                    alleles_tenx_hap1.append(tenx_rec_sample.alleles[0])
                    alleles_tenx_hap2.append(tenx_rec_sample.alleles[1])

        overlap_with_trio = len(alleles_trio_hap1)
        if overlap_with_trio > 0:
            tenx_phase_sets[ps_id].rephased = True
            tenx_phase_sets[ps_id].trio_overlap = overlap_with_trio
            tenx_phase_sets[ps_id].hd_1 = hamming_dist(alleles_trio_hap1, alleles_tenx_hap1) + hamming_dist(alleles_trio_hap2, alleles_tenx_hap2) + 1 #  add pseudocounts for log-ratio
            tenx_phase_sets[ps_id].hd_2 = hamming_dist(alleles_trio_hap1, alleles_tenx_hap2) + hamming_dist(alleles_trio_hap2, alleles_tenx_hap1) + 1
            tenx_phase_sets[ps_id].set_log2ratio()

    return tenx_phase_sets


def write_rephased_tenx_vcf(tenx_vcf, tenx_records, tenx_phase_sets, threshold, workdir):
    """ Writes new 10X VCF file and switches genotypes if logratios above /
    below threshold """
    basename = os.path.basename(tenx_vcf)
    if basename.endswith('.vcf'):
        offset = -4
    elif basename.endswith('.vcf.gz'):
        offset = -7
    else:
        return
    tenx_rephased_vcf = workdir + '/' + basename[:offset] + '.filtered.het.rephased.vcf'

    vcf_in = VariantFile(tenx_vcf)
    vcf_out = VariantFile(tenx_rephased_vcf, 'w', header=vcf_in.header)

    for ps_id in tenx_phase_sets:
        if tenx_phase_sets[ps_id].rephased:
            chrom = tenx_phase_sets[ps_id].chrom
            if tenx_phase_sets[ps_id].log2ratio <= -threshold:
                for pos in tenx_phase_sets[ps_id].positions:
                    tenx_records[chrom + ':' + str(pos)].samples[0]['PS'] = 1
                    vcf_out.write(tenx_records[chrom + ':' + str(pos)])
            elif tenx_phase_sets[ps_id].log2ratio >= threshold:
                for pos in tenx_phase_sets[ps_id].positions:
                    tenx_records[chrom + ':' + str(pos)].samples[0]['PS'] = 1
                    GT_swapped = (tenx_records[chrom + ':' + str(pos)].samples[0]['GT'][1], tenx_records[chrom + ':' + str(pos)].samples[0]['GT'][0])
                    tenx_records[chrom + ':' + str(pos)].samples[0]['GT'] = GT_swapped
                    tenx_records[chrom + ':' + str(pos)].samples[0].phased = True
                    vcf_out.write(tenx_records[chrom + ':' + str(pos)])
    return tenx_rephased_vcf


def postprocess_tenx_rephased_vcf(tenx_rephased_vcf):
    """ Sorts, compresses and indexes 10X rephased VCF """
    tenx_rephased_sorted_vcf = tenx_rephased_vcf[:-4] + '.sorted.vcf'
    tenx_rephased_sorted_zipped_vcf = tenx_rephased_sorted_vcf + '.gz'

    command_sort = ['bcftools', 'sort', tenx_rephased_vcf, '>', tenx_rephased_sorted_vcf]
    command_zip = ['bgzip', tenx_rephased_sorted_vcf]
    command_index = ['tabix', tenx_rephased_sorted_zipped_vcf]
    command_rm_unsorted = ['rm', tenx_rephased_vcf]

    run(' '.join(command_sort), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_zip), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_index), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_rm_unsorted), shell=True, check=True, executable='/bin/bash')
