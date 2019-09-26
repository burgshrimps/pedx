import logging
from pysam import VariantFile, AlignmentFile
import os
import numpy as np
import matplotlib.pyplot as plt

def phase_structural_variants(sv_vcf, long_reads_bam, workdir):
    sv_vcf_basename = os.path.basename(sv_vcf)
    if sv_vcf_basename.endswith('.vcf'):
        offset = -4
    elif sv_vcf_basename.endswith('.vcf.gz'):
        offset = -7
    else:
        return
    sv_filtered_phased_vcf = workdir + '/' + sv_vcf_basename[:offset] + '.filtered.phased.vcf'
    vcf_in = VariantFile(sv_vcf)
    vcf_out = VariantFile(sv_filtered_phased_vcf, 'w', header=vcf_in.header)
    bam_in = AlignmentFile(long_reads_bam)

    chr_to_include = ['1',
                      '2',
                      '3',
                      '4',
                      '5'
                      '6',
                      '7',
                      '8',
                      '9',
                      '10',
                      '11',
                      '12',
                      '13',
                      '14',
                      '15',
                      '16',
                      '17',
                      '18',
                      '19',
                      '20',
                      '21',
                      '22',
                      'X',
                      'Y']

    num_phased_svs = 0
    for rec in vcf_in.fetch():
        sv_chrom = rec.chrom
        if sv_chrom in chr_to_include:
            if rec.filter.keys()[0] == 'PASS':
                rec_sample = rec.samples[0]
                sv_pos = rec.pos
                sv_read_ids = rec.info['READS']
                sv_support = rec.info['SUPPORT']

                begin_pos = sv_pos - 1
                if 'END' in rec.info:
                    end_pos = rec.info['END']
                else:
                    end_pos = sv_pos
                print(sv_chrom, begin_pos, end_pos)

                hap1_counter = 0
                hap2_counter = 0
                for read in bam_in.fetch(sv_chrom, begin_pos, end_pos):
                    if read.query_name in sv_read_ids:
                        if read.has_tag('HP'):
                            read_hp = read.get_tag('HP')
                            hap1_counter += read_hp == 1
                            hap2_counter += read_hp == 2

                threshold_read_count = max(int(0.85 * sv_support), 5)
                threshold_het = 0.8
                threshold_hom = 0.2

                if (hap1_counter + hap2_counter) >= threshold_read_count:
                    allele_frequency_hap1 = hap1_counter / float(hap1_counter + hap2_counter)
                    allele_frequency_hap2 = hap2_counter / float(hap1_counter + hap2_counter)

                    if allele_frequency_hap1 >= threshold_hom and allele_frequency_hap1 < threshold_het:
                        rec.samples[0]['GT'] = (1, 1)
                        rec.samples[0].phased = True
                    elif allele_frequency_hap1 >= threshold_het:
                        rec.samples[0]['GT'] = (1, 0)
                        rec.samples[0].phased = True
                    elif allele_frequency_hap2 >= threshold_het:
                        rec.samples[0]['GT'] = (0, 1)
                        rec.samples[0].phased = True

                    vcf_out.write(rec)
