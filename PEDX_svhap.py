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
    phasing_stat_f = open(workdir + '/' + 'phasing_stat.txt')

    
    chr_to_include = ['1',
                      '2',
                      '3',
                      '4',
                      '5',
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
    
    """
    chr_to_include = ['chr1',
                      'chr2',
                      'chr3',
                      'chr4',
                      'chr5',
                      'chr6',
                      'chr7',
                      'chr8',
                      'chr9',
                      'chr10',
                      'chr11',
                      'chr12',
                      'chr13',
                      'chr14',
                      'chr15',
                      'chr16',
                      'chr17',
                      'chr18',
                      'chr19',
                      'chr20',
                      'chr21',
                      'chr22',
                      'chrX',
                      'chrY']
    """

    phasing_stat = {'INS' : {'Total':0, 'Phased HOM':0, 'Phased HET':0},
                    'DEL' : {'Total':0, 'Phased HOM':0, 'Phased HET':0},
                    'INV' : {'Total':0, 'Phased HOM':0, 'Phased HET':0},
                    'BND' : {'Total':0, 'Phased HOM':0, 'Phased HET':0},
                    'DUP:TANDEM' : {'Total':0, 'Phased HOM':0, 'Phased HET':0},
                    'DUP_INT' : {'Total':0, 'Phased HOM':0, 'Phased HET':0}}

    prev_chrom = ''
    for rec in vcf_in.fetch():
        sv_chrom = rec.chrom
        if sv_chrom in chr_to_include:
            if sv_chrom != prev_chrom:
                logging.info('Processing {0}'.format(sv_chrom))
            prev_chrom = sv_chrom
            if rec.filter.keys()[0] == 'PASS':
                sv_pos = rec.pos
                sv_read_ids = rec.info['READS']
                sv_support = rec.info['SUPPORT']
                sv_type = rec.info['SVTYPE']

                phasing_stat[sv_type]['Total'] += 1

                begin_pos = sv_pos - 1
                if 'END' in rec.info:
                    end_pos = rec.info['END']
                else:
                    end_pos = sv_pos

                hap1_counter = 0
                hap2_counter = 0
                try:
                    read_iterator = bam_in.fetch(sv_chrom, begin_pos-2000, end_pos+2000)
                except ValueError:
                    read_iterator = bam_in.fetch(sv_chrom, begin_pos, end_pos)
                for read in read_iterator:
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
                        phasing_stat[sv_type]['Phased HOM'] += 1
                    elif allele_frequency_hap1 >= threshold_het:
                        rec.samples[0]['GT'] = (1, 0)
                        rec.samples[0].phased = True
                        phasing_stat[sv_type]['Phased HET'] += 1
                    elif allele_frequency_hap2 >= threshold_het:
                        rec.samples[0]['GT'] = (0, 1)
                        rec.samples[0].phased = True
                        phasing_stat[sv_type]['Phased HET'] += 1

                    vcf_out.write(rec)
    
    phasing_stat_f.write('\tTotal\tPhased HOM\tPhased HET\n')
    for sv in phasing_stat:
        phasing_stat_f.write('{0}:\t{1}\t{2}\t{3}\n'.format(sv, phasing_stat[sv]['Total'], phasing_stat[sv]['Phased HOM'], phasing_stat[sv]['Phased HET']))
    phasing_stat_f.close()