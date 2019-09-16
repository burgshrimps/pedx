import logging
from pysam import VariantFile, AlignmentFile
import os
import numpy as np

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

    num_phased_svs = 0
    for rec in vcf_in.fetch():
        print(rec.chrom + ':' + str(rec.pos))
        if rec.filter.keys()[0] == 'PASS':
            rec_sample = rec.samples[0]
            sv_chrom = rec.chrom
            sv_pos = rec.pos
            sv_read_ids = rec.info['READS']
            sv_support = rec.info['SUPPORT']

            begin_pos = sv_pos - 1
            if 'END' in rec.info:
                end_pos = rec.info['END']
            else:
                end_pos = sv_pos

            maternal_counter = 0
            paternal_counter = 0
            for read in bam_in.fetch(sv_chrom, begin_pos, end_pos):
                if read.query_name in sv_read_ids:
                    if read.has_tag('HP'):
                        read_hp = read.get_tag('HP')
                        maternal_counter += read_hp == 1
                        paternal_counter += read_hp == 2
            threshold = int(0.85 * sv_support)
            if (maternal_counter + paternal_counter) >= threshold:
                if rec_sample['GT'][0] == rec_sample['GT'][1]:
                    rec.samples[0].phased = True
                else:
                    if np.log2((maternal_counter + 1) / float(paternal_counter + 1)) > 1:
                        rec.samples[0]['GT'][0] = 1
                        rec.samples[0]['GT'][1] = 0
                        rec.samples[0].phased = True
                    elif np.log2((maternal_counter + 1) / float(paternal_counter + 1)) < -1:
                        rec.samples[0]['GT'][0] = 0
                        rec.samples[0]['GT'][1] = 1
                        rec.samples[0].phased = True
                vcf_out.write(rec)
