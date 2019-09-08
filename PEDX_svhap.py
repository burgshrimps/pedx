import logging
from pysam import VariantFile, AlignmentFile
import os

def phase_structural_variants(sv_vcf, long_reads_bam, workdir):
    sv_vcf_basename = os.path.basename(sv_vcf)
    if sv_vcf_basename.endswith('.vcf'):
        offset = -4
    elif sv_vcf_basename.endswith('.vcf.gz'):
        offset = -7
    else:
        return
    sv_filtered_phased_vcf = workdir + sv_vcf_basename[:offset] + '.filtered.phased.vcf'
    vcf_in = VariantFile(sv_vcf)
    vcf_out = VariantFile(sv_filtered_phased_vcf, 'w', header=vcf_in.header)
    bam_in = AlignmentFile(long_reads_bam)

    num_phased_svs = 0
    for rec in vcf_in.fetch():
        if rec.filter.keys()[0] == 'PASS':
            sv_chrom = rec.chrom
            sv_pos = rec.pos 
            sv_read_ids = rec.info['READS']

            maternal_counter = 0
            paternal_counter = 0
            for read in bam_in.fetch(sv_chrom, sv_pos-1, sv_pos):
                if read.query_name in sv_read_ids:
                    if read.has_tag('HP'):
                        read_hp = read.get_tag('HP')
                        maternal_counter += read_hp == 1
                        paternal_counter += read_hp == 2
            if maternal_counter > 2 or paternal_counter > 2:
                num_phased_svs += 1
    print(num_phased_svs)




