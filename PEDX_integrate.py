import logging
from subprocess import run
import os
from pysam import VariantFile


def get_filtered_phased_het_trio_variants(trio_vcf, trio_filtered_het_phased_vcf, sample_name):
    vcf_in = VariantFile(trio_vcf)
    vcf_in.subset_samples([sample_name])
    vcf_out = VariantFile(trio_filtered_het_phased_vcf, 'w', header=vcf_in.header)
    
    for rec in vcf_in.fetch():
    	if rec.filter.keys()[0] == 'PASS':
    	    rec_sample = rec.samples[0]
    	    if rec_sample.phased and rec_sample['GT'][0] != rec_sample['GT'][1]:
    	        rec.samples[0].update({'PS':1})
    	        vcf_out.write(rec)
    return 0
    	        

def filter_trio_vcf(trio_vcf, workdir, sample_name):
    """ Uses command line tools to filter trio VCF file and add PS tag """
    trio_vcf_basename = os.path.basename(trio_vcf)
    if trio_vcf_basename.endswith('.vcf'):
        offset = -4
    elif trio_vcf_basename.endswith('.vcf.gz'):
        offset = -7
    else:
        return
    tmp_header = workdir + '/tmp_header.vcf'
    tmp_variants = workdir + '/tmp_variants.vcf'
    tmp_reheadered = workdir + '/tmp_reheadered.vcf'
    trio_filtered_het_phased_vcf = workdir + '/' + trio_vcf_basename[:offset] + '.filtered.het.phased.pstag.vcf'
    trio_filtered_het_phased_zipped_vcf = trio_filtered_het_phased_vcf + '.gz'
    
    command_get_header = ['bcftools', 'view', '-h', trio_vcf, '>', tmp_header]
    command_modify_header = 'sed -i \'5i##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ID of Phase Set for Variant\">\' ' + str(tmp_header)
    command_get_variants = ['bcftools', 'view', '-H', trio_vcf, '>', tmp_variants]
    command_reheader = ['cat', tmp_header, tmp_variants, '>', tmp_reheadered]
    command_zip = ['bgzip', trio_filtered_het_phased_vcf]
    command_index = ['tabix', trio_filtered_het_phased_zipped_vcf]
    command_clean = ['rm', workdir + '/tmp*']
    
    logging.info(' -> Adding PS FORMAT to header')
    run(' '.join(command_get_header), shell=True, check=True, executable='/bin/bash')
    run(command_modify_header, shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_get_variants), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_reheader), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Write filtered, phased and heterozygous variants to {0}'.format(trio_filtered_het_phased_vcf))
    get_filtered_phased_het_trio_variants(tmp_reheadered, trio_filtered_het_phased_vcf, sample_name)
    
    logging.info(' -> Compress VCF file')
    run(' '.join(command_zip), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Index VCF file')
    run(' '.join(command_index), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Clean temporary files')
    run(' '.join(command_clean), shell=True, check=True, executable='/bin/bash')
    
    return trio_filtered_het_phased_zipped_vcf


def merge_trio_10X_vcf(tenx_rephased, trio_filtered, workdir):
    """ Merge filtered trio VCF and rephased 10x VCF """
    tenx_trio_merged_vcf = workdir + '/10X_and_trio_merged.vcf'
    tenx_trio_merged_sorted_vcf = tenx_trio_merged_vcf[:-4] + '.sorted.vcf'
    tenx_trio_merged_sorted_zipped_vcf = tenx_trio_merged_sorted_vcf + '.gz'
    
    command_merge = ['bcftools', 'concat', '-a', '-d', 'all', tenx_rephased, trio_filtered, '>', tenx_trio_merged_vcf]
    command_sort = ['bcftools', 'sort', tenx_trio_merged_vcf, '>', tenx_trio_merged_sorted_vcf]
    command_zip = ['bgzip', tenx_trio_merged_sorted_vcf]
    command_index = ['tabix', tenx_trio_merged_sorted_zipped_vcf]
    command_rm = ['rm', tenx_trio_merged_vcf]
    
    logging.info(' -> Merge 10X and trio VCF files to {0}'.format(tenx_trio_merged_vcf))
    run(' '.join(command_merge), shell=True, check=False, executable='/bin/bash')
    
    logging.info(' -> Sort merged VCF file')
    run(' '.join(command_sort), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Compress VCF file')
    run(' '.join(command_zip), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Index VCF file')
    run(' '.join(command_index), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Remove intermediate VCF file')
    run(' '.join(command_rm), shell=True, check=True, executable='/bin/bash')
