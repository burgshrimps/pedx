import logging
from subprocess import run
import os


def filter_trio_vcf(trio_vcf, workdir, sample_name):
    """ Uses command line tools to filter trio VCF file and add PS tag """
    tmp_header_vcf = workdir + '/tmp_header.vcf'
    tmp_variants_vcf = workdir + '/tmp_variants.vcf'
    trio_vcf_basename = os.path.basename(trio_vcf)
    if trio_vcf_basename.endswith('.vcf'):
        offset = -4
    elif trio_vcf_basename.endswith('.vcf.gz'):
        offset = -7
    else:
        return
    trio_filtered_het_phased_vcf = workdir + '/' + trio_vcf_basename[:offset] + '.filtered.het.phased.pstag.vcf'
    trio_filtered_het_phased_zipped_vcf = trio_filtered_het_phased_vcf + '.gz'

    command_get_header = ['bcftools', 'view', '-h', '-s', sample_name, trio_vcf, '>', tmp_header_vcf]
    command_modify_header = ['sed', '-i', '4i##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">', tmp_header_vcf]
    command_get_variants = ['bcftools', 'view', '-H', '-s', sample_name, trio_vcf, '|', 'grep', '1|0\|0|1', '|', 'awk', 'BEGIN{OFS="\t"} {$9=$9":PS"; $10=$10":1"} { if ($7 == "PASS") print $0}', '>', tmp_variants_vcf]
    command_cat_header_variants = ['cat', tmp_header_vcf, tmp_variants_vcf, '>', trio_filtered_het_phased_vcf]
    command_zip = ['bgzip', trio_filtered_het_phased_vcf]
    command_index = ['tabix', trio_filtered_het_phased_zipped_vcf]
    command_rm = ['rm', tmp_header_vcf, tmp_variants_vcf]

    run(' '.join(command_get_header), shell=True, check=True, executable='/bin/bash')

    f = open(tmp_header_vcf, 'r')
    header_lines = f.readlines()
    f.close()
    phase_set_header_line = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">'
    header_lines.insert(3, phase_set_header_line)
    f = open(tmp_header_vcf, 'w')
    header_lines = ''.join(header_lines)
    f.write(header_lines)
    f.close()

    """
    run(' '.join(command_modify_header), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_get_variants), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_cat_header_variants), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_zip), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_index), shell=True, check=True, executable='/bin/bash')
    run(' '.join(command_rm), shell=True, check=True, executable='/bin/bash')
    """
