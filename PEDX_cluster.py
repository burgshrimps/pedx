import logging
from subprocess import run
import os


def cluster_reads_by_hap(workdir, in_bam, in_vcf):
    """ Add HP tag (1: maternal, 2: paternal) to alignments in BAM file based on VCF file """
    bam_basename = os.path.basename(in_bam)
    bam_haplotagged = workdir + '/' + bam_basename[:-4] + '.haplotagged.bam'
    
    command_haplotag = ['whatshap', 'haplotag', '-o', bam_haplotagged, in_vcf, in_bam, '&>>', workdir + '/haplotag.log']
    command_idx_bam = ['samtools', 'index', bam_haplotagged]
    
    logging.info(' -> Run haplotag')
    run(' '.join(command_haplotag), shell=True, check=True, executable='/bin/bash')
    
    logging.info(' -> Index haplotagged BAM file')
    run(' '.join(command_idx_bam), shell=True, check=True, executable='/bin/bash')
