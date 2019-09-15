import argparse
import sys
import os
from subprocess import run, CalledProcessError


class ToolMissingError(Exception): pass


def check_prerequisits():
    devnull = open(os.devnull, 'w')
    try:
        run(['bcftools', '--version'], stdout=devnull, stderr=devnull, check=True)
        run(['bgzip', '--version'], stdout=devnull, stderr=devnull, check=True)
        run(['tabix', '--version'], stdout=devnull, stderr=devnull, check=True)
        run(['whatshap', '--version'], stdout=devnull, stderr=devnull, check=True)
    except FileNotFoundError as e:
        raise ToolMissingError('{0} was not found'.format(e.filename)) from e
    except CalledProcessError as e:
        raise ToolMissingError('{0} has failed.'.format(' '.join(e.cmd))) from e     

     
def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='pedX')

    subparsers = parser.add_subparsers(help='modes', dest='sub')

    parser_rephase = subparsers.add_parser('rephase',
                                           help='Check if genotype of 10X phase sets needs \
                                                 needs to be switched in order to comply \
                                                 with trio phasing.')
    parser_rephase.add_argument('workdir',
                               type=str,
                               help='Working and output directory.')
    parser_rephase.add_argument('tenx',
                               type=str,
                               help='Phased variants by 10X in .vcf or .vcf.gz format.')
    parser_rephase.add_argument('trio',
                               type=str,
                               help='Phased variants by trio in .vcf or .vcf.gz format.')
    parser_rephase.add_argument('--trio_sample_idx',
                                 metavar='INT',
                                 type=str,
                                 default=0,
                                 help='Sample index of child in the trio VCF file.')
    parser_rephase.add_argument('--log2ratio_threshold',
                                 dest='thr',
                                 metavar='INT',
                                 type=str,
                                 default=1,
                                 help='Minimum log2ratio between different genotype \
                                       switched in order to end up in resulting \
                                       VCF file.')

    parser_integrate = subparsers.add_parser('integrate',
                                           help='Filter and integrate Trio VCF into \
                                                 rephased 10X VCF.')
    parser_integrate.add_argument('workdir',
                               type=str,
                               help='Working and output directory.')
    parser_integrate.add_argument('tenx_rephased',
                               type=str,
                               help='Phased 10X variants rephased by pedx rephase in .vcf or .vcf.gz format.')
    parser_integrate.add_argument('trio',
                               type=str,
                               help='Phased variants by trio in .vcf or .vcf.gz format.')
    parser_integrate.add_argument('sample_name',
                               type=str,
                               help='Name of sample in VCF files.')
                               
    parser_cluster = subparsers.add_parser('cluster',
                                            help='Add HP tag (1: maternal, 2: paternal) to alignments in BAM based on VCF.')
    parser_cluster.add_argument('workdir', 
                                 type=str,
                                 help='Input and output directory.')
    parser_cluster.add_argument('bam',
                                 type=str,
                                 help='Alignments to be tagged.')
    parser_cluster.add_argument('vcf',
                                 type=str,
                                 help='Phased variants used for clustering.')

    parser_svhap = subparsers.add_parser('svhap',
                                           help='Phase structural variants called with SVIM.')
    parser_svhap.add_argument('workdir',
                               type=str,
                               help='Working and output directory.')
    parser_svhap.add_argument('sv',
                               type=str,
                               help='Structural variants called by SVIM in .vcf or .vcf.gz format.')
    parser_svhap.add_argument('bam',
                               type=str,
                               help='Haplotagged long reads used to generate SV.vcf file.')
    
    return parser.parse_args(arguments)
