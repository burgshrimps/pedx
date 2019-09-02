import argparse
import sys

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='pedX')

    subparsers = parser.add_subparsers(help='modes', dest='sub')

    parser_index = subparsers.add_parser('index',
                                         help='Create index of input vcf files \
                                               and addtitional index of 10X phase sets.')
    parser_index.add_argument('tenx',
                               type=str,
                               help='Phased variants by 10X in .vcf or .vcf.gz format.')
    parser_index.add_argument('trio',
                               type=str,
                               help='Phased variants by trio in .vcf or .vcf.gz format.')

    parser_rephase = subparsers.add_parser('rephase',
                                           help='Check if genotype of 10X phase sets needs \
                                                 needs to be switched in order to comply \
                                                 with trio phasing.')
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

    return parser.parse_args(arguments)
