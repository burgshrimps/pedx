import logging
import os
import sys

from PEDX_input_parsing import parse_arguments, check_prerequisits
from PEDX_index import create_index_positions
from PEDX_rephase import compare_genotypes, write_rephased_tenx_vcf, postprocess_tenx_rephased_vcf
from PEDX_integrate import filter_trio_vcf


def main():
    # Fetch command line argumnets
    options = parse_arguments()

    if not options.sub:
        print('Choose between the modes ("index", "rephase", "cluster", "svcall" or "svhap").')
        return

    # Set up logging
    logFormatter = logging.Formatter('%(asctime)s %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info('############### Start PEDX ###############')
    logging.info('CMD: python3 {0}'.format(' '.join(sys.argv)))
    for arg in vars(options):
        logging.info('PARAMETER: {0}, VALUE: {1}'.format(arg, getattr(options, arg)))

    logging.info('# Check prerequisits')
    check_prerequisits()

    if options.sub == 'rephase':
        logging.info('MODE: rephase')
        logging.info('WORKDIR: {0}'.format(options.workdir))
        logging.info('10X VCF: {0}'.format(options.tenx))
        logging.info('TRIO VCF: {0}'.format(options.trio))
        logging.info('TRIO SAMPLE IDX: {0}'.format(options.trio_sample_idx))
        logging.info('THRESHOLD: {0}'.format(options.thr))
        logging.info('# Create index for trio VCF')
        trio_indexed_records = create_index_positions(options.trio)
        logging.info('# Create index for 10X VCF')
        tenx_indexed_records, tenx_phase_sets = create_index_positions(options.tenx, create_indexed_phase_sets=True)
        logging.info('# Compare genotypes')
        tenx_phase_sets_with_logratios = compare_genotypes(trio_indexed_records, tenx_indexed_records, tenx_phase_sets, options.trio_sample_idx)
        logging.info('# Write rephased 10X VCF')
        tenx_rephased_vcf = write_rephased_tenx_vcf(options.tenx, tenx_indexed_records, tenx_phase_sets_with_logratios, options.thr, options.workdir)
        logging.info('# Postprocess rephased 10X VCF')
        postprocess_tenx_rephased_vcf(tenx_rephased_vcf)

    elif options.sub == 'integrate':
        logging.info('MODE: integrate')
        logging.info('WORKDIR: {0}'.format(options.workdir))
        logging.info('REPHASED 10X VCF: {0}'.format(options.tenx_rephased))
        logging.info('TRIO VCF: {0}'.format(options.trio))
        logging.info('SAMPLE NAME: {0}'.format(options.sample_name))
        logging.info('# Filter trio VCF')
        filter_trio_vcf(options.trio, options.workdir, options.sample_name)


main()
