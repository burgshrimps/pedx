import logging
import os
import sys

from PEDX_input_parsing import parse_arguments
from PEDX_index import create_index_positions
from PEDX_rephase import compare_genotypes


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


    if options.sub == 'index':
        logging.info('MODE: index')
        logging.info('10X VCF: {0}'.format(options.tenx))
        logging.info('TRIO VCF: {0}'.format(options.trio))

        logging.info('########## CREATE 10X POSITIONS INDEX ##########')
        create_index_positions(options.tenx, phase_set=True)
        logging.info('########## CREATE TRIO POSITIONS INDEX ##########')
        create_index_positions(options.trio)
        logging.info('Wrote index files to same directories as input VCF files.')
    elif options.sub == 'rephase':
        logging.info('MODE: rephase')
        logging.info('10X VCF: {0}'.format(options.tenx))
        logging.info('TRIO VCF: {0}'.format(options.trio))
        logging.info('TRIO SAMPLE IDX: {0}'.format(options.trio_sample_idx))
        logging.info('########## COMPARE GENOTYPES ##########')
        compare_genotypes(options.tenx, options.trio, options.trio_sample_idx)


main()
