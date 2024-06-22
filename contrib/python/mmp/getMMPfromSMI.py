#!/usr/bin/env python

###################################################################
# Summary: Script to generate matched moleular pairs from input smiles
#
# Example Usage:
#         python getMMPfromSMI.py --smi_file /home/my_username/input.smi
#         python getMMPfromSMI.py -i input.smi -o output.csv
#         python getMMPfromSMI.py -i input.smi -o output.csv -c double -f remove_nonrings -r rdkit
#
#
#############################################################

import argparse
import sys
import logging

from mmp.mmp_objects import MMPObjectClass

#
# Main
#
def main():

    parser = argparse.ArgumentParser(description='Generate matched molecular pairs from input smiles')
    parser.add_argument("-i", required=True, help="Input smiles file")
    parser.add_argument("-o", required=True, help="Output pairs file CSV")
    parser.add_argument("-c", default='BOTH',
                        choices=['SINGLE', 'DOUBLE', 'BOTH', 'single', 'double', 'both'],
                        help="SINGLE|DOUBLE|BOTH (default = BOTH)")
    parser.add_argument("-f", default='NONE',
                        help="REMOVE_NONRINGS|REMOVE_RINGS|NONE")
#    parser.add_option("-e", "--prescan_smi", dest="prescan_smi", action="store_true", \
#                      help="Invoke if you want to avoid error check of smiles file via pre-scan", default=True)
    parser.add_argument('-b', type=float, default=0.50001,
                        help='Dicer Threshold for Fragmentation, default = 0.50001')
    parser.add_argument("-l", default=None,
                        help="Name of file to place log info in")
    parser.add_argument("-d", default=None,
                        help="CRITICAL|ERROR|WARNING|INFO|DEBUG")
    parser.add_argument("-t", action="store_true", default=False,
                        help="Invoke this flag (no value needed) if you want a second header line with datatypes")

    args = parser.parse_args()
    types_header = args.t

    smi_fi = args.i
    out_fi = args.o

    # check cut type is valid
    cut_type = args.c

    cut_type_options = ['SINGLE', 'DOUBLE', 'BOTH']
    if not any(opt in cut_type.upper() for opt in cut_type_options):
        parser.print_help()
        parser.error('Need to specify a valid cut type consult --help\n')
        sys.exit()

    # check filter type is valid
    filter_type = args.f
    filter_type_options = ['REMOVE_NONRINGS', 'REMOVE_RINGS', 'NONE']
    if not any(opt in filter_type.upper() for opt in filter_type_options):
        parser.print_help()
        parser.error('Need to specify a valid filter type consult --help\n')
        sys.exit()

    try:
        dicer_threshold = float(args.b)
        if dicer_threshold < 0.1 or dicer_threshold > 0.9:
            sys.exit('Please specify a value for the dicer threshold between 0.1 and 0.9')
    except:
        sys.exit('Please specify a value for the dicer threshold between 0.1 and 0.9')

    # logging

    log_file = args.l
    log_level = args.d

    mmplogger = logging.getLogger('lillymol_file_logger')

    if log_file is not None:

        fh = logging.FileHandler(log_file)

        if log_level is None:
            mmplogger.setLevel(logging.INFO)
        else:
            numeric_level = getattr(logging, log_level.upper(), None)
            if not isinstance(numeric_level, int):
                raise ValueError('Invalid log level: %s' % log_level)
            mmplogger.setLevel(numeric_level)

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(module)s: %(message)s',
                                      datefmt='%m/%d/%Y %I:%M:%S %p')

        fh.setFormatter(formatter)
        mmplogger.addHandler(fh)
 
    else:
        logging.disable(logging.CRITICAL)

    #
    mmplogger.info('Instantiating MMP Objects')
    my_mmp_object = MMPObjectClass(mmplogger)

    # pre-check the smiles for odd stuff
    my_mmp_object.scan_input_smiles(smi_fi)

    # execute dicer cmd and parse results to mmp_object
    mmplogger.info('Parse Dicer Fragments to Pairs')
    my_mmp_object.build_from_dicer(smi_fi, cut_type, filter_type, threshold=dicer_threshold)

    mmplogger.info('Write out final pairs')
    if types_header:
        my_mmp_object.print_to_file(out_fi, cut_type, inc_types_header=types_header)
    else:
        my_mmp_object.print_to_file(out_fi, cut_type)
    mmplogger.info('Complete.')


if __name__ == '__main__':

    main()
    sys.exit()
