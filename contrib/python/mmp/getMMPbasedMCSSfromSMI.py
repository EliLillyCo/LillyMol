#!/usr/bin/env python

###################################################################
# Summary: Script to generate largest MCSS based MMP for each pair
#          using input smiles
#
# Example Usage:
#         python getMCSSbasedMMPfromSMI.py -i input.smi -o output.arff
#                                          -c double -f remove_nonrings
#
#############################################################
import argparse
import os
import sys
import logging

from mmp.mmp_mcss_objects import MMPbasedMCSSObjectClass


#
# Main
#
def main():
    #
    parser = argparse.ArgumentParser(description='Generate MMP based MCSS from SMI input')
    parser.add_argument("-i", "--smi_file",
                        required=True,
                        help="Input smiles file")
    parser.add_argument("-o", "--out_file",
                        required=True,
                        help="Output pairs file CSV")
    parser.add_argument("-x", "--mdc_atm_hard",
                        type=int,
                        help="Max Double Cut Atom cutoff (Hard)"
                             "This implements a filter that will remove *any* double cuts (and therefore MCSS) where one "
                             "half has num_atoms <= mdc_atm_hard.  This is a hard cutoff and may therefore mean you will "
                             "fail to find any MCSS between two different compounds if the only valid MCSS is a double "
                             "cut with num_atoms <= mdc_atm_hard.  This is not a recommended option, use mdc_atm_soft ",
                        default=None)
    parser.add_argument("-s", "--mdc_atm_soft",
                        type=int,
                        help="Max Double Cut Atom cutoff (Soft)"
                             "* must be used with mdc_atm_soft_threshold *"
                             "When a double cut has a greater number of atoms than a single cut, the double cut will be "
                             "discarded in preference to the smaller number of atom single cut only when (a) and (b): "
                             "(a) either part of the double cut context has num_atoms <= mdc_atm_soft "
                             "(b) total double cut atom <= (single cut atoms + mdc_atm_soft_threshold) ",
                        default=None)
    parser.add_argument("-t", "--mdc_atm_soft_threshold",
                        type=int,
                        help="Threshold value for mdc_atm_soft"
                             "* must be used with mdc_atm_soft *"
                             "See help text for mdc_atm_soft",
                        default=None)
    parser.add_argument("-c", "--cut_type",
                        default="BOTH",
                        choices=['SINGLE', 'DOUBLE', 'BOTH'],
                        help="SINGLE|DOUBLE|BOTH (default = BOTH)")
    parser.add_argument("-f", "--filter",
                        default="NONE",
                        choices=['REMOVE_NONRINGS', 'REMOVE_RINGS', 'NONE'],
                        help="REMOVE_NONRINGS|REMOVE_RINGS|NONE, default = NONE")
    parser.add_argument("-l", "--log_file", dest="log_file", help="Name of file to place log info in")
    parser.add_argument("-d", "--log_level", dest="log_level", help="CRITICAL|ERROR|WARNING|INFO|DEBUG")

    opts = parser.parse_args()
    smi_fi = opts.smi_file
    out_fi = opts.out_file
    # minimum double cut (mdc) atoms == The minimum number of atoms allowed in any double cut fragment
    madc_hard = opts.mdc_atm_hard
    madc_soft = opts.mdc_atm_soft
    madc_soft_threshold = opts.mdc_atm_soft_threshold
    cut_type = opts.cut_type.upper()
    filter_type = opts.filter

    # should be a group?
    if (madc_soft is None) and (madc_soft_threshold is not None) or\
       (madc_soft is not None) and (madc_soft_threshold is None):
        parser.print_help()
        parser.error('Please specify both -s and -t together, you cannot use one or the other, consult --help\n')
        sys.exit()

    # get log details if they are defined or set to None
    log_file = opts.log_file
    log_level = opts.log_level

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
    my_mmp_mcss_object = MMPbasedMCSSObjectClass(mmplogger)

    # pre-check the smiles for for odd stuff
    # if prescan_smi:
    my_mmp_mcss_object.scan_input_smiles(smi_fi)

    mmplogger.info('Parse Dicer Fragments to Pairs')
    my_mmp_mcss_object.build_from_dicer(smi_fi, cut_type, filter_type)
    mmplogger.info('Get NAtoms for all frags')
    my_mmp_mcss_object.enumerate_fragment_properties()
    mmplogger.info('Write out final pairs')

    if madc_soft is not None:
        my_mmp_mcss_object.get_largest_mcs_pairs(out_fi, cut_type, mdc_atm_soft=madc_soft,
                                                 mdc_atm_soft_threshold=madc_soft_threshold)

    elif madc_hard is not None:
        my_mmp_mcss_object.get_largest_mcs_pairs(out_fi, cut_type, mdc_atm_hard=madc_hard)

    else:
        my_mmp_mcss_object.get_largest_mcs_pairs(out_fi, cut_type)

    mmplogger.info('Complete.')


if __name__ == '__main__':
    main()
    sys.exit()
