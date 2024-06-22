#!/usr/bin/env python
##############################################################
"""
Summary: Input SMI to script.  SMI will be fragmented using DICER
and the resulting R-groups/fragments will be used to query the MMPs
for MMP based improvements/transformations for command line specified
ADME end points (JAL)

"""
#############################################################
import logging
import argparse
import sys
import os

from mmp.mmp_enum_mols_from_pairs import MMPEnumerateNewMols


def main():

    def restricted_dicer_float(x):
        """check maxff switch on dicer is valid, where 0.3 is smaller
        and 0.5 is default and 0.7 is huge/avoided"""
        x = float(x)
        if x < 0.01 or x > 0.99:
            raise argparse.ArgumentTypeError("%r not in range [0.01, 0.99]" % (x,))
        return x

    parser = argparse.ArgumentParser(description='Enumerate new molecules from input SMI and a matched pairs file. ')

    parser.add_argument("-i", nargs=1, required=True,
                        help="Specify input smiles file")
    parser.add_argument("-p", nargs=1, required=True,
                       help="Specify input pairs file")
    parser.add_argument("--frag_left_col", nargs=1, required=True,
                        help="Specify name of frag left column in input pairs file")
    parser.add_argument("--frag_right_col", nargs=1, required=True,
                        help="Specify name of frag right column in input pairs file")
    parser.add_argument("-o", nargs=1, required=True,
                        help="Specify output CSV filename")
    parser.add_argument('-H', action='store_false', default=True,
                        help='Invoke flag to include all H transformations (usually lots!)')
    parser.add_argument('-b', type=restricted_dicer_float, default=0.3,
                        help='Dicer Threshold for Fragmentation, default=0.3 for small side chain replace or set to '
                             '0.50001 for larger fragment removal including core replacements')
    parser.add_argument('-l', nargs=1, default=False, help='optional log filename')
    parser.add_argument('-L', nargs=1, default='INFO', help='debug log level CRITICAL|ERROR|WARNING|INFO|DEBUG')

    args = parser.parse_args()

    # create logger
    mmplogger = logging.getLogger('lillymol_file_logger')

    log_level = args.L[0]

    if args.l:
        log_file = args.l[0]

    else:
        log_file = None
        log_level = None

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
                                      datefmt='%Y/%m/%d %I:%M:%S %p')

        fh.setFormatter(formatter)

        mmplogger.addHandler(fh)

    else:
        logging.disable(logging.CRITICAL)

    #######################
    # main stuff
    #######################

    # print('Instantiating MMP ADME DB Objects'
    mmplogger.info('Instantiating Enumeration Object')
    mmp_newmol_object = MMPEnumerateNewMols(mmplogger)

    # read in the input smiles then fragment
    mmp_newmol_object.scan_input_smiles(args.i[0], injest=True)
    # fragment
    mmp_newmol_object.fragment_reference_smi('BOTH', args.b, exclude_h_subst=args.H)

    # get transformations for enumeration
    temp_dict = mmp_newmol_object.pairs_file_to_dict(args.p[0],
                                                     frag_l_col=args.frag_left_col[0],
                                                     frag_r_col=args.frag_right_col[0],
                                                     exclude_h_subst=args.H)

    mmp_newmol_object.add_transformation_group('pairs_from_file', temp_dict)

    # enumerate
    mmp_newmol_object.write_rxn_files()
    mmp_newmol_object.write_reactants_mol_frag_dict()
    mmp_newmol_object.do_reactions()
    mmp_newmol_object.write_products_to_csv(args.o[0])


if __name__ == '__main__':
    main()
