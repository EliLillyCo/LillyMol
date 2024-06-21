#!/usr/bin/env python
###################################################################
# Summary: Script to generate matched moleular pairs from input smiles
#
#############################################################
import argparse
import os
import sys
import logging

from mmp.mmp_pairs_objects import MMPPairsObjectClass, validate_agg_method


def create_parser():

    parser = argparse.ArgumentParser(description='Generate matched molecular pairs and summarised Transforms from CSV')
    parser.add_argument("-i",
                        required=True,
                        help="Input CSV filename")
    parser.add_argument("-o",
                        help="Output filename to write pairs and stats data")
    parser.add_argument("-s", "--smiles",
                        required=True,
                        help="Column name for Molecule SMILES")
    parser.add_argument("-n", "--molid",
                        required=True,
                        help="Column name for Molecule ID")
    parser.add_argument("-d", "--diff_calc",
                        default='DIFF',
                        choices=['DIFF', 'RATIO'],
                        help="Difference calc DIFF|RATIO (DIFF = FRAG_R - FRAG_L / RATIO = FRAG_R / FRAG_L), default = DIFF")
    parser.add_argument("-c", "--cut_type",
                        default="BOTH",
                        choices=['SINGLE', 'DOUBLE', 'BOTH'],
                        help="SINGLE|DOUBLE|BOTH (default = BOTH)")
    parser.add_argument("-f", "--filter",
                        default="NONE",
                        choices=['REMOVE_NONRINGS', 'REMOVE_RINGS', 'NONE'],
                        help="REMOVE_NONRINGS|REMOVE_RINGS|NONE, default = NONE",)
    parser.add_argument("-l", "--log_file",
                        help="Name of file to place log info in")
    parser.add_argument("-L", "--log_level",
                        help="CRITICAL|ERROR|WARNING|INFO|DEBUG")
    parser.add_argument("-k", "--skip_col_check",
                        help="Skip auto detect of numeric columns. Code will attempt to get diff data on all/any column.",
                        action="store_true",
                        default=False)
    parser.add_argument("-t", "--types_header",
                        help="Invoke this flag (no value needed) if you want a second header line with datatypes",
                        action="store_true",
                        default=False)
    parser.add_argument("-a", "--act_data_col",
                        help="Optional flag to specify column to calculate differences on, else will do all columns")
    parser.add_argument("-A", "--agg_method",
                        default=False,
                        help="Invoke flag with value MEAN_DIFF|MEAN_DIFF_INVLOG|CATEGORICAL|DIFF50|DIFFXX to get summarised"
                             " frag pairs, also needs -a. MEAN simply calculates the mean of the input data (after "
                             " aggregation). MEAN_DIFF calculates MEAN, confidence limits are added, input should already "
                             " be in the Log scale (e.g. LogD). MEAN_DIFF_INVLOG also inverts the resulting mean out of log"
                             " scale to get FOLD_CHANGE (e.g.: Solubility, Clearance). CATEGORICAL expects binary 0/1 two "
                             " category data (e.g.: Permeation). DIFFXX where 0 < X < 100 implements an index function, for"
                             " use with CYP Pc Inhibition & metabolism data.")
    parser.add_argument("-g", "--agg_type",
                        default="FRAG",
                        choices=['FRAG', 'ATTACH'],
                        help="Invoke flag with FRAG|ATTACH to exclude|include atom attachment points in aggregation, "
                             "default = FRAG")
    parser.add_argument("-p", "--prop_data",
                        default='NONE',
                        choices=['BASIC', 'ALL'],
                        help="Invoke flag with BASIC|ALL to get MWT/CLP (BASIC) or full iwdecr property data (ALL) diffs")
    parser.add_argument("-q", "--add_quality",
                        dest="add_quality",
                        help="Invoke flag if you want to add Quality metric to summarised pairs",
                        action="store_true", default=False)
    parser.add_argument("-x", "--del_id_dupes",
                        help="Invoke flag if you want to keep and count duplicate positional isomer fragL->fragR changes"
                             "for same ID",
                        action="store_false",
                        default=True)
    parser.add_argument("-y", "--get_low_n_pairs",
                        help="Invoke flag if you want to get a summarised data point for pairs with N<3 with DIFFxx and "
                             "MEAN_DIFF Pc calcs. Other method always return values for low N pairs. This flag impacts speed, "
                             "particularly for DIFFxx due to calculation of confidence intervals.  At a rough estimate the "
                             "execution time of the script will increase by <=x Pc where x is equal to the number of summarised "
                             "pairs with n<3. When 50 Pc of summarised pairs have n<3 script is <=50 Pc slower.",
                        action="store_true",
                        default=False)

    args = parser.parse_args()

    # rely on module level function to validate this option
    if args.agg_method:
        args.agg_method = args.agg_method.upper()  # validate this below using module level function
        validate_agg_method(args.agg_method)

    # check -A and -a flags
    # check this behaviour early on rather than half way through the code as the pairs
    # aggregations code on the imported object will replicate this behaviour
    if args.agg_method and args.act_data_col is None:
            parser.print_help()
            parser.error('Can only invoke the -A when you specify the activity column using -a (because it can be slow!)')
            sys.exit()

    if (args.add_quality is True and args.prop_data.upper() == 'NONE') or \
            (args.add_quality is True and args.agg_type is False):
        parser.print_help()
        parser.error(
            'Can only generate Quality metric when using -A (aggregation) and -p (properties), consult --help\n')
        sys.exit()

    if args.agg_type:
        args.agg_type = args.agg_type.upper()

    return parser

#
# Main
#
def main():

    parser = create_parser()
    args = parser.parse_args()

    if str(args.prop_data).upper() == 'BASIC':
        add_prop_diff = 1
    elif str(args.prop_data).upper() == 'ALL':
        add_prop_diff = 2
    else:
        add_prop_diff = 0

    if args.diff_calc.upper() == 'RATIO':
        fold_diff = True
    else:
        fold_diff = False

    mmplogger = logging.getLogger('lillymol_file_logger')

    if args.log_file is not None:

        fh = logging.FileHandler(args.log_file)

        if args.log_level is None:
            mmplogger.setLevel(logging.INFO)
        else:
            numeric_level = getattr(logging, args.log_level.upper(), None)
            if not isinstance(numeric_level, int):
                raise ValueError('Invalid log level: %s' % args.log_level)
            mmplogger.setLevel(numeric_level)

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(module)s: %(message)s',
                                      datefmt='%m/%d/%Y %I:%M:%S %p')

        fh.setFormatter(formatter)
        mmplogger.addHandler(fh)

    else:
        logging.disable(logging.CRITICAL)
    # instantiate objects
    mmplogger.info('Instantiating the class')
    my_mmp_data_object = MMPPairsObjectClass(mmplogger)

    # parse CSV
    mmplogger.info('Will now parse CSV into the object')
    
    # my_mmp_data_object.csv_to_pair_and_data_objects(args.i, args.smiles, args.molid, args.cut_type.upper(), args.filter)
    if args.skip_col_check:
        my_mmp_data_object.csv_sniffer(args.i, args.smiles, args.molid, skip_num_check=True)
    else:
        my_mmp_data_object.csv_sniffer(args.i, args.smiles, args.molid)
    
    if args.act_data_col is None and add_prop_diff == 0:
        my_mmp_data_object.csv_to_data_objects(args.i, args.smiles, args.molid)
    
    elif args.act_data_col is None and add_prop_diff > 0:
        my_mmp_data_object.csv_to_data_objects(args.i, args.smiles, args.molid, add_prop_diff=add_prop_diff)
    
    elif args.act_data_col is not None and add_prop_diff == 0:
        my_mmp_data_object.csv_to_data_objects(args.i, args.smiles, args.molid, act_data_col=args.act_data_col)
    
    else:
        my_mmp_data_object.csv_to_data_objects(args.i, args.smiles, args.molid, act_data_col=args.act_data_col,
                                               add_prop_diff=add_prop_diff)

    tmp_dicer_file = my_mmp_data_object.write_mol_smi_dict_tofile()
    my_mmp_data_object.build_from_dicer(tmp_dicer_file, args.cut_type.upper(), args.filter)

    # get pairs
    mmplogger.info('Now getting pairs and writing to file: %s' % args.o)
    if args.types_header is True:
        my_mmp_data_object.get_pairs_and_diffs(args.o, args.cut_type.upper(),
                                               fold_diff=fold_diff,
                                               inc_types_header=args.types_header,
                                               add_qmetric=args.add_quality)
    
    else:
        my_mmp_data_object.get_pairs_and_diffs(args.o, args.cut_type.upper(),
                                               fold_diff=fold_diff,
                                               add_qmetric=args.add_quality)
    
    mmplogger.info('Completed generation of raw pairs.')

    # Now do pair aggregation if requested, consider
    #  Essential: agg_type / agg_method
    #  Optional and interdependent: prop_data (default False) / args.act_data_col (default False)
    #
    if args.agg_method:
   
        #########################################################
        # free up some memory
        # - Running on HPC environment, with 66,514 smi -> 14,514,719 pairs I get various MemoryError returns
        #   from python/pandas specifically as we do the aggregation so need to free up some space
        #########################################################
        # This should empty the object out, we'll then read back from CSV
        mmplogger.info('Cleaning out all objects from memory')
        my_mmp_data_object.clean_out_data()

        mmplogger.info('Now getting summarised/aggregated pairs')
        
        mmplogger.debug('Generating pandas dataframe from the CSV file we just wrote')
        # seems like this is very slow
        # my_mmp_data_object.pairsdataobj_to_pd(cut_type, args.act_data_col)
        # try this instead
        my_mmp_data_object.pd_read_csv(args.o)

        mmplogger.debug('Done reading CSV, now aggregating')
        
        if args.act_data_col is None:

            if add_prop_diff > 0:
                my_mmp_data_object.pd_aggregate_pairs_to_csv(args.o+'.sum', args.agg_type,
                                                             prop_data=True,
                                                             agg_method=args.agg_method,
                                                             remove_id_dupes=args.del_id_dupes,
                                                             inc_low_n_vals=args.get_low_n_pairs,
                                                             add_qmetric=args.add_quality)
            else:
                my_mmp_data_object.pd_aggregate_pairs_to_csv(args.o+'.sum', args.agg_type,
                                                             agg_method=args.agg_method,
                                                             remove_id_dupes=args.del_id_dupes,
                                                             inc_low_n_vals=args.get_low_n_pairs,
                                                             add_qmetric=args.add_quality)
        else:
            if add_prop_diff > 0:
                my_mmp_data_object.pd_aggregate_pairs_to_csv(args.o+'.sum', args.agg_type,
                                                             prop_data=True,
                                                             agg_method=args.agg_method,
                                                             act_col=args.act_data_col,
                                                             remove_id_dupes=args.del_id_dupes,
                                                             inc_low_n_vals=args.get_low_n_pairs,
                                                             add_qmetric=args.add_quality)
            else:
                my_mmp_data_object.pd_aggregate_pairs_to_csv(args.o+'.sum', args.agg_type,
                                                             agg_method=args.agg_method,
                                                             act_col=args.act_data_col,
                                                             remove_id_dupes=args.del_id_dupes,
                                                             inc_low_n_vals=args.get_low_n_pairs,
                                                             add_qmetric=args.add_quality)
        
        mmplogger.info('Completed generation of summarised pairs.')


#
#
if __name__ == '__main__':

    main()
    sys.exit()
