#!/usr/bin/env python
##############################################################################
#
# Summary: 
# Generate all Matched Series from input SMI file and prints to output
#
# Example usage:
#  python compare_mmp_sum_files.py
#   -i input1.pairs.sum input2.pairs.sum 
#   -o output_diffs.csv
#
##############################################################################
import argparse 
import sys
import logging

from mmp.mmp_series_object import MMPSeriesObjectClass


def main():

    def restricted_float(x):
        x = float(x)
        if x < 0.1 or x > 0.99:
            raise argparse.ArgumentTypeError("%r not in range [0.1, 0.99]" % (x,))
        return x

    parser = argparse.ArgumentParser(description='Generate Matched Series data from an input CSV. '
                                                 'Script prints all valid series but does not exhaustively enumerate '
                                                 'all sub-series of a series as this is too combinatorially explosive. '
                                                 'Use output with script getMMPSeriesSuggestfromCSV to utilise data. '
                                                 '-m / -f flags sets as per Keefer and Chang MedChemComm 2017 where a '
                                                 'series must have >=5 compounds, activity spread >= 0.5 and skew <= 3')
    parser.add_argument('-i', nargs=1, required=True, help='Input CSV file of smiles and activity data')
    parser.add_argument('-o', nargs=1, required=True, help='Output filename to write series to')
    parser.add_argument('-s', nargs=1, required=True, help='Column name for Molecule SMILES')
    parser.add_argument('-n', nargs=1, required=True, help='Column name for Molecule ID')
    parser.add_argument('-a', nargs=1, required=True, help='Column name for Activity Data')
    parser.add_argument('-m', type=int, default=5,
                        help="Invoke to set a min limit on the reported series length (default >= 5). "
                             "The minimum allowed value is 3 as a value of 2 will simple find matched pairs")
    parser.add_argument('-t', type=restricted_float, default=0.50001,
                        help="Set the percentage of the molecule to be retained as context. Use smaller values (0.3) "
                             "for small side chain replacements only, default 0.50001")
    parser.add_argument('-c', type=str, default='BOTH', choices=['SINGLE', 'DOUBLE', 'BOTH'])
    parser.add_argument('-l', nargs=1, help='optional log filename', default=False)
    parser.add_argument('-L', nargs=1, choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'], default='INFO')

    args = parser.parse_args()

    input_file = args.i[0]
    smi_col = args.s[0]
    id_col = args.n[0]
    act_col = args.a[0]

    if args.m <= 2:
        parser.print_help()
        sys.exit('Please ensure your min series length (-m) value is greater than 2')

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
                                      datefmt='%m/%d/%Y %I:%M:%S %p')

        fh.setFormatter(formatter)
        
        mmplogger.addHandler(fh)

    else:
        logging.disable(logging.CRITICAL)

    #######################
    # main stuff
    #######################

    mmplogger.info('Instantiating MMP Series Objects')
    mmp_series_object = MMPSeriesObjectClass(mmplogger)

    # setup
    mmp_series_object.setup_mmp_data_for_mms(input_file, smi_col, id_col, act_col, args.m, args.t, cut_type=args.c)

    # write raw series
    mmp_series_object.write_raw_series_to_file(csv_out=args.o[0])


if __name__ == "__main__":

    main()
