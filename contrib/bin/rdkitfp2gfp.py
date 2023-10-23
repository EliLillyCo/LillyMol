#!/usr/bin/env python

"""
    rdkitfp2gfp.py

    C3/Eli Lilly and Co

    Takes as input an rdkit generated fp file; returns an bit vector file

    Example:
    rdkitfp2gfp.sh -p 1000 -n 1 -pfn sample.gfp sample.fps

"""

from rdkit import DataStructs

import argparse
import string


def main():
    """Module mainline (for standalone execution)"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-pfn", "--productFileName", help="output file to place products (bit vector file)") 
    parser.add_argument("rdkitfps", help="input file containing hashed fingerprints generated using rdkit") 
    parser.add_argument("-gfp", "--gfpOutputFormat", help="write output in gfp format", action="store_true", default=0)
    parser.add_argument("-dlm", "--delimiter", help="delimiter/separator in input file", default='\t')
    parser.add_argument("-p", "--progressReport", help="report progress every number of lines", type=int, default=10000)
    parser.add_argument("-nid", "--numericID", help="assign id starting from this number", type=int, default=-99999999)
    parser.add_argument("-v", "--verbose", help="verbose mode", action="store_true", default=0)
    parser.add_argument("-d", "--debug", help="enable debug mode", action="store_true", default=0)

    # get the arguments
    args = parser.parse_args()

    # inform...
    if args.debug:
        print (args)
        print ("About to transform fingerprints in {}".format(args.rdkitfps))

    #
    transform(args.rdkitfps, args.productFileName, args.gfpOutputFormat,
              args.delimiter, args.numericID, args.progressReport,
              verbose=args.verbose, debug=args.debug)

    return


# read-transform-write numlines by numlines
def transform(infile, outfile, outformat, delimiter, numericID, numlines, verbose, debug):

    ins = open(infile)
    outs = open(outfile, 'w')

    # inform
    if debug or verbose:
        print ("About to read rdkit fingerprint file {} and transform to bivector file ... {}".format(infile, outfile))

    # run 
    nid = numericID
    ctr = 0
    # vlen = 0
    for line in ins:
        line.replace(r'\r', '')
        # if comment line ignore
        if line[0] == '#':
            if debug:
                print (line)
            continue

        # split line
        chunky_line = string.split(line, delimiter)

        # get rid of the newline character if it exists
        if chunky_line[-1][-1] == '\n' or chunky_line[-1][-1] == '^M':
            chunky_line[-1] = chunky_line[-1][:-1]

        # get the fp needed
        fp = DataStructs.CreateFromFPSText(chunky_line[0])
        bs = fp.ToBitString()

        if nid != -99999999:
            pcn = nid
            nid += 1
        else:
            try:
                pcn = chunky_line[1]
            except:
                #raise IOError("Line " + str(ctr) + " missing id element...")
                print ("Line " + str(ctr) + " missing id element...")
                continue

        # write
        if outformat:
            outs.write("$SMI<C>\n")
            outs.write("PCN<"+str(pcn)+">\n")
            outs.write("FPRDK<"+bs+">\n")
            outs.write("|\n")
        else:
            outs.write(bs + " " + str(pcn) + "\n")

        if debug:
            print (len(bs))

        if verbose and ctr % numlines == 0:
            print ("Processed {} fingerprints...".format(ctr))

        ctr += 1

    # inform
    if debug or verbose:
        print ("\nRead/transformed/wrote {} fingerprints.".format(ctr))

    ins.close()
    
    return


#
#
#
if __name__ == "__main__":
    main()
