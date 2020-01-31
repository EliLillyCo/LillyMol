###################################################################
"""
Summary: Enumerate new molecules from matched pairs

About: Expects a smiles file of input that it will fragment using dicer into the
separate parts: context + frag_l.  It also expects a matched pairs summary file with
columns frag_l and frag_r as minimum.  Script will iterate over all possible
fragmentation patterns of every input SMI (there is usually more than one per smi) and
by searching the matched pairs summary file for matching frag_l it will then use trxn
to switch frag_l on original molecule with a matched pair frag_r, fully
enumerated new molecules from the matched pairs. (JAL)

"""########################################################

import logging
import sys
import csv
# import re
import tempfile
import os
import unittest

from mmp.mmp_dicer_functions import execute_dicer

if 'LILLYMOL_HOME' in os.environ:
    import pybase.pyopmo as pymo
else:
    import pybase.pymo as pymo


class MMPEnumerateNewMols(object):
    """Object and methods for querying a set of MMPs matching input mol for transfomation
    Instantiation of the object requires a valid python logger object to be
    passed in as a parameter, even if the logger is switched off.

    Example usage:
        mmplogger = logging.getLogger('lillymol_file_logger')
        logging.disable(logging.CRITICAL)
        mmp_admedb_object = MMPQueryADMEDB(mmplogger)

    """
    def __init__(self, logger_object):

        self.logger = logger_object
        if len(logging.Logger.manager.loggerDict) < 1:
            # exit with system status 1 and custom error
            sys.exit("Invalid or no logger object passed to MMPObjectClass.  Please create \
                    and pass a logger and set to use logging.disable if you don't want logging")

        # files for single
        self.rxn_file_ctx_sgl = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_frags_sgl = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_rxn_map_sgl = tempfile.NamedTemporaryFile(delete=False, suffix='.rxn', encoding='utf-8', mode='wt')
        self.rxn_file_makethese_sgl = tempfile.NamedTemporaryFile(delete=False, suffix='.txt', encoding='utf-8', mode='wt')
        self.rxn_file_output_sgl = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        # files for double
        self.rxn_file_ctx_dbl1 = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_ctx_dbl2 = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_frags1_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_frags2_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.smi', encoding='utf-8', mode='wt')
        self.rxn_file_rxn_map1_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.rxn', encoding='utf-8', mode='wt')
        self.rxn_file_rxn_map2_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.rxn', encoding='utf-8', mode='wt')
        self.rxn_file_rxn_map3_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.rxn', encoding='utf-8', mode='wt')
        self.rxn_file_rxn_map4_dbl = tempfile.NamedTemporaryFile(delete=False, suffix='.rxn', encoding='utf-8', mode='wt')
        self.rxn_file_makethese_dbl1 = tempfile.NamedTemporaryFile(delete=False, suffix='.txt', encoding='utf-8', mode='wt')
        self.rxn_file_makethese_dbl2 = tempfile.NamedTemporaryFile(delete=False, suffix='.txt', encoding='utf-8', mode='wt')
        self.rxn_file_output_dbl1 = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        self.rxn_file_output_dbl2 = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        #
        self.do_single = False
        self.do_double = False
        #
        self.single_exit_code1 = None
        self.double_exit_code1 = None
        self.double_exit_code2 = None
        #
        self.smi_file = None
        self.csv_file = None
        # storage
        self.mol_smi_dict = {}
        self.mol_fragments_dict_single = {}
        self.mol_fragments_dict_double = {}
        self.unique_frags = set()
        self.transformation_groups = {}
        self.new_mols = {}

        # this is the reaction specification that trxn needs to combine isotopically labelled mmp fragmentation points
        self.single_rxn = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n )\n"
        self.single_rxn += " (1 Sidechain\n  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n  (A I join (0 0))\n )\n)"

        # double cut needs two types of reaction to deal with the special case [12 isotopic label
        self.double_rxn_1 = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[1*].[2*]\")\n )\n (1 Sidechain\n"
        self.double_rxn_1 += "  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n  (A I join (0 0))\n )\n)"

        self.double_rxn_2 = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[1*].[2*]\")\n  (A I isotope (0 0))\n"
        self.double_rxn_2 += "  (A I isotope (1 0))\n )\n (2 Sidechain\n  (A C smarts \"[!0*]\")\n"
        self.double_rxn_2 += "  (A I isotope (0 0))\n  (A I join (1 0))\n )\n)"

        self.double_rxn_3 = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[12*]\")\n )\n (1 Sidechain\n"
        self.double_rxn_3 += "  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n  (A I join (0 0))\n )\n)"

        self.double_rxn_4 = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[12*]\")\n  (A I isotope (0 0))\n"
        self.double_rxn_4 += " )\n (2 Sidechain\n  (A C smarts \"[!0*]\")\n"
        self.double_rxn_4 += "  (A I isotope (0 0))\n  (A I join (0 0))\n )\n)"

    def clean_out_data_mmpenumerate(self):
        """Method to clean out all objects in class"""

        #
        self.smi_file = None
        self.csv_file = None
        # storage
        self.mol_smi_dict = {}
        self.mol_fragments_dict_single = {}
        self.mol_fragments_dict_double = {}
        self.unique_frags = set()
        self.transformation_groups = {}
        self.new_mols = {}

    def file_exists_check(self, file_to_check):
        """        """
        try:
            if os.path.isfile(file_to_check) > 0:
                return True
            else:
                return False
        except OSError:
            return False

    def set_do_single(self):
        """Method to set self.do_single.  Used when enumerating from Pairs file or methods that don't directly
        confirm single/double fragmentation"""
        self.do_single = True

    def set_do_double(self):
        """Method to set self.do_single.  Used when enumerating from Pairs file or methods that don't directly
        confirm single/double fragmentation"""
        self.do_double = True

    def scan_input_smiles(self, smi_fi, injest=False):
        """
        TODO : Method ripped from mmp_object, only ingest line differs, could extract to utils function
        Dicer currently finishes with an exit code of zero, even if it fails and bombs half way through. This means the
        MMP code can complete with results but has missed many pairs due to failed fragmentation, this injects a
        pre-check into the process to remove/drop smiles likely to fail dicer BEFORE they get parsed to dicer or just
        crash out with an error code.
        """

        self.smi_file = smi_fi

        with open(smi_fi, "rt") as input_smifi:

            self.logger.info('Beginning pre-scan of input smiles file')

            line_num = 0
            for line in input_smifi:

                line_num += 1
                line_list = line.split()

                # check it's valid smiles file format which is "smiles\sid"
                if len(line_list) != 2:
                    # random exit code of 4
                    sys.exit("Failed to parse smiles from line %s of file %s" % (line_num, smi_fi))

                # Is the smiles_id numeric?
                try:
                    smiles_id = int(line_list[1])
                except:
                    sys.exit("Failed smiles parser as non-numeric id on line %s of file %s" % (line_num, smi_fi))

                smiles = line_list[0]

                if len(smiles) < 2:
                    sys.exit("Failed smiles parser as tiny smi size on line %s of file %s" % (line_num, smi_fi))

                if "." in smiles:
                    sys.exit("Failed smiles parser as salt on line %s of file %s (pls clean using fileconv)" % (
                        line_num, smi_fi))

                if injest:
                    # add to the mol_smi_dict as mol_id => smiles
                    self.mol_smi_dict[smiles_id] = line_list[0]

            self.logger.info('Completed pre-scan of input smiles file with no errors')

        input_smifi.close()

    def fragment_reference_smi(self, cut_type, dicer_threshold, exclude_h_subst=True):
        """ Method to fragment using dicer, all smiles we've read into self.original_mols
        :dicer_threshold = The dicer threshold value to use during fragmentation (maxff)
        """

        # using key word "DOUBLE" will execute single and double cut in dicer so need to weed out
        # the double cuts or we'll get the same as "BOTH"
        parse_single = False
        parse_double = False
        if cut_type == "SINGLE" or cut_type == "BOTH":
            parse_single = True
        if cut_type == "DOUBLE" or cut_type == "BOTH":
            parse_double = True

        if self.smi_file is None:
            sys.exit('Please use method scan_input_smiles to set input smiles file')

        # execute dicer to get fragments of key_smi
        for num_cuts, mol_id, ctx_orig, frag, fattach, cattach in execute_dicer(self.smi_file, cut_type,
                                                                                dicer_threshold, self.logger):

            if num_cuts == 1:
                if parse_single:
                    if mol_id in self.mol_fragments_dict_single:
                        self.mol_fragments_dict_single[mol_id][ctx_orig] = frag
                    else:
                        self.mol_fragments_dict_single[mol_id] = {}
                        # it's possible we'll get duplicate fragmentation patterns e.g.: CF3 but just overwrite
                        self.mol_fragments_dict_single[mol_id][ctx_orig] = frag

            elif num_cuts == 2:
                if parse_double:
                    if mol_id in self.mol_fragments_dict_double:
                        self.mol_fragments_dict_double[mol_id][ctx_orig] = frag
                    else:
                        self.mol_fragments_dict_double[mol_id] = {}
                        # it's possible we'll get duplicate fragmentation patterns e.g.: CF3 but just overwrite
                        self.mol_fragments_dict_double[mol_id][ctx_orig] = frag

            else:
                sys.exit("Error: Dicer gave a cut type other than 1 and 2!! (method call fragment_reference_smi)")

        # now get a unique set of frags for querying the DB
        for mol_id in self.mol_fragments_dict_single:
            for ctx, frag in self.mol_fragments_dict_single[mol_id].items():
                self.unique_frags.add(frag)

        for mol_id in self.mol_fragments_dict_double:
            for ctx, frag in self.mol_fragments_dict_double[mol_id].items():
                self.unique_frags.add(frag)

        ###################################################################
        # This single line prevents all H -> substitutions from occurring #
        # if it gets removed, then edit the unit tests to include [1H]    #
        ###################################################################
        if exclude_h_subst is True:
            if '[1H]' in self.unique_frags:
                self.unique_frags.remove('[1H]')

        self.logger.debug("Found %s unique fragments for search" % len(self.unique_frags))

    def pairs_file_to_dict(self, pairs_file,
                           frag_l_col='FRAG_L', frag_r_col='FRAG_R',
                           exclude_h_subst=True, num_pairs_limit=999999):
        """
        :param pairs_file file to read
        :param frag_l_col='FRAG_L' column name of left hand fragment column
        :param frag_r_col='FRAG_R' column name of right hand fragment column
        :param exclude_h_subst Default False means we will not enumerate H -> R changes
        :param num_pairs_limit Parses only the first n lines of the file, use case is get Top 100 from
        most common Lilly transforms
        """
        # I will return this as frag_L: {frag_r: [data], frag_r: [data]}
        temp_dict = {}

        sniffer = csv.Sniffer()
        if sniffer.has_header(pairs_file) == 0:
            sys.exit("Sniffer module did not detect a header in your CSV file, please check!")

        pairs_count = 0
        with open(pairs_file, 'rt') as csvfile:

            csv_dialect = sniffer.sniff(csvfile.read(65536))
            csvfile.seek(0)
            reader = csv.reader(csvfile, csv_dialect, skipinitialspace=True)

            headers = next(reader)
            self.logger.info("Got headers from CSV: %s" % headers)

            csv_items_per_line = len(headers)
            self.logger.info("Expected number of items per line: %s" % csv_items_per_line)

            if frag_l_col not in headers:
                self.logger.warn('Cannot find column FRAG_L in CSV file')
                sys.exit('Cannot find column FRAG_L in CSV file')
            else:
                header_fragl_num = headers.index(frag_l_col)

            if frag_r_col not in headers:
                self.logger.warn('Cannot find column FRAG_R in CSV file')
                sys.exit('Cannot find column FRAG_R in CSV file')
            else:
                header_fragr_num = headers.index(frag_r_col)

            self.logger.info('Column headers look great, proceed with parse')

            for row in reader:

                pairs_count += 1
                if pairs_count <= num_pairs_limit:

                    frag_l = row[header_fragl_num]
                    frag_r = row[header_fragr_num]

                    if frag_l in temp_dict:
                        temp_dict[frag_l][frag_r] = row
                    else:
                        temp_dict[frag_l] = {}
                        temp_dict[frag_l][frag_r] = row

                else:
                    self.logger.info('Early Termination parsing %d pairs from file as n limit set' % pairs_count)
                    break

        self.logger.info('Done reading from csv')

        if exclude_h_subst is True:
            if '[1H]' in temp_dict:
                del(temp_dict['[1H]'])

        return temp_dict

    def add_transformation_group(self, name, pairs_dict):
        """ """

        if name in self.transformation_groups:
            sys.exit('Error, you tried to add the same transformation group twice')
        else:
            self.transformation_groups[name] = pairs_dict

    # write the rxn files we need to disk
    def write_rxn_files(self, rxn_type='both'):
        """pre-defined rxn is written to file ready for trxn"""

        if rxn_type is 'single' or rxn_type is 'both':
            self.rxn_file_rxn_map_sgl.write(self.single_rxn)
            self.rxn_file_rxn_map_sgl.close()

        if rxn_type is 'double' or rxn_type is 'both':
            self.rxn_file_rxn_map1_dbl.write(self.double_rxn_1)
            self.rxn_file_rxn_map1_dbl.close()
            self.rxn_file_rxn_map2_dbl.write(self.double_rxn_2)
            self.rxn_file_rxn_map2_dbl.close()

            self.rxn_file_rxn_map3_dbl.write(self.double_rxn_3)
            self.rxn_file_rxn_map3_dbl.close()
            self.rxn_file_rxn_map4_dbl.write(self.double_rxn_4)
            self.rxn_file_rxn_map4_dbl.close()

    def write_reactants_mol_frag_dict(self):
        """Take the various context we got from dicer for each molecule
        (found in self.mol_fragments_dict_single[mol_id]) then see if the
        matching frag_l has alternates for a given transformation group (which
        will be found in self.transformation_groups[name]).  If so, write
        out for enumeration
        """

        random_id = 0
        confirm_write = False

        for mol_id in self.mol_fragments_dict_single:

            for ctx, frag_l in self.mol_fragments_dict_single[mol_id].items():

                for name in self.transformation_groups:

                    if frag_l in self.transformation_groups[name]:

                        random_id += 1
                        # delimit the id & frag_l using '___' as it's not valid smiles
                        # we will split this back again on results read to get original id and frag_l
                        mol_seq_id = str(mol_id) + "___" + frag_l + "___" + str(random_id)

                        # write 1 - context to add new frag to
                        self.rxn_file_ctx_sgl.write(ctx + " " + mol_seq_id + "\n")

                        for frag_r in self.transformation_groups[name][frag_l]:

                            # at this point only are we completely sure we have a rxn for a given single or double
                            # cut, therefore set the vars so we know we should read back / expect results
                            confirm_write = True
                            # write 2 - new fragment to be added to ctx
                            self.rxn_file_frags_sgl.write(frag_r + " " + frag_r + "\n")
                            # write 3 - add to the full 'make' list
                            self.rxn_file_makethese_sgl.write(mol_seq_id + " " + frag_r + "\n")

        self.do_single = confirm_write
        self.rxn_file_ctx_sgl.close()
        self.rxn_file_frags_sgl.close()
        self.rxn_file_makethese_sgl.close()

        # now repeat for double but this time we write to three files not two
        confirm_write = False

        # need this to help find and convert smiles labelled with [12 to duplicate but singly labelled smi
        # regex = r"(.*)\[12(.{1,4})\](.*)"

        for mol_id in self.mol_fragments_dict_double:

            for ctx, frag_l in self.mol_fragments_dict_double[mol_id].items():

                ctx1 = ctx.split('.')[0]
                ctx2 = ctx.split('.')[1]

                for name in self.transformation_groups:

                    if frag_l in self.transformation_groups[name]:

                        # write 1 - ctx has to be split, easier to make these the frags
                        self.rxn_file_frags1_dbl.write(ctx1 + " " + ctx1 + "\n")
                        self.rxn_file_frags2_dbl.write(ctx2 + " " + ctx2 + "\n")

                        for frag_r in self.transformation_groups[name][frag_l]:

                            # delimit the id & frag_l using '___' as it's not valid smiles
                            # we will split this back again on results read to get original id and frag_l
                            mol_seq_id = str(mol_id) + "___" + frag_l + "___" + frag_r

                            # print frag_r
                            # at this point only are we completely sure we have a rxn for a given single or double
                            # cut, therefore set the vars so we know we should read back / expect results
                            confirm_write = True

                            if "[12" in frag_r:
                                # match = re.search(regex, frag_r)
                                # print match.group(1)
                                # frag_r_before = frag_r
                                # frag_r = match.group(1) + "[1" + match.group(2) + "].[2" + match.group(2) + "]" +
                                #     match.group(3)
                                # print ("Edited frag from %s to %s" % (frag_r_before, frag_r))
                                self.rxn_file_ctx_dbl2.write(frag_r + " " + mol_seq_id + "\n")
                                self.rxn_file_makethese_dbl2.write(mol_seq_id + " " + ctx1 + " " + ctx2 + "\n")

                            else:
                                # print ("Edited frag to %s" % (frag_r))
                                self.rxn_file_ctx_dbl1.write(frag_r + " " + mol_seq_id + "\n")
                                self.rxn_file_makethese_dbl1.write(mol_seq_id + " " + ctx1 + " " + ctx2 + "\n")

        self.do_double = confirm_write
        # Yes I know I used 'with' but I still had to close these files to avoid odd unittest bug
        self.rxn_file_ctx_dbl1.close()
        self.rxn_file_ctx_dbl2.close()
        self.rxn_file_frags1_dbl.close()
        self.rxn_file_frags2_dbl.close()
        self.rxn_file_makethese_dbl1.close()
        self.rxn_file_makethese_dbl2.close()

    def write_reactants_simple_dict(self, simple_dict):
        """Take the various context + fragment pairs we have in a pandas dataframe and write out to file
        for enumeration.  Writes to same source as write_reactants_mol_frag_dict so do not use this
        method at the same time as any other method as files will get overwritten.The assumption is that
        all ctx / frag pairs in the dataframe already have attachment points within the smiles
        """
        random_id = 0

        if len(simple_dict) < 1:
            raise Exception("Error, no items in dict so cant enumerate anything")

        for (ctx_smi, frag_smi) in list(simple_dict.keys()):

            # delimit the id & frag_l using '___' as it's not valid smiles
            # we will split this back again on results read to get original id and frag_l
            random_id += 1
            mol_seq_id1 = str(random_id) + "___" + ctx_smi
            random_id += 1
            mol_seq_id2 = str(random_id) + "___" + frag_smi

            # write 1 - context to add new frag to
            self.rxn_file_ctx_sgl.write(ctx_smi + " " + mol_seq_id1 + "\n")
            # write 2 - new fragment to be added to ctx
            self.rxn_file_frags_sgl.write(frag_smi + " " + mol_seq_id2 + "\n")
            # write 3 - add to the full 'make' list
            self.rxn_file_makethese_sgl.write(mol_seq_id1 + " " + mol_seq_id2 + "\n")

        self.do_single = True
        self.rxn_file_ctx_sgl.close()
        self.rxn_file_frags_sgl.close()
        self.rxn_file_makethese_sgl.close()

    def do_reactions(self):
        """
        Added function to pymo pymo.make_these_mols as Alternative to trxn, need this for the MMP code for generating
        new mols from MMP's, trxn version would be this (two components, single cut fragments):
            trxn.sh  -S - -r oneConnection.rxn partOne.smi partTwo.smi
        For connecting two bonds (three component, two cut fragments) it would be this:
            trxn.sh  -S -rxn -S - -r twoConnection.rxn partThree.smi  partOne.smi partTwo.smi
        BUT, if have a long list of different contexts (partOne) and don't want exhaustive enumeration, specify rxn's:
            make_these_molecules.sh -R oneConnection.rxn -M m2Make.txt -S - partOne.smi partTwo.smi
        In this case, you can put all context fragments SMILES (context1a, context 1b, ...)  in one reagent file, and
        all fragments SMILES (frag1, frag2, ...) in the second reagent file.  If have something like (context1a frag1\n
        context1a frag2\ncontext1b frag3\n...) in your m2Make.txt file, you will create the molecules you wanted
        """
        # pymo will strip file extension from output filename, trxn requires an stem not full name/suffix but will then
        # write to an output file {stem}.smi which makes temp file handling more complicated. This is exhaustive
        # exit_status = pymo.trxn([self.rxn_file_ctx_sgl.name, self.rxn_file_frags_sgl.name],
        #                        self.rxn_file_rxn_map_sgl.name, outfile=self.rxn_file_output.name)
        # this method / code is not exhaustive and work on instructions in self.rxn_file_makethese.name to make only
        # the listed combinations of reagents

        if self.do_single:
            self.single_exit_code1 = pymo.make_these_molecules(
                [self.rxn_file_ctx_sgl.name, self.rxn_file_frags_sgl.name], self.rxn_file_makethese_sgl.name,
                [self.rxn_file_rxn_map_sgl.name], outfile=self.rxn_file_output_sgl.name, debug=False)

            # print ("Ran make_these_molecules with exit code %s and output %s.smi" % (
            #        self.single_exit_code1, self.rxn_file_output_sgl.name))
            self.logger.info("Ran make_these_molecules with exit code %s and output %s.smi" % (
                self.single_exit_code1, self.rxn_file_output_sgl.name))

        if self.do_double:
            self.double_exit_code1 = pymo.make_these_molecules(
                [self.rxn_file_ctx_dbl1.name, self.rxn_file_frags1_dbl.name, self.rxn_file_frags2_dbl.name],
                self.rxn_file_makethese_dbl1.name, [self.rxn_file_rxn_map1_dbl.name, self.rxn_file_rxn_map2_dbl.name],
                outfile=self.rxn_file_output_dbl1.name, debug=False)

            # print ("Ran make_these_molecules with exit code %s and output %s.smi" % (self.double_exit_code1,
            #        self.rxn_file_output_dbl1.name))
            self.logger.info("Ran make_these_molecules with exit code %s and output %s.smi" % (
                self.double_exit_code1, self.rxn_file_output_dbl1.name))

            self.double_exit_code2 = pymo.make_these_molecules(
                [self.rxn_file_ctx_dbl2.name, self.rxn_file_frags1_dbl.name, self.rxn_file_frags2_dbl.name],
                self.rxn_file_makethese_dbl2.name, [self.rxn_file_rxn_map3_dbl.name, self.rxn_file_rxn_map4_dbl.name],
                outfile=self.rxn_file_output_dbl2.name, debug=False)

            # print ("Ran make_these_molecules with exit code %s and output %s.smi" % (
            #        self.double_exit_code2, self.rxn_file_output_dbl2.name))
            self.logger.info("Ran make_these_molecules with exit code %s and output %s.smi" % (
                self.double_exit_code2, self.rxn_file_output_dbl2.name))

    def yield_products_complex_id(self):
        """Reads output from do_reactions function, output file self.rxn_file_output.name plus .smi extension"""

        self.logger.info("Reading enumerated mols from output files")

        def yield_products_submethod(filename, cut_type):
            with open(filename + ".smi", "rt") as outfile:

                for line in outfile:

                    # Return will be new_mol_smi, id_including_the_new_frag, context (1 or two smi)
                    # do not know what fragment we have added
                    #
                    line = line.rstrip("\n")
                    line_list = line.split(" ")
                    # example output:
                    # [{smiles}, '1234___{frag_smi1}}', '+', '{frag_smi2}']
                    try:
                        new_mol = line_list[0]
                        orig_mol_id = int(line_list[1].split("___")[0])
                        frag_removed = line_list[1].split("___")[1]
                        if cut_type == 1:
                            frag_added = line_list[3]
                        elif cut_type == 2:
                            frag_added = line_list[1].split("___")[2]
                    except:
                        self.logger.debug("Skipped line as cannot parse id properly")
                        continue

                    self.logger.debug("I now have %s %s %s" % (new_mol, orig_mol_id, frag_added))
                    yield new_mol, orig_mol_id, frag_removed, frag_added

        if self.do_single:

            if self.single_exit_code1 == 0 and self.file_exists_check(self.rxn_file_output_sgl.name):
                for new_mol, orig_mol_id, frag_removed, frag_added in \
                        yield_products_submethod(self.rxn_file_output_sgl.name, 1):
                    yield new_mol, orig_mol_id, 'single', frag_removed, frag_added
            else:
                self.logger.debug("No results for Single Cut Enumeration")

        #
        if self.do_double:

            if self.double_exit_code1 == 0 and self.file_exists_check(self.rxn_file_output_dbl1.name):
                for new_mol, orig_mol_id, frag_removed, frag_added in \
                        yield_products_submethod(self.rxn_file_output_dbl1.name, 2):
                    yield new_mol, orig_mol_id, 'double', frag_removed, frag_added
            else:
                self.logger.debug("No results for Double Cut Enumeration (1 of 2)")

            if self.double_exit_code2 == 0 and self.file_exists_check(self.rxn_file_output_dbl2.name):
                for new_mol, orig_mol_id, frag_removed, frag_added in \
                        yield_products_submethod(self.rxn_file_output_dbl2.name, 2):
                    yield new_mol, orig_mol_id, 'double', frag_removed, frag_added
            else:
                self.logger.debug("No results for Double Cut Enumeration (2 of 2)")

    def yield_products_simple_dict_input(self):
        """Reads output from do_reactions function, output file self.rxn_file_output.name plus .smi extension"""

        self.logger.info("Reading enumerated mols from output files")

        def yield_products_simple_submethod(filename):
            with open(filename + ".smi", "rt") as outfile:

                for line in outfile:

                    line = line.rstrip("\n")
                    line_list = line.split(" ")

                    try:
                        new_mol = line_list[0]
                        context = line_list[1].split("___")[1]
                        frag = line_list[3].split("___")[1]

                    except:
                        self.logger.debug("Skipped line as cannot parse it")
                        continue

                    self.logger.debug("I now have %s %s which made %s" % (context, frag, new_mol))
                    yield context, frag, new_mol

        if self.do_single:

            if self.single_exit_code1 == 0 and self.file_exists_check(self.rxn_file_output_sgl.name):
                for context, frag, new_mol in yield_products_simple_submethod(self.rxn_file_output_sgl.name):
                    # print context, frag, new_mol
                    yield 'single', context, frag, new_mol

            else:
                self.logger.debug("No results for Single Cut Enumeration")

        if self.do_double:
            pass

    def write_products_to_csv(self, csv_filename):
        """
        Method to write out enumerated products to a CSV file
        :csv_filename name of CSV file to write output to
        """

        with open(csv_filename, 'w') as csv_out:

            headers = "NEW_MOL_SMI,ORIG_MOL_ID,CUT_TYPE,FRAG_REMOVED,FRAG_ADDED"
            csv_out.write(headers + "\n")

            for new_mol, orig_mol_id, cut_type, frag_removed, frag_added in self.yield_products_complex_id():

                write_string = ','.join([new_mol, str(orig_mol_id), cut_type, frag_removed, frag_added])
                csv_out.write(write_string + "\n")


#
# unittest everything
#
class _TestMMPEnumerateNewMols(unittest.TestCase):
    """Test class to test the object and methods"""

    @classmethod
    def setUpClass(self):
        #
        self.maxDiff = None

        # setup test data locations
        self.temp_file_input_csv = tempfile.NamedTemporaryFile(delete=False,
                                                               encoding='utf-8', mode='wt')
        self.temp_file_input_smi = tempfile.NamedTemporaryFile(delete=False, suffix='.smi',
                                                               encoding='utf-8', mode='wt')

        # setup a logger object
        self.mmplogger = logging.getLogger('mmpobjectclass_testlogger')
        logging.disable(logging.CRITICAL)

        #################################
        #
        # Now create test data
        #
        #################################

        self.test_input_smi_data = {
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL3552650/
            'CCN1CC2(CN(Cc3ccncc3)CCN(C2)C(=O)C(C)C)CC1=O': 3552650,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL2426642/
            'CC(C)(C(c1ccncc1)c2ccc3c(cnn3c4ccc(F)cc4)c2)C(=O)Nc5nccs5': 2426642
        }

        self.test_input_single_csv_header = \
            'FRAG_L,FRAG_R,MOL_CLP_DIFF,MOL_CLP_STDEV,NUM_VALS,PERCENT_THAT_INCREASE,PERCENT_THAT_DECREASE,PERCENT_WITH_NOCHANGE,MOL_MWT_DIFF'
        self.test_input_single_csv_data = {
            'Fc1c[n][1cH][n]c1,[1nH]1[n]ccc1,-0.202,,1,0.000,0.000,100.000,-30.000': None,
            '[1nH]1[n]ccc1,Fc1c[n][1cH][n]c1,0.202,,1,0.000,0.000,100.000,30.000': None,
            'O=[1CH]C(C)C,O=[1CH]CC,-0.309,0.000,10,10.000,10.000,80.000,-14.100': None,
            'O=[1CH]C(C)C,[1CH3]C(C)C,1.599,,1,0.000,0.000,100.000,-14.000': None,
            'O=[1CH]C(C)C,[1H],-0.134,0.335,7,0.000,42.860,57.140,-70.100': None,
            'O=[1CH]C(C)C,O=[1CH]c1ccc(C#N)cc1,0.355,,1,0.000,0.000,100.000,59.000': None,
            'O=[1CH]CC,O=[1CH]C(C)C,0.309,0.000,10,10.000,10.000,80.000,14.100': None,
            'O=[1CH]CC,[1CH3]C(C)C,1.908,,1,0.000,0.000,100.000,0.100': None,
            'O=[1CH]CC,[1H],-0.118,0.497,4,0.000,25.000,75.000,-56.000': None,
            'O=[1CH]CC,O=[1CH]c1ccc(C#N)cc1,0.664,,1,0.000,0.000,100.000,73.100': None,
            '[1CH3]C(C)C,O=[1CH]C(C)C,-1.599,,1,0.000,0.000,100.000,14.000': None,
            '[1CH3]C(C)C,O=[1CH]CC,-1.908,,1,0.000,0.000,100.000,-0.100': None,
            '[1CH3]c1cc[n]cc1,[1CH3]c1c[n]ccc1,0.000,0.000,12,0.000,0.000,100.000,0.000': None,
            '[n]1c[1cH]ccc1,[1H],-0.388,0.434,17,11.760,0.000,88.240,-77.100': None,
            '[n]1c[1cH]ccc1,CCOc1[1cH]cccc1,1.945,,1,0.000,0.000,100.000,43.100': None,
            '[1nH]1[n]ccc1,FC(c1[n][1nH]c(C2COC2)c1)(F)F,0.068,,1,0.000,0.000,100.000,124.100': None,
            '[1FH],[1BrH],0.630,0.186,9,0.000,0.000,100.000,60.900': None,
            # junk data but a good double cut example for testing
            "[12CH4],[12NH3],0.630,0.186,9,0.000,0.000,100.000,60.900": None,
            "[2CH3][1CH3],O=[2CH][1CH2]C,0.630,0.186,9,0.000,0.000,100.000,60.900": None,
            # test methyl replacement so change to something significant
            "[1CH4],[1IH],0.630,0.186,9,0.000,0.000,100.000,60.900": None,
            # included to allow H subst
            '[1H],O=[1CH]CC,-0.118,0.497,4,0.000,25.000,75.000,-56.000': None,
        }

        self.test_output_products_single = {
            'C(I)N1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(CC1=CC=NC=C1)C2': [3552650, '[1CH4]', '[1IH]'],
            'O=C1N(CC)CC2(CN(C(=O)C(C)C)CCN(C2)CC2=CN=CC=C2)C1': [3552650, '[1CH3]c1cc[n]cc1', '[1CH3]c1c[n]ccc1'],
            'O=C1N(CC2(CN(CC3=CC=NC=C3)CCN(C2)C(=O)CC)C1)CC': [3552650, 'O=[1CH]C(C)C', 'O=[1CH]CC'],
            'O=C1N(CC2(CN(CC3=CC=NC=C3)CCN(C2)CC(C)C)C1)CC': [3552650, 'O=[1CH]C(C)C', '[1CH3]C(C)C'],
            'O=C1N(CC2(CN(CC3=CC=NC=C3)CCN([H])C2)C1)CC': [3552650, 'O=[1CH]C(C)C', '[1H]'],
            'O=C1N(CC2(CN(CC3=CC=NC=C3)CCN(C2)C(=O)C2=CC=C(C#N)C=C2)C1)CC': [3552650, 'O=[1CH]C(C)C',
                                                                             'O=[1CH]c1ccc(C#N)cc1'],
            'O=C(C(C)I)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1': [3552650, '[1CH4]', '[1IH]'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(I)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1': [2426642, '[1CH4]', '[1IH]'],
            'O=C(NC1=NC=CS1)C(C)(C)C(C1=CC2=C(N(N=C2)C2=CC=C(Br)C=C2)C=C1)C1=CC=NC=C1': [2426642, '[1FH]', '[1BrH]']
        }

        self.test_output_products_single_incH = self.test_output_products_single.copy()
        self.test_output_products_single_incH.update({
            'C(CN1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(CC1=CC=NC=C1)C2)C(=O)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(CN(C(=O)C2)C(C)C(=O)CC)CN(CC2=CC=NC=C2)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(C(N(C(=O)C2)CC)C(=O)CC)CN(CC2=CC=NC=C2)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(C(N(CC3=CC=NC=C3)CC1)C(=O)CC)CN(C(=O)C2)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(C(C2=CC=NC=C2)C(=O)CC)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(CC2=C(C=NC=C2)C(=O)CC)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC(=C2)C(=O)CC)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC(N(CC2=CC=NC=C2)CC2(CN(C(=O)C2)CC)C1)C(=O)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1C(CN(CC2=CC=NC=C2)CC2(CN(C(=O)C2)CC)C1)C(=O)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1C(C2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1)C(=O)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)(C)C(=O)CC)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'C(C(C)C(=O)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1)C(=O)CC': [3552650, '[1H]', 'O=[1CH]CC'],
            'O=C(C(C)C)N1CC2(C(C(=O)N(C2)CC)C(=O)CC)CN(CC2=CC=NC=C2)CC1': [3552650, '[1H]', 'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C(=O)NC2=NC=CS2)CC(=O)CC)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                             'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)(C2=CC=NC=C2)C(=O)CC)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=C(C=NC=C2)C(=O)CC)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC(=NC=C2)C(=O)CC)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=C(C=C23)C(=O)CC)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC(=C23)C(=O)CC)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2C3=CC=C(C=C3C(=N2)C(=O)CC)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC(=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1)C(=O)CC': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=C(C=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1)C(=O)CC': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=C2C=CC(=C3C(=O)CC)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                             'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)N(C2=NC=CS2)C(=O)CC)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC(=CS2)C(=O)CC)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC'],
            'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=C(S2)C(=O)CC)C2=CC=NC=C2)C=C1': [2426642, '[1H]',
                                                                                               'O=[1CH]CC']
        })

        self.test_output_products_double = {
            'N(C)N1CC2(CN(C(=O)C(C)C)CCN(C2)CC2=CC=NC=C2)CC1=O': [3552650, '[12CH4]', '[12NH3]'],
            'N(N1CCN(CC2(CN(C(=O)C2)CC)C1)C(=O)C(C)C)C1=CC=NC=C1': [3552650, '[12CH4]', '[12NH3]']
        }

        # smi file - basic test
        for smi, smi_id in self.test_input_smi_data.items():
            self.temp_file_input_smi.write(smi+" "+str(smi_id)+"\n")
        self.temp_file_input_smi.close()

        # csv file - basic test
        self.temp_file_input_csv.write(self.test_input_single_csv_header+"\n")
        for data in list(self.test_input_single_csv_data.keys()):
            self.temp_file_input_csv.write(data+"\n")
        self.temp_file_input_csv.close()

        # container for results data
        self.test_dataset_testresults = {}


    @classmethod
    def tearDownClass(self):
        """Cleanup for end of all tests"""

        # remove(self.temp_file_input_smi.name)
        os.remove(self.temp_file_input_csv.name)

    def setUp(self):
        """Setup object for clean reuse in further tests"""
        # create empty mmp object
        self.test_mmp_pairs_object = MMPEnumerateNewMols(self.mmplogger)

    def tearDown(self):
        """Tear down object for clean reuse in further tests"""
        #
        self.test_mmp_pairs_object.clean_out_data_mmpenumerate()
        self.test_dataset_testresults.clear()

    def test_scan_input_smiles(self):
        """Test build_graph_from_pairs"""
        #
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)

        # returned chembl compounds
        self.assertEqual(self.test_mmp_pairs_object.mol_smi_dict, {2426642: 'CC(C)(C(c1ccncc1)c2ccc3c(cnn3c4ccc(F)cc4)c2)C(=O)Nc5nccs5',
                                                                   3552650: 'CCN1CC2(CN(Cc3ccncc3)CCN(C2)C(=O)C(C)C)CC1=O'})

    def test_fragment_reference_smi(self):
        """ """
        #
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        # new stuff
        self.test_mmp_pairs_object.fragment_reference_smi('SINGLE', 0.3)

        #print(self.test_mmp_pairs_object.unique_frags)
        # removed '[1H]' as we don't allow H substitution
        self.assertEqual(self.test_mmp_pairs_object.unique_frags,
                         {'s1[1cH][n]cc1', '[n]1cc[1cH]cc1', '[1CH4]', 'O=[1CH]Nc1scc[n]1', '[1CH3]C',
                          '[1NH2]c1scc[n]1', '[1FH]', '[1CH3]c1cc[n]cc1', 'C[1CH2]C', 'Fc1cc[1cH]cc1', 'O=[1CH]C(C)C'})

    def test_fragment_reference_smi_incH(self):
        """ """
        #
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        # new stuff
        self.test_mmp_pairs_object.fragment_reference_smi('SINGLE', 0.3, exclude_h_subst=False)

        #print(self.test_mmp_pairs_object.unique_frags)
        # removed '[1H]' as we don't allow H substitution
        self.assertEqual(self.test_mmp_pairs_object.unique_frags,
                         {'[1H]', 'O=[1CH]Nc1scc[n]1', 's1[1cH][n]cc1', '[1NH2]c1scc[n]1', 'C[1CH2]C', '[n]1cc[1cH]cc1',
                          'O=[1CH]C(C)C', '[1CH3]C', '[1CH4]', 'Fc1cc[1cH]cc1', '[1FH]', '[1CH3]c1cc[n]cc1'})

    def test_pairs_file_to_dict(self):
        """ """
        #
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name)

        self.assertEqual(len(result), 11)

    def test_pairs_file_to_dict_numpairslimit(self):
        """ """
        #
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name, num_pairs_limit=10)

        # result is not a 1:1 mapping of pairs to items but should be less than previous test
        self.assertEqual(len(result), 4)

    def test_add_transformation_group(self):
        """ """
        #
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name)
        self.test_mmp_pairs_object.add_transformation_group('test', result)

        self.assertEqual(len(self.test_mmp_pairs_object.transformation_groups), 1)

    def test_everything_at_once_single(self):
        """ """
        # input smi
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        self.test_mmp_pairs_object.fragment_reference_smi('SINGLE', 0.3)
        # get pairs
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name)
        self.test_mmp_pairs_object.add_transformation_group('test', result)
        # write rxn files
        self.test_mmp_pairs_object.write_rxn_files()
        self.test_mmp_pairs_object.write_reactants_mol_frag_dict()
        self.test_mmp_pairs_object.do_reactions()
        # new for test
        for new_mol, orig_mol_id, cut_type, frag_removed, frag_added in self.test_mmp_pairs_object.yield_products_complex_id():
            self.test_dataset_testresults[new_mol] = [orig_mol_id, frag_removed, frag_added]

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_output_products_single)

    def test_everything_at_once_single_simple_dict(self):
        """ """

        # fragmented forms of original test data - chembl molecules
        test_data = {
            ('[1CH3]N1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(Cc1cc[n]cc1)C2', '[1CH4]'): None,
            ('O=C1N(CC)CC2(CN(C(=O)C(C)C)CC[1NH]C2)C1', '[1CH3]c1cc[n]cc1'): None,
            ('[1CH3]N1CC2(CN(C(=O)C2)CC)CN(C(=O)C(C)C)CC1','[n]1cc[1cH]cc1'): None,
            ('O=C1N(CC2(CN(Cc3cc[n]cc3)CC[1NH]C2)C1)CC', 'O=[1CH]C(C)C'): None,
            ('CCN1C(=O)CC2(C1)CN(Cc1cc[n]cc1)CCN([1CH]=O)C2', 'C[1CH2]C'): None,
            ('O=C([1CH2]C)N1CC2(CN(C(=O)C2)CC)CN(Cc2cc[n]cc2)CC1', '[1CH4]'): None,
            ('O=C(C(C)C)N1CC2(C[1NH]C(=O)C2)CN(Cc2cc[n]cc2)CC1', '[1CH3]C'): None,
            ('[1CH3]CN1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(Cc1cc[n]cc1)C2', '[1H]'): None,
            ('O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN([1CH2]c2cc[n]cc2)CC1', '[1H]'): None,
            ('Fc1ccc([n]2[n]cc3c2ccc(C([1CH](C(=O)Nc2scc[n]2)C)c2cc[n]cc2)c3)cc1', '[1CH4]'): None,
            ('Fc1ccc([n]2[n]cc3c2ccc([1CH2]C(C(=O)Nc2scc[n]2)(C)C)c3)cc1', '[n]1cc[1cH]cc1'): None,
            ('O=C(Nc1scc[n]1)C(C(c1cc2c([1nH][n]c2)cc1)c1cc[n]cc1)(C)C', 'Fc1cc[1cH]cc1'): None,
        }

        # write rxn files
        self.test_mmp_pairs_object.write_rxn_files()
        self.test_mmp_pairs_object.write_reactants_simple_dict(test_data)
        self.test_mmp_pairs_object.do_reactions()
        #
        for cut_type, context, frag, new_mol in self.test_mmp_pairs_object.yield_products_simple_dict_input():
            self.test_dataset_testresults[context, frag, new_mol] = None

        #print(self.test_dataset_testresults)
        self.assertEqual({('[1CH3]N1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(Cc1cc[n]cc1)C2', '[1CH4]',
                           'C(C)N1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(CC1=CC=NC=C1)C2'): None, (
                          'O=C1N(CC)CC2(CN(C(=O)C(C)C)CC[1NH]C2)C1', '[1CH3]c1cc[n]cc1',
                          'O=C1N(CC)CC2(CN(C(=O)C(C)C)CCN(C2)CC2=CC=NC=C2)C1'): None, (
                          '[1CH3]N1CC2(CN(C(=O)C2)CC)CN(C(=O)C(C)C)CC1', '[n]1cc[1cH]cc1',
                          'C(N1CC2(CN(C(=O)C2)CC)CN(C(=O)C(C)C)CC1)C1=CC=NC=C1'): None, (
                          'O=C1N(CC2(CN(Cc3cc[n]cc3)CC[1NH]C2)C1)CC', 'O=[1CH]C(C)C',
                          'O=C1N(CC2(CN(CC3=CC=NC=C3)CCN(C2)C(=O)C(C)C)C1)CC'): None, (
                          'CCN1C(=O)CC2(C1)CN(Cc1cc[n]cc1)CCN([1CH]=O)C2', 'C[1CH2]C',
                          'CCN1C(=O)CC2(C1)CN(CC1=CC=NC=C1)CCN(C(=O)C(C)C)C2'): None, (
                          'O=C([1CH2]C)N1CC2(CN(C(=O)C2)CC)CN(Cc2cc[n]cc2)CC1', '[1CH4]',
                          'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1'): None, (
                          'O=C(C(C)C)N1CC2(C[1NH]C(=O)C2)CN(Cc2cc[n]cc2)CC1', '[1CH3]C',
                          'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(CC2=CC=NC=C2)CC1'): None, (
                          '[1CH3]CN1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(Cc1cc[n]cc1)C2', '[1H]',
                          'C([H])CN1C(=O)CC2(C1)CN(C(=O)C(C)C)CCN(CC1=CC=NC=C1)C2'): None, (
                          'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN([1CH2]c2cc[n]cc2)CC1', '[1H]',
                          'O=C(C(C)C)N1CC2(CN(C(=O)C2)CC)CN(C([H])C2=CC=NC=C2)CC1'): None, (
                          'Fc1ccc([n]2[n]cc3c2ccc(C([1CH](C(=O)Nc2scc[n]2)C)c2cc[n]cc2)c3)cc1', '[1CH4]',
                          'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1'): None, (
                          'Fc1ccc([n]2[n]cc3c2ccc([1CH2]C(C(=O)Nc2scc[n]2)(C)C)c3)cc1', '[n]1cc[1cH]cc1',
                          'FC1=CC=C(N2N=CC3=CC(=CC=C23)C(C(C)(C)C(=O)NC2=NC=CS2)C2=CC=NC=C2)C=C1'): None, (
                          'O=C(Nc1scc[n]1)C(C(c1cc2c([1nH][n]c2)cc1)c1cc[n]cc1)(C)C', 'Fc1cc[1cH]cc1',
                          'O=C(NC1=NC=CS1)C(C)(C)C(C1=CC2=C(N(N=C2)C2=CC=C(F)C=C2)C=C1)C1=CC=NC=C1'): None},
                          self.test_dataset_testresults)

    def test_everything_at_once_single_incH(self):
        """ """
        # input smi
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        self.test_mmp_pairs_object.fragment_reference_smi('SINGLE', 0.3, exclude_h_subst=False)
        # get pairs
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name, exclude_h_subst=False)
        self.test_mmp_pairs_object.add_transformation_group('test', result)
        # write rxn files
        self.test_mmp_pairs_object.write_rxn_files()
        self.test_mmp_pairs_object.write_reactants_mol_frag_dict()

        self.test_mmp_pairs_object.do_reactions()
        # new for test
        for new_mol, orig_mol_id, cut_type, frag_removed, frag_added in \
                self.test_mmp_pairs_object.yield_products_complex_id():
            self.test_dataset_testresults[new_mol] = [orig_mol_id, frag_removed, frag_added]

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_output_products_single_incH)

    def test_everything_at_once_double(self):
        """ """
        # input smi
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        self.test_mmp_pairs_object.fragment_reference_smi('DOUBLE', 0.3)
        # get pairs
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name)
        self.test_mmp_pairs_object.add_transformation_group('test', result)
        # write rxn files
        self.test_mmp_pairs_object.write_rxn_files()
        self.test_mmp_pairs_object.write_reactants_mol_frag_dict()
        self.test_mmp_pairs_object.do_reactions()

        # pp.pprint(self.test_mmp_pairs_object.mol_fragments_dict_double)
        # new for test
        for new_mol, orig_mol_id, cut_type, frag_removed, frag_added in \
                self.test_mmp_pairs_object.yield_products_complex_id():
            self.test_dataset_testresults[new_mol] = [orig_mol_id, frag_removed, frag_added]

        #print(self.test_dataset_testresults)
        #
        self.assertEqual(self.test_dataset_testresults, self.test_output_products_double)

    def test_everything_at_once_both(self):
        """ """
        # input smi
        self.test_mmp_pairs_object.scan_input_smiles(self.temp_file_input_smi.name, injest=True)
        self.test_mmp_pairs_object.fragment_reference_smi('BOTH', 0.3)
        # get pairs
        result = self.test_mmp_pairs_object.pairs_file_to_dict(self.temp_file_input_csv.name)
        self.test_mmp_pairs_object.add_transformation_group('test', result)
        # write rxn files
        self.test_mmp_pairs_object.write_rxn_files()
        self.test_mmp_pairs_object.write_reactants_mol_frag_dict()
        self.test_mmp_pairs_object.do_reactions()

        # introspect
        #print(self.test_mmp_pairs_object.mol_fragments_dict_single)
        #print(self.test_mmp_pairs_object.mol_fragments_dict_double)

        # pp.pprint(self.test_mmp_pairs_object.mol_fragments_dict_double)
        # new for test
        for new_mol, orig_mol_id, cut_type, frag_removed, frag_added in \
                self.test_mmp_pairs_object.yield_products_complex_id():
            self.test_dataset_testresults[new_mol] = [orig_mol_id, frag_removed, frag_added]

        # get the results by merging single and double
        results = self.test_output_products_single.copy()
        results.update(self.test_output_products_double)
        #print(results)
        #
        self.assertEqual(self.test_dataset_testresults, results)

if __name__ == '__main__':
    unittest.main()
