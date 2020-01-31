###################################################################
""" Summary: Class and Methods for deriving MCSS based MMP's

About:  Derive a matched pair based MCSS from a pair molecules

To do:  - extend the method enumerate_fragment_properties to also
          enumerate self.mol_smi_dict as this would allow the addition 
          of a flag '-p' that prints out whole molecule props alongside
          MCSS and therefore compute %molecule that the MCSS covers
        - could use other descriptors from IW code to get MCSS via 
          bond count not #Atoms or 
        - Should move iterators in process_mcss_list_to_string to be numeric
          and store numeric ID's in self.largest_mcs_mmp_double/single
        - could allow further switched to change behaviour of tie break
          where single/double  or double alone give tie break MCSS
          [connected substructures versus disconnected or both/either]

        - Extension to triple cut would allow improved search/match e.g.:
          N1(C(c2c(cc3c(c2)OCO3)CC1)c4cc(c(c(c4)OC)O)OC)C(=O)OC CHEMBL311765
          N1(C(c2c(cc(cc2)O)CC1)c3ccc(cc3)OCCN4CCCC4)C(=O)OCC CHEMBL94080
"""
###################################################################
import logging
import csv
import os
import sys
import unittest
import tempfile

from builtins import range

from mmp.mmp_data_objects import MMPDataObjectClass

if 'LILLYMOL_HOME' in os.environ:
    import pybase.pyopmo as pymo
else:
    import pybase.pymo as pymo


class MMPbasedMCSSObjectClass(MMPDataObjectClass):

    def __init__(self, logger_object):
        """
        Example usage:
            mmplogger = logging.getLogger('lillymol_file_logger')
            logging.disable(logging.CRITICAL)
            my_mmp_mcss_object = MMPbasedMCSSObjectClass(mmplogger)

        """

        MMPDataObjectClass.__init__(self, logger_object)

        self.logger = logger_object
        if len(logging.Logger.manager.loggerDict) < 1:
            # exit with system status 1 and custom error
            sys.exit("Invalid or no logger object passed to MMPObjectClass.  Please create \
                    and pass a logger and set to use logging.disable if you don't want logging")

        # this is used for storing the largest MCS MMP for given pair
        self.largest_mcs_mmp_result = {}
        self.ref_smi_props = {}

    def clean_out_data_mcss_obj(self):
        """Method to clean out all objects in class"""

        self.clean_out_data()
        self.mcs_mmp.clear()

    def enumerate_fragment_properties(self):
        """Writes out the ref_smi_dict to disk, calculates natoms, returns data to self.ref_smi_props
        Some complexities in method such as double cut fragments (iw_descr only calcs largest frag)"""

        frag_smi_file = tempfile.NamedTemporaryFile(delete=False, suffix='.smi')
        frag_smi_props_out = tempfile.NamedTemporaryFile(delete=False)

        with open(frag_smi_file.name, "w") as f:
            for item in self.refsmi_dict:
                if isinstance(item, int):
                    # can't see an easy way to do this except string compare, [1H] causes iw_descr to crash out
                    if self.refsmi_dict[item] != '[1H]':
                        f.write(self.refsmi_dict[item]+" "+str(item)+"\n")

        # run pymo.iwdescr
        self.logger.info("Running pymo.iwdescr on %s smi with in:%s, out:%s" %
                         (len(self.refsmi_dict), frag_smi_file.name, frag_smi_props_out.name))
        exit_status = pymo.iwdescr(frag_smi_file.name, frag_smi_props_out.name, params_dict={'-l': '', '-v': ''},
                                   loggero=self.logger)
        self.logger.debug("Ran iwdescr with exit status %s" % exit_status)

        with open(frag_smi_props_out.name, "r") as csv_file:

            reader = csv.reader(csv_file, delimiter=' ')
            i = -1
            for row in reader:
               
                i += 1
                # if header row, append headers
                if i == 0:
                    if row[1] != 'w_natoms':
                        self.logger.warn("When this was written, NATOMs was in array position 1 (zero indexed) with "
                                         "column title w_natoms.  Now it's not, it's: %s" % row[1])
                        sys.exit("When this was written, NATOMs was in array position 1 (zero indexed) with column "
                                 "title w_natom.  Now it's not, it's: %s" % row[1])
                    continue
                
                # we trust there is only one entry per id
                # print row[0], row[1]
                self.ref_smi_props[int(row[0])] = int(row[1])
        
        frag_smi_props_out.close()

        self.logger.debug("Completed load of %s mol props from dict of %s from file %s" %
                          (len(self.ref_smi_props), len(self.refsmi_dict)/2, frag_smi_props_out.name))

    def get_largest_mcs_pairs(self, out_file, cut_type, mdc_atm_soft=None, mdc_atm_soft_threshold=None,
                              mdc_atm_hard=None):
        """Method to print out a single smi - smi pair from the input CSV with data differences. Selection of the
        exact matched pair for a given smi - smi combination is based on the largest Maximum Common Substructure
        which equates to the MMP with the smallest MWT/#Atoms difference across all MMP's for that smi/smi combo

        out_file:
          The user specified output file
        
        cut_type:
          Specifies the type of fragmentation required.  Allowed values are SINGLE,
          DOUBLE or BOTH.  Currently this class does not support anything greater than
          double cut fragmentation

        mdc_atm_hard:
          max double cut atom cutoff (hard)
          Never consider double cut context fragments where one half has num_atoms <= mdc_atm_hard
          i.e.: this is a hard cutoff filter implemented during dicer parsing

        mdc_atm_soft:
          max double cut atom cutoff (soft)
          * must be used with mdc_atm_soft_threshold
          When double cut is greater than single, if one part of double context has num_atoms <= mdc_atm_soft and
          total double cut atom <= single cut atoms + mdc_atm_soft_threshold then discard

        mdc_atm_soft_threshold:
          max double cut atom cutoff threshold (soft)
          * must be used with mdc_atm_soft
          This gets added to single cut num atoms each comparison that's done, if and when mdc_atm_soft is set
          see details of mdc_atm_soft

        Example usage:

          # give me a CSV named my_output.pairs of all MCS based pairs:
          my_mmp_object.get_largest_mcs_pairs('myoutput.csv', 'BOTH', 'DICER')

          # give me a CSV of only the DOUBLE cut MCS based pairs with RDKit attachment points:
          my_mmp_object.get_largest_mcs_pairs('myoutput.csv', 'DOUBLE', 'RDKIT')
        """

        if (mdc_atm_soft is not None and mdc_atm_soft_threshold is None) or\
           (mdc_atm_soft is None and mdc_atm_soft_threshold is not None):
            sys.exit("Error, mdc_atm_soft and mdc_atm_soft_threshold must be specified together.")

        def process_mcss_list_to_string(prefix, input_list):
            """sub method to build a printable string from input list of specific structure"""

            out_string = ''
            num_of_entries = len(input_list)

            if num_of_entries > 4:
                
                for i_ in range(0, num_of_entries, 4):
                    out_string = out_string + prefix + "_" + str((i_/4)+1) + "," + str(molid_L) + "," + str(molid_R)
                    out_string = out_string + "," + str(sum(input_list[0 + i_])) + "," + str(input_list[1 + i_]) + ","
                    out_string = out_string + str(input_list[2 + i_]) + "," + str(input_list[3 + i_])
                    out_string += "\n"

            else:
                if len(input_list[1]) > 1:
                    ctx_smi = self.refsmi_dict[input_list[1][0]] + "." + self.refsmi_dict[input_list[1][1]]
                else:
                    ctx_smi = self.refsmi_dict[input_list[1][0]]
                out_string = prefix + "," + str(molid_L) + "," + str(molid_R) + ","
                out_string = out_string + str(sum(input_list[0])) + "," + ctx_smi + ","
                out_string = out_string + str(self.refsmi_dict[input_list[2]]) + "," \
                             + str(self.refsmi_dict[input_list[3]])
                out_string += "\n"

            return out_string

        def disambiguate_double_list(input_list):
            """sub method to untangle double cut tie break cases"""
            
            num_of_entries = len(input_list)

            filtered_list = []

            # The tie code should have only saved the example with the largest 'smallest fragment' size
            # so now we just take the first example where atom numbering [1 before [2
            # Theoretically, if two different examples of a double cut fragmentation pattern exist with the same number
            # of atoms *in both parts* of the context, then there is another tie break here.  e.g.:
            # num_atoms in context = (2,10) should always appear not (1,11) but can't disentangle many (1,11)
            # Decided not to handle this and instead just take the first one with the ordered numbering
            for i_ in range(0, num_of_entries, 4):

                # only use if the isomeric label is the right way round, [1 before [2
                if '[1' in self.refsmi_dict[input_list[1 + i_][0]]:
                    filtered_list = input_list[(0 + i_): (4 + i_)]
                else:
                    continue
                
            return filtered_list

        def remove_atom_num_dupes(input_list):
            """sub method to get only 1 example of simple isomeric numbering flip"""

            # only use if the isomeric label is the right way round, [1 before [2
            if '[1' in self.refsmi_dict[input_list[1][0]]:
                # take the first 4 items
                output_list = input_list[:4]

            else:
                # just take the last 4 items
                output_list = input_list[-4:]

            return output_list

        self.logger.info('Opening output file for write: %s' % out_file)

        # check cut_type, convert to int
        if cut_type.upper() == 'DOUBLE':
            # confusing but faster later
            cut_type_id = 3
        elif cut_type.upper() == 'BOTH':
            # confusing but faster later
            cut_type_id = 2
        elif cut_type.upper() == 'SINGLE':
            cut_type_id = 1
        else:
            self.logger.warn('cut_type specification is incorrect, using single cut: %s' % cut_type.upper())
            cut_type_id = 1

        # fail if both single_pairs_dict and double_pairs_dict are empty
        if (len(self.single_pairs_dict) == 0) and (len(self.double_pairs_dict) == 0):
            self.logger.debug('No data found in single_pairs_dict and/or double_pairs_dict, expect no results')
            # sys.exit("Error: no data found in single_pairs_dict and/or double_pairs_dict, nothing to find and write")
        
        #
        # Here we build data structures of type:
        #  self.largest_mcs_mmp_result[(molid_L, molid_R)] = [(#atoms, #atoms or None),
        #     (context_id, context_id or None), frag_Left_id, frag_Right_id]
        #

        # single - this is easy as we only keep/store the one with the greatest number of atoms
        if cut_type_id <= 2:
            
            for molid_L, molid_R, ctx_id, frag_L_id, frag_R_id in \
                    self.iterator_single_pairs_dict_numeric(inc_attachpt=False):
            
                if (molid_L, molid_R) in self.largest_mcs_mmp_result:

                    if self.largest_mcs_mmp_result[(molid_L, molid_R)][0][0] <= self.ref_smi_props[ctx_id]:

                        if self.largest_mcs_mmp_result[(molid_L, molid_R)][0][0] == self.ref_smi_props[ctx_id]:

                            self.largest_mcs_mmp_result[(molid_L, molid_R)].extend(
                                [(self.ref_smi_props[ctx_id], ), (ctx_id, ), frag_L_id, frag_R_id])
                        
                        else:
                            self.largest_mcs_mmp_result[(molid_L, molid_R)] = [
                                (self.ref_smi_props[ctx_id], ), (ctx_id, ), frag_L_id, frag_R_id]
        
                else:
                    self.largest_mcs_mmp_result[(molid_L, molid_R)] = [
                        (self.ref_smi_props[ctx_id], ), (ctx_id, ), frag_L_id, frag_R_id]

        # now build the final results on the fly
        # double - for each one we compare against what we already have in self.largest_mcs_mmp_result

        ctx_natoms = None

        if cut_type_id >= 2:

            for molid_L, molid_R, ctx1_id, ctx2_id, frag_L_id, frag_R_id in \
                    self.iterator_double_pairs_dict_numeric(inc_attachpt=False):

                #
                if ctx1_id in self.ref_smi_props:
                    ctx_natoms = (self.ref_smi_props[ctx1_id], )
                else:
                    ctx1_smi = self.refsmi_dict[ctx1_id]
                    ctx1_smi = ctx1_smi.replace("[1", "[9")
                    ctx1_smi = ctx1_smi.replace("[2", "[1")
                    ctx1_smi = ctx1_smi.replace("[9", "[2")
                    try:
                        ctx_natoms = (self.ref_smi_props[self.refsmi_dict[ctx1_smi]], )
                    except:
                        print("ERR >>>")
                        print(("{} {} {} {} {} {}".format(molid_L, molid_R, ctx1_id, ctx2_id, frag_L_id, frag_R_id)))
                        print(("{} {} {}".format(ctx1_id, ctx1_smi, self.refsmi_dict[ctx1_smi])))
                        print("")

                if ctx2_id in self.ref_smi_props:
                    ctx_natoms = ctx_natoms + (self.ref_smi_props[ctx2_id], )
                else:
                    ctx2_smi = self.refsmi_dict[ctx2_id]
                    ctx2_smi = ctx2_smi.replace("[1", "[9")
                    ctx2_smi = ctx2_smi.replace("[2", "[1")
                    ctx2_smi = ctx2_smi.replace("[9", "[2")
                    ctx_natoms = ctx_natoms + (self.ref_smi_props[self.refsmi_dict[ctx2_smi]], )

                # If the indicator flag check_all_context is set to true we need to pre-filter all ctx fragments
                # to ensure they are greater than or equal to the specified limit for mdc_atm_hard (maximum double
                # cut atoms hard limit).  This is a crude filter and could remove valid double cut MCSS.
                if mdc_atm_hard is not None:
                    if ctx_natoms[0] <= mdc_atm_hard:
                        continue

                    elif ctx_natoms[1] <= mdc_atm_hard:
                        continue

                #
                # Main
                # have we seen this smi - smi pair before?
                if (molid_L, molid_R) in self.largest_mcs_mmp_result:

                    # get the number of atoms in the context
                    num_atoms_existing = self.largest_mcs_mmp_result[(molid_L, molid_R)][0]
                    if len(num_atoms_existing) > 1:
                        total_num_atoms_existing = sum(num_atoms_existing)
                    else:
                        total_num_atoms_existing = num_atoms_existing[0]
                    total_num_atoms_new = sum(ctx_natoms)

                    if total_num_atoms_new > total_num_atoms_existing:

                        # if it is a double and we have a min fragment setting
                        if mdc_atm_soft is not None:

                            # if it falls below the threshold at which we apply this min frag setting
                            if total_num_atoms_new <= (total_num_atoms_existing + mdc_atm_soft_threshold):
                                # only keep if both frag sizes are legal
                                if '[1' in self.refsmi_dict[ctx1_id]:
                                    if (ctx_natoms[0] > mdc_atm_soft) and (ctx_natoms[1] > mdc_atm_soft):
                                        self.largest_mcs_mmp_result[(molid_L, molid_R)] = \
                                            [ctx_natoms, (ctx1_id, ctx2_id), frag_L_id, frag_R_id]

                            # above threshold so keep anyway
                            else:
                                if '[1' in self.refsmi_dict[ctx1_id]:
                                    self.largest_mcs_mmp_result[(molid_L, molid_R)] = \
                                        [ctx_natoms, (ctx1_id, ctx2_id), frag_L_id, frag_R_id]

                        else:
                            if '[1' in self.refsmi_dict[ctx1_id]:
                                self.largest_mcs_mmp_result[(molid_L, molid_R)] = \
                                    [ctx_natoms, (ctx1_id, ctx2_id), frag_L_id, frag_R_id]

                    # tie-break
                    elif total_num_atoms_new == total_num_atoms_existing:

                        # single always wins over double, so only consider this if existing is double
                        # double cut tie breaks get disambiguated later using custom function
                        if len(num_atoms_existing) == 1:
                            continue

                        else:
                            # consider the size of the 'smallest fragment' and add if same, replace if bigger,
                            # drop if smaller
                            if min(ctx_natoms) > min(num_atoms_existing):
                                if '[1' in self.refsmi_dict[ctx1_id]:
                                    self.largest_mcs_mmp_result[(molid_L, molid_R)] = \
                                        [ctx_natoms, (ctx1_id, ctx2_id), frag_L_id, frag_R_id]

                            elif min(ctx_natoms) == min(num_atoms_existing):
                                self.largest_mcs_mmp_result[(molid_L, molid_R)].extend(
                                    [ctx_natoms, (ctx1_id, ctx2_id), frag_L_id, frag_R_id])

                            else:
                                # don't store as we have a better context with a larger 'smallest fragment'
                                continue

                    # double cut context must be smaller than what we already have so discard this new one
                    else:
                        continue

                else:
                    # new result, case where we only have a double cut MCSS so add it!
                    if '[1' in self.refsmi_dict[ctx1_id]:
                        self.largest_mcs_mmp_result[(molid_L, molid_R)] = [ctx_natoms, (ctx1_id, ctx2_id),
                                                                           frag_L_id, frag_R_id]

        with open(out_file, "w") as final_out:

            final_out.write('CUT_TYPE,MOL_ID_L,MOL_ID_R,NATOMS,MCSS,FRAG_L,FRAG_R\n')
            
            # do single cut first as these take precedence above a double
            for (molid_L, molid_R) in self.largest_mcs_mmp_result:

                list_length = len(self.largest_mcs_mmp_result[(molid_L, molid_R)])
                # the list self.largest_mcs_mmp_result[(molid_L, molid_R)] contains an ordered list of items
                # the first 4 are (1) a tuple of the num_atoms (2) fragment (3&4) context in two parts
                # Therefore if the list is greater than 8 items it means we have more than one double
                # cut that we need to consider, possibly as a double cut tie break.  We do not consider the
                # case where there are 8 items as we know this will be two identical fragmentation patterns
                # with differing isomeric numbering on the atom attachment points therefore we use >8 not >=8
                if list_length > 8:
                    if len(self.largest_mcs_mmp_result[(molid_L, molid_R)][0]) == 1:
                        # disambiguate single cut list
                        final_out.write(process_mcss_list_to_string('SINGLE', self.largest_mcs_mmp_result[
                                                                                  (molid_L, molid_R)][0:4]))

                    else:
                        # print("Double won (a): ", molid_L, molid_R, self.largest_mcs_mmp_result[(molid_L, molid_R)])
                        new_list = disambiguate_double_list(self.largest_mcs_mmp_result[(molid_L, molid_R)])
                        final_out.write(process_mcss_list_to_string('DOUBLE', new_list))

                elif list_length == 4:
                    # print("Single won (a): ", molid_L, molid_R, self.largest_mcs_mmp_result[(molid_L, molid_R)])
                    final_out.write(process_mcss_list_to_string('SINGLE', self.largest_mcs_mmp_result[
                        (molid_L, molid_R)]))

                else:
                    # print("Double wins (b): ", molid_L, molid_R, self.largest_mcs_mmp_result[(molid_L, molid_R)])
                    # need to remove atom numbering dupes then print
                    new_list = remove_atom_num_dupes(self.largest_mcs_mmp_result[(molid_L, molid_R)])
                    final_out.write(process_mcss_list_to_string('DOUBLE', new_list))


class _TestMMPbasedMCSSObjectClass(unittest.TestCase):
    """Test class for MMPDataObjectClass(object) written to use pythons unittest

    Example usage:

     python mmp_mcss_objects.py

     coverage run mmp_mcss_objects.py
     coverage report mmp_mcss_objects.py

    """

    def setUp(self):

        """Instantiate temp file names, test data objects that get written to temp files
        a silent logger object (needed to instantiate class) and the mmp object we'll test"""

        self.maxDiff = None

        # setup test data location use tempfile.NamedTemporaryFile(delete=False) to persist data on disk
        self.temp_file_input_smi_01 = tempfile.NamedTemporaryFile(delete=False, suffix=".smi",
                                                                  encoding='utf-8', mode='wt')
        self.temp_file_input_smi_03 = tempfile.NamedTemporaryFile(delete=False, suffix=".smi",
                                                                  encoding='utf-8', mode='wt')
        self.temp_file_output_pairs = tempfile.NamedTemporaryFile(delete=False)

        # setup a logger object
        self.mmplogger = logging.getLogger('mmpobjectclass_testlogger')
        # logging.disable(logging.CRITICAL)

        # create empty mmp object
        self.test_mmp_mcss_object = MMPbasedMCSSObjectClass(self.mmplogger)

        # data set for use in testing input
        self.test_dataset_goldeninput_smi_01 = {
            # The following represent synthetic data, analogues of CHEMBL1382609
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1382609/
            # 1. substituents are added to the pyrazole ring to generate side chain MMPs
            #    H on CHEMBL1382609 between two methyls is changed to Br, F, C, I to
            #    visually see the change in the smiles string (avoiding Cl as already present)
            #    e.g.: N1C(=C(Br)C(=N1)C)C
            # 2. core ring system is modified (phenyl to pyridine) to see ring switch MMP's
            #    Presence/Absence of Pyridine-N and N-positional isomerism in Cl-Ph ring
            #    e.g.: C2=NC(=CS2)C2=CC=C(Cl)C=C2 + addition of N ->
            #          C2=NC(=CS2)C2=CN=C(Cl)C=C2 + move N around ring ->
            #          C2=NC(=CS2)C2=NC=C(Cl)C=C2
            # for 1,2 single wins
            '001': 'N1(C2=NC(=CS2)C2=CC=C(Cl)C=C2)C(=C(Br)C(=N1)C)C',
            '002': 'N1(C2=NC(=CS2)C2=CC=C(Cl)C=C2)C(=C(F)C(=N1)C)C',
            # for 2,5 double wins tie
            '003': 'N1(C2=NC(=CS2)C2=CN=C(Cl)C=C2)C(=C(F)C(=N1)C)C',
            # The following represent synthetic data, analogues of CHEMBL1341352
            # for 1341352 and it's synthetic unsubstituted analogue there is no double
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1341352/
            '1341352': 'Cc1cc(nn1CC(=O)NCc2ccccc2)C(F)(F)F',
            '004': 'c1cc(nn1CC(=O)NCc2ccccc2)',
            # more double cut only
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL6211
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL6232
            '6211': 'O=C(OCC1N(C(=O)c2cc(c(OC)c(c2)OC)OC)CCN(C1)C(=O)c1cc(c(OC)c(OC)c1)OC)CCCCCCC',
            '6232': 'O=C(N1C(CN(C(=O)c2cc(c(OC)c(c2)OC)OC)CC1)COC(=O)CC(C)(C)C)c1cc(c(OC)c(OC)c1)OC'
        }

        self.test_dataset_goldeninput_smi_03 = {
            # repeat of above
            '001': 'N1(C2=NC(=CS2)C2=CC=C(Cl)C=C2)C(=C(Br)C(=N1)C)C',
            '002': 'N1(C2=NC(=CS2)C2=CC=C(Cl)C=C2)C(=C(F)C(=N1)C)C',
        }

        # all smiles are output from above input as either a repeat smiles or a fragment of them
        self.test_dataset_golden_output_01 = {'CUT_TYPE,MOL_ID_L,MOL_ID_R,NATOMS,MCSS,FRAG_L,FRAG_R': None,
                                              'SINGLE,1,2,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1BrH],[1FH]': None,
                                              'SINGLE,2,1,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1FH],[1BrH]': None,
                                              'DOUBLE,2,3,14,[1ClH].Fc1c([n](c2sc[2cH][n]2)[n]c1C)C,[1cH]1cc[2cH]cc1,[n]1[1cH]cc[2cH]c1': None,
                                              'DOUBLE,3,2,14,[1ClH].Fc1c([n](c2sc[2cH][n]2)[n]c1C)C,[n]1[1cH]cc[2cH]c1,[1cH]1cc[2cH]cc1': None,
                                              'SINGLE,1341352,4,11,O=C(NCc1ccccc1)[1CH3],Cc1[1nH][n]c(C(F)(F)F)c1,[1nH]1[n]ccc1': None,
                                              'SINGLE,4,1341352,11,O=C(NCc1ccccc1)[1CH3],[1nH]1[n]ccc1,Cc1[1nH][n]c(C(F)(F)F)c1': None,
                                              'DOUBLE,6211,6232,40,[1CH4].[2CH3]C(=O)OCC1N(C(=O)c2cc(c(OC)c(c2)OC)OC)CCN(C1)C(=O)c1cc(c(OC)c(OC)c1)OC,[1CH3]CCC[2CH3],C[12CH2]C': None,
                                              'DOUBLE,6232,6211,40,[1CH4].[2CH3]C(=O)OCC1N(C(=O)c2cc(c(OC)c(c2)OC)OC)CCN(C1)C(=O)c1cc(c(OC)c(OC)c1)OC,C[12CH2]C,[1CH3]CCC[2CH3]': None}

        self.test_dataset_golden_output_02 = {'CUT_TYPE,MOL_ID_L,MOL_ID_R,NATOMS,MCSS,FRAG_L,FRAG_R': None,
                                              'SINGLE,1,2,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1BrH],[1FH]': None,
                                              'SINGLE,2,1,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1FH],[1BrH]': None,
                                              'SINGLE,2,3,13,Fc1c([n](c2sc[1cH][n]2)[n]c1C)C,Clc1cc[1cH]cc1,Clc1[n]c[1cH]cc1': None,
                                              'SINGLE,3,2,13,Fc1c([n](c2sc[1cH][n]2)[n]c1C)C,Clc1[n]c[1cH]cc1,Clc1cc[1cH]cc1': None,
                                              'SINGLE,1341352,4,11,O=C(NCc1ccccc1)[1CH3],Cc1[1nH][n]c(C(F)(F)F)c1,[1nH]1[n]ccc1': None,
                                              'SINGLE,4,1341352,11,O=C(NCc1ccccc1)[1CH3],[1nH]1[n]ccc1,Cc1[1nH][n]c(C(F)(F)F)c1': None,
                                              'SINGLE,6211,6232,39,[1CH3]C(=O)OCC1N(C(=O)c2cc(c(OC)c(c2)OC)OC)CCN(C1)C(=O)c1cc(c(OC)c(OC)c1)OC,[1CH3]CCCCC,C[1CH](C)C': None,
                                              'SINGLE,6232,6211,39,[1CH3]C(=O)OCC1N(C(=O)c2cc(c(OC)c(c2)OC)OC)CCN(C1)C(=O)c1cc(c(OC)c(OC)c1)OC,C[1CH](C)C,[1CH3]CCCCC': None}

        self.test_dataset_golden_output_03 = {'CUT_TYPE,MOL_ID_L,MOL_ID_R,NATOMS,MCSS,FRAG_L,FRAG_R': None,
                                              'SINGLE,1,2,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1BrH],[1FH]': None,
                                              'SINGLE,2,1,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1FH],[1BrH]': None,
                                              'DOUBLE,1,2,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1BrH],[1FH]': None,
                                              'DOUBLE,2,1,19,Clc1ccc(c2csc([n]3[n]c([1cH]c3C)C)[n]2)cc1,[1FH],[1BrH]': None}

        # write test data to temp file (smi)
        for smi_id, smi in list(self.test_dataset_goldeninput_smi_01.items()):
            self.temp_file_input_smi_01.write(smi + " " + smi_id + "\n")
        self.temp_file_input_smi_01.close()

        # write test data to temp file (smi)
        for smi_id, smi in list(self.test_dataset_goldeninput_smi_03.items()):
            self.temp_file_input_smi_03.write(smi + " " + smi_id + "\n")
        self.temp_file_input_smi_03.close()

        # container for results data
        self.test_dataset_testresults = {}

    def tearDown(self):

        """Tear down object for clean reuse in further tests"""
        # clean out the object
        self.test_mmp_mcss_object.clean_out_data()
        # clean out the temp data store
        self.test_dataset_testresults.clear()

        os.remove(self.temp_file_input_smi_01.name)

    def test_get_largest_mcs_pairs_with_diff(self):
        """Test method to get largest MCS MMP for given smi - smi pair"""

        # 6. full build then write of pairs to file, but only for a single named column
        self.test_mmp_mcss_object.build_from_dicer(self.temp_file_input_smi_01.name, 'BOTH', 'NONE')
        self.test_mmp_mcss_object.enumerate_fragment_properties()
        self.test_mmp_mcss_object.get_largest_mcs_pairs(self.temp_file_output_pairs.name, 'BOTH')

        # now read it back into temp object and check it's what we wrote out!
        test_results_filehandle = open(self.temp_file_output_pairs.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None
        test_results_filehandle.close()

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_golden_output_01, self.test_dataset_testresults)

    def test_get_largest_mcs_pairs_mdc_atm_hard(self):
        """Test method to get largest MCS MMP for given smi - smi pair"""

        # 6. full build then write of pairs to file, but only for a single named column
        self.test_mmp_mcss_object.build_from_dicer(self.temp_file_input_smi_01.name, 'BOTH', 'NONE')
        self.test_mmp_mcss_object.enumerate_fragment_properties()
        self.test_mmp_mcss_object.get_largest_mcs_pairs(self.temp_file_output_pairs.name, 'BOTH', mdc_atm_hard=4)

        # now read it back into temp object and check it's what we wrote out!
        test_results_filehandle = open(self.temp_file_output_pairs.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None
        test_results_filehandle.close()

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_golden_output_02, self.test_dataset_testresults)

    def test_get_largest_mcs_pairs_mdc_atm_soft(self):
        """Test method to get largest MCS MMP for given smi - smi pair"""

        # 6. full build then write of pairs to file, but only for a single named column
        self.test_mmp_mcss_object.build_from_dicer(self.temp_file_input_smi_03.name, 'BOTH', 'NONE')
        self.test_mmp_mcss_object.enumerate_fragment_properties()

        #
        self.test_mmp_mcss_object.get_largest_mcs_pairs(self.temp_file_output_pairs.name, 'BOTH')
        # now read it back into temp object and check it's what we wrote out!
        test_results_filehandle = open(self.temp_file_output_pairs.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None
        test_results_filehandle.close()

        self.test_mmp_mcss_object.get_largest_mcs_pairs(self.temp_file_output_pairs.name, 'BOTH', mdc_atm_soft=3,
                                                        mdc_atm_soft_threshold=4)
        # now read it back into temp object and check it's what we wrote out!
        test_results_filehandle = open(self.temp_file_output_pairs.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None
        test_results_filehandle.close()

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_golden_output_03, self.test_dataset_testresults)


if __name__ == '__main__':
    unittest.main()
