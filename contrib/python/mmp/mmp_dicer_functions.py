###################################################################
""" Summary: Builds dicer commend line for MMP generation JAL

Notes:
 -G iso01 (force deprioritisation of isomeric numbering in canonicalisation)
 -B nosmi (removes need for grep -v "B=" | grep -v "FRAGID")
 -c (discard chirality)
 -B lostchiral (check each fragment for lost chirality)
 -m 0 (include H cuts)
 -k 2 (single and double cuts, 1 = single only)
 -i ICTE will skip connection table errors rather than fail silently half way through with exit code 0
 
 -B atype=sb (write smiles and complementary smiles with atom type labels on break points i.e.: attachment points)
 -B MAXFF= (discard fragments that comprise more than <f> fraction of atoms in parent)
 -C auto (write smiles and complementary smiles with auto label on break points)
 -s "ClC" -s "BrC" -s "FC" (smarts for cut points - breaks bond between 1st and 2nd matched atom)
 -B addq (run the -q queries in addition to the default rules)
 -B bscb (allow Carbon-Carbon single bonds to break)
 -M 15 (discard fragments with more than 15 atoms)
 -X 5000 (maximum number of fragments per molecule to produce)
 -A I (enable input of aromatic structures)
 -A D (use Daylight aromaticity definitions)

 Generates fragment and complimentary frag - can still get the same cansmi with different isomeric numbering
  [1xx]xxxxx[2xx] or [2xx]xxxxx[1xx]
 For complimentary part (context) - Guarantee first fragment is labelled 1 always, second fragment is 2:
  [1xx].xxxxx[2xx]xxxx

"""
##########################################################################
import sys
import os
import logging
from subprocess import Popen, PIPE

import unittest
import tempfile


def build_mmpdicer_cmd(smi_fi, cut_type, threshold, return_threshold=False):

    try:
        home_dir = os.environ['C3TK_HOME']

    except KeyError:
        home_dir = os.environ['LILLYMOL_HOME']

    except KeyError:
        sys.exit("Failed to identiy home dir, please set C3TK_HOME or LILLYMOL_HOME")

    build_dir = 'Linux'

    try:
        build_dir = os.environ['BUILD_DIR']
    except KeyError:
        pass

    root_dir = home_dir + '/bin/' + build_dir

    dicer_binary = root_dir + '/dicer'

    try:
        # default was 0.50001
        # smaller value like 0.3 will result in less fragmentation (smaller fragments produced)
        maxff = float(threshold)
    except:
        sys.exit('Cannot cast threshold value to float: %s' % threshold)

    if threshold < 0.01 or threshold > 0.999:
        sys.exit('Invalid input to dicer command line (build_mmpdicer_cmd) please try >= 0.1 and <= 0.9')

    # full command line, see comment above
    dicer_cmd = '-B atype=sb -B MAXFF=' + str(maxff) + \
                ' -C auto -s "ClC" -s "BrC" -s "FC" -B addq -B bscb -m 0 -M 15 -X 5000 -i smi -A I -A D' \
                ' -B nosmi -G iso01 -c -i ICTE'

    if cut_type.upper() == 'DOUBLE':
        dicer_cut_opt = '-k 2'

    elif cut_type.upper() == 'BOTH':
        dicer_cut_opt = '-k 2'
        
    elif cut_type.upper() == 'SINGLE':
        dicer_cut_opt = '-k 1'
        
    else:
        dicer_cut_opt = '-k 1'

    if return_threshold:
        dicer_cmd += " -B appnatoms"

    dicer_full_cmd = dicer_binary + " " + dicer_cmd + " " + dicer_cut_opt + " " + smi_fi

    return dicer_full_cmd


def invoke_dicer_subprocess(dicer_cmd, bufsize):

    return Popen(dicer_cmd, stdout=PIPE, stderr=PIPE, shell=True, bufsize=-1)


def subprocess_dicer_communicate(dicer_cmd, bufsize):
    
    proc = invoke_dicer_subprocess(dicer_cmd, bufsize)

    while proc.returncode is None:
        
        (stdout, stderr) = proc.communicate()
        for line in stdout.decode('utf-8').splitlines():
            
            yield line

    if proc.returncode != 0:

        sys.exit("Please check your input SMI, the dicer code failed with exit code %s" % str(proc.returncode))


def execute_dicer(smi_fi, cut_type, threshold, logger_object, return_threshold=False):

        logger = logger_object

        if len(logging.Logger.manager.loggerDict) < 1:
            # exit with system status 1 and custom error
            sys.exit("Invalid or no logger object passed to MMPObjectClass.  Please create \
                    and pass a logger and set to use logging.disable if you don't want logging")

        dicer_cmd = build_mmpdicer_cmd(smi_fi, cut_type, threshold, return_threshold=return_threshold)

        logger.info('Attempting execution of dicer:\n %s' % dicer_cmd)

        for line in subprocess_dicer_communicate(dicer_cmd, 1):
            
            # logger.debug('Parsing Dicer output...')

            line = line.strip('\r')
            line = line.strip('\n')
            line_list = line.split()
            
            #
            # expect to see all lines in following format:
            # [1CH4] --some id-- AT=[1:C2] --some smiles-- COMP --some id-- AT=[1:C3]
            if return_threshold:
                try:
                    mol_id = int(line_list[1])
                    ctx_orig = str(line_list[5].lstrip('0'))
                    frag = str(line_list[0])
                    # these next two items are the fragment and context attachment points
                    fattach = str(line_list[4]).replace("AT=", "")
                    cattach = str(line_list[8]).replace("AT=", "")
                    cut_threshold_frag = float(line_list[3])
                    cut_num_atoms_frag = float(line_list[2])
                    cut_threshold_ctx = float(line_list[10])
                    cut_num_atoms_ctx = float(line_list[9])
                except:
                    logger.warn('Error parsing Dicer output (with threshold), unexpected line format: %s' % line)
                    try:
                        mol_id = int(line_list[1])
                    except:
                        logger.warn('Looks like your Input ID is not an integer')
                    continue
            else:
                try:
                    mol_id = int(line_list[1])
                    ctx_orig = str(line_list[3].lstrip('0'))
                    frag = str(line_list[0])
                    # these next two items are the fragment and context attachment points
                    fattach = str(line_list[2]).replace("AT=", "")
                    cattach = str(line_list[6]).replace("AT=", "")
                except:
                    logger.warn('Error parsing Dicer output, unexpected line format: %s' % line)
                    try:
                        mol_id = int(line_list[1])
                    except:
                        logger.warn('Looks like your Input ID is not an integer')
                    continue
            
            # now separate out single from double cuts
            # we can do this by counting the number of dot disconnected smiles in the context smi string
            num_cuts = ctx_orig.count('.') + 1

            if return_threshold:
                yield num_cuts, mol_id, ctx_orig, frag, fattach, cattach, \
                      cut_num_atoms_frag, cut_threshold_frag, cut_num_atoms_ctx, cut_threshold_ctx
            else:
                yield num_cuts, mol_id, ctx_orig, frag, fattach, cattach


class _TestMMPDicerFunctions(unittest.TestCase):

    """Test class for Dicer execution
    """
    
    def setUp(self):

        self.maxDiff = None

        self.temp_file_input_smi = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_output_data = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')

        self.mmplogger = logging.getLogger('mmpobjectclass_testlogger')
        logging.disable(logging.CRITICAL)

        self.test_dataset_input_smi_01 = {
                # All the following test data is from CHEMBL
                # CHEMBL2105127 https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL2105127/
                '2105127': 'CC1(CC(=O)N(CN2CCOCC2)C1=O)c3ccccc3',
                # CHEMBL3989502 https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL3989502/
                '3989502': 'CN1[C@@H]2CC[C@H]1C[C@@H](C2)OC3c4ccccc4N(C)S(=O)(=O)c5ccccc35',
                }
        
        self.test_dataset_goldenoutput_fragmentedsmols = {
                 (1, 2105127, 'O=C1N(CN2CCOCC2)C(=O)C[1CH]1c1ccccc1', '[1CH4]', '[1:C3]', '[1:C3]'): None,
                 (1, 2105127, 'CC1(c2ccccc2)CC(=O)[1NH]C1=O', '[1CH3]N1CCOCC1', '[1:NAM]', '[1:C3]'): None,
                 (1, 2105127, '[1CH3]N1C(=O)C(c2ccccc2)(CC1=O)C', 'O1CC[1NH]CC1', '[1:C3]', '[1:N3]'): None,
                 (1, 2105127, 'C[1CH]1CC(=O)N(CN2CCOCC2)C1=O', '[1cH]1ccccc1', '[1:C3]', '[1:CAR]'): None,
                 (2, 2105127, '[1CH3]N1CCOCC1.[2cH]1ccccc1', 'C[2CH]1CC(=O)[1NH]C1=O', '[2:CAR|1:C3]', '[1:NAM|2:C3]'): None,
                 (2, 2105127, 'CC1(c2ccccc2)CC(=O)[1NH]C1=O.O1CC[2NH]CC1', '[12CH4]', '[2:N3|1:NAM]', '[1:C3|2:C3]'): None,
                 (2, 2105127, 'O1CC[1NH]CC1.[2cH]1ccccc1', '[1CH3]N1C(=O)[2CH](CC1=O)C', '[1:N3|2:CAR]', '[1:C3|2:C3]'): None,
                 (1, 2105127, '[1CH3]C1(c2ccccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2ccccc2)[1CH2]C(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2ccccc2)CC(=O)N([1CH2]N2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2ccccc2)CC(=O)N(CN2[1CH2]COCC2)C1=O', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2ccccc2)CC(=O)N(CN2CCO[1CH2]C2)C1=O', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2[1cH]cccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2c[1cH]ccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 2105127, 'CC1(c2cc[1cH]cc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1S(=O)(=O)c2c(C(OC3CC4[1NH]C(C3)CC4)c3c1cccc3)cccc2', '[1CH4]', '[1:N3]', '[1:C3]'): None,
                 (2, 3989502, '[1CH4].[2OH]C1c2c(N(S(=O)(=O)c3c1cccc3)C)cccc2', '[1NH]1C2C[2CH2]CC1CC2', '[1:C3|2:O3]', '[1:N3|2:C3]'): None,
                 (2, 3989502, '[1CH4].CN1S(=O)(=O)c2c([2CH2]c3c1cccc3)cccc2', '[2OH]C1CC2[1NH]C(C1)CC2', '[2:C3|1:C3]', '[1:N3|2:O3]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c([1NH]S(=O)(=O)c5c3cccc5)cccc4)CC1CC2', '[1CH4]', '[1:N3]', '[1:C3]'): None,
                 (1, 3989502, '[1OH]C1c2c(N(S(=O)(=O)c3c1cccc3)C)cccc2', 'CN1C2C[1CH2]CC1CC2', '[1:O3]', '[1:C3]'): None,
                 (2, 3989502, 'CN1C2C[1CH2]CC1CC2.CN1S(=O)(=O)c2c([2CH2]c3c1cccc3)cccc2', '[12OH2]', '[2:C3|1:C3]', '[1:O3|2:O3]'): None,
                 (1, 3989502, 'CN1S(=O)(=O)c2c([1CH2]c3c1cccc3)cccc2', '[1OH]C1CC2N(C)C(C1)CC2', '[1:C3]', '[1:O3]'): None,
                 (1, 3989502, '[1CH3]N1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1[1CH]2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2C[1CH2]C1CC(OC1c3c(N(S(=O)(=O)c4c1cccc4)C)cccc3)C2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)[1CH2]C1CC2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2C[1CH](OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(O[1CH]3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4[1cH]cccc4N(S(=O)(=O)c4c3cccc4)C)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cc[1cH]c4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)c[1cH]cc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(S(=O)(=O)N(C)c5[1cH]cccc35)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, '[1CH3]N1S(=O)(=O)c2c(C(OC3CC4N(C)C(C3)CC4)c3c1cccc3)cccc2', '[1H]', '[1:C3]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5[1cH]cccc35)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cc[1cH]c5)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3c[1cH]cc5)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None,
                 (1, 3989502, 'CN1C2CC(OC3c4[1cH]cccc4S(=O)(=O)N(c4c3cccc4)C)CC1CC2', '[1H]', '[1:CAR]', '[1:H]'): None}

        self.test_dataset_goldenoutput_fragmentedsmols_withthreshold = {
            (1, 2105127, 'O=C1N(CN2CCOCC2)C(=O)C[1CH]1c1ccccc1', '[1CH4]', '[1:C3]', '[1:C3]', 1.0, 0.04761905, 20.0,
              0.952381): None, (
             1, 2105127, 'CC1(c2ccccc2)CC(=O)[1NH]C1=O', '[1CH3]N1CCOCC1', '[1:NAM]', '[1:C3]', 7.0, 0.3333333, 14.0,
             0.6666667): None, (
             1, 2105127, '[1CH3]N1C(=O)C(c2ccccc2)(CC1=O)C', 'O1CC[1NH]CC1', '[1:C3]', '[1:N3]', 6.0, 0.2857143, 15.0,
             0.7142857): None, (
             1, 2105127, 'C[1CH]1CC(=O)N(CN2CCOCC2)C1=O', '[1cH]1ccccc1', '[1:C3]', '[1:CAR]', 6.0, 0.2857143, 15.0,
             0.7142857): None, (
             2, 2105127, '[1CH3]N1CCOCC1.[2cH]1ccccc1', 'C[2CH]1CC(=O)[1NH]C1=O', '[2:CAR|1:C3]', '[1:NAM|2:C3]', 8.0,
             0.3809524, 13.0, 0.6190476): None, (
             2, 2105127, 'CC1(c2ccccc2)CC(=O)[1NH]C1=O.O1CC[2NH]CC1', '[12CH4]', '[2:N3|1:NAM]', '[1:C3|2:C3]', 1.0,
             0.04761905, 20.0, 0.952381): None, (
             2, 2105127, 'O1CC[1NH]CC1.[2cH]1ccccc1', '[1CH3]N1C(=O)[2CH](CC1=O)C', '[1:N3|2:CAR]', '[1:C3|2:C3]', 9.0,
             0.4285714, 12.0, 0.5714286): None, (
             1, 2105127, '[1CH3]C1(c2ccccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2ccccc2)[1CH2]C(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2ccccc2)CC(=O)N([1CH2]N2CCOCC2)C1=O', '[1H]', '[1:C3]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2ccccc2)CC(=O)N(CN2[1CH2]COCC2)C1=O', '[1H]', '[1:C3]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2ccccc2)CC(=O)N(CN2CCO[1CH2]C2)C1=O', '[1H]', '[1:C3]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2[1cH]cccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2c[1cH]ccc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 2105127, 'CC1(c2cc[1cH]cc2)CC(=O)N(CN2CCOCC2)C1=O', '[1H]', '[1:CAR]', '[1:H]', 1.0, 0.04761905, 21.0,
             1.0): None, (
             1, 3989502, 'CN1S(=O)(=O)c2c(C(OC3CC4[1NH]C(C3)CC4)c3c1cccc3)cccc2', '[1CH4]', '[1:N3]', '[1:C3]', 1.0,
             0.03571429, 27.0, 0.9642857): None, (
             2, 3989502, '[1CH4].[2OH]C1c2c(N(S(=O)(=O)c3c1cccc3)C)cccc2', '[1NH]1C2C[2CH2]CC1CC2', '[1:C3|2:O3]',
             '[1:N3|2:C3]', 8.0, 0.2857143, 20.0, 0.7142857): None, (
             2, 3989502, '[1CH4].CN1S(=O)(=O)c2c([2CH2]c3c1cccc3)cccc2', '[2OH]C1CC2[1NH]C(C1)CC2', '[2:C3|1:C3]',
             '[1:N3|2:O3]', 9.0, 0.3214286, 19.0, 0.6785714): None, (
             1, 3989502, 'CN1C2CC(OC3c4c([1NH]S(=O)(=O)c5c3cccc5)cccc4)CC1CC2', '[1CH4]', '[1:N3]', '[1:C3]', 1.0,
             0.03571429, 27.0, 0.9642857): None, (
             1, 3989502, '[1OH]C1c2c(N(S(=O)(=O)c3c1cccc3)C)cccc2', 'CN1C2C[1CH2]CC1CC2', '[1:O3]', '[1:C3]', 9.0,
             0.3214286, 19.0, 0.6785714): None, (
             2, 3989502, 'CN1C2C[1CH2]CC1CC2.CN1S(=O)(=O)c2c([2CH2]c3c1cccc3)cccc2', '[12OH2]', '[2:C3|1:C3]',
             '[1:O3|2:O3]', 1.0, 0.03571429, 27.0, 0.9642857): None, (
             1, 3989502, 'CN1S(=O)(=O)c2c([1CH2]c3c1cccc3)cccc2', '[1OH]C1CC2N(C)C(C1)CC2', '[1:C3]', '[1:O3]', 10.0,
             0.3571429, 18.0, 0.6428571): None, (
             1, 3989502, '[1CH3]N1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1[1CH]2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2C[1CH2]C1CC(OC1c3c(N(S(=O)(=O)c4c1cccc4)C)cccc3)C2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)[1CH2]C1CC2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2C[1CH](OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(O[1CH]3c4c(N(S(=O)(=O)c5c3cccc5)C)cccc4)CC1CC2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4[1cH]cccc4N(S(=O)(=O)c4c3cccc4)C)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)cc[1cH]c4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cccc5)C)c[1cH]cc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(S(=O)(=O)N(C)c5[1cH]cccc35)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, '[1CH3]N1S(=O)(=O)c2c(C(OC3CC4N(C)C(C3)CC4)c3c1cccc3)cccc2', '[1H]', '[1:C3]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5[1cH]cccc35)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3cc[1cH]c5)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4c(N(S(=O)(=O)c5c3c[1cH]cc5)C)cccc4)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None, (
             1, 3989502, 'CN1C2CC(OC3c4[1cH]cccc4S(=O)(=O)N(c4c3cccc4)C)CC1CC2', '[1H]', '[1:CAR]', '[1:H]', 1.0,
             0.03571429, 28.0, 1.0): None
        }

        for smi_id, smi in list(self.test_dataset_input_smi_01.items()):
            self.temp_file_input_smi.write(smi+" "+smi_id+"\n")
        self.temp_file_input_smi.close()

        # container for results data
        self.test_dataset_testresults = {}

    def test_execute_dicer(self):
        """Test the dicer execution"""

        for cut_type, mol_id, ctx_orig, frag, fattach, cattach in execute_dicer(self.temp_file_input_smi.name,
                                                                                'BOTH',
                                                                                0.50001,
                                                                                self.mmplogger):
            self.test_dataset_testresults[cut_type, mol_id, ctx_orig, frag, fattach, cattach] = None

        # print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_fragmentedsmols)

    def test_execute_dicer_with_threshold(self):
        """Test dicer execution with return threshold"""

        for cut_type, mol_id, ctx_orig, frag, fattach, cattach, na_frag, thr_frag, na_ctx, thr_ctx in \
                execute_dicer(self.temp_file_input_smi.name,
                              'BOTH',
                              0.50001,
                              self.mmplogger,
                              return_threshold=True):
            self.test_dataset_testresults[cut_type, mol_id, ctx_orig, frag, fattach, cattach,
                                          na_frag, thr_frag, na_ctx, thr_ctx] = None

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_fragmentedsmols_withthreshold)

    def test_execute_dicer_with_threshold2(self):
        """Test dicer execution with return threshold"""

        for cut_type, mol_id, ctx_orig, frag, fattach, cattach, na_frag, thr_frag, na_ctx, thr_ctx in \
                execute_dicer(self.temp_file_input_smi.name, 'BOTH', 0.7, self.mmplogger, return_threshold=True):
            self.test_dataset_testresults[cut_type, mol_id, ctx_orig, frag, fattach, cattach,
                                          na_frag, thr_frag, na_ctx, thr_ctx] = None

        # more as we allow larger dicer fragments, default 37
        self.assertEqual(len(self.test_dataset_testresults), 41)

    def test_execute_dicer_with_threshold3(self):
        """Test dicer execution with return threshold"""

        for cut_type, mol_id, ctx_orig, frag, fattach, cattach, na_frag, thr_frag, na_ctx, thr_ctx in \
                execute_dicer(self.temp_file_input_smi.name, 'BOTH', 0.3, self.mmplogger, return_threshold=True):
            self.test_dataset_testresults[cut_type, mol_id, ctx_orig, frag, fattach, cattach,
                                          na_frag, thr_frag, na_ctx, thr_ctx] = None

        # less as we restrict dicer to smaller fragments, default 37
        self.assertEqual(len(self.test_dataset_testresults), 31)

if __name__ == '__main__':
    unittest.main()
