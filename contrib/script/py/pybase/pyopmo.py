"""
Contains Python "bindings" to molecule object and Python functions
utilizing molecule object key functionalities
The intent is to provide a Pythonic interface to mo utilities and, thus,
enable easy/consistent use of mo from within python programs.

Testing:
$ python <path-to-location>/pymo.py
"""

import os
import subprocess as sp
import shlex
import shutil
import sys
import argparse
import logging
import unittest
from tempfile import NamedTemporaryFile, mkdtemp

log = logging.getLogger('lilly.' + __name__)
log.addHandler(logging.NullHandler())

build_dir = 'Linux'

try:
      build_dir = os.environ['BUILD_DIR']
except EnvironmentError:
     pass

home_dir = os.environ['LILLYMOL_HOME']
root_dir = home_dir + '/bin/' + build_dir

# dictionary of commands that will be turned into functions
# function_name: [script_name, debug_message, default_params_dict]
mo_tool_map = {
    'dicer': [root_dir + '/dicer',
              'Recursively cuts molecules based on queries',
              {}],
    'fileconv': [root_dir + '/fileconv',
                 'file/molecule conversion utilities',
                 {}],
    'iwdescr': [root_dir + '/iwdescr',
                'compute physicochemical descriptors using iwdescr',
                {'-l': ''}],
    'make_these_molecules': [root_dir + '/make_these_molecules',
                             'Makes molecules from isotopically labelled '
                             'reagents according to make file, not '
                             'combinatorial join',
                             {}],
    'preferred_smiles': [root_dir + '/preferred_smiles',
                         '',
                         {}],
    'alogp': [root_dir + '/abraham',
            '',
            {'-l': '',
             '-F': home_dir + '/contrib/data/queries/abraham/Abraham',
             '-P': home_dir + '/contrib/data/queries/abraham/Alpha2H',
             '-g': 'all' }]
}

def __function_generator(fct_name, 
                         script_name, 
                         debug_msg,
                         default_params_dict,
                         expect_zero=True):
    """
    A generator for functions which runs one of LillyMol scripts with a
    predefined set of parameters.

    :param str fct_name: the function name of the newly generated function
    :param script_name: your LillyMol script path (from mo_tool above)
    :param debug_msg: quick message about what the script does (from mo_tool above)
    :param default_params_dict: default parameters
    :param expect_zero: whether to expect zero as a return value from the script
    :return: a function which you can call to run the script
    :rtype: callable
    """

    def funct(infile, outfile=None, params_dict=default_params_dict,
              params_input_string=None, loggero=None, pretend=False,
              altret=False, inpipe=False):
        # Use param_input_string when a dictionary can not be used
        # e.g. tsubstructure.sh -A M -A D

        log.debug('funct = %s: ' % script_name + debug_msg)

        params_string = ' '

        if params_input_string:
            params_string += params_input_string + ' '

        for k, v in list(params_dict.items()):
            if type(v) is list:
                for vv in v:
                    params_string += k + ' ' + vv + ' '
            else:
                params_string += k + ' ' + str(v) + ' '
        params_string = params_string[:-1]

        if type(infile) is list:
            infile_s = ' '.join(infile)
        elif infile is None:
            infile_s = ''
        else:
            infile_s = infile

        cmd_2_execute = script_name + params_string + ' '

        if not inpipe:
            cmd_2_execute += infile_s

        out_fh = None

        if altret:
            out_fh = sp.PIPE
        elif outfile:
            out_fh = open(outfile, 'w')

        if pretend:
            log.warning('Just pretending')
            exit_status = 0
        else:
            cmd = shlex.split(cmd_2_execute)
            log.info('Executing: {}'.format(cmd_2_execute))
            #print('Executing: {}'.format(cmd_2_execute))

            # FIXME: this gets really ugly now...
            if inpipe:
                my_proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=out_fh,
                                   stderr=sp.PIPE, shell=False)

                # FIXME: Python2
                out, err = my_proc.communicate(infile_s.encode('utf-8'))
            else:
                my_proc = sp.Popen(cmd, stdout=out_fh, stderr=sp.PIPE,
                                   shell=False)
                out, err = my_proc.communicate()

            exit_status = my_proc.returncode

            # NOTE: Python2 returns strings so need to check
            # FIXME: this needs to be simplified as soon as everything is
            #        Python3
            if type(err) == bytes:
                err = err.decode('utf-8')
            if type(out) == bytes:
                out = out.decode('utf-8')

        if outfile:
            out_fh.close()

        if exit_status and expect_zero:
            log.error("%s failed:\n%s" % (script_name, cmd_2_execute))

            if err:
                log.error(err)
        else:
            log.debug("Done: " + debug_msg)

        # very ugly workaround for this mis-designed function
        if not altret:
            return exit_status
        else:
            return exit_status, out, err

    funct.__name__ = fct_name
    funct.__doc__ = debug_msg

    return funct


for name, params in list(mo_tool_map.items()):
    nparams = len(params)

    if not (3 <= nparams <= 5):
        raise IndexError('mo_tool_map: "{}" has {:d} parameter(s) but should '
                         'have 3-5'.format(name, nparams))

    locals()[name] = __function_generator(name, *params)


def make_these_molecules(rgnt_list, make_these_file, reaction_file_list, outfile=None, params_dict={}, debug=True,
                         loggero=None):
    """
    Alternative to trxn, used in MMP code for generating new mols from MMP's, trxn version would be alternative:
    For connecting one bond (two components, single cut fragments):
        trxn.sh  -S - -r oneConnection.rxn partOne.smi partTwo.smi
    For connecting two bonds (three component, two cut fragments):
        trxn.sh  -S -rxn -S - -r twoConnection.rxn partThree.smi  partOne.smi partTwo.smi
    BUT, if we have a long list of different contexts (partOne) and don't want exhaustive enumeration, specify rxn's:
        make_these_molecules.sh -R oneConnection.rxn -M m2Make.txt -S - partOne.smi partTwo.smi
    In this case, you can put all context fragments SMILES (context1a, context 1b, ...)  in one reagent file, and
    all fragments SMILES (frag1, frag2, ...) in the second reagent file.  If you have something like (context1a frag1\n
    context1a frag2\ncontext1b frag3\n...) in your m2Make.txt file, you will create the molecules you wanted
    """

    log.debug("Generating virtual compounds using rxn and reagents supplied plus specified combinations file")

    # prep reagents file string
    rgnt_string = " ".join(rgnt_list)

    log.debug("These are the reagent files...." + str(rgnt_string))

    # prep params string
    params_string = " "
    for k, v in list(params_dict.items()):
        params_string += k + " " + v + " "
    params_string = params_string[:-1]

    # set outfile
    # improved a bit to handle files with '.' in main name, other than in the extension
    if outfile:
        if outfile[-4:] == ".smi" or outfile[-4:] == ".txt":
            params_string += " -S " + os.path.splitext(outfile)[0]
        else:
            params_string += " -S " + outfile

    reaction_file = ""
    for rxn_file in reaction_file_list:  # todo: if single string, this is split in characters
        reaction_file += ' -R ' + rxn_file

    cmd_line = (mo_tool_map['make_these_molecules'][0] + reaction_file +
                ' -M ' + make_these_file + " " + params_string + " " +
                rgnt_string)

    log.debug("Executing: %s" % cmd_line)
    #if debug:
        #print(cmd_line)

    my_proc = sp.Popen(shlex.split(cmd_line), stdout=None, stderr=sp.PIPE,
                       shell=False)

    for line in my_proc.stderr.readlines():
        log.debug(line.rstrip())

    exit_status = my_proc.wait()

    log.debug("Done generating compounds")

    return exit_status


#####################################################
class _TestPymo(unittest.TestCase):
    """Test class for pymo module

    Example usage:

     python pymo.py (to execute all tests)

     python pymo.py -c (for verbose console logging)
     python pymo.py -f mylog.log (for logging to file mylog.log)
     python pymo.py -c -f mylog.log (for both)

     python pymo.py _Test_pymo.test_fetch_smiles # (to execute only the specified test)

     coverage run pymo.py (to run test code coverage analysis)
     coverage report pymo.py (to view the result of the test code coverage analysis)
    """

    def setUp(self):
        """setup test data location, unittest config and logger"""

        # location of test data
        self.test_data_location = root_dir + '/contrib/script/py/mmp/testdata/'

        # temp output file and dir
        self.temp_inp_file = NamedTemporaryFile(encoding='utf-8', mode='wt', suffix='.smi', delete=False)
        self.temp_out_file = NamedTemporaryFile(encoding='utf-8', mode='wt', delete=False)
        self.temp_out_dir = mkdtemp()

        test_smiles = {
            # basic test set - all the below id's and structures are CHEMBL
            '3105327': 'Cc1ccc2c(ccn2c3nc(cs3)c4cc(ccc4F)C(F)(F)F)c1',
            '1526778': 'CC(=O)c1c(C)n(c(C)c1C(=O)C)c2nc(c(C)s2)c3ccc(C)c(C)c3',
            '1494678': 'CC(=O)c1c(C)n(c(C)c1C(=O)C)c2nc(c(C)s2)c3ccccc3',
            '472166': 'OC(CCn1ccnc1)(c2ccccc2)c3ccccc3',
            '69798': 'Cc1nccn1CCC(O)(c2ccccc2)c3ccccc3',
            '367346': 'Cc1sc(N)nc1c2cccc(Cl)c2',
            '366881': 'Cc1sc(N)nc1c2ccc(Cl)c(Cl)c2',
            '1477460': 'COc1ccc(cc1)c2nc(sc2C)n3c(C)c(C(=O)C)c(C(=O)C)c3C',
            '1441050': 'COc1ccc(cc1OC)c2nc(sc2C)n3c(C)c(C(=O)C)c(C(=O)C)c3C'
        }

        # write test data to temp file 01
        for smi_id, smi in test_smiles.items():
            string = smi+' '+smi_id+'\n'
            self.temp_inp_file.write(string)
        self.temp_inp_file.close()

    def tearDown(self):
        """cleanup test data and settings"""

        # Clean up the directory
        # os.removedirs(self.temp_out_dir)
        shutil.rmtree(self.temp_out_dir)

    def test_fileconv(self):
        log.debug("Testing fileconv")
        exit_status = fileconv(os.path.join(self.test_data_location, self.temp_inp_file.name),
                               params_dict={'-v': '',
                                            '-c': '10',
                                            '-C': '20',
                                            '-S': os.path.join(self.temp_out_dir, 'atomcountfilter'),
                                            '-o': 'smi'}
                              )
        log.debug("fileconv return code was: %s" % exit_status)
        self.assertEqual(exit_status, 0)

    def test_make_these_molecules(self):
        log.debug("Testing _make_these_molecules")

        # test data from mmp_enum_mols_from_pairs.py
        context = NamedTemporaryFile(encoding='utf-8', mode='wt', suffix='.smi', delete=False)
        examp_context = "[1CH3]CCCCCC partOne_1\n[1CH3]CCCCCCC partOne_2\n[1CH3]CCCCCCCC partOne_3\n[1CH3]CCCCCCCCC partOne_4\n[1CH3]CCCCCCCCCC partOne_5"
        context.write(examp_context)
        context.close()

        frags = NamedTemporaryFile(encoding='utf-8', mode='wt', suffix='.smi', delete=False)
        examp_frags = "[1OH]CCCCCC partTwo_1\n[1OH]CCCCCCC partTwo_2\n[1OH]CCCCCCCC partTwo_3\n[1OH]CCCCCCCCC partTwo_4\n[1OH]CCCCCCCCCC partTwo_5\n[1OH]CCCCCCCCCCC partTwo_6\n"
        frags.write(examp_frags)
        frags.close()

        make_instr = NamedTemporaryFile(encoding='utf-8', mode='wt', delete=False)
        examp_make_instr = "partOne_2 partTwo_3\npartOne_4 partTwo_5\n"
        make_instr.write(examp_make_instr)
        make_instr.close()

        rxn = NamedTemporaryFile(encoding='utf-8', mode='wt', suffix='.rxn', delete=False)
        # this is the reaction specification that trxn needs to combine isotopically labelled mmp fragmentation points
        single_rxn = "(0 Reaction\n (0 Scaffold\n  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n )\n"
        single_rxn += " (1 Sidechain\n  (A C smarts \"[!0*]\")\n  (A I isotope (0 0))\n  (A I join (0 0))\n )\n)"
        rxn.write(single_rxn)
        rxn.close()

        exit_status = make_these_molecules([context.name, frags.name],
                                           make_instr.name,
                                           [rxn.name])

        log.debug("make_these_molecules return code was: %s" % exit_status)
        self.assertEqual(exit_status, 0)

    def test_dicer(self):
        log.debug("Testing dicer")
        log.debug("Testing %s" % dicer.__name__)
        exit_status = dicer(os.path.join(self.test_data_location, self.temp_inp_file.name),
                               os.path.join(self.temp_out_dir, 'dicer.out'),
                               params_dict={}, )
        log.debug("dicer return code was: %s" % exit_status)
        self.assertEqual(exit_status, 0)

    def test_alogp(self):
        log.debug("Testing CMI")
        exit_status = alogp(os.path.join(self.test_data_location, self.temp_inp_file.name),
                               os.path.join(self.temp_out_dir, 'alogp.out')
                          )
        log.debug("CMI return code was: %s" % exit_status)
        self.assertEqual(exit_status, 0)

if __name__ == "__main__":

    # optional command line flags
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log_file',
                        help='Name of file to place debug log info in')
    parser.add_argument('-c', '--console',
                        help='Switch on console logging',
                        default=False,
                        action='store_true')

    args = parser.parse_args()
    #print(args)
    logger_file = args.log_file
    console_on = args.console

    pymotest_logger = logging.getLogger("pymo.testlogger")
    pymotest_logger.setLevel(logging.DEBUG)

    log_formatter = logging.Formatter("%(asctime)s [%(funcName)-12.12s] "
                                      "[%(levelname)-5.5s] %(message)s",
                                      datefmt="%Y-%m-%d %H:%M:%S")

    if console_on:
        print("Switched on console")
        h1 = logging.StreamHandler(stream=sys.stdout)
        h1.setLevel(logging.DEBUG)
        h1.setFormatter(log_formatter)
        pymotest_logger.addHandler(h1)
    else:
        print("console off")

    if logger_file is not None:
        print(("Switched on logging to file: {}".format(logger_file)))
        fileHandler = logging.FileHandler(filename=logger_file)
        fileHandler.setFormatter(log_formatter)
        fileHandler.setLevel(logging.DEBUG)
        pymotest_logger.addHandler(fileHandler)
    else:
        print("file logging off")

    if console_on is False and logger_file is None:
        pymotest_logger.setLevel(logging.CRITICAL)

    unittest.main()
