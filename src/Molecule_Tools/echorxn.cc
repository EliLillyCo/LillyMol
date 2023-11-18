#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/rxn_file.h"
#include "Molecule_Lib/molecule_to_query.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString output_filename;

static int
echorxn (IWReaction & rxn,
         const char * fname)
{
  std::ofstream output;

  if (output_filename.length())
  {
    output.open (output_filename.null_terminated_chars());
  }
  else
  {
    IWString newfname = fname;

    if (newfname.ends_with(".rxn"))
      newfname.chop(4);

    newfname += "_echo.rxn";

    output.open(newfname.null_terminated_chars());
  }

  if (! output.good())
  {
    cerr << "Looks like I cannot open the output file\n";
    return 0;
  }

  return rxn.write_msi(output);
}

static int
echorxn (const char * fname, int query_files_in_reaction_directory)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open reaction file '" << fname << "'\n";
    return 0;
  }

  IWReaction rxn;

  if (query_files_in_reaction_directory)
    rxn.set_query_files_in_current_directory(0);

  Sidechain_Match_Conditions sidechain_match_conditions;
  if (! rxn.do_read(fname, sidechain_match_conditions))
  {
    cerr << "Cannot read reaction from '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
  {
    cerr << "Reaction from '" << fname << "'\n";
    rxn.debug_print(cerr);
  }

  if (! rxn.ok())
  {
    cerr << "Yipes, ok fails for newly created reaction!!!\n";
  }

  return echorxn(rxn, fname);
}

static int
do_echo_isis_reaction_file(RXN_File & ISIS_rxn,
                            int remove_fragments)
{
  if (remove_fragments)
    ISIS_rxn.set_remove_product_fragments(1);
 
  RXN_File_Create_Reaction_Options rxnfcro;
  Molecule_to_Query_Specifications mqs;

  IWReaction rxn;

  if (! ISIS_rxn.create_reaction(rxn, rxnfcro, mqs))
  {
    cerr << "OOPs, could not convert intermediate representation to reaction object\n";
    return 5;
  }

  if (verbose > 1)
    rxn.debug_print(cerr);

  IWString output_name;
  if (output_filename.length() > 0)
  {
    output_name = output_filename;
    if (! output_name.ends_with(".rxn"))
      output_name << ".rxn";
  }
  else
  {
    output_name = ISIS_rxn.fname();

    if (output_name.ends_with(".rxn"))
      output_name.chop(4);

    output_name << "_echo.rxn";
  }

  std::ofstream output(output_name.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "Cannot open output file '" << output_name << "'\n";
    return 0;
  }

  return rxn.write_msi(output);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage: " << prog_name << " <options> <reaction1> <reaction2> ...\n";

  cerr << "  -h             query files in same directory as reaction files\n";
  cerr << "  -D <fname>     ISIS reaction file\n";
  cerr << "  -S <fname>     specify file name stem for echo'd reaction files\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
echorxn(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vD:S:h");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  int query_files_in_reaction_directory = 0;
  if (cl.option_present('h'))
  {
    query_files_in_reaction_directory = 1;
    if (verbose)
      cerr << "Query files in same directory as reaction file\n";
  }

  if (cl.option_present('S'))
  {
    if (1 == cl.number_elements())
      ;
    else if (0 == cl.number_elements() && cl.option_present('D'))
      ;
    else
    {
      cerr << "When used with the -S option, only one reaction file is allowed\n";
      usage(4);
    }

    cl.value('S', output_filename);
    if (verbose)
      cerr << "Echo'd reaction will be written to '" << output_filename << "'\n";
  }

  if (cl.option_present('D'))
  {
    set_auto_create_new_elements(1);

    RXN_File ISIS_rxn;
    if (! parse_isis_rxn_file_options(cl, 'D', ISIS_rxn))
    {
      cerr << "Cannot parse -D option\n";
      usage(4);
    }

    int remove_fragments = cl.option_present('l');

    do_echo_isis_reaction_file(ISIS_rxn, remove_fragments);

    return 0;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(5);
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! echorxn(cl[i], query_files_in_reaction_directory))
      return i + 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = echorxn(argc, argv);

  return rc;
}
