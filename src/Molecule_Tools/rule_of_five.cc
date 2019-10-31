/*
  Computes Pfizer rule of five.
  (Advanced Drug Delivery Reviews, 23 (1997) 3-25)
  If input is from clogp, we can do clogp as well
*/

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <limits>
using namespace std;

#include "cmdline.h"
#include "accumulator.h"
#include "misc.h"

#include "qry_wstats.h"
#include "output.h"
#include "istream_and_type.h"
#include "iwstandard.h"
#include "aromatic.h"
#include "nvrtspsa.h"
#include "misc2.h"

static Chemical_Standardisation chemical_standardisation;

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int reduce_to_largest_fragment = 1;    // always

static int append_violations_to_name = 0;

static int reading_input_from_biobyte_clogp = 0;

static int include_psa = 0;

static extending_resizable_array<int> feature_count;

/*
  Some times we just want all molecules written with an annotation of their score
*/

static Molecule_Output_Object stream_for_all_molecules;

static int function_as_filter = 0;

/*
  Statistics on the number of donors and acceptors
*/

static Accumulator_Int<int> acceptor_count_accumulator;
static extending_resizable_array<int> acceptor_count;
static Accumulator_Int<int> donor_count_accumulator;
static extending_resizable_array<int> donor_count;
static Accumulator<float> clogp_accumulator;

/*
  We can also do the Pfizer rule of 5 stuff
*/

static Molecule_Output_Object Pfizer_druglike;
static Molecule_Output_Object Pfizer_non_druglike;

static molecular_weight_t upper_molecular_weight_cutoff = 500.0;

static float upper_clogp_cutoff = 5.0;

/*
  How many violations before the molecule is considered non-druglike
*/

static int violation_cutoff = 2;

static int violation_count[5] = {0, 0, 0, 0, 0,};

static int druglike_molecules_found = 0;
static int non_druglike_molecules_found = 0;

static IWString ro5_tag("RO5<");
static IWString ro5_violation_tag("VIOLATES<");
static IWString ro5_violation_count_tag("VIOLATIONS<");

static int primary_amines_contribute = 1;

/*
  It will help if we can group all the attributes together
*/

class RO5_Attributes
{
  private:
    int _natoms;
    molecular_weight_t _amw;
    int _acount;
    int _dcount;
    float _clogp;

    int _number_violations;

    IWString _string_form;

    float _nvrtspsa;

//  private functions

    int _initialise (Molecule &, int *);
    int _identify_violations ();

  public:
    RO5_Attributes();

    int initialise (Molecule &);

    void set_clogp (float c);

    int write_descriptors (const IWString &, IWString_and_File_Descriptor &);

    int number_violations () const { return _number_violations;}
    const IWString & violations_in_string_form () const { return _string_form;}
};

RO5_Attributes::RO5_Attributes()
{
  _natoms = 0;
  _clogp = static_cast<float>(-100.0);
  _amw = static_cast<molecular_weight_t>(0.0);

  _number_violations = 0;

  return;
}

int
RO5_Attributes::write_descriptors (const IWString & mname,
                                   IWString_and_File_Descriptor & output)
{
  if (0 == mname.length())
    output << "RO5UNK" << molecules_read;
  else
    append_first_token_of_name(mname, output);

  output << ' ' << _dcount << ' ' << _acount << ' ' << _natoms << ' ';

  if (_amw > static_cast<molecular_weight_t>(0.0))
    output << _amw;
  else
    output << '.';

  if (_clogp > -100.0)
    output << ' ' << _clogp;
  else
    output << " .";

  output << ' ' << _number_violations;

  if (include_psa)
    output << ' ' << _nvrtspsa;

  output << '\n';

  return output.good();
}

int
RO5_Attributes::_identify_violations ()
{
  assert (_natoms > 0);

  _number_violations = 0;
  _string_form.resize(0);

  if (_clogp > upper_clogp_cutoff)
  {
    if (append_violations_to_name)
      _string_form << " CP " << _clogp;
    _number_violations++;
    if (verbose > 1)
      cerr << " Clogp out of range " << _clogp << endl;
  }

  if (_amw > upper_molecular_weight_cutoff)
  {
    if (append_violations_to_name)
      _string_form << " AMW " << _amw;
    _number_violations++;
    if (verbose > 1)
      cerr << " Molecular weight = " << _amw << endl;
  }

  if (_dcount > 5)
  {
    if (append_violations_to_name)
      _string_form << " donor " << _dcount;
    _number_violations++;
    if (verbose > 1)
      cerr << " Contains " << _dcount << " donors\n";
  }

  if (_acount > 10)
  {
    if (append_violations_to_name)
      _string_form << " acceptor " << _acount;
    _number_violations++;
    if (verbose > 1)
      cerr << " Contains " << _acount << " acceptors\n";
  }

  violation_count[_number_violations]++;

  return _number_violations;
}

static int
assign_acceptors_simplistic (Molecule & m, int * isotope)
{
  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);

    if (8 == ai->atomic_number() || 7 == ai->atomic_number())
    {
      isotope[i] = 1;
      rc++;
    }
  }

  return rc;
}

/*
  In the simplistic form, a donor is an O or N with an attached hydrogen
*/

static int
assign_donors_simplistic (Molecule & m, int * isotope)
{
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);

    if ((8 == ai->atomic_number() || 7 == ai->atomic_number()) && m.hcount(i))
    {
      if (0 == isotope[i])
        isotope[i] = 3;
      else
        isotope[i] = 2;

      rc++;
    }
  }

  return rc;
}

static int
assign_donors (Molecule & m, int * isotope)
{
  return assign_donors_simplistic(m, isotope);
}

static int
assign_acceptors (Molecule & m, int * isotope)
{
  return assign_acceptors_simplistic(m, isotope);
}

int
RO5_Attributes::initialise (Molecule & m)
{
  _natoms = m.natoms();

  if (0 == _natoms)
  {
    cerr << "RO5_Attributes::initialise: cannot process empty molecule '" << m.name() << "'\n";
    return 0;
  }

  _amw = m.molecular_weight();

  int * isotope = new_int(_natoms); std::unique_ptr<int[]> free_isotope(isotope);

  return _initialise(m, isotope);
}

void
RO5_Attributes::set_clogp (float c)
{
  _clogp = c;

  if (c > upper_clogp_cutoff)
    _number_violations++;

  return;
}

int
RO5_Attributes::_initialise (Molecule & m,
                             int * isotope)
{
  _acount = assign_acceptors(m, isotope);

  acceptor_count_accumulator.extra(_acount);
  acceptor_count[_acount]++;

  _dcount = assign_donors(m, isotope);

  donor_count_accumulator.extra(_dcount);
  donor_count[_dcount]++;

  int matched_atoms = _acount + _dcount;

  feature_count[matched_atoms]++;

  if (verbose > 1)
    cerr << molecules_read << " '" << m.molecule_name() << "' " << _acount << " acceptors " << _dcount << " donors\n";

  if (include_psa)
    _nvrtspsa = novartis_polar_surface_area(m);

  _identify_violations();

  return 1;
}

static int
do_append_violations_to_name (Molecule & m,
                              int number_violations,
                              const RO5_Attributes & ro5)
{
  IWString tmp = m.molecule_name();
  tmp << " PR5:" << number_violations << ro5.violations_in_string_form();

  m.set_name(tmp);

  return 1;
}

/*
  Return 1 for a pass, and 0 for a failure
*/

static int
check_pass_fail (Molecule & m,
                 int number_violations)
{
  int rc;
  if (number_violations >= violation_cutoff)
  {
    if (Pfizer_non_druglike.active())
      Pfizer_non_druglike.write(m);

    non_druglike_molecules_found++;
    rc = 0;
  }
  else
  {
    if (Pfizer_druglike.active())
      Pfizer_druglike.write(m);

    druglike_molecules_found++;
    rc = 1;
  }

  return rc;
}

static int
rule_of_five (Molecule & m,
              RO5_Attributes & ro5,
              IWString_and_File_Descriptor & output)
{
  int number_violations = ro5.number_violations();

  if (append_violations_to_name)
    do_append_violations_to_name(m, number_violations, ro5);

  if (stream_for_all_molecules.active())
    stream_for_all_molecules.write(m);

  if (number_violations >= violation_cutoff)
  {
    non_druglike_molecules_found++;
    if (Pfizer_non_druglike.active())
      Pfizer_non_druglike.write(m);
  }
  else
  {
    druglike_molecules_found++;
    if (Pfizer_druglike.active())
      Pfizer_druglike.write(m);
  }

  if (function_as_filter)
  {
    if (number_violations < violation_cutoff)
      output << m.smiles() << ' ' << m.name() << '\n';

    return 1;
  }

  return ro5.write_descriptors(m.name(), output);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment && m.number_fragments() > 1)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    (void) chemical_standardisation.process(m);

  return;
}

static int
compute_rule_of_five_attributes (Molecule & m, 
                                 RO5_Attributes & ro5)
{
  molecules_read++;

  preprocess(m);

// Now start the stuff we came here for

  if (! ro5.initialise(m))
  {
    cerr << "Fatal error processing '" << m.name() << "'\n";
    return 0;
  }

  return 1;
}

static int
rule_of_five (data_source_and_type<Molecule> & input,
              IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    unique_ptr<Molecule> free_m(m);

    RO5_Attributes ro5;
    if (! compute_rule_of_five_attributes(*m, ro5))
      return 0;

    if (! rule_of_five(*m, ro5, output))
      return 0;
  }

  return 1;
}

/*
  Input looks like   
  3.719  0 SMILES ID
*/

static int
rule_of_five_record (const const_IWSubstring & buffer,
                     IWString_and_File_Descriptor & output)
{
  const_IWSubstring token;
  int i = 0;

  if (! buffer.nextword(token, i))
  {
    cerr << "Cannot extract logp\n";
    return 0;
  }

  float clogp;
  if (! token.numeric_value(clogp))
  {
    cerr << "Invalid logp '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "Cannot extract errlvl\n";
    return 0;
  }

// Get smiles and id here so we have the id for reporting errors

  const_IWSubstring smiles, id;

  if (! buffer.nextword(smiles, i) || ! buffer.nextword(id, i))
  {
    cerr << "No smiles/id data\n";
    return 0;
  }

  int errlvl;
  if ('0' == token)
    errlvl = 0;
  else if (! token.numeric_value(errlvl))
  {
    cerr << "Invalid clogp errlvl '" << token << "'\n";
    return 0;
  }
  else if (errlvl < 0)
    cerr << "Complete failure for clogp computation on '" << id << "', beware erroneous 0.0 value\n";
  else if (errlvl > 59)
    cerr << "Warning, possibly invalid logp '" << id << "', errlvl " << errlvl << " clogp " << clogp << " molecule '" << id << "'\n";

  Molecule m;

  if (! m.build_from_smiles(smiles))
  {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return 0;
  }

  m.set_name(id);

  RO5_Attributes ro5;
  if (! compute_rule_of_five_attributes(m, ro5))
    return 0;

  ro5.set_clogp(clogp);
//cerr << "Set logp to " << clogp << endl;

  return rule_of_five(m, ro5, output);
}

static int
rule_of_five (iwstring_data_source & input,
              IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! rule_of_five_record(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return output.good();
}

static int
rule_of_five (const char * fname,
              IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return rule_of_five(input, output);
}

static int
rule_of_five (const char * fname,
              int input_type,
              IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  return rule_of_five(input, output);
}

static int
open_output_object (const Command_Line & cl, Molecule_Output_Object & mstream,
                    char cflag,
                    const char * text)
{
  assert (cl.option_present(cflag));
  assert (0 == mstream.number_elements());

  if (cl.option_present('o'))
  {
    if (! mstream.determine_output_types(cl))
    {
      cerr << "Cannot determin output type(s)\n";
      return 46;
    }
  }
  else
    mstream.add_output_type(SMI);

  const_IWSubstring fname;
  (void) cl.value(cflag, fname);

  if (! mstream.new_stem(fname, 1))     // 2nd arg means keep existing suffix
  {
    cerr << "Cannot set '" << cflag << " stream name to '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << text << " written to '" << fname << "'\n";

  return 1;
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << endl;

  cerr << "Pfizer fule of five computation\n";
  cerr << "Usage: " << prog_name << " <options> <input_file> <input_file> ...\n";
  cerr << "  -P <fname>     output molecules OK with Pfizer rule of 5\n";
  cerr << "  -p <fname>     output molecules failing Pfizer rule of 5\n";
  cerr << "  -W <number>    specify upper molecular weight cutoff\n";
  cerr << "  -z <number>    specify number of violations for rejection (default " << violation_cutoff << ")\n";
  cerr << "  -z APPEND      append violations to molecule name\n";
  cerr << "  -O <value>     specify upper clogp cutoff value (default " << upper_clogp_cutoff << ")\n";
  cerr << "  -L <name>      specify name for writing labelled molecules\n";
  cerr << "  -S <name>      specify name for writing all molecules\n";
  cerr << "  -B             input from Biobyte clogp\n";
  cerr << "  -N             include Novartis PSA value in ouput\n";
  cerr << "  -m             NH2 groups contribute 2 donors (in accord with Lipinski)\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type\n";
  cerr << "  -E ...         standard Element options\n";
  cerr << "  -A ...         standard aromaticity options\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

int
rule_of_five (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:o:S:g:t:P:p:lW:z:O:BNfm");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl))
    usage(2);

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl))
      usage(3);
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      usage(6);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will strip multi-fragment molecules to largest fragment\n";
  }

#ifdef NOT_USED_QUQUQU
// They must be getting some kind of output

  if (cl.option_present('P') || cl.option_present('p'))   // rule of 5 pass/fail individual files
    ;
  else if (cl.option_present('S'))   // all molecules
    ;
  else if (cl.option_present('n'))   // output suppressed
    ;
  else
  {
    Pfizer_druglike.add_output_type(SMI);
    if (! Pfizer_druglike.new_stem("-"))
    {
      cerr << "Cannot direct passing molecules to stdout\n";
      return 8;
    }

    if (verbose)
      cerr << "Writing druglike molecules only to stdout\n";
  }
#endif

  if (cl.option_present('O'))
  {
    if (! cl.value('O', upper_clogp_cutoff))
    {
      cerr << "Cannot parse -O option\n";
      usage(17);
    }

    if (verbose)
      cerr << "Molecules with clogp > " << upper_clogp_cutoff << " will fail the rule of 5\n";
  }

  if (cl.option_present('m'))
  {
    primary_amines_contribute = 2;

    if (verbose)
      cerr << "Primary amines will count as two donors\n";
  }

  int input_type = 0;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++))
    {
      if ("APPEND" == z)
      {
        append_violations_to_name = 1;
        if (verbose)
          cerr << "Rule of five violations appended to the name\n";
      }
      else if (! z.numeric_value(violation_cutoff) || violation_cutoff < 1 || violation_cutoff > 5)
      {
        cerr << "The -z option must be followed by a positive number '" << z << "' invalid\n";
        usage(8);
      }
      else if (verbose)
      {
        cerr << "Any molecule with " << violation_cutoff << " or more violations will be non-druglike\n";
      }
    }
  }

  if (cl.option_present('W'))
  {
    if (! cl.value('W', upper_molecular_weight_cutoff) || upper_molecular_weight_cutoff < 12.0)
    {
      cerr << "The -W option must be followed by a positive number\n";
      usage(8);
    }

    if (verbose)
      cerr << "Upper molecular weight cutoff for rule of 5 is " << upper_molecular_weight_cutoff << endl;
  }

  if (cl.option_present('B'))
  {
    reading_input_from_biobyte_clogp = 1;

    if (verbose)
      cerr << "Input from Biotype clogp\n";
  }

  if (cl.option_present('f'))
  {
    function_as_filter = 1;

    if (verbose)
      cerr << "Will work as a filter\n";
  }

  if (cl.option_present('N'))
  {
    include_psa = 1;

    if (verbose)
      cerr << "Will also include Novartis PSA values\n";
  }

  if (function_as_filter && include_psa)
  {
    cerr << "Sorry, the -f and -N options are mutually incompatible, ignoring -N option\n";
    include_psa = 0;
  }

  if (cl.option_present('S'))
  {
    if (! append_violations_to_name)
    {
      append_violations_to_name = 1;
      cerr << "Will automatically append violation information to the output\n";
    }

    if (! open_output_object(cl, stream_for_all_molecules, 'S', "All molecules"))
      return 6;
  }

  if (cl.option_present('P'))
  {
    if (! open_output_object(cl, Pfizer_druglike, 'P', "Pfizer drug-like"))
      return 6;
  }

  if (cl.option_present('p'))
  {
    if (! open_output_object(cl, Pfizer_non_druglike, 'p', "Pfizer non-drug-like"))
      return 6;
  }

  if (reading_input_from_biobyte_clogp)
    ;
  else if (function_as_filter)
    ;
  else if (0 == input_type && ! all_files_recognised_by_suffix(cl))
  {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  acceptor_count.resize(100);
  donor_count.resize(100);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  if (! function_as_filter)
  {
    output << "Name ro5_nhoh ro5_no natoms amw clogp violations";

    if (include_psa)
      output << " nvrtspsa";
    output << '\n';
  }

  int rc = 0;
  if (reading_input_from_biobyte_clogp)
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! rule_of_five(cl[i], output))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! rule_of_five(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    if (0 == molecules_read)
      return rc;

    cerr << druglike_molecules_found << " passed Pfizer rule of five, " << non_druglike_molecules_found << " did not\n";
    for (int i = 0; i < 5; i++)
    {
      cerr << violation_count[i] << " molecules had " << i << " violations\n";
    }

    cerr << "Molecules had between " << acceptor_count_accumulator.minval() << " and " << acceptor_count_accumulator.maxval() << " acceptors";
    if (acceptor_count_accumulator.n())
      cerr << ", ave " << acceptor_count_accumulator.average();
    cerr << endl;

    for (int i = 0; i < acceptor_count.number_elements(); i++)
    {
      if (acceptor_count[i])
        cerr << acceptor_count[i] << " molecules had " << i << " acceptors\n";
    }

    cerr << "Molecules had between " << donor_count_accumulator.minval() << " and " << donor_count_accumulator.maxval() << " donors";
    if (donor_count_accumulator.n())
      cerr << ", ave " << donor_count_accumulator.average();
    cerr << endl;

    for (int i = 0; i < donor_count.number_elements(); i++)
    {
      if (donor_count[i])
        cerr << donor_count[i] << " molecules had " << i << " donors\n";
    }

    if (clogp_accumulator.n())
    {
      cerr << "clogp values between " << clogp_accumulator.minval();
      if (clogp_accumulator.n() > 1)
        cerr << " and " << clogp_accumulator.maxval() << " ave " << clogp_accumulator.average();
      cerr << endl;
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rule_of_five(argc, argv);

  return rc;
}
