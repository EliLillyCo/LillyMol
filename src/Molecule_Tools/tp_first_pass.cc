/*
  First phase for 3rd party processing.

  Can do the following:   (although driven by command line options)

    Remove any molecule in which no fragment contains more than 11 atoms.
    Remove any molecule with covalently bound non-organic atoms
    Remove any molecule containing atoms types we don't want
    Remove any molecule which does not contain C and (O or N)
    Remove hydrogens.
    Apply all chemical standardisations

    Sept 99
    We want to remove molecules in which every bond is a ring bond. I'm
    generalising that to use a ring_bond_ratio. If the ratio falls below
    a threshold, the molecule is rejected.
*/

#include <iostream>
#include <iostream>
#include <memory>
#include <limits>
using std::cerr;
using std::endl;
#include <assert.h>

#include "cmdline.h"
#include "iwstring.h"
#include "accumulator.h"
#include "misc.h"

#include "molecule.h"
#include "smiles.h"
#include "aromatic.h"
#include "path.h"
#include "misc2.h"
#include "iwstandard.h"

#include "istream_and_type.h"
#include "output.h"

/*
  The programme can operate in the mode of just scanning its input
  and not writing anything
*/

static int molecules_read = 0;
static int molecules_written = 0;

static Chemical_Standardisation chemical_standardisation;

static int skip_molecules_with_abnormal_valences = 0;
static int molecules_with_abnormal_valences = 0;

static int exclude_isotopes = 0;
static int molecules_containing_isotopes = 0;

static int convert_isotopes = 0;
static int isotopes_converted = 0;

static int molecules_containing_colvalent_non_organics = 0;
static int molecules_containing_non_allowed_atom_types = 0;

static int exclude_molecules_containing_non_periodic_table_elements = 1;
static int molecules_containing_non_periodic_table_elements = 0;

static int run_all_checks = 0;

/*
  One could imagine a poorly defined counterion entered as [X]. We can
  optionally allow non-periodic table elements if they aren't connected
  to anything
*/

static int allow_non_periodic_table_elements_if_not_connected = 0;

static double ring_bond_ratio = 2.0;
static int molecules_with_bad_ring_bond_ratios = 0;

static Accumulator_Int<int> natoms_accumulator;

// Note that these are applied differently from fileconv

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = 0;

int molecules_below_atom_count_cutoff = 0;
int molecules_above_atom_count_cutoff = 0;

static int exclude_molecules_with_no_interesting_atoms = 1;
static int molecules_with_no_interesting_atoms = 0;

static double min_fraction_interesting_atoms = 0.0;
static int molecules_with_not_enough_interesting_atoms = 0;

static int reject_if_fragements_with_this_many_atoms = std::numeric_limits<int>::max();
static int mixtures_rejected = 0;

static int append_rejection_reason_to_name = 0;

static int write_rejection_reason_like_tsubstructure = 0;

/*
  We can make the output easier to parse if we have a fixed string before 
  the rejection reason
*/

static IWString prepend_before_reason;

/*
  We don't want strange characters in molecule names
*/

static char translate_non_printing_chars = '\0';

static int remove_non_printing_chars = 0;

/*
  With the -B option, one specify the number of connection table errors
  allowed before programme exit.
*/

static int connection_table_errors_allowed = 0;

static IWString output_file_stem;

/*
  We can append any arbitrary text to the name of each molecule
*/

static IWString text_to_append;

static int lower_ring_count_cutoff = 0;
static int molecules_with_too_few_rings = 0;
static int upper_ring_count_cutoff = 0;
static int molecules_with_too_many_rings = 0;

static int upper_ring_size_cutoff = 0;
static int molecules_with_ring_sizes_out_of_range = 0;

/*
  We can optionally assign each molecule written a number R(number).
*/

#include "numass.h"

static Number_Assigner number_assigner;

#include "rmele.h"

static Elements_to_Remove elements_to_remove;

/*
  We can transform element types.
*/

#include "etrans.h"

static Element_Transformations element_transformations;

static Molecule_Output_Object rejections_output_object;

const char *prog_name;

static int verbose = 0;

static void
usage (int rc = 1)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "usage: " << prog_name << " -i <input type> -o <output type> file1 file2...\n";
  cerr << "The following options are recognised\n";
  cerr << "  -c <number>    exclude molecules with atom count below <number>\n";
  cerr << "  -C <number>    exclude molecules with atom count above <number>\n";
  cerr << "  -r <number>    omit molecules with fewer than <number> rings\n";
  cerr << "  -R <number>    omit molecules with more than <number> rings\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -Z <ringsize>  upper ring size limit\n";
  cerr << "  -I <0,1>       <exclude,include> molecules containing isotopes\n";
  cerr << "  -I change      change all isotopic atoms to their natural form\n";
  cerr << "  -s discard     discard all chiral information in the input\n";
  cerr << "  -s good        ignore erroneous chiral input\n";
  cerr << "  -s 1           include chiral info in output (default)\n";
  cerr << "  -s 0           exclude chiral info from output\n";
//cerr << "  -n <number>    assign sequential numbers R(%d) starting with <number>\n";
  cerr << "  -b <ratio>     skip molecules with ring bond ratio's >= than <ratio>\n";
  cerr << "  -w             run all checks - normally discards molecules once problem found\n";
  cerr << "  -V             skip any molecule with abnormal valences\n";
  cerr << "  -k             allow molecules having no \"interesting\" atoms to pass\n";
  cerr << "  -f <fraction>  minimum fraction of interesting atoms required\n";
  cerr << "  -y             allow non periodic table elements if they are not connected\n";
  cerr << "  -L <fname>     write rejected molecules to <fname>\n";
  cerr << "  -a             append rejection reason to name in reject log\n";
  cerr << "  -u             write the rejection reason like tsubstructure\n";
  cerr << "  -K <text>      prepend 'text' before rejection reason\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -P <text>      append <text> to the name of every molecule\n";
  cerr << "  -p remove      remove non printing characters in molecule names\n";
  cerr << "  -p <char>      convert non printing characters in molecule names to '<char>'\n";
  cerr << "  -x <natoms>    reject if molecule contains > 1 fragment with at least <natoms> atoms\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol> (use -E '*' for auto create)\n";
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  (void) display_standard_etrans_options (cerr, 't');
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -v             verbose output\n";

  exit (rc);
}

/*
  An unfortunate problem with the variable molecules_containing_isotopes.
  It only gets incremented for this molecule in the upper level functions
  if we are excluding isotopes.
*/

static int
do_convert_isotopes (Molecule & m)
{

  int rc = m.transform_to_non_isotopic_form ();

  if (rc)
    molecules_containing_isotopes++;

  return rc;
}

static int
interesting_atoms (Molecule & m)
{
  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;

  int matoms = m.natoms ();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number (i);
    if (6 == z)
    {
      carbon = 1;
      if (nitrogen || oxygen)
        return 1;
    }
    else if (7 == z)
    {
      nitrogen = 1;
      if (carbon)
        return 1;
    }
    else if (8 == z)
    {
      oxygen = 1;
      if (carbon)
        return 1;
    }
  }
  
  return 0;
}

static int
count_interesting_atoms (Molecule & m,
                         const int * include_atom,
                         const int flag)
{
  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;

  const auto matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (NULL != include_atom && flag != include_atom[i])
      continue;

    atomic_number_t z = m.atomic_number (i);
    if (6 == z)
      carbon++;
    else if (7 == z)
      nitrogen++;
    
    else if (8 == z)
      oxygen++;
  }

  if (0 == oxygen && 0 == nitrogen)
    return 0;

  if (0 == carbon)
    return 0;
  
  return nitrogen + oxygen;
}

/*
  We can optionally append the reason for rejection to the molecule name
*/

static IWString rejection_reason;

static int
exclude_for_no_interesting_atoms (Molecule & m)
{
  int nf = m.number_fragments ();
  if (1 == nf)
  {
    if (interesting_atoms (m))
      return 0;    // do not reject this molecule
    else
    {
      rejection_reason = "no_interesting_atoms";
      return 1;    // yes, reject this molecule
    }
  }

  resizable_array_p<Molecule> components;
  m.create_components (components);

  int largest_frag = 0;
  int interesting_atoms_in_largest_frag = 0;

  for (int i = 0; i < nf; i++)
  {
    Molecule * c = components[i];

    int catoms = c->natoms ();

    assert (catoms == m.atoms_in_fragment (i));

    if (catoms < largest_frag)
      continue;

    interesting_atoms_in_largest_frag = interesting_atoms (*c);
    largest_frag = catoms;
  }

  if (interesting_atoms_in_largest_frag)
    return 0;     // do not exclude it
  else
  {
    rejection_reason = "no_interesting_atoms";
    return 1;     // yes, exclude this molecule
  }
}

static int
exclude_for_too_few_interesting_atoms (Molecule & m)
{
  const int nf = m.number_fragments();

  if (1 == nf)
  {
    const int interesting_atoms = count_interesting_atoms(m, NULL, 0);
    if (static_cast<double>(interesting_atoms) / static_cast<double>(m.natoms()) >= min_fraction_interesting_atoms)
      return 0;                // not rejected
    rejection_reason = "too few interesting atoms";
    return 1;          // yes, exclude this molecule
  }

  const int matoms = m.natoms();

  int * f = new int[matoms]; std::unique_ptr<int[]> free_f(f);

  m.fragment_membership(f);

  int atoms_in_largest_fragment = 0;
  int interesting_atoms_in_largest_fragment = 0;

  for (int i = 0; i < nf; ++i)
  {
    const auto interesting_atoms = count_interesting_atoms(m, f, i+1);
    const auto aif = m.atoms_in_fragment(i);

    if (aif > atoms_in_largest_fragment)
    {
      atoms_in_largest_fragment = aif;
      interesting_atoms_in_largest_fragment = interesting_atoms;
    }
    else if (aif == atoms_in_largest_fragment && interesting_atoms > interesting_atoms_in_largest_fragment)
      interesting_atoms_in_largest_fragment = interesting_atoms;
  }

  if (static_cast<double>(interesting_atoms_in_largest_fragment) / static_cast<double>(atoms_in_largest_fragment) >= min_fraction_interesting_atoms)
    return 0;            // do not reject this molecule

  rejection_reason = "Too_few_interesting_atoms";

  return 1;        // yes, reject this molecule
}

static int ok_elements[HIGHEST_ATOMIC_NUMBER + 1];

static int
exclude_for_atom_type (const Molecule & m)
{

  int matoms = m.natoms ();
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi (i);

    if (1 == a->atomic_number () && 1 != a->ncon ())
    {
      if (verbose > 1)
        cerr << "Contains two valent hydrogen\n";
      molecules_containing_non_allowed_atom_types++;
      rejection_reason = "Two_valent_Hydrogen";
      return 1;    // yes, exclude this molecule
    }

    const Element * e = a->element ();

    atomic_number_t z = e->atomic_number ();

    if (e->organic ())
      continue;


    if (exclude_molecules_containing_non_periodic_table_elements && ! e->is_in_periodic_table ())
    {
      if (allow_non_periodic_table_elements_if_not_connected && 0 == a->ncon ())
      {
        if (verbose > 1)
          cerr << "Allowing singly connected non-periodic table element '" << e->symbol () << "'\n";
        continue;
      }

      if (verbose > 1)
        cerr << "Contains non periodic table atom '" << e->symbol () << "'\n";
      molecules_containing_non_periodic_table_elements++;
      rejection_reason = "non_periodic_table_atom";
      return 1;    // yes, exclude this molecule
    }

    if (! ok_elements[z])
    {
      if (verbose > 1)
        cerr << "Contains non-allowed atom '" << e->symbol () << "'\n";
      molecules_containing_non_allowed_atom_types++;
      rejection_reason = "non_allowed_atom";
      return 1;    // yes, exclude this molecule
    }

    if (a->ncon () > 0)
    {
      if (verbose > 1)
        cerr << "Contains covalently bound non-organic '" << e->symbol () << "'\n";
      molecules_containing_colvalent_non_organics++;
      rejection_reason = "covalent_non-organic";
      return 1;    // yes, exclude this molecule
    }
  }

  return 0;    // do not exclude
}

static int
process_fragments (Molecule & m)
{
  int atoms_in_largest_fragment = 0;

  int number_large_fragments = 0;

  int nf = m.number_fragments ();
  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment (i);

    if (aif >= reject_if_fragements_with_this_many_atoms)
      number_large_fragments++;

    if (aif > atoms_in_largest_fragment)
      atoms_in_largest_fragment = aif;
  }

  if (number_large_fragments > 1)
  {
    if (verbose > 1)
      cerr << "Mixture " << number_large_fragments << " large fragments\n";
    rejection_reason = "mixture";
    mixtures_rejected++;
    return 0;
  }

  if (upper_atom_count_cutoff > 0 && atoms_in_largest_fragment > upper_atom_count_cutoff)
  {
    if (verbose > 1)
      cerr << "Too many atoms " << atoms_in_largest_fragment << " in largest fragment\n";
    rejection_reason = "too_many_atoms";
    molecules_above_atom_count_cutoff++;
    return 0;
  }

  if (atoms_in_largest_fragment < lower_atom_count_cutoff)
  {
    molecules_below_atom_count_cutoff++;
    rejection_reason = "not_enough_atoms";
    if (verbose > 1)
      cerr << "Too few atoms " << atoms_in_largest_fragment << " in " << nf << " fragments\n";
    return 0;
  }

  return 1;
}

static int
reject_for_ring_size_condition (Molecule & m,
                                int upper_ring_size_cutoff)
{
  for (int i = m.nrings() - 1; i >= 0; i--)
  {
    const Ring * ri = m.ringi(i);

    if (ri->number_elements() > upper_ring_size_cutoff)
      return 1;
  }

  return 0;
}

static int
reject_for_ring_bond_ratio (Molecule & m)
{
  if (ring_bond_ratio > 1.0)
    return 0;

  (void) m.ring_membership ();

  int nb = m.nedges ();

  int ring_bonds = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi (i);
    if (b->nrings ())
      ring_bonds++;
  }

  double ratio = static_cast<double>(ring_bonds) / static_cast<double>(nb);

//cerr << "My ratio " << ratio << " target " << ring_bond_ratio << endl;

  if (ratio >= ring_bond_ratio)
    return 1;     // return 1 means reject this molecule

  return 0;       // return 0 means don't reject this molecule
}

static int 
_apply_all_filters (Molecule & m)
{
  assert (m.ok ());

  rejection_reason.resize_keep_storage(0);

// Must do chemical standaridsation first because it may have
// explicit hydrogens, which messes up the atom count things

  (void) chemical_standardisation.process(m);

  int matoms = m.natoms ();

  if (0 == matoms)
    return 0;

  if (lower_atom_count_cutoff > 0 || upper_atom_count_cutoff > 0 ||
      std::numeric_limits<int>::max() != reject_if_fragements_with_this_many_atoms)
  {
    if (process_fragments (m))   // great, everything OK
      ;
    else if (run_all_checks)
      ;
    else
      return 0;
  }

  elements_to_remove.process (m);

  element_transformations.process (m);

  int keep = 0;

  if (exclude_for_atom_type (m))
  {
    if (verbose > 1)
      cerr << "Found bad atom types\n";
  }
  else if (exclude_molecules_with_no_interesting_atoms && exclude_for_no_interesting_atoms (m))
  {
    if (verbose > 1)
      cerr << "No interesting atoms\n";

    molecules_with_no_interesting_atoms++;
  }
  else if (min_fraction_interesting_atoms > 0.0 && exclude_for_too_few_interesting_atoms(m))
  {
    if (verbose > 1)
      cerr << "Too few interesting atoms\n";
    molecules_with_not_enough_interesting_atoms++;
  }
  else if (lower_ring_count_cutoff && 
           m.nrings () < lower_ring_count_cutoff)
  {
    molecules_with_too_few_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings () << " rings, which is below cutoff\n";
    rejection_reason << "not_enough_rings";
  }
  else if (upper_ring_count_cutoff &&
           m.nrings () > upper_ring_count_cutoff)
  {
    molecules_with_too_many_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings () << " rings, which is above cutoff\n";
    rejection_reason << "too_many_rings";
  }
  else if (exclude_isotopes && m.number_isotopic_atoms ())
  {
    molecules_containing_isotopes++;
    if (verbose > 1)
      cerr << "Molecule contains isotopes\n";
    rejection_reason << "isotopes";
  }
  else if (skip_molecules_with_abnormal_valences && ! m.valence_ok ())
  {
    molecules_with_abnormal_valences++;
    if (verbose > 1)
      cerr << "Molecule contains abnormal valence(s)\n";
    rejection_reason << "abnormal_valence";
  }
  else if (reject_for_ring_bond_ratio (m))
  {
    molecules_with_bad_ring_bond_ratios++;
    if (verbose > 1)
      cerr << "Ring bond ratio out of range\n";
    rejection_reason << "ring_bond_ratio";
  }
  else if (upper_ring_size_cutoff > 0 && reject_for_ring_size_condition (m, upper_ring_size_cutoff))
  {
    molecules_with_ring_sizes_out_of_range++;
    if (verbose > 1)
      cerr << "Ring size out of range\n";

    rejection_reason << "ring size";
  }
  else
    keep = 1;     // molecule is good!

  if (! keep)
    return 0;

  if (convert_isotopes)
    do_convert_isotopes (m);

  return 1;
}

static int
do_append_rejection_reason (Molecule & m,
                            const IWString & rejection_reason)
{
  IWString tmp = m.molecule_name ();
  tmp += ' ';

  if (write_rejection_reason_like_tsubstructure)
  {
    tmp << "(1 matches to '" << rejection_reason << "')";
  }
  else
  {
    if (prepend_before_reason.length ())
      tmp << prepend_before_reason;
    tmp += rejection_reason;

  }
  m.set_name (tmp);

  return 1;
}

static int
apply_all_filters (Molecule & m, int molecule_number)
{
  int rc = _apply_all_filters (m);

  if (run_all_checks && rejection_reason.length() > 0)
    ;
  else if (rc)      // molecule is OK.
    return rc;

  if (append_rejection_reason_to_name && rejection_reason.length() > 0)
    do_append_rejection_reason (m, rejection_reason);
  if (rejections_output_object.good ())
    rejections_output_object.write (m);

  return 0;
}

/*
  We return the number of non-printing characters in the name
*/

static int
non_printing_characters_in_name (Molecule & m)
{
  IWString mname = m.name ();

  int rc = 0;
  for (int i = mname.length () - 1; i >= 0; i--)
  {
    char c = mname[i];

    if (isprint (c))
      continue;

    if ('\0' != translate_non_printing_chars)
      mname[i] = translate_non_printing_chars;
    else if (remove_non_printing_chars)
    {
      mname.remove_item (i);
      i++;
    }

    rc++;
  }

  if (rc && ('\0' != translate_non_printing_chars || remove_non_printing_chars))
    m.set_name (mname);

  if (verbose > 1 && rc)
    cerr << "Removed/changed " << rc << " non printing chars in name\n";

  return rc;
}

static int
tp_first_pass (data_source_and_type<Molecule> & input,
          Molecule_Output_Object & output_object)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule ()))
  {

    std::unique_ptr<Molecule> free_m (m);

//  if (verbose > 1)
//    cerr << "Molecule " << input.molecules_read () << " finishes at line " << input.lines_read () << endl;

    if (verbose)
      natoms_accumulator.extra (m->natoms ());

    (void) non_printing_characters_in_name (*m);

    if (! apply_all_filters (*m, input.molecules_read ()))
      continue;

//  if (number_assigner.active ())
//    number_assigner.process (*m);

    if (text_to_append.length ())
      m->append_to_name (text_to_append);

    if (! output_object.write (m))
      return 0;
    
    molecules_written++;
  }

  return 0;
}

/*
*/

int
tp_first_pass (const char *fname,
               int input_type,
               Molecule_Output_Object & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name (fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good ())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1)
    input.set_verbose (1);

  if (connection_table_errors_allowed)
    input.set_connection_table_errors_allowed (connection_table_errors_allowed);

// Avoid name collisions before they occur

  if (output_file_stem.length ())
  {
    if (output.would_use_name (output_file_stem.null_terminated_chars (), fname))
    {
      cerr << "tp_first_pass: input '" << fname << "' and output stem '" << 
               output_file_stem << "' must not collide\n";
      return 4;
    }
  }
  else if (output.would_use_name (fname))
  {
    cerr << "tp_first_pass: input '" << fname << "' and output must be distinct\n";
    return 3;
  }

// Set up the output object for this file stem
// If there is a new stem for output files, make sure we get it.
// Make sure we deal properly with the case of multiple input files
// and a single output file (via the -S option)

  static int first_call = 1;

  int rc = 0;
  if (output_file_stem.length ())
  {
    if (first_call)
      rc = output.new_stem (output_file_stem);
    else
      rc = output.ok ();

    first_call = 0;
  }
  else
    rc = output.new_stem (fname);

  if (0 == rc)
  {
    cerr << "Output object could not open file\n";
    return 2;
  }

  rc = tp_first_pass (input, output);

  molecules_read += input.molecules_read ();

  return rc;
}

static int
tp_first_pass (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "aI:g:t:n:L:S:d:A:K:X:c:C:E:vVi:o:r:R:B:P:p:b:kyue:x:Z:wf:");

  verbose = cl.option_count ('v');

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  if (! process_elements (cl))
    usage (2);

  if (cl.option_present ('g'))
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose > 1, 'g'))
    {
      usage (6);
    }
  }

  if (! process_standard_aromaticity_options (cl, verbose))
  {
    usage (5);
  }

  if (cl.option_present ('t'))
  {
    if (! element_transformations.construct_from_command_line (cl, verbose, 't'))
      usage (8);
  }

  if (cl.option_present ('k'))
  {
    exclude_molecules_with_no_interesting_atoms = 0;

    if (verbose)
      cerr << "Will allow molecules with no interesting atoms to pass\n";
  }

  if (cl.option_present('f'))
  {
    if (! exclude_molecules_with_no_interesting_atoms)
    {
      cerr << "The -f and -k options cannot be used together\n";
      usage(1);
    }

    if (! cl.value('f', min_fraction_interesting_atoms) || min_fraction_interesting_atoms < 0.0 || min_fraction_interesting_atoms > 1.0)
    {
      cerr << "The minimum fraction of interesting atoms needed (-f) must be a valid fraction\n";
      usage(1);
    }

    if (verbose)
      cerr << "Molecules must have a minimum fraction " << min_fraction_interesting_atoms << " of interesting atoms\n";
  }

  if (cl.option_present ('y'))
  {
    allow_non_periodic_table_elements_if_not_connected = 1;

    if (verbose)
      cerr << "Non periodic table elements allowed if not connected\n";
  }

  int input_type = 0;
  if (cl.option_present ('i'))
  {
    if (! process_input_type (cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;

  if (0 == input_type && ! all_files_recognised_by_suffix (cl))
    return 4;

  Molecule_Output_Object output;

  if (! cl.option_present ('o'))
  {
    output.add_output_type (SMI);
  }
  else if (! output.determine_output_types (cl))
  {
    cerr << "Cannot determine output types\n";
    usage (8);
  }

  if (cl.option_present ('S'))
  {
    cl.value ('S', output_file_stem);

    if (output.would_overwrite_input_files(cl, output_file_stem))
    {
      cerr << "Cannot overwrite input file(s)\n";
      return 4;
    }

    if (verbose)
      cerr << "New files will be created with stem '" << output_file_stem << "'\n";
  }

// By default we echo any chiral info present in the input

  set_include_chiral_info_in_smiles (1);

  if (cl.option_present ('L'))
  {
    if (! cl.option_present ('o'))
      rejections_output_object.add_output_type (SMI);
    else if (! rejections_output_object.determine_output_types (cl))
    {
      cerr << "Cannot discern output types for rejections file\n";
      return 2;
    }

    IWString reject_log_file_name;
    cl.value ('L', reject_log_file_name);

    if (! rejections_output_object.new_stem (reject_log_file_name))
    {
      cerr << "Rejections file cannot use stem '" << reject_log_file_name << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Rejected structures will be written to '" << reject_log_file_name << "'\n";

    if (cl.option_present ('a'))
    {
      append_rejection_reason_to_name = 1;
      if (verbose)
        cerr << "The reason for rejection will be appended to the molecule name\n";
    }

    if (cl.option_present ('u'))
    {
      append_rejection_reason_to_name = 1;
      write_rejection_reason_like_tsubstructure = 1;
      if (verbose)
        cerr << "Rejection reasons written like tsubstructure\n";
    }
  }

  if (cl.option_present ('K'))
  {
    prepend_before_reason = cl.string_value ('K');

    append_rejection_reason_to_name = 1;

    if (verbose)
      cerr << "Will prepend '" << prepend_before_reason << "' before rejection reasons\n";

    if (! prepend_before_reason.ends_with (' '))
      prepend_before_reason << ' ';
  }

  if (cl.option_present ('X'))
  {
    if (! elements_to_remove.construct_from_command_line (cl, verbose, 'X'))
    {
      cerr << "Cannot discern elements to remove from -X switch\n";
      usage (18);
    }
  }

  if (! cl.option_present ('c'))
    ;
  else if (cl.value ('c', lower_atom_count_cutoff) && lower_atom_count_cutoff > 0)
  {
    if (verbose)
      cerr << "Will exclude molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }
  else
  {
    cerr << "Cannot discern lower atom count cutoff from " << cl.option_value ('c') << "'\n";
    usage (48);
  }

  if (! cl.option_present ('C'))
    ;
  else if (cl.value ('C', upper_atom_count_cutoff))
  {
    if (upper_atom_count_cutoff < lower_atom_count_cutoff)
    {
      cerr << "Upper atom count cutoff " << upper_atom_count_cutoff << 
              " must be greater than lower atom count cutoff " << lower_atom_count_cutoff << endl;
      usage (49);
    }
    if (verbose)
      cerr << "Will exclude molecules with more than " << upper_atom_count_cutoff << " atoms\n";
  }
  else
  {
    cerr << "Cannot discern upper atom count cutoff from '" << cl.option_value ('C') << "'\n";
    usage (50);
  }

//if (! number_assigner.initialise (cl, 'n', verbose))
//{
//  cerr << "Cannot process -n option\n";
//  usage (51);
//}

  if (cl.option_present ('r'))
  {
    if (! cl.value ('r', lower_ring_count_cutoff) ||
          lower_ring_count_cutoff < 1)
    {
      cerr << "-r option needs a whole number > 0\n";
      usage (52);
    }

    if (verbose)
      cerr << "Molecules containing fewer than " << lower_ring_count_cutoff <<
              " will be ignored\n";
  }

  if (cl.option_present ('R'))
  {
    if (! cl.value ('R', upper_ring_count_cutoff) ||
          upper_ring_count_cutoff < lower_ring_count_cutoff)
    {
      cerr << "-R option needs a whole number > " << lower_ring_count_cutoff << endl;
      usage (53);
    }

    if (verbose)
      cerr << "Molecules containing more than " << upper_ring_count_cutoff <<
              " rings will be ignored\n";
  }

  if (cl.option_present('w'))
  {
    run_all_checks = 1;

    if (verbose)
      cerr << "Will run all checks\n";
  }

  if (cl.option_present('Z'))
  {
    if (! cl.value('Z', upper_ring_size_cutoff) || upper_ring_size_cutoff < 3)
    {
      cerr << "The upper ring size cutoff value (-Z) must be a valid ring size\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard molecules with ring sizes > " << upper_ring_size_cutoff << endl;
  }

  if (cl.option_present ('V'))
  {
    skip_molecules_with_abnormal_valences = 1;
    if (verbose)
      cerr << "Molecules containing abnormal valences will be skipped\n";
  }

  if (cl.option_present ('I'))
  {
    int i = 0;
    IWString tmp;
    while (cl.value ('I', tmp, i++))
    {
      if ('0' == tmp)
      {
        exclude_isotopes = 1;
        if (verbose)
          cerr << "Molecules containing isotopes will be excluded\n";
      }
      else if ('1' == tmp)
      {
        exclude_isotopes = 0;
        if (verbose)
          cerr << "No action taken on molecules containing isotopes\n";
      }
      else if ("convert" == tmp)
      {
        convert_isotopes = 1;
        if (verbose)
          cerr << "All isotopic atoms will be converted to their non isotopic form\n";
      }
      else
      {
        cerr << "Unrecognised -I qualifier '" << tmp << "'\n";
        usage (55);
      }
    }
  }

  if (cl.option_present ('B'))
  {
    if (1 != cl.option_count ('B'))
    {
      cerr << "Only one -B option is allowed\n";
      usage (57);
    }

    if (! cl.value ('B', connection_table_errors_allowed) ||
          connection_table_errors_allowed < 0)
    {
      cerr << "The -B option requires a non negative integer argument\n";
      usage (34);
    }

    if (verbose)
      cerr << connection_table_errors_allowed << " connection table errors allowed\n";
  }

  if (cl.option_present ('P'))
  {
    cl.value ('P', text_to_append);
    if (verbose)
      cerr << "Will append '" << text_to_append << "' to the name of each molecule written\n";

    if (! text_to_append.starts_with (' '))
      text_to_append.insert_at_beginning (' ');
    
  }

  if (cl.option_present ('p'))
  {
    const_IWSubstring p;
    cl.value ('p', p);

    if ("remove" == p)
    {
      remove_non_printing_chars = 1;
      if (verbose)
        cerr << "Non printing characters removed from molecule names\n";
    }
    else if (1 == p.length ())
    {
      translate_non_printing_chars = p[0];
      if (verbose)
        cerr << "Non printing characters translated to '" << translate_non_printing_chars << "'\n";
    }
    else
    {
      cerr << "Invalid/unrecognised -p qualifier '" << p << "'\n";
      usage (6);
    }
  }

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', ring_bond_ratio) || ring_bond_ratio < 0.0 || ring_bond_ratio > 1.0)
    {
      cerr << "The lower ring bond ratio (-b) option must be followed by a valid fraction\n";
      usage (14);
    }

    if (verbose)
      cerr << "Molecules with ring bond ratio's of " << ring_bond_ratio << " or less are rejected\n";

  }

// initialise the array of allowable elements

  set_vector(ok_elements, HIGHEST_ATOMIC_NUMBER + 1, 0);

  ok_elements[6] = 1;
  ok_elements[7] = 1;
  ok_elements[8] = 1;
  ok_elements[9] = 1;
  ok_elements[15] = 1;
  ok_elements[16] = 1;
  ok_elements[17] = 1;
  ok_elements[35] = 1;
  ok_elements[53] = 1;
  ok_elements[3]  = 1;     // Li
  ok_elements[11] = 1;    // Na
  ok_elements[12] = 1;    // Mg
  ok_elements[19] = 1;    // K
  ok_elements[20] = 1;    // Ca

  if (cl.option_present('e'))
  {
    int i = 0;
    IWString e;

    while (cl.value('e', e, i++))
    {
      const Element * o = get_element_from_symbol_no_case_conversion(e);
      if (NULL == o)
      {
        cerr << "Sorry, non periodic table element '" << e << "', cannot be OK\n";
        return 5;
      }

      atomic_number_t z = o->atomic_number();

      assert (z >= 0 && z <= HIGHEST_ATOMIC_NUMBER);

      ok_elements[z] = 1;

      if (verbose)
        cerr << "Element " << e << " atomic number " << z << " allowed\n";

      if (5 == z || 14 == z)
      {
        Element * x = const_cast<Element *>(o);
        x->set_organic(1);
      }
    }
  }

  for (int i = 1; i < HIGHEST_ATOMIC_NUMBER; i++)
  {
    if (ok_elements[i])
      continue;

    const Element * e = get_element_from_atomic_number (i);
    if (! e->organic ())
      continue;

    ok_elements[i] = 1;
    if (verbose)
      cerr << "Element " << e->symbol () << " atomic number " << i << " allowed\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << prog_name << ": insufficient arguments " << argc << "\n";
    usage (56);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    const char *fname = cl[i];

    if (verbose)
      cerr << prog_name << " processing '" << fname << "'\n";

    rc += tp_first_pass (fname, input_type, output);
  }

  if (verbose)
  {
    if (cl.number_elements () > 1)
      cerr << molecules_read << " molecules read, ";
    cerr << molecules_written << " molecules written\n";

    if (0 == molecules_read)
      return rc;

    cerr << "Molecules had between " << natoms_accumulator.minval () << " and " <<
            natoms_accumulator.maxval () << " atoms\n";

    if (molecules_containing_non_allowed_atom_types)
      cerr << "Skipped " << molecules_containing_non_allowed_atom_types <<
        " molecules containing non allowed atoms\n";

    if (molecules_containing_non_periodic_table_elements)
      cerr << "Skipped " << molecules_containing_non_periodic_table_elements <<
        " molecules containing non periodic table atoms\n";

    if (molecules_containing_colvalent_non_organics)
      cerr << "Skipped " << molecules_containing_colvalent_non_organics <<
        " molecules containing covalently bonded non organics\n";

    if (molecules_with_no_interesting_atoms)
      cerr << "Skipped " << molecules_with_no_interesting_atoms <<
              " molecules with no interesting atoms\n";

    if (molecules_below_atom_count_cutoff)
      cerr << "Skipped " << molecules_below_atom_count_cutoff << 
              " molecules with atom count less than " << lower_atom_count_cutoff << endl;
    if (molecules_above_atom_count_cutoff)
      cerr << "Skipped " << molecules_above_atom_count_cutoff << 
              " molecules with atom count greater than " << upper_atom_count_cutoff << endl;

    if (molecules_with_too_few_rings)
      cerr << "Skipped " << molecules_with_too_few_rings << 
              " molecules having fewer than " << lower_ring_count_cutoff << " rings\n";
    if (molecules_with_too_many_rings)
      cerr << "Skipped " << molecules_with_too_many_rings << 
              " molecules having more than " << upper_ring_count_cutoff << " rings\n";
    if (molecules_with_ring_sizes_out_of_range)
      cerr << "Skipped " << molecules_with_ring_sizes_out_of_range << " molecules with rings containing more than " << upper_ring_size_cutoff << " atoms\n";
    
    if (molecules_containing_isotopes)
      cerr << molecules_containing_isotopes << " molecules contained isotopic atoms\n";

    if (isotopes_converted)
      cerr << isotopes_converted << " isotopes converted to natural form\n";

    elements_to_remove.report (cerr);

    if (element_transformations.number_elements ())
      element_transformations.debug_print (cerr);

    if (molecules_with_abnormal_valences)
      cerr << molecules_with_abnormal_valences << " molecules containing abnormal valences\n";

    if (chemical_standardisation.active ())
      chemical_standardisation.report (cerr);

    if (molecules_with_bad_ring_bond_ratios)
      cerr << molecules_with_bad_ring_bond_ratios << " molecules with ring bond ratios out of range\n";
  }

  return rc;
}

int
main (int argc, char **argv)
{
  prog_name = argv[0];

  int rc = tp_first_pass (argc, argv);

  return rc;
}


/* arch-tag: 4d4389d5-b769-4e7b-ab7b-465769722c88

*/
