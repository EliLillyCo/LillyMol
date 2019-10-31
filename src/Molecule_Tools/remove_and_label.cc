/*
  Removes singly connected atoms and places isotopic labels on
  the adjacent atoms.
  The reason for this is that I want to be able to remove
  things like X and Y from molecules and identify the
  attachment points. Unfortunately, I currently don't have
  a means of specifying element X in a query. but once I change
  the substructure stuff to use element symbol hash values, then
  this becomes un-necessary
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"

#include "istream_and_type.h"
#include "mdl_molecule.h"
#include "ematch.h"
#include "molecule.h"
#include "output.h"
#include "aromatic.h"
#include "molecule_to_query.h"
#include "iwstandard.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

/*
  We can create a command line for the -j option to substitutions
*/

static IWString_and_File_Descriptor stream_for_substitutions_J_option;

/*
  Some programmes need to know which elements are removed at each point
*/

static std::ofstream stream_for_elements_being_replaced;

static int remove_chirality = 0;

static IWString fname_for_query;

static Molecule_Output_Object stream_for_R1R2;

/*
  Significant complexity with only allowing substitutions at the R group
  attachment point.
  If the query produced will be used to search molecules that might, or
  might not, have explicit hydrogens, we need to build the query
  differently.

  If we know that explicit hydrogens will never be present, we can
  specify ncon for all the non R-group atoms in the molecule.
  But if explicit hydrogens might be present, we need to specify
  that information via the hydrogen count
*/

static int only_allow_substitution_at_r_group_sites = 0;

#define ONLY_SUB_RGROUP_NCON 1
#define ONLY_SUB_RGROUP_HCOUNT 2

static int extra_substitutions_allowed_at_substitution_points = 0;

static int extra_bonds_allowed_at_substitution_points = 0;   // 0 means unspecified, obviously must be at least one

static int preserve_ring_membership_of_attachment_points = 0;

static int replace_lost_atom_with_dummy = 0;

static const Element * star_element = NULL;

/*
  Kludge. Only allow certain non-organic elements
*/

static IW_Regular_Expression rx_for_allowable_non_organics;

class Element_and_Isotope
{
  private:
    Element_Matcher _ematch;
    int             _isotope;
    int             _replace_with_dummy_atom;

  public:
    Element_and_Isotope ();

    int build (const const_IWSubstring & buffer);

    int isotope () const { return _isotope;}

    int replace_with_dummy_atom() const { return _replace_with_dummy_atom;}

    int match (const Element * e) { return _ematch.matches (e);}
};

Element_and_Isotope::Element_and_Isotope ()
{
  _isotope = 0;

  _replace_with_dummy_atom = 0;

  return;
}

int
Element_and_Isotope::build (const const_IWSubstring & buffer)
{
  _isotope = 0;

  if (! buffer.contains('='))
  {
    if (! _ematch.construct_from_string(buffer))
    {
      cerr << "Element_and_Isotope::build:invalid element match specification '" << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  const_IWSubstring e, i;
  
  if (! buffer.split(e, '=', i) || 0 == e.length() || 0 == i.length())
  {
    cerr << "Element_and_Isotope::build:must be of the form 'ele=iso'\n";
    return 0;
  }

  if (! _ematch.construct_from_string(e))
  {
    cerr << "Element_and_Isotope::build:invalid element match '" << e << "'\n";
    return 0;
  }

  if (! i.numeric_value(_isotope) || _isotope < 0)
  {
    cerr << "Element_and_Isotope::build:invalid isotopic specification '" << i << "'\n";
    return 0;
  }

  return 1;
}

static Element_and_Isotope * element_and_isotope = NULL;
static int nelei = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Removes atoms and isotopically labels the adjacent atoms\n";
  cerr << "  -R <s=nn>     specify element(s) to be removed and label\n";
  cerr << "  -S <stem>     output file name stem\n";
  cerr << "  -j <fname>    options needed by substitutions\n";
  cerr << "  -J <fname>    rewrite molecule with substitution points changed to R1, R2, etc\n";
  cerr << "  -K <fname>    write connection point and element lost\n";
  cerr << "  -Q <fname>    create query file\n";
  cerr << "  -r            only allow substitution at R group sites\n";
  cerr << "  -h            the only allow substitution directive implemented as hcount rather than ncon\n";
  cerr << "  -p            preserve ring membership of attachment points in query file\n";
  cerr << "  -e <econ>     number of extra connections allowed at substitution points\n";
  cerr << "  -H <rx>       only allow non-organic elements that match <rx>\n";
  cerr << "  -d            replace lost elements with dummy atom (default *)\n";
  cerr << "  -c            remove chirality from input molecules\n";
  cerr << "  -o <type>     output type\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
create_query (MDL_Molecule & m,
              const resizable_array<int> & attachment_points)
{
  if (verbose)
    cerr << "Creating query from '" << m.smiles() << "', will write to '" << fname_for_query << "'\n";

  Molecule_to_Query_Specifications mqs;

//mqs.set_query_must_match_both_explicit_and_implicit_hydrogens(1);

  int matoms = m.natoms();

  if (0 == only_allow_substitution_at_r_group_sites)
    ;
  else if (ONLY_SUB_RGROUP_NCON == only_allow_substitution_at_r_group_sites)
  {
    for (int i = 0; i < matoms; i++)
    {
      MDL_Atom_Data * mad = m.mdl_atom(i);

      mad->set_substitution(m.ncon(i));

      if (NULL == star_element)
        ;
      else if (star_element == m.elementi(i))
        mad->set_substitution(0);
    }

//  Allow any substitution at the attachment points

    for (int i = 0; i < attachment_points.number_elements(); i++)
    {
      int j = attachment_points[i];
  
      MDL_Atom_Data * mad = m.mdl_atom(j);
  
      if (extra_substitutions_allowed_at_substitution_points > 0)
      {
        mad->set_max_ncon(m.ncon(j) + extra_substitutions_allowed_at_substitution_points);
        mad->set_substitution(0);
      }
      else
        mad->set_substitution(0);
  
      if (extra_bonds_allowed_at_substitution_points)
        mad->set_nbonds(m.nbonds(j) + extra_bonds_allowed_at_substitution_points);
  
      if (preserve_ring_membership_of_attachment_points)
        mad->set_ring_bond(m.nrings(j) + 1);
    }
  }
  else     // convey the information about connections via hcount
  {
    for (auto i = 0; i < matoms; ++i)
    {
      MDL_Atom_Data * mad = m.mdl_atom(i);

      mad->set_hcount(m.hcount(i) + 1);    // ISIS convention for how hydrogens are specified
      mad->set_substitution(0);

      if (NULL != star_element && star_element == m.elementi(i))
        mad->set_substitution(0);
    }

    for (auto i = 0; i < attachment_points.number_elements(); ++i)
    {
      auto j = attachment_points[i];

      MDL_Atom_Data * mad = m.mdl_atom(j);

      mad->set_hcount(0);    // which means not specified
    }
  }

// And at any * atoms

  if (NULL != star_element)
  {
    for (int i = 0; i < matoms; i++)
    {
      if (star_element == m.elementi(i))
      {
        MDL_Atom_Data * mad = m.mdl_atom(i);
        mad->set_substitution(0);
      }
    }
  }

  Substructure_Query q;
  if (! q.create_from_molecule(m, mqs))
  {
    cerr << "Cannot create query from '" << m.name() << "'\n";
    return 0;
  }

  return q.write_msi(fname_for_query);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (remove_chirality)
    m.remove_all_chiral_centres();

  return;
}

static int
add_r_group(Molecule & m,
            int & ndx,
            atom_number_t a1,
            atom_number_t a2 = INVALID_ATOM_NUMBER)
{
  IWString s;

  ndx++;

  s << 'R'<< ndx;

  const Element * e = get_element_from_symbol_no_case_conversion(s);

  if (NULL == e)
    e = create_element_with_symbol(s);

  if (NULL == e)
  {
    cerr << "Yipes, cannot fetch element for '" << s << "'\n";
    return 0;
  }

  m.add(e);
  m.add_bond(m.natoms() - 1, a1, SINGLE_BOND);

  if (INVALID_ATOM_NUMBER != a2)
    m.add_bond(m.natoms() - 1, a2, SINGLE_BOND);

  return 1;
}

/*
  We loop over the Element_and_Isotope array as the outer loop.
  That way, if someone has an ordering of elements on the command
  line, that ordering will be reflected in the final output
*/

static int
remove_and_label (MDL_Molecule & m)
{
  cerr << "remove_and_label processing '" << m.name() << "'\n";

  if (! m.arrays_allocated())   // ensure the MDL_File_Data stuff is present
    m.build(m);

  m.remove_explicit_hydrogens();

  int matoms = m.natoms();

  Set_of_Atoms atoms_to_be_removed;
  resizable_array<const Element *> elements_being_removed;
  Set_of_Atoms attachment_points;

  Molecule molecule_for_R1R2(m);

// We need to keep track of whether each item found is singly or doubly connected

  resizable_array<int> single_or_double;

  int ndx = 0;

  for (int i = 0; i < nelei; i++)
  {
    Element_and_Isotope & eli = element_and_isotope[i];

    for (int j = 0; j < matoms; j++)
    {
      const Element * e = m.elementi(j);

      if (! eli.match(e))
        continue;

      int jcon = m.ncon(j);

      if (jcon > 2)
      {
        cerr << "Ignoring non-singly connected '" << e->symbol() << "', atom " << j << " in '" << m.name() << "'\n";
        continue;
      }

      if (! rx_for_allowable_non_organics.active())  // no checking to be done
        ;
      else if (rx_for_allowable_non_organics.matches(e->symbol()))   // great, OK match
        ;
      else
      {
        cerr << "Non-allowed atomic symbol in match atom '" << e->symbol() << "', must match '" << rx_for_allowable_non_organics.source() << "'\n";
        return 0;
      }

      atoms_to_be_removed.add(j);
      elements_being_removed.add(e);

      if (1 == jcon)
      {
        const Atom * a = m.atomi(j);

        const Bond * b = a->item(0);

        if (0 == extra_bonds_allowed_at_substitution_points)   // no need to check anything
          ;
        else if (b->is_single_bond())
          ;
        else if (1 == extra_bonds_allowed_at_substitution_points && b->is_double_bond())
        {
          cerr << "Removing doubly bonded atom, but only 1 extra bond allowed, impossible\n";
          return 0;
        }
        else if (b->is_triple_bond() && extra_bonds_allowed_at_substitution_points < 3)
        {
          cerr << "Removing triply bonded atom, but only " << extra_bonds_allowed_at_substitution_points << " extra bonds allowed, impossible\n";
          return 0;
        }

        atom_number_t k = b->other(j);

        m.set_isotope(k, eli.isotope());

        attachment_points.add(k);

        if (verbose > 2)
          cerr << "Atom " << j << " (" << m.atomic_symbol(j) << ") will be removed, attached to atom " << k << " type " << m.atomic_symbol(k) << endl;

        if (stream_for_R1R2.active())
          add_r_group(molecule_for_R1R2, ndx, k);

        single_or_double.add(1);
      }
      else if (2 == jcon)
      {
        atom_number_t k1 = m.other(j, 0);
        m.set_isotope(k1, eli.isotope());
        attachment_points.add(k1);
        atom_number_t k2 = m.other(j, 1);
        m.set_isotope(k2, eli.isotope());
        attachment_points.add(k2);

        if (verbose > 2)
          cerr << "Atom " << j << " (" << m.atomic_symbol(j) << ") will be removed, attached to atoms " << k1 << " type " << m.atomic_symbol(k1) << " and " << k2 << " type " << m.atomic_symbol(k2) << endl;

        if (stream_for_R1R2.active())
          add_r_group(molecule_for_R1R2, ndx, k1, k2);

        single_or_double.add(2);
      }
    }
  }

  int nr = atoms_to_be_removed.number_elements();

  if (0 == nr)
  {
    cerr << m.name() << " no atoms to be removed\n";
    return 0;
  }

  int * atoms_to_be_removed_array = new_int(matoms); std::unique_ptr<int[]> free_atoms_to_be_removed_array(atoms_to_be_removed_array);

  atoms_to_be_removed.set_vector(atoms_to_be_removed_array, 1);

  if (! stream_for_substitutions_J_option.is_open())
  {
    m.remove_atoms(atoms_to_be_removed_array);
    return 1;
  }

// We need to write out the atom numbers of the attachment points. Adjust atom numbers

  int na = attachment_points.number_elements();

  cerr << "attachment points " << attachment_points << " Removing atoms " << atoms_to_be_removed << endl;

  if (replace_lost_atom_with_dummy)   // no need to adjust the atom numbers
  {
    for (int i = 0; i < na; i++)
    {
      atom_number_t ai = attachment_points[i];
      m.set_element(ai, star_element);
    }
  }
  else
  {
    for (int i = 0; i < na; i++)
    {
      atom_number_t ai = attachment_points[i];
  
      int number_atoms_being_removed_less_than_ai = 0;
  
      for (int j = 0; j < nr; j++)
      {
        if (atoms_to_be_removed[j] < ai)
          number_atoms_being_removed_less_than_ai++;
      }
  
      if (number_atoms_being_removed_less_than_ai)
        attachment_points[i] -= number_atoms_being_removed_less_than_ai;
    }
  }

  if (stream_for_R1R2.active())
  {
    molecule_for_R1R2.remove_atoms(atoms_to_be_removed);
    stream_for_R1R2.write(molecule_for_R1R2);
  }

  assert (single_or_double.number_elements() == nr);

  if (! replace_lost_atom_with_dummy)
    m.remove_atoms(atoms_to_be_removed_array);

  stream_for_substitutions_J_option << " -j ";

  ndx = 0;    // index into the attachment_points array

  for (int i = 0; i < nr; i++)
  {
    if (i > 0)
      stream_for_substitutions_J_option << ',';

    if (1 == single_or_double[i])
    {
      stream_for_substitutions_J_option << attachment_points[ndx];
      ndx++;
    }
    else
    {
      stream_for_substitutions_J_option << attachment_points[ndx] << '-' << attachment_points[ndx + 1];
      ndx += 2;
    }
  }

  if (stream_for_elements_being_replaced.rdbuf()->is_open())
  {
    for (int i = 0; i < nr; i++)
    {
      atom_number_t j = attachment_points[i];
      const Element * e = elements_being_removed[i];
      stream_for_elements_being_replaced << j << ' ' << e->symbol() << '\n';
    }
  }

// We need to warn if the resulting molecule is symmetric at the attachment points

  matoms = m.natoms();

  for (int i = 0; i < na; i++)
  {
    atom_number_t j = attachment_points[i];

    int s = m.symmetry_class(j);

    int number_with_same_symmetry_class = 0;

    for (int k = 0; k < na; k++)
    {
      atom_number_t l = attachment_points[k];

      if (s == m.symmetry_class(l))
        number_with_same_symmetry_class++;
    }

    if (1 == number_with_same_symmetry_class)
      continue;

    cerr << "Symmetry warning, atom " << j << " in '" << m.name() << "' " << number_with_same_symmetry_class << " symmetry equivalent atoms\n";
  }

  if (fname_for_query.length())
  {
    m.transform_to_non_isotopic_form();
    create_query(m, attachment_points);
  }

  return nr;
}

static int
remove_and_label (data_source_and_type<MDL_Molecule> & input,
                  Molecule_Output_Object & output)
{
  MDL_Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<MDL_Molecule> free_m(m);

    preprocess(*m);

    remove_and_label(*m);

    output.write(*m);
  }

  return 1;
}

static int
remove_and_label (const char * fname, int input_type, 
                  Molecule_Output_Object & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<MDL_Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return remove_and_label(input, output);
}

static int
remove_and_label (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lR:rS:o:j:J:Q:K:H:e:b:cdph");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  set_auto_create_new_elements(1);

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
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
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  nelei = cl.option_count('R');

  if (0 == nelei)
  {
    cerr << "Must specify one or more elements to remove via the -R option\n";
    usage(3);
  }

  if (nelei > 0)
  {
    element_and_isotope = new Element_and_Isotope[nelei];

    int i = 0;
    const_IWSubstring r;
    while (cl.value('R', r, i))
    {
      if (! element_and_isotope[i].build(r))
      {
        cerr << "Invalid element replacement specification '" << r << "'\n";
        usage(3);
      }

      i++;
    }
  }

  if (cl.option_present('H'))
  {
    const_IWSubstring h = cl.string_value('H');

    if (! rx_for_allowable_non_organics.set_pattern(h))
    {
      cerr << "Cannot initialise regular expression for allowable non-organics '" << h << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Only non organic elements matching '" << rx_for_allowable_non_organics.source() << "' will be allowed\n";
  }

  if (cl.option_present('r'))
  { 
    only_allow_substitution_at_r_group_sites = ONLY_SUB_RGROUP_NCON;
    if (verbose)
      cerr << "Will only allow substitution at R group site(s)\n";
  }

  if (cl.option_present('h'))
  {
    only_allow_substitution_at_r_group_sites = ONLY_SUB_RGROUP_HCOUNT;

    if (verbose)
      cerr << "Query will only allow substitution at R group sites (directed by Hcount)\n";
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', extra_substitutions_allowed_at_substitution_points) || extra_substitutions_allowed_at_substitution_points < 1)
    {
      cerr << "The extra connections at substitution points option (-e) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Must have " << extra_substitutions_allowed_at_substitution_points << " extra substitutions at attachment points\n";

    if (! only_allow_substitution_at_r_group_sites)
    {
      cerr << "Only allowing substitution at R group sites (-r)\n";
      only_allow_substitution_at_r_group_sites = 1;
    }
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', extra_bonds_allowed_at_substitution_points) || extra_bonds_allowed_at_substitution_points < 1)
    {
      cerr << "The extra bonds allowed at substitution points (-b) option must be a whole number >= 1\n";
      usage(4);
    }

    if (verbose)
      cerr << "A maximum of " <<extra_bonds_allowed_at_substitution_points << " extra bonds allowed at substitution points\n";
  }

  if (cl.option_present('p'))
  {
    preserve_ring_membership_of_attachment_points = 1;

    if (verbose)
      cerr << "Ring membership of attachment points will be preserved\n";
  }

  if (cl.option_present('c'))
  {
    remove_chirality = 1;
    if (verbose)
      cerr << "Chirality will be removed from all input molecules\n";
  }

  if (cl.option_present('d'))
  {
    replace_lost_atom_with_dummy = 1;

    star_element = get_element_from_symbol_no_case_conversion("*");
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  Molecule_Output_Object output;

  if (cl.option_present('o'))
  {
    if (! output.determine_output_types(cl, 'o'))
    {
      cerr << "Cannot determine output type(s), -o option\n";
      return 3;
    }
  }
  else
    output.add_output_type(SDF);

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s), '" << s << "'\n";
      return 4;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot open output stream(s), stem '" << s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Output written to '" << s << "'\n";
  }
  else
  {
    output.new_stem("-");
  }

  if (cl.option_present('j'))
  {
    const char * j = cl.option_value('j');

    if (! stream_for_substitutions_J_option.open(j))
    {
      cerr << "Cannot open stream for substitutions '" << j << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Information for substitutions written to '" << j << "'\n";
  }

  if (cl.option_present('K'))
  {
    const char * k = cl.option_value('K');

    stream_for_elements_being_replaced.open(k, std::ios::out);

    if (! stream_for_elements_being_replaced.good())
    {
      cerr << "Cannot open stream for element replacements '" << k << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Information for elements being replaced written to '" << k << "'\n";
  }

  if (cl.option_present('J'))
  {
    const const_IWSubstring j = cl.string_value('J');

    stream_for_R1R2.add_output_type(SMI);

    if (stream_for_R1R2.would_overwrite_input_files(cl, j))
    {
      cerr << "Cannot overwrite input file(s) with -J option '" << j << "'\n";
      return 4;
    }

    if (! stream_for_R1R2.new_stem(j))
    {
      cerr << "Cannot initialise -J file '" << j << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Substitution pattern R1, R2, etc, writen to '" << j << ".smi'\n";
  }

  if (cl.option_present('Q'))
  {
    cl.value('Q', fname_for_query);

    if (! fname_for_query.ends_with(".qry"))
      fname_for_query << ".qry";

    if (verbose)
      cerr << "Will write query to '" << fname_for_query << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! remove_and_label(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (cl.option_present('j'))
    stream_for_substitutions_J_option << '\n';

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = remove_and_label(argc, argv);

  return rc;
}
