/*
  Tester for symmetry determinations.
*/

#include <iomanip>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/minmaxspc.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwminmax.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/symmetry.h"

using std::cerr;

static Element_Transformations etrans;

const char * prog_name = nullptr;

static int molecules_read = 0;
static int molecules_passing = 0;
static int molecules_failing = 0;
static int do_pass_fail = 0;

/*
  We can filter on 
   1. the number of symmetric atoms
   2. the fraction of symmetric atoms
   3. the number of non-symmetric atoms
   4. the fraction of non-symmetric atoms
*/

static Min_Max_Specifier<int> symmat;
static Min_Max_Specifier<float> fsymmat;
static Min_Max_Specifier<int> nsymmat;
static Min_Max_Specifier<float> fnsymmat;

static int verbose = 0;

static int debug_print_symmetry_info = 0;

/*
  We can function as a generic tester for symmetry perception, or we
  can do the symmetry demerit for third party lists.
*/

static int function_as_tester = 1;

static IWString tag_for_symmetry_scaling;

static Molecule_Output_Object stream_for_labelled_atoms;

static Molecule_Output_Object stream_for_molecules_that_pass;
static Molecule_Output_Object stream_for_molecules_that_fail;

static int append_rejection_reason_to_name = 0;

/*
  For third party lists, we produce a scaling factor. The scale
  factor already in the file will be multiplied by this value.
  The lowest scale factor we can emit is 0.1

  If the molecule is entirely three-fold or higher symmetry, the
  scale factor is 0.1
  If the molecule is entirely two-fold symmetry, the scale factor
  is 0.2
*/

static int too_many_classes_count = 0;
static int molecules_with_errors = 0;

/*
  We can ignore uninteresting symmetry - cf3, benzene and such. We can also
  ignore symmetry groupings that are below a given size
*/

static int remove_uninteresting_symmetry = 0;

/*
  Sometimes we aren't interested in localised symmetry, but longer range
*/

static int ignore_if_closer_than = 0;

/*
  Generally, widely separated symmetry is more interesting than short
  range symmetry. We can examine symmetry groupings that stretch beyond
  a given bond separation
*/

static int max_separation_threshold = 4;

/*
  One kind of demerit is to demerit anything with more than
  a certain number of symmetry related atoms
*/

static int too_many_symmetric_atoms = 0;

/*
  Another possibility is to compute symmetry related descriptors
*/

static std::ofstream stream_for_symmetry_descriptors;
static int write_symmetry_descriptors_to_stdout = 0;

static int include_centre_atoms_with_symmetry_groupings = 0;

static int ignore_molecules_with_fewer_atoms = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Usage " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << "Computes, tests or filters symmetry related values\n";
  cerr << "  -T <tag>       do third party symmetry demerit processing\n";
  cerr << "  -x <number>    demerit anything with more than <number> symmetric atoms\n";
  cerr << "  -S <file>      compute symmetry descriptors and write to <file>\n";
  cerr << "  -L <fname>     file for labelled atoms - by symmetry group\n";
  cerr << "  -u <number>    remove uninteresting symmetry. If number > 1, symmetry\n";
  cerr << "                 groupings with <= number atoms will be discarded\n";
  cerr << "  -d <dist>      ignore symmetry when atoms <dist> bonds apart or shorter\n";
  cerr << "  -w <dist>      symmetric atoms > <dist> bonds apart considered well separated\n";
  cerr << "  -W ...         symmetry rejection condition, enter '-W help'\n";
  cerr << "  -m <stem>      write molecules that pass symmetry conditions to <stem>\n";
  cerr << "  -n <stem>      write molecules that fail symmetry conditions to <stem>\n";
  cerr << "  -a             append rejection reason to name in -n file\n";
  cerr << "  -c             include centre atoms in symmetry groupings\n";
  cerr << "  -t ...         element transformation options, enter '-t help'\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
// clang-format on

  exit(rc);
}

static void
do_remove_uninteresting_symmetry (Molecule & m,
                                  Symmetry_Info & symmetry_info)
{
  symmetry_info.remove_trivial_cf3(m);

  symmetry_info.remove_benzene(m);

  if (1 == remove_uninteresting_symmetry)
    return;

  for (int i = symmetry_info.symmetry_groupings() - 1; i >= 0; i--)
  {
    const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(i);
    if(s->number_elements() <= remove_uninteresting_symmetry)
      symmetry_info.remove_symmetry_grouping(i);
  }

  return;
}

/*
  We have the neighbours of an atom.  Identify those neighbours that
  are singly connected and all the same
*/

static int
identify_singly_attached_equivalent_neighbours (Molecule & m,
                                                Set_of_Atoms & s)
{
  int symmetry_class = -1;

  for (int i = s.number_elements() - 1; i >= 0; i--)
  {
    atom_number_t j = s[i];

    const Atom * aj = m.atomi(j);

    if (1 != aj->ncon())     // not singly connected, omit
      s.remove_item(i);
    else if (symmetry_class < 0)    // first singly connected atom we've seen
      symmetry_class = m.symmetry_class(j);
    else if (symmetry_class != m.symmetry_class(j))    // different symmetry classes, done
      return 0;
  }

  return s.number_elements() > 1;
}

/*
  We need to identify two atoms with the same symmetry grouping
*/

static int
short_range_symmetry_only (Molecule & m,
                           int s,
                           int ignore_if_closer_than)
{
  int matoms = m.natoms();

  atom_number_t a1 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; i++)
  {
    if (s != m.symmetry_class(i))
      continue;

    if (INVALID_ATOM_NUMBER == a1)
    {
      a1 = i;
      continue;
    }

    if (m.fragment_membership(i) != m.fragment_membership(a1))
      continue;

    int d = m.bonds_between(a1, i);

    return d <= ignore_if_closer_than;
  }

  return 0;
}

static int
identify_uninteresting_symmetry (Molecule & m,
                                 int * is_uninteresting)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (4 != a->ncon())
      continue;

    Set_of_Atoms s;
    a->connections(i, s);

    if (! identify_singly_attached_equivalent_neighbours(m, s))
      continue;

    s.set_vector(is_uninteresting, 1);
    rc++;
  }

  return rc;
}

template <typename T>
void
do_append(Molecule & m,
          const char * zreason,
          T v)
{
  IWString tmp(m.name());

  tmp.append_with_spacer(zreason);

  if (! tmp.ends_with(' '))
    tmp << ' ';

  tmp << v;

  m.set_name(tmp);

  return;
}

template void do_append(Molecule &, const char *, int);
template void do_append(Molecule &, const char *, float);

/*
  Rather than using symmetry groupings, we focus instead on actually
  symmetric atoms
*/


static int
_check_pass_fail(Molecule & m)
{
  const int matoms = m.natoms();

  int * is_uninteresting = new_int(matoms); std::unique_ptr<int[]> free_is_uninteresting(is_uninteresting);

  identify_uninteresting_symmetry(m, is_uninteresting);

  int * scount = new_int(matoms + 1); std::unique_ptr<int[]> free_scount(scount);

  const int * s = m.symmetry_classes();

  for (int i = 0; i < matoms; i++)
  {
    if(is_uninteresting[i])
      continue;

    const int si = s[i];

    assert (si >= 0 && si <= matoms);

    scount[si]++;
  }

  int symmetry_groupings = 0;
  int symmetric_atoms = 0;

  for (int i = 0; i <= matoms; i++)
  {
    if (scount[i] <= 1)
      continue;

    if (remove_uninteresting_symmetry > 0 && scount[i] <= remove_uninteresting_symmetry)
      continue;

    if (ignore_if_closer_than > 0 && short_range_symmetry_only(m, i, ignore_if_closer_than))
      continue;

    symmetry_groupings++;
    symmetric_atoms += scount[i];
  }

  if (symmetric_atoms > matoms)
  {
    cerr << m.name() << " too many symmetric atoms " << symmetric_atoms << " but " << matoms << " atoms in molecule\n";
  }

  if (verbose > 2)
    cerr << m.name() << " natoms = " << matoms << " found " << symmetric_atoms << " symmetric atoms in " << symmetry_groupings << " groups\n";

  if (symmat.is_set())
  {
    if (! symmat.matches(symmetric_atoms))
    {
      if (append_rejection_reason_to_name)
        do_append(m, "symmat", symmetric_atoms);
      return 0;
    }
  }

  if(fsymmat.is_set())
  {
    const float tmp = iwmisc::Fraction<float>(symmetric_atoms, matoms);

    if (! fsymmat.matches(tmp))
    {
      if (append_rejection_reason_to_name)
        do_append(m, "fsymat", tmp);
      return 0;
    }
  }

  int non_symmetric_atoms = matoms - symmetric_atoms;

  if (! nsymmat.matches(non_symmetric_atoms))
  {
    if (append_rejection_reason_to_name)
      do_append(m, "nsymmat", non_symmetric_atoms);
    return 0;
  }

  if(fnsymmat.is_set())
  {
    const float tmp = iwmisc::Fraction<float>(non_symmetric_atoms, matoms);
    if (! fnsymmat.matches(tmp))
    {
      if (append_rejection_reason_to_name)
        do_append(m, "fnsymmat", tmp);
      return 0;
    }
  }

  return 1;
}

static int
check_pass_fail(Molecule & m)
{
  if (m.natoms() < ignore_molecules_with_fewer_atoms)
  {
//  cerr << "too few atoms " << m.natoms() << '\n';
    if (stream_for_molecules_that_pass.active())
      stream_for_molecules_that_pass.write(m);

    return 1;
  }

  const int pf = _check_pass_fail(m);

  if(pf)
  {
    molecules_passing++;
    if(stream_for_molecules_that_pass.active())
      stream_for_molecules_that_pass.write(m);
  }
  else
  {
    molecules_failing++;
    if(stream_for_molecules_that_fail.active())
      stream_for_molecules_that_fail.write(m);
  }

  return 1;
}

static int
connected_to_two_equivalent_atoms (const Molecule & m,
                                   atom_number_t zatom,
                                   const int * is_in_symmetry_group)
{
  const Atom * a = m.atomi(zatom);

  resizable_array<int> symmetries_found;
  int symmetric_neighbours_encountered = 0;

  for (int i = a->ncon() - 1; i >= 0; i--)
  {
    const atom_number_t j = a->other(zatom, i);

    if (0 == is_in_symmetry_group[j])
      continue;

    symmetric_neighbours_encountered++;
    symmetries_found.add_if_not_already_present(is_in_symmetry_group[j]);
  }

  if (symmetric_neighbours_encountered == symmetries_found.number_elements())
    return 0;

  return 1;
}

static int
do_include_centre_atoms_with_symmetry_groupings(Molecule & m,
                                const Symmetry_Info & symmetry_info,
                                int * is_in_symmetry_group)
{
  const int matoms = m.natoms();

  set_vector(is_in_symmetry_group, matoms, 0);

  const int groups = symmetry_info.symmetry_groupings();

//cerr << "Has " << groups << " groups\n";

  for (int i = 0; i < groups; i++)
  {
    const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(i);

//  cerr << "Group " << i << " has " << s->number_elements() << " atoms\n";

    if (s->number_elements() > 1)
      s->set_vector(is_in_symmetry_group, i + 1);
  }

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
//  cerr << " atom " << i << " value " << is_in_symmetry_group[i] << '\n';

    if (is_in_symmetry_group[i])
      continue;

    if (connected_to_two_equivalent_atoms (m, i, is_in_symmetry_group))
      rc++;
  }

  return rc;
}

static void
compute_atoms_in_wide_ranging_symmetry_groups(Molecule & m,
                                              int & max_sym_sep,
                                              int & atoms_in_wide_ranging_symmetry_groups)
{
  max_sym_sep = 0;
  atoms_in_wide_ranging_symmetry_groups = 0;

  const int * s = m.symmetry_classes();

  const int matoms = m.natoms();

  resizable_array<int> symmetry_classes_with_long_separations;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == s[i])
      continue;

    for (int j = i + 1; j < matoms; j++)
    {
      if (s[j] != s[i])
        continue;

      if (m.fragment_membership(i) != m.fragment_membership(j))
        continue;

      const int d = m.bonds_between(i, j);

      if (d > max_sym_sep)
        max_sym_sep = d;

      if (d > max_separation_threshold)
        symmetry_classes_with_long_separations.add_if_not_already_present(s[i]);
    }
  }

  if (0 == symmetry_classes_with_long_separations.number_elements())
    return;

  for (int i = 0; i < matoms; i++)
  {
    if (symmetry_classes_with_long_separations.contains(s[i]))
      atoms_in_wide_ranging_symmetry_groups++;
  }

  return;
}

static int
compute_symmetry_descriptors(Molecule & m,
                             std::ostream & os)
{
  write_space_suppressed_string(m.name(), os);

  Symmetry_Info symmetry_info;

  symmetry_info.initialise(m);

  if (debug_print_symmetry_info)
    symmetry_info.debug_print(cerr);

  const int matoms = m.natoms();

  os << ' ' << matoms;

  int groups = symmetry_info.symmetry_groupings();
  if (0 == groups)      // no symmetry in the molecule
  {
    os << " 0 0 0 0 0 0 0 0 0 0\n";

    return os.good();
  }

  double dmatoms = static_cast<double>(matoms);

  int largest_symmetric_grouping = 0;
  int highest_degree = 0;

  double si = 0.0;

  if (verbose > 2)
    cerr << m.name() << " contains " << groups << " symmetry groups\n";

  int * is_in_symmetry_group = new_int(matoms); std::unique_ptr<int[]> free_is_in_symmetry_group(is_in_symmetry_group);

  for (int i = 0; i < groups; i++)
  {
    const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(i);

    s->set_vector(is_in_symmetry_group, 1);

    const int ns = s->number_elements();

    if (ns > largest_symmetric_grouping)
      largest_symmetric_grouping = ns;

    if(s->degree() > highest_degree)
      highest_degree = s->degree();

    const double tmp = static_cast<double>(ns) / dmatoms;

    si -= tmp * log10(tmp);
  }

  int atoms_in_symmetry_groupings = count_non_zero_occurrences_in_array(is_in_symmetry_group, matoms);

  if (atoms_in_symmetry_groupings > matoms)
    cerr << m.name() << " too many symmetric atoms " << atoms_in_symmetry_groupings << " but " << matoms << " atoms in molecule\n";

  if (include_centre_atoms_with_symmetry_groupings)
    atoms_in_symmetry_groupings += do_include_centre_atoms_with_symmetry_groupings(m, symmetry_info, is_in_symmetry_group);

  int max_sym_sep = 0;
  int atoms_in_wide_ranging_symmetry_groups = 0;
  compute_atoms_in_wide_ranging_symmetry_groups(m, max_sym_sep, atoms_in_wide_ranging_symmetry_groups);

  os << ' ' << groups;
  os << ' ' << m.number_symmetry_classes();
  os << ' ' << atoms_in_symmetry_groupings;
  os << ' ' << static_cast<double>(atoms_in_symmetry_groupings) / dmatoms;
  os << ' ' << largest_symmetric_grouping;
  os << ' ' << highest_degree;
  os << ' ' << si;
  os << ' ' << max_sym_sep;
  os << ' ' << atoms_in_wide_ranging_symmetry_groups;
  os << ' ' << iwmisc::Fraction<float>(atoms_in_wide_ranging_symmetry_groups, matoms);
  os << '\n';

  return os.good();
}

static int
do_third_party_symmetry_demerit (Molecule & m,
                                 Symmetry_Info & symmetry_info,
                                 float & symmetry_scale_factor)
{
  symmetry_scale_factor = 1.0;

  int matoms = m.natoms();

  int groups = symmetry_info.symmetry_groupings();
  // If no groups, no scaling applied.
  if (0 == groups)
    return 1;

  if (1 == groups)
  {
    const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(0);

    int gsize = s->number_elements();
    if (gsize == matoms)
    {
      if(s->degree() >= 3)
        symmetry_scale_factor = 0.1;
      else
        symmetry_scale_factor = 0.2;

      return 1;
    }

    if (gsize < 10)
      symmetry_scale_factor = 1.0 - 0.3 * float(gsize) / float(matoms);
    else if (gsize < 15)
      symmetry_scale_factor = 1.0 - 0.5 * float(gsize) / float(matoms);
    else
      symmetry_scale_factor = 1.0 - 0.8 * float(gsize) / float(matoms);

    if(s->degree() > 2)
      symmetry_scale_factor *= 0.5;

    return 1;
  }

// We have multiple groups

  int atoms_in_symmetry_groupings = 0;

  for (int i = 0; i < groups; i++)
  {
    const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(i);

    int ns = s->number_elements();
    if (1 == ns)       // not symmetric
      continue;

    atoms_in_symmetry_groupings += s->number_elements();

    float tmp = 0.6f * iwmisc::Fraction<float>(s->number_elements(), matoms);

    if(s->degree() > 2)
      tmp *= 0.5f;

    symmetry_scale_factor -= tmp;
  }

  return 1;
}

static void
preprocess(Molecule & m)
{
  m.reduce_to_largest_fragment();

  if(etrans.number_elements())
   (void) etrans.process(m);

  return;
}

static int
do_third_party_symmetry_demerit (Molecule & m,
                                 std::ostream & os)
{
  IWString smiles = m.smiles();

  preprocess(m);

  Symmetry_Info symmetry_info;

  symmetry_info.initialise(m);

  if(verbose)
    symmetry_info.debug_print(cerr);

  if(remove_uninteresting_symmetry)
    do_remove_uninteresting_symmetry(m, symmetry_info);

  float symmetry_scale_factor = 1.0;

 (void) do_third_party_symmetry_demerit(m, symmetry_info, symmetry_scale_factor);

  os << smiles << ' ' << m.name() << ' ' << tag_for_symmetry_scaling << ':' << symmetry_scale_factor << '\n';

  return os.good();
}

static int
test_symmetry_class (Molecule & m,
                     const Set_of_Atoms * s)
{
  int ns = s->number_elements();

  if (1 == ns)     // no checking to do
    return 1;

  atom_number_t a0 = s->item(0);
  atomic_number_t z = m.atomic_number(a0);
  int con = m.ncon(a0);
  formal_charge_t fc = m.formal_charge(a0);
  int is_ring = m.is_ring_atom(a0);
  int mfrag = m.fragment_membership(a0);

  for (int i = 1; i < ns; i++)
  {
    atom_number_t a = s->item(i);
    if (m.atomic_number(a) != z)
      return 0;
    if (m.ncon(a) != con)
      return 0;
    if (m.formal_charge(a) != fc)
      return 0;
    if (m.is_ring_atom(a) != is_ring)
      return 0;
    if (m.fragment_membership(a) != mfrag)
     return 0;
  }

  return 1;
}

static int
tsymmetry (Molecule & m)
{
  if(verbose)
    cerr << "Processing " << m.name() << '\n';

  preprocess(m);

  int matoms = m.natoms();

  Symmetry_Info symmetry_info;

  symmetry_info.initialise(m);

  if(remove_uninteresting_symmetry)
    do_remove_uninteresting_symmetry(m, symmetry_info);

  if(verbose)
    symmetry_info.debug_print(cerr);

  int classes = symmetry_info.number_elements();

  if (verbose || classes > matoms)
  {
    cerr << " Molecule with " << m.natoms() << " atoms has " << classes << " symmetry class";
    if (classes > 1)
      cerr << "es";
    cerr << '\n';
  }

  int errors_encountered_this_molecule = 0;

  if (classes > matoms)
  {
    cerr << "Impossible\n";
    too_many_classes_count++;
    errors_encountered_this_molecule++;
  }

// Scan each symmetry class

  for (int i = 0; i < classes; i++)
  {
    const Set_of_Atoms * s = symmetry_info[i];

    assert(s->number_elements() <= matoms);

    if (0 == test_symmetry_class(m, s))
    {
      cerr << "Yipes, error for symmetry class " <<(*s) << '\n';
      errors_encountered_this_molecule++;
    }
  }


  if(errors_encountered_this_molecule)
  {
    molecules_with_errors++;
    m.debug_print(cerr);
    return 1;
  }

  if(stream_for_labelled_atoms.active() && symmetry_info.symmetry_groupings())
  {
    IWString mname = m.name();

    float symmetry_scale_factor = 1.0;

   (void) do_third_party_symmetry_demerit(m, symmetry_info, symmetry_scale_factor);

    IWString tmp = mname;

    tmp << " SYMMSCALE: " << symmetry_scale_factor;
    m.set_name(tmp);

    stream_for_labelled_atoms.write(m);

    m.set_name(mname);

    for (int i = 0; i < symmetry_info.symmetry_groupings(); i++)
    {
      if (i > 0)
        m.transform_to_non_isotopic_form();

      const Symmetric_Atoms * s = symmetry_info.symmetry_grouping(i);
      for (int j = 0; j < s->number_elements(); j++)
      {
        atom_number_t k = s->item(j);
        m.set_isotope(k, i + 1);
      }

      stream_for_labelled_atoms.write(m);
    }
  }

  return 0;
}

static int
tsymmetry(data_source_and_type<Molecule> & input)
{
  Molecule * m;
  int rc = 0;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1)
      cerr << molecules_read << ' ' << m->name() << '\n';

    preprocess(*m);

    if (0 == m->natoms())    // ignore
      cerr << "Ignoring molecule with 0 atoms '" << m->name() << "'\n";
    else if(function_as_tester)
      rc += tsymmetry(*m);
    else if(do_pass_fail)
      rc += check_pass_fail(*m);
    else if(stream_for_symmetry_descriptors.rdbuf()->is_open())
     (void) compute_symmetry_descriptors(*m, stream_for_symmetry_descriptors);
    else if(write_symmetry_descriptors_to_stdout)
     (void) compute_symmetry_descriptors(*m, std::cout);
    else if(tag_for_symmetry_scaling.length())
     (void) do_third_party_symmetry_demerit(*m, std::cout);
    else
    {
      cerr << "What am I supposed to be doing?\n";
      return 0;
    }
  }

  return rc;
}

static int
tsymmetry (const char * fname, FileType input_type)
{
  assert (nullptr != fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return tsymmetry(input);
}

static int
handle_file_opening (Command_Line & cl,
                     char flag,
                     Molecule_Output_Object & output)
{
  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SMI);
  else if (! output.determine_output_types(cl, 'o'))
  {
    cerr << "Cannot determine output type(s)\n";
    usage(8);
  }

  const_IWSubstring stem = cl.string_value(flag);

  if (output.would_overwrite_input_files(cl, stem))
  {
    cerr << "Cannot overwrite input files\n";
    return 0;
  }

  if (! output.new_stem(stem))
  {
    cerr << "Cannot open output '" << stem << "'\n";
    return 0;
  }

  return 1;
}

static int
tsymmetry (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "S:x:T:vi:E:A:t:L:o:u:W:m:n:ad:czw:");
  if(cl.unrecognised_options_encountered())
    usage(2);

  verbose = cl.option_count('v');

  if (! process_elements(cl))
    usage(3);

  if (! process_standard_smiles_options(cl, verbose))
  {
    usage(4);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    usage(5);
  }

  if (cl.option_present('t'))
  {
    if (! process_element_transformations(cl, etrans, verbose))
      usage(6);
  }

  if (cl.option_present('T'))
  {
    function_as_tester = 0;

    tag_for_symmetry_scaling = cl.string_value('T');

    if(verbose)
      cerr << "Will compute symmetry demerit for third party acquisitions, tag '" << tag_for_symmetry_scaling << "'\n";
  }

  if (cl.option_present('u'))
  {
    if (! cl.value('u', remove_uninteresting_symmetry) || remove_uninteresting_symmetry <= 0)
    {
      cerr << "The -u option must be followed by a whole positive number\n";
      usage(17);
    }

    if(verbose)
    {
      cerr << "Uninteresting symmetry items discarded";
      if (remove_uninteresting_symmetry > 1)
        cerr << ", will ignore groups with less than " << remove_uninteresting_symmetry << " atoms";
      cerr << '\n';
    }
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', ignore_if_closer_than) || ignore_if_closer_than < 1)
    {
      cerr << "The ignore closer than (-d option) specifiecation must be a value whole number\n";
      usage(4);
    }

    if(verbose)
      cerr << "Will ignore symmetry groupings where the atoms are " << ignore_if_closer_than << " or bonds apart\n";
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', max_separation_threshold) || max_separation_threshold < 1)
    {
      cerr << "The max separation threshold (-w option) specifiecation must be a value whole number\n";
      usage(4);
    }

    if(verbose)
      cerr << "Symmetric atoms > " << max_separation_threshold << " bonds will be considered well separated\n";
  }

  if (cl.option_present('x'))
  {
    if(function_as_tester)
    {
      cerr << "The -x option only makes sense with the -T option\n";
      usage(63);
    }

    if (! cl.value('x', too_many_symmetric_atoms) || too_many_symmetric_atoms < 0)
    {
      cerr << "The -x switch requires a non-negative whole number\n";
      usage(12);
    }

    if(verbose)
      cerr << "Molecules with more than " << too_many_symmetric_atoms << " will be demerited\n";
  }

  if (cl.option_present('c'))
  {
    include_centre_atoms_with_symmetry_groupings = 1;

    if (verbose)
      cerr << "Will include centre atoms with symmetry groups\n";
  }

  if (cl.option_present('z'))
  {
    debug_print_symmetry_info = 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(9);
  }

  if (cl.option_present('S'))
  {
    if (cl.option_present('n') || cl.option_present('m') || cl.option_present('W'))
    {
      cerr << "Sorry, the -S option is not compatible with the -n and -m options\n";
      usage(4);
    }

    function_as_tester = 0;

    IWString fname = cl.string_value('S');

    IWString header("id natoms ngroups nclass symat fsymat lrgsymgr hgsymdgr si maxsep wellsep fwellsep\n");

    if ('-' == fname)
    {
      write_symmetry_descriptors_to_stdout = 1;
      std::cout << header;
    }
    else
    { 
      stream_for_symmetry_descriptors.open(fname.null_terminated_chars());
      if (! stream_for_symmetry_descriptors.good())
      {
        cerr << "Sorry, cannot open '" << fname << "' for output\n";
        return 24;
      }

      if(verbose)
        cerr << "Symmetry related descriptors will be written to '" << fname << "'\n";

      stream_for_symmetry_descriptors << header;
    }

  }

  if (cl.option_present('L'))
  {
    if (! handle_file_opening(cl, 'L', stream_for_labelled_atoms))
    {
      cerr << "Cannot initialise labelled atom stream(-L)\n";
      return 3;
    }

    if(verbose)
      cerr << "Labelled atoms written to '" << stream_for_labelled_atoms.stem() << "'\n";
  }

  if (cl.option_present('W'))
  {
    function_as_tester = 0;

    if (cl.option_present('m'))
    {
      if (! handle_file_opening(cl, 'm', stream_for_molecules_that_pass))
      {
        cerr << "Cannot open stream for passing molecules(-m)\n";
        return 5;
      }

      if(verbose)
        cerr << "Molecules that pass the symmetry conditions written to '" << stream_for_molecules_that_pass.stem() << "'\n";
    }

    if (cl.option_present('n'))
    {
      if (! handle_file_opening(cl, 'n', stream_for_molecules_that_fail))
      {
        cerr << "Cannot open stream for failing molecules(-n)\n";
        return 5;
      }

      if(verbose)
        cerr << "Molecules that fail the symmetry conditions written to '" << stream_for_molecules_that_fail.stem() << "'\n";

      if (cl.option_present('a'))
      {
        append_rejection_reason_to_name = 1;
        if (verbose)
          cerr << "Rejection reason appended to molecule name\n";
      }
    }

    int i = 0;
    const_IWSubstring w;
    while (cl.value('W', w, i++))
    {
      if (w.starts_with("symmat"))
      {
        w.remove_leading_chars(6);
        if (! symmat.initialise(w))
        {
          cerr << "Cannot initialise condition for 'symmat', '" << w << "'\n";
          return 3;
        }
      }
      else if (w.starts_with("fsymmat"))
      {
        w.remove_leading_chars(7);
        if (! fsymmat.initialise(w))
        {
          cerr << "Cannot initialise condition for 'fsymmat', '" << w << "'\n";
          return 3;
        }
      }
      else if (w.starts_with("nsymmat"))
      {
        w.remove_leading_chars(7);
        if (! nsymmat.initialise(w))
        {
          cerr << "Cannot initialise condition for 'nsymmat', '" << w << "'\n";
          return 3;
        }
      }
      else if (w.starts_with("fnsymmat"))
      {
        w.remove_leading_chars(8);
        if (! fnsymmat.initialise(w))
        {
          cerr << "Cannot initialise condition for 'fnsymmat', '" << w << "'\n";
          return 3;
        }
      }
      else if (w.starts_with("minatm="))
      {
        w.remove_leading_chars(7);
        if (! w.numeric_value(ignore_molecules_with_fewer_atoms) || ignore_molecules_with_fewer_atoms < 1)
        {
          cerr << "The ignore molecules below atom count 'minatm=' directive must have a whole +ve number\n";
          return 1;
        }

        if (verbose)
          cerr << "Will not consider molecules with fewer than " << ignore_molecules_with_fewer_atoms << " atoms\n";
      }
      else if ("help" == w)
      {
        // clang-format off
        cerr << " -W symmat<op><number>      number   of atoms in symmetry relationships\n";
        cerr << " -W fsymmat<op><ratio>      fraction of atoms in symmetry relationships\n";
        cerr << " -W nsymmat<op><number>     number   of atoms NOT in symmetry relationships\n";
        cerr << " -W nfsymmat<op><ratio>     fraction of atoms NOT in symmetry relationships\n";
        cerr << " where <op> is one of '=', '>=' or '<='\n";
        cerr << " -W minatm=<natoms>         just PASS molecules with fewer than <natoms> atoms - avoid problems with benzene\n";
        cerr << "                            note that these molecules are written to the -m stream\n";
        // clang-format on

        return 0;
      }
      else
      {
        cerr << "Unrecognised -W qualifier '" << w << "'\n";
        usage(6);
      }
    }

    do_pass_fail = 1;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    rc += tsymmetry(cl[i], input_type);
  }

  if(verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    if(do_pass_fail)
      cerr << molecules_passing << " molecules passed conditions, " << molecules_failing << " molecules failed\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tsymmetry(argc, argv);

  return rc;
}
