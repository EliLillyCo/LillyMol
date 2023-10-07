/*
  Computes ring substitution fingerprint
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/path_around_ring.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int molecules_read = 0;

static int molecules_with_no_rings = 0;

static IWString fingerprint_tag ("NCRS<");

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int positional_information_only = 1;

static int simple_atom_types = 0;

static int full_atom_types = 0;

static int verbose = 0;

static int max_path_length = 24;

static int flatten_fingerprints_to_01 = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int work_as_filter = 0;

static Accumulator_Int<int> acc_nset;

/*
  We can lower the distances between molecules by including presence
  or absence of the features
*/

static int place_single_feature_bits = 0;

/*
  Make sure the ring atom types have gaps in them because we
  incrememt by 1 if the ring atom is aromatic
*/

#define RSTYPE_DOUBLE_BOND_OUTSIDE_RING 1
#define RSTYPE_SPIRO 2
#define RSTYPE_RING_JOIN 3
#define RSTYPE_TWO_SUBSTITUENTS 4
#define RSTYPE_SUBSTITUTED 5
#define RSTYPE_ANOTHER_RING 7

#define RSTYPE_SATURATED_CARBON 9
#define RSTYPE_UNSATURATED_CARBON 11
#define RSTYPE_SATURATED_HETEROATOM 13
#define RSTYPE_UNSATURATED_HETEROATOM 15

#define RSTYPE_METHYL 17
#define RSTYPE_TERMINAL_N 19
#define RSTYPE_NITRO 21
#define RSTYPE_SATURATED_NITROGEN 23
#define RSTYPE_SP2_NITROGEN 25
#define RSTYPE_HYDROXY 27
#define RSTYPE_ETHER 29
#define RSTYPE_SULPH 31
#define RSTYPE_FLUORINE 33
#define RSTYPE_HALOGEN 35

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "  -J <tag>       tag for fingerprints\n";
  cerr << "  -P <bonds>    maximum path length around edge of ring system\n";
//cerr << "  -b            break after finding adjacent match\n";
  cerr << "  -M ...        options for specifying just what kind of information to generate\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -A ...        aromaticity specifications, enter '-A help' for info\n";
  cerr << "  -E ...        element specifications, enter '-E help' for info\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
ring_substitution(const IWString & mname,
                  const resizable_array<int> & abstract_path,
                  int * tmp,
                  Sparse_Fingerprint_Creator & sfpc)
{
  int n = abstract_path.number_elements();

//#define DEBUG_RING_SUBSTITUTION
#ifdef DEBUG_RING_SUBSTITUTION
  cerr << mname << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " abstract_path " << abstract_path[i] << endl;
  }
#endif

  assert (n > 2);

  copy_vector(tmp, abstract_path.rawdata(), n);
  copy_vector(tmp + n, abstract_path.rawdata(), n);

  int jstop = n;
  if (jstop > max_path_length)
    jstop = max_path_length;

  for (int i = 0; i < n; i++)
  {
    int ta = tmp[i];

    if (0 == ta)
      continue;

    for (int j = 1; j < jstop; j++)
    {
      int k = j + i;    // index of next atom

      int tb = tmp[k];
      if (0 == tb)
        continue;

      unsigned int b = 40 + 30 + j * 30 * 30 * max_path_length + (ta * tb);   // 30 is larger than all atom types assigned

#ifdef DEBUG_RING_SUBSTITUTION
      cerr << " i = " << i << " ar " << ta << " j = " << j << " k = " << k << " ar " << tb << " b = " << b << endl;
#endif
      sfpc.hit_bit(b);
    }
  }

  return 1;
}

static int
ring_substitution(Molecule & m,
                  const int * atype,
                  const Set_of_Atoms & par,
                  int * tmp,
                  Sparse_Fingerprint_Creator & sfpc)
{
  resizable_array<int> abstract_path;
 
  const int n = par.number_elements();

  int number_features = 0;
  int the_feature = -1;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = par[i];

    abstract_path.add(atype[j]);

    if (0 == atype[j])
      continue;

    if (0 == number_features)
      the_feature = atype[j];

    number_features++;

    if (place_single_feature_bits)
      sfpc.hit_bit(atype[j]);
  }

//cerr << "Abstract ring contains " << number_features << " features\n";

  if (1 == number_features)
  {
    sfpc.hit_bit(1383 * n + the_feature);
    return 1;
  }

  assert (abstract_path.number_elements() == n);

  return ring_substitution(m.name(), abstract_path, tmp, sfpc);
}

static int
ring_substitution(Molecule & m,
                  const int * atype,
                  Sparse_Fingerprint_Creator & sfpc)
{
  const int nr = m.nrings ();

  if (0 == nr)
  {
    molecules_with_no_rings++;
    return 1;
  }

  const int matoms = m.natoms();

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int * tmp = new int[matoms + matoms]; std::unique_ptr<int[]> free_tmp(tmp);

  m.compute_aromaticity_if_needed();

// Do all the non-fused rings first. Includes spiro fused

  int rings_processed = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

//  cerr << "Ring " << i << " has " << ri->largest_number_of_bonds_shared_with_another_ring() << " shared bonds\n";

    if (ri->largest_number_of_bonds_shared_with_another_ring() > 0)
      continue;

    ring_substitution(m, atype, *ri, tmp, sfpc);

    ring_already_done[i] = 1;
    rings_processed++;
  }

  if (nr == rings_processed)
    return 1;

// Now any fused rings

  int * in_ring_system = new_int(matoms); std::unique_ptr<int[]> free_in_ring_system(in_ring_system);

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    set_vector(in_ring_system, matoms, 0);

    const Ring * ri = m.ringi(i);

    ri->set_vector(in_ring_system, 1);

    rings_processed++;

    for (int j = i + 1; j < nr; j++)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);

      if (rj->fused_system_identifier() != ri->fused_system_identifier())
        continue;

      rj->set_vector(in_ring_system, 1);
      ring_already_done[j] = 1;

      rings_processed++;
    }

    Set_of_Atoms s;
    if (! path_around_edge_of_ring_system(m, in_ring_system, 1, s))   // strongly fused
    {
      sfpc.hit_bit(count_occurrences_of_item_in_array(1, matoms, in_ring_system));   // unprocessed ring, hit bit according to size
      continue;
    }

    ring_substitution(m, atype, s, tmp, sfpc);

    if (nr == rings_processed)
      break;
  }

  return 1;
}

static int
determine_substitution_type(Molecule & m,
                            atom_number_t zatom,
                            const Atom & a)
{
  int acon = a.ncon();
//cerr << "determine_substitution_type for atom type " << m.smarts_equivalent_for_atom(zatom) << '\n';

  if (2 == acon)    // unsubstituted
    return 0;

  if (m.nrings(zatom) > 1)
    return RSTYPE_RING_JOIN;

  if (4 == acon)
    return RSTYPE_TWO_SUBSTITUENTS;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a[i];

    if (b->nrings())
      continue;

    if (b->is_double_bond())
      return RSTYPE_DOUBLE_BOND_OUTSIDE_RING;

    if (positional_information_only)
      return RSTYPE_SUBSTITUTED;

    atom_number_t o = b->other(zatom);

    if (m.is_ring_atom(o))
      return RSTYPE_ANOTHER_RING;

    atomic_number_t zo = m.atomic_number(o);

    int ocon = m.ncon(o);

    if (simple_atom_types)
    {
      int unsaturation = m.nbonds(o) - ocon;
      if (6 == zo)
      {
        if (0 == unsaturation)
          return RSTYPE_SATURATED_CARBON;
        else
          return RSTYPE_UNSATURATED_CARBON;
      }
      else if (1 == ocon && (9 == zo || 17 == zo || 35 == zo || 53 == zo))
        return RSTYPE_HALOGEN;
      else if (0 == unsaturation)
        return RSTYPE_SATURATED_HETEROATOM;
      else
        return RSTYPE_UNSATURATED_HETEROATOM;

      continue;
    }

//  Full atom typing

    if (6 == zo)
    {
      if (1 == ocon)
        return RSTYPE_METHYL;

      int nbonds = m.nbonds(o);
      if (ocon == nbonds)
        return RSTYPE_SATURATED_CARBON;
      else
        return RSTYPE_UNSATURATED_CARBON;
    }
    else if (7 == zo)
    {
      if (1 == ocon)
        return RSTYPE_TERMINAL_N;

      int nbonds = m.nbonds(o);

      if (3 == ocon && 5 == nbonds)
        return RSTYPE_NITRO;
      else if (ocon == nbonds)
        return RSTYPE_SATURATED_NITROGEN;
      else
        return RSTYPE_SP2_NITROGEN;
    }
    else if (8 == zo)
    {
      if (1 == ocon)
        return RSTYPE_HYDROXY;
      else
        return RSTYPE_ETHER;
    }
    else if (16 == zo)
    {
      if (1 == ocon)
        return RSTYPE_HYDROXY;
      else if (2 == ocon)
        return RSTYPE_ETHER;
      else
        return RSTYPE_SULPH;   // some other kind of state, who knows
    }
    else if (9 == zo)
      return RSTYPE_FLUORINE;
    else if (17 == zo || 35 == zo || 53 == zo)
      return RSTYPE_HALOGEN;

    return RSTYPE_SUBSTITUTED;    // of not classified above
  }

  assert (nullptr == "Could not classify non-ring bond???");

  return 0;
}

static int
is_spiro_fused(Molecule & m,
               atom_number_t zatom,
               const Atom & a)
{
  int acon = a.ncon();

  if (4 != acon)
    return 0;

  if (2 != m.nrings(zatom))   // too hard and rare otherwise
    return 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a[i];

    if (0 == b->nrings())
      return 0;
  }

  return 1;
}

static int
has_double_bond_outside_ring(Molecule & m,
                             atom_number_t zatom,
                             const Atom & a)
{
  int acon = a.ncon();

  if (3 != acon)
    return 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a[i];
//  cerr << " from atom " << zatom << " have bond to " << b->other(zatom) << " ring " << b->nrings() << " single? " << b->is_single_bond() << '\n';

    if (b->nrings())
      continue;

    if (b->is_double_bond())
      return 1;
  }

  return 0;
}

static int
assign_atom_types(Molecule & m,
                  int * atype)
{
  m.ring_membership();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    int nr = m.nrings(i);

    if (0 == nr)
      continue;

    const Atom * a = m.atomi(i);

    if (has_double_bond_outside_ring(m, i, *a))
    {
      atype[i] = RSTYPE_DOUBLE_BOND_OUTSIDE_RING;
      if (m.is_aromatic(i))
        atype[i]++;
      continue;
    }

    if (is_spiro_fused(m, i, *a))
    {
      atype[i] = RSTYPE_SPIRO;
      continue;
    }

    atype[i] = determine_substitution_type(m, i, *a);

    if (0 == atype[i])
      ;
    else if (m.is_aromatic(i))
      atype[i]++;
  }

  return 1;
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment ();

  if (chemical_standardisation.active ())
    chemical_standardisation.process (m);

  return;
}

static int
ring_substitution(Molecule & m,
                  IWString_and_File_Descriptor & output)
{
  preprocess(m);

  int * atype = new_int(m.natoms()); std::unique_ptr<int[]> free_atype(atype);

  assign_atom_types(m, atype);

#ifdef DEBUG_RING_SUBSTITUTION
  for (int i = 0; i < m.natoms(); i++)
  {
    cerr << "Atom " << i << " assigned type " << atype[i] << endl;
  }
#endif

  Sparse_Fingerprint_Creator sfpc;

  if (! ring_substitution(m, atype, sfpc))
    return 0;

  if (! work_as_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  if (flatten_fingerprints_to_01) {
    sfpc.flatten_to_01();
  }

  IWString tmp;
  sfpc.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);

  output << tmp << '\n';

  if (verbose > 2) {
    sfpc.debug_print(cerr);
  }

  if (! work_as_filter) {
    output << "|\n";
  }

  if (verbose) {
    acc_nset.extra(sfpc.nset());
  }

  if (verbose > 2) {
    cerr << molecules_read << " '" << m.name() << "' hits " << sfpc.nbits() << " bits\n";
  }

  return 1;
}

static int
ring_substitution(data_source_and_type<Molecule> & input,
                  IWString_and_File_Descriptor & output)
{
  Molecule * m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (! ring_substitution(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
ring_substitution(const char * fname,
                  FileType input_type,
                  IWString_and_File_Descriptor & output)
{
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 2)
    input.set_verbose(1);

  return ring_substitution(input, output);
}

static int
ring_substitution_filter_record(const const_IWSubstring & buffer,
                                IWString_and_File_Descriptor & output)
{
  const_IWSubstring tmp(buffer);
  tmp.remove_leading_chars(smiles_tag.length());
  tmp.chop(1);

  Molecule m;
  if (! m.build_from_smiles(tmp)) {
    cerr << "Invalid smiles '" << tmp << "'\n";
    return 0;
  }

  return ring_substitution(m, output);
}

static int
ring_substitution_filter(iwstring_data_source & input,
                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << "\n";

    if (! buffer.starts_with(smiles_tag))
      continue;

    if (! ring_substitution_filter_record(buffer, output)) {
      cerr << "Fatal error, line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static void
display_dash_M_options(std::ostream & os)
{
  os << "  -M posi    set bits to indicate positions on rings only\n";
  os << "  -M simp    set bits to indicate positions and simple atom types\n";
  os << "  -M full    set bits to indicate positions and complete atom types\n";
  os << "  -M sfb     set bits to indicate presence or absence of single features\n";
  os << "  -M 01      remove count information from fingerprints, just presence of absence\n";

  exit(2);
}

static int
ring_substitution_filter(const char * fname,
                         IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return 0;
  }

  return ring_substitution_filter(input, output);
}

static int
ring_substitution(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:J:lg:M:fbP:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose > 1))
  {
    cerr << "Cannot process -A option\n";
    usage(11);
  }

  if (! process_elements(cl, verbose > 1, 'E'))
  {
    cerr << "Cannot initialise elements\n";
    usage(8);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('b'))
  {
//  break_at_subsitutents = 1;
    if (verbose)
      cerr << "Will break after encountering closest substituent\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (cl.option_present('f'))
    ;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('P'))
  {
    if (! cl.value('P', max_path_length) || max_path_length < 2)
    {
      cerr << "The maximum path length (-m) must be a whole +ve number\n";
      usage(7);
    }

    if (verbose)
      cerr << "Max path length " << max_path_length << endl;
  }

  if (cl.option_present('M'))
  {
    int nset = 0;

    positional_information_only = 0;

    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (m.starts_with("posi"))
      {
        positional_information_only = 1;
        nset++;
        if (verbose)
          cerr << "Only positional information will be recorded\n";
      }
      else if (m.starts_with("simp"))
      {
        simple_atom_types = 1;
        nset++;
        if (verbose)
          cerr << "Simple atom types for substitutent atoms\n";
      }
      else if ("full" == m)
      {
        full_atom_types = 1;
        nset++;
        if (verbose)
          cerr << "Full atom typing\n";
      }
      else if ("sfb" == m)
      {
        place_single_feature_bits = 1;
        if (verbose)
          cerr << "Will set bits according to presence or absence of features\n";
      }
      else if ("01" == m)
      {
        flatten_fingerprints_to_01 = 1;
        if (verbose)
          cerr << "Fingerprints will be flattened to just 01 forms\n";
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_M_options(cerr);
        usage(4);
      }
    }

    if (1 != nset)
    {
      cerr << "Sorry, must specify one, and exactly one of 'posi', 'simple' or 'full' for substituent atom types\n";
      usage(6);
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f'))
  {
    work_as_filter = 1;
    ring_substitution_filter(cl[0], output);
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! ring_substitution(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Processed " << molecules_read << " molecules, " << molecules_with_no_rings << " had no rings\n";
    if (acc_nset.n() > 1)
      cerr << "Fingerprints have between " << acc_nset.minval() << " and " << acc_nset.maxval() << " bits set, ave " << static_cast<float>(acc_nset.average()) << endl;
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ring_substitution(argc, argv);

  return rc;
}
