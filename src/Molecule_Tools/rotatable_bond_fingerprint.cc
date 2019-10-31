/*
  Atom pair fingerprint where the number of bonds is the number
  of rotatable bonds
*/

#include <iostream>
#include <memory>
#include <limits>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"
#include "sparse_fp_creator.h"

#include "istream_and_type.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "rotbond_common.h"
#include "atom_typing.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static Atom_Typing_Specification atom_typing_specification;

static int min_distance = 0;
static int max_distance = std::numeric_limits<int>::max();

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

static IWString tag("NCRB<");

static int function_as_filter = 0;

static int include_atom_pair_bits = 0;

static int only_form_bits_for_all_rotatable_bonds_between = 0;

static extending_resizable_array<int> global_chain_length;

static IWString_and_File_Descriptor stream_for_rbcount;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Generates atom pair fingerprints, distance is rotatable bond count\n";
  cerr << " -c <number>    shortest bond distance to process\n";
  cerr << " -C <number>    longest  bond distance to process\n";
  cerr << " -P <type>      atom type specification, enter '-P help' for info\n";
  cerr << " -J ...         fingerprint specification\n";
  cerr << " -p             include atom pair bits\n";
  cerr << " -k             only produce a bit if all bonds between two atoms are rotatable\n";
  cerr << " -D <fname>     write detailed bond separation info to <fname>\n";
  cerr << " -l             reduce to largest fragment\n";
  cerr << " -f             function as a tdt filter\n";
  cerr << " -i <type>      input specification\n";
  cerr << " -g ...         chemical standardisation options\n";
  cerr << " -E ...         standard element specifications\n";
  cerr << " -A ...         standard aromaticity specifications\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
no_rotatable_bonds(Molecule & m,
                   IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;
  sfc.hit_bit(2776704, 1);

  sfc.hit_bit(49682238, m.nedges());

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << "\n";

  return 1;
}

static int
write_rotbond_data(Molecule & m,
                   const atom_number_t a1,
                   const atom_number_t a2,
                   const int d,
                   const int rb,
                   IWString_and_File_Descriptor & output)
{
  Molecule mcopy(m);
  mcopy.set_atom_map_number(a1, rb);
  mcopy.set_atom_map_number(a2, rb);
  output << mcopy.smiles() << ' ' << m.name() << ' ' << d << '\n';
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static void
form_bit_not_bonded(unsigned int atype1,
                    const  int distance,
                    unsigned int atype2, 
                    Sparse_Fingerprint_Creator & sfc)
{
  if (atype1 > atype2)
    std::swap(atype1, atype2);

  sfc.hit_bit(975535 * atype1 + distance * 7927 + atype2);

  return;
}

static unsigned int
bond_constant(const Bond * b)
{
  if (b->is_aromatic())
    return 4;

  if (b->is_single_bond())
    return 1;

  if (b->is_double_bond())
    return 2;

  if (b->is_triple_bond())
    return 3;

  return 5;
}

static void
form_bit_rotatable(unsigned int atype1,
                   const int distance,
                   const int rb,
                   unsigned int atype2,
                   Sparse_Fingerprint_Creator & sfc)
{
  if (atype1 > atype2)
    std::swap(atype1, atype2);

  sfc.hit_bit(333161 * atype1 + rb * 19912 + atype2);

  if (rb < distance / 2)
  {
    sfc.hit_bit(1, 1);
    sfc.hit_bit(atype1 + atype2, 1);
  }
  else
  {
    sfc.hit_bit(2, 1);
    sfc.hit_bit(2 * atype1 + atype2, 1);
  }

  sfc.hit_bit(3, 10 * rb / distance + 1);
  sfc.hit_bit(3 * atype1 + atype2, rb/distance + 1);

  return;
}

static void
form_bit_bonded(unsigned int atype1,
                const Bond * b,
                unsigned int atype2,
                Sparse_Fingerprint_Creator & sfc)
{
  if (atype1 > atype2)
    std::swap(atype1, atype2);

  auto bc = bond_constant(b);

  if (b->nrings())
    bc = bc * 11;

  sfc.hit_bit(41202 * atype1 + bc * 2121 + atype2);

  return;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
rotatable_bond_fingerprint(Molecule & m,
                           const int * bonds_between,
                           const int * rotatable_bonds_between,
                           IWString_and_File_Descriptor & output)
{
  const auto matoms = m.natoms();

  unsigned int * atype = new unsigned int[matoms]; std::unique_ptr<unsigned int[]> free_atype(atype);

  if (! atom_typing_specification.assign_atom_types(m, atype)) 
  {
    cerr << "Cannot assign atom types '" << m.name() << "'\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;

  extending_resizable_array<int> chain_length;

  for (auto i = 0; i < matoms; ++i)
  {
    for (auto j = i + 1; j < matoms; ++j)
    {
      const auto d = bonds_between[i * matoms + j];

      if (d > max_distance)
        continue;

      if (! include_atom_pair_bits)
        ;
      else if (1 == d)
        form_bit_bonded(atype[i], m.bond_between_atoms(i, j), atype[j], sfc);
      else
        form_bit_not_bonded(atype[i], d, atype[j], sfc);

      const auto rb = rotatable_bonds_between[i * matoms + j];

      if (0 == rb)
        continue;

      if (stream_for_rbcount.is_open())
        write_rotbond_data(m, i, j, d, rb, stream_for_rbcount);

      if (rb == d)
        chain_length[rb]++;

      if (only_form_bits_for_all_rotatable_bonds_between && rb < d)
        continue;

      form_bit_rotatable(atype[i], d, rb, atype[j], sfc);
    }
  }

  int longest_chain = 0;
  int tot = 0;
  int nchains = 0;
  for (int i = 1; i < chain_length.number_elements(); ++i)
  {
    if (chain_length[i] > 0)
    {
      sfc.hit_bit(i, chain_length[i]);

      global_chain_length[i]++;

      tot += i * chain_length[i];
      nchains += 1;

      longest_chain = i;
    }
  }

  sfc.hit_bit(898231222, longest_chain);

  sfc.hit_bit(7231, tot / nchains + 1);

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << "\n";

  return 1;
}

static int
rotatable_bond_fingerprint(Molecule & m,
                           IWString_and_File_Descriptor & output)
{
  const auto matoms = m.natoms();

  (void) m.compute_aromaticity_if_needed();   // the atom typing probably needs this

  int * rotatable = new_int(matoms * matoms); std::unique_ptr<int[]> free_rotatable(rotatable);
  
  const auto nedges = m.nedges();

  int rotatable_bonds = 0;

  const Atom ** atom = new const Atom *[matoms]; std::unique_ptr<const Atom *[]> free_atom(atom);

  m.atoms(atom);

  for (auto i = 0; i < nedges; ++i)
  {
    const auto * b = m.bondi(i);

    if (! b->is_single_bond())
      continue;

    if (b->nrings())    // ring bonds are never flexible here
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (1 == m.ncon(a1))
      continue;

    if (1 == m.ncon(a2))
      continue;

    if (triple_bond_at_either_end(m, b))
      continue;

    if (is_non_rotatable_amide(m, a1, a2))    // O=C-[NH]- is non-rotatable
      continue;

    if (is_non_rotatable_sulphonamide(m, a1, a2))    // O=S(=O)-[NH]- is non-rotatable
      continue;

    rotatable[a1 * matoms + a2] = 1;
    rotatable[a2 * matoms + a1] = 1;

    rotatable_bonds++;
  }

  if (0 == rotatable_bonds)
    return no_rotatable_bonds(m, output);

  int maxcon = 4;
  for (auto i = 0; i < matoms; ++i)
  {
    const auto c = atom[i]->ncon();
    if (c > maxcon)
      maxcon = c;
  }

  maxcon++;

  const auto distance_not_set = matoms + matoms + matoms;

  int * bonds_between = new_int(matoms*matoms, distance_not_set); std::unique_ptr<int[]> free_bonds_between(bonds_between);
  int * rotatable_bonds_between = new_int(matoms*matoms, -1); std::unique_ptr<int[]> free_rotatable_bonds_between(rotatable_bonds_between);

  int * connection_table = new int[matoms * maxcon]; std::unique_ptr<int[]> free_ct(connection_table);

  for (auto i = 0; i < matoms; ++i)
  {
    connection_table[i * maxcon] = atom[i]->connections(i, connection_table + i * maxcon + 1);
  }


/*int * all_needed = new int[matoms* matoms + matoms*matoms + matoms]; std::unique_ptr<int[]> free_all(all_needed);
  int * bonds_between = all_needed;
  int * rotatable_bonds_between = all_needed + matoms*matoms;
  int * connections = all_needed + matoms*matoms + matoms* matoms;

  set_vector(bonds_between, matoms*matoms, distance_not_set);
  set_vector(rotatable_bonds_between, matoms*matoms, -1);*/

  for (auto i = 0; i < matoms; ++i)
  {
    bonds_between[i * matoms + i] = 0;

    const auto icon = connection_table[i * maxcon];
    const int * cti = connection_table + i * maxcon + 1;

    for (auto j = 0; j < icon; ++j)
    {
      const auto k = cti[j];

      const auto r = rotatable[i * matoms + k];

      bonds_between[i * matoms + k] = 1;
      bonds_between[k * matoms + i] = 1;

      if (r)
      {
        rotatable_bonds_between[i * matoms + k] = 1;
        rotatable_bonds_between[k * matoms + i] = 1;
      }
      else
      {
        rotatable_bonds_between[i * matoms + k] = 0;
        rotatable_bonds_between[k * matoms + i] = 0;
      }
    }
  }

  for (auto d = 2; d < matoms; ++d)
  {
    bool keep_going = false;

    for (auto i = 0; i < matoms; ++i)
    {
      int * distance_to_i  = bonds_between + i * matoms;

      for (auto j = 0; j < matoms; ++j)
      {
        if (j == i)
          continue;

        if (distance_not_set == distance_to_i[j])
          continue;

        const auto jcon = connection_table[j * maxcon];
        const int * ctj = connection_table + j * maxcon + 1;

        for (auto k = 0; k < jcon; ++k)
        {
          const auto l = ctj[k];

          if (distance_to_i[j] + 1 < distance_to_i[l])
          {
            distance_to_i[l] = distance_to_i[j] + 1;
            bonds_between[l * matoms + i] = distance_to_i[l];
            rotatable_bonds_between[i * matoms + l] = rotatable_bonds_between[i * matoms + j] + rotatable[j * matoms + l];
            rotatable_bonds_between[l * matoms + i] = rotatable_bonds_between[i * matoms + l];
            keep_going = true;
          }
        }
      }
    }

    if (! keep_going)
      break;
  }

#ifdef PRINT_DISTANCE_MATRIX
  cout << "      ";
  for (auto i = 0; i < matoms; ++i)
  {
    cout << ' ' << i;
  }
  cout << endl;
  for (auto i = 0; i < matoms; ++i)
  {
    cout << "Atom " << i;
    for (auto j = 0; j < matoms; ++j)
    {
      if (j == i)
        continue;

      cout << ' ' << bonds_between[i * matoms + j];
    }
    cout << endl;
  }
#endif

//#define CHECK_VS_MOLECULE
#ifdef CHECK_VS_MOLECULE
  m.recompute_distance_matrix();
  for (auto i = 0; i < matoms; ++i)
  {
    for (auto j = 0; j < matoms; ++j)
    {
      if (i == j)
        continue;

      if (bonds_between[i * matoms + j] != m.bonds_between(i, j))
      {
        cerr << "Distance matrix inconsistency, computed " << bonds_between[i * matoms + j] << " molecule yields " << m.bonds_between(i, j) << ", atoms " << i << " and " << j << endl;
        return 0;
      }
    }
  }
#endif

  if (stream_for_rbcount.is_open())
  {
    IWString tmp(m.name());
    tmp << ' ' << rotatable_bonds;
    m.set_name(tmp);
  }

  return rotatable_bond_fingerprint(m, bonds_between, rotatable_bonds_between, output);
}

static int
rotatable_bond_fingerprint(data_source_and_type<Molecule> & input,
                           IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    output << smiles_tag << m->smiles() << ">\n";
    output << identifier_tag << m->name() << ">\n";

    if (! rotatable_bond_fingerprint(*m, output))
      return 0;

    output << "|\n";

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
rotatable_bond_fingerprint(const char * fname, int input_type, 
                           IWString_and_File_Descriptor & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return rotatable_bond_fingerprint(input, output);
}

static int
rotatable_bond_fingerprint_filter_molecule (const_IWSubstring buffer,     // we have our own copy
                                            IWString_and_File_Descriptor & output)
{
  buffer.chop();

  buffer.remove_leading_chars(smiles_tag.length());

  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Cannot parse smiles '" << buffer << "'\n";
    return 0;
  }

  preprocess(m);

  return rotatable_bond_fingerprint(m, output);
}

static int
rotatable_bond_fingerprint_filter(iwstring_data_source & input,
                                  IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if (! buffer.starts_with(smiles_tag))
      continue;

    if (! rotatable_bond_fingerprint_filter_molecule(buffer, output))
      return 0;
  }

  return 1;
}

static int
rotatable_bond_fingerprint_filter(const char * fname,
                                  IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open for filter processing\n";
    return 0;
  }

  return rotatable_bond_fingerprint_filter(input, output);
}

static int
rotatable_bond_fingerprint(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lJ:P:c:C:fpkD:");

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
  else
    set_global_aromaticity_type(Daylight);

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

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('k'))
  {
    only_form_bits_for_all_rotatable_bonds_between = 1;
    if (verbose)
      cerr << "Will only fingerprint atom pairs with all rotatable bonds between\n";
  }

  int input_type = 0;

  if (cl.option_present('f'))
  {
    function_as_filter = 1;
  }
  else if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('J'))
  {
    cl.value('J', tag);

    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "Fingerprints written with tag '" << tag << "'\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', min_distance)  || min_distance < 1)
    {
      cerr << "The minimum distance option (-c) must be a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only fingerprint distances " << min_distance << " or shorter\n";
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', max_distance)  || max_distance < min_distance)
    {
      cerr << "The maximum distance option (-C) must be a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only fingerprint distances " << max_distance << " or longer\n";
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 2;
    }
  }
  else
  {
    atom_typing_specification.set_atom_type(IWATTYPE_COMPLEX);
  }

  if (cl.option_present('p'))
  {
    include_atom_pair_bits = 1;

    if (verbose)
      cerr << "Will also include atom pair bits\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('D'))
  {
    IWString fname = cl.string_value('D');
    if (! fname.ends_with(".smi"))
      fname << ".smi";

    if (! stream_for_rbcount.open(fname.null_terminated_chars()))
    {
      cerr << "Cannot open stream for rotatable bond info '" << fname << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Detailed rotatable bond info written to '" << fname << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (function_as_filter)
  {
    if (! rotatable_bond_fingerprint_filter(cl[0], output))
      rc = 2;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! rotatable_bond_fingerprint(cl[i], input_type, output))
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

    for (int i = 1; i < global_chain_length.number_elements(); ++i)
    {
      if (global_chain_length[i] > 0)
        cerr << global_chain_length[i] << " molecules had a chain of length " << i << endl;
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rotatable_bond_fingerprint(argc, argv);

  return rc;
}
