/*
  For testing 3D related things. Make arbitrary changes to the molecules
*/

#include <iostream>
#include <limits>
#include <memory>
#include <random>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int number_transformations_per_molecule = 1;

static int max_atoms_translated_per_molecule = std::numeric_limits<int>::max();

static int arbitrary_rotations = 0;
static int arbitrary_molecular_translations = 0;
static int arbitrary_atomic_translations = 0;
static int invert_molecule = 0;

static distance_t max_molecular_translation = 10.0;
static distance_t max_atomic_translation = 0.05;

static distance_t bump_check = static_cast<distance_t>(0.0);

static resizable_array_p<Substructure_Query> only_translate_these_atoms;

static int write_unchanged_molecules = 1;

static std::random_device rd;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Applies random geometric transformations to (mostly rigid) molecules\n";
  cerr << "  -e            invert coordinates\n";
  cerr << "  -r            perform arbitrary rotations\n";
  cerr << "  -T <max>      perform arbitrary molecular translations, max offset <max>\n";
  cerr << "  -t <max>      perform arbitrary atomic translations, max offset <max>\n";
  cerr << "  -a <number>   for atomic translations, max number of atoms per molecule moved\n";
  cerr << "  -s <smarts>   for atomic translations, only move atoms that match <smarts>\n";
  cerr << "  -n <number>   number transformations per molecule\n";
  cerr << "  -b <dist>     bump check, reject molecule if any two atoms closer than <dist>\n";
  cerr << "  -c            if writing smiles output, include coordinates with smiles\n";
  cerr << "  -z            write molecules even if they are unchanged\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input  specification\n";
  cerr << "  -o <type>     output specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
ok_bump_check (const Molecule & m,
               const distance_t bump_check)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    const auto ai = m.atomi(i);

    for (auto j = i + 1; j < matoms; ++j)
    {
      if (m.atomi(j)->distance(*ai) < bump_check)
      {
        cerr << "Bad distance " << m.atomi(j)->distance(*ai) << endl;
        return 0;
      }
    }
  }

  return 1;
}

/*
  VERY inefficient
*/

static int
only_match_queries_match(Molecule_to_Match & target,
                         const atom_number_t zatom)
{
  for (int i = 0; i < only_translate_these_atoms.number_elements(); ++i)
  {
    Substructure_Results sresults;
    const int nhits = only_translate_these_atoms[i]->substructure_search(target, sresults);

    for (int j = 0; j < nhits; ++j)
    {
      if (sresults.embedding(j)->contains(zatom))
        return 1;
    }
  }

  return 0;
}

static void
do_translate_atom (Molecule & m,
                   const atom_number_t zatom,
                   std::uniform_real_distribution<float> & u)
{
  const Atom * a = m.atomi(zatom);

  coord_t t = u(rd);
  if (u(rd) < 0.5)
    t = -t;

  m.setx(zatom, a->x() + t);

  t = u(rd);
  if (u(rd) < 0.5)
    t = -t;

  m.sety(zatom, a->y() + t);

  t = u(rd);
  if (u(rd) < 0.5)
    t = -t;

  m.setz(zatom, a->z() + t);

  return;
}

static int
do_arbitrary_atomic_translations_max (Molecule & m,
                                      const int max_atoms_translated_per_molecule)
{
  Molecule_to_Match target(&m);    // may not get used

  const int matoms = m.natoms();

  std::uniform_int_distribution<int> um(0, matoms-1);

  int istart = um(rd);

  std::uniform_real_distribution<float> u(0.0f, max_atomic_translation);

  int processed = 0;
  int attempts = 0;

  for (int i = istart; processed < max_atoms_translated_per_molecule && attempts <= matoms; i = (i+1) % matoms, ++attempts)
  {
    if (only_translate_these_atoms.number_elements() && ! only_match_queries_match(target, i))
      continue;

    do_translate_atom(m, i, u);
    processed++;
  }

  return processed;
}

static int
do_arbitrary_atomic_translations (Molecule & m)
{
  Molecule_to_Match target(&m);    // may not get used

  const int matoms = m.natoms();

  if (max_atoms_translated_per_molecule < matoms)
    return do_arbitrary_atomic_translations_max(m, max_atoms_translated_per_molecule);

  std::uniform_real_distribution<float> u(0.0f, max_atomic_translation);

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (only_translate_these_atoms.number_elements() && ! only_match_queries_match(target, i))
      continue;

    do_translate_atom(m, i, u);
    rc++;
  }

  return rc;
}

static int
do_invert_molecule(Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    coord_t x = a->x();
    coord_t y = a->y();
    coord_t z = a->z();

    m.setxyz(i, -x, -y, -z);
  }

  int nc = m.chiral_centres();

  for (int i = 0; i < nc; i++)
  {
    const Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    m.invert_chirality_on_atom(c->a());
  }

  return 1;
}

static int
do_arbitrary_rotation(Molecule & m)
{
  std::uniform_int_distribution<int> u(0, m.natoms() - 1);

  int i = u(rd);

  Coordinates origin;
  m.get_coords(i, origin);
  m.translate_atoms(-origin);

  std::uniform_real_distribution<float> u01(0.0f, 1.0f);

  coord_t x = u01(rd);
  coord_t y = u01(rd);
  coord_t z = u01(rd);

  Coordinates rotation_axis;
  rotation_axis.setxyz(x, y, z);
  rotation_axis.normalise();

  std::uniform_real_distribution<float> uhalfpi(0.0f, M_PI * 0.50);
//angle_t theta = DEG2RAD * intbtwij(0, 180);
  const angle_t theta = uhalfpi(rd);

  m.rotate_atoms(rotation_axis, theta);

  m.translate_atoms(origin);

  return 1;
}

static int
do_arbitrary_molecular_translation (Molecule & m)
{
  std::uniform_real_distribution<float> u(0.0f, 1.0f);

  coord_t x = u(rd) * max_molecular_translation;
  if (u(rd) < 0.5)
    x = -x;

  coord_t y = u(rd) * max_molecular_translation;
  if (u(rd) < 0.5)
    y = -y;

  coord_t z = u(rd) * max_molecular_translation;
  if (u(rd) < 0.5)
    z = -z;

  if (verbose > 2)
    cerr << "Translation by (" << x << ',' << y << ',' << z << ")\n";

  m.translate_atoms(x, y, z);

  return 1;
}

/*
*/

static int
arbitrary_geometric_changes (Molecule & m,
                             Molecule_Output_Object & output)
{
  if (bump_check > 0.0 && (! ok_bump_check(m, bump_check)))
  {
    cerr << "Warning, initial molecule '" << m.name() << "' fails bump check " << bump_check << endl;
  }

  Coordinates * c = nullptr;

  if (number_transformations_per_molecule > 1)   // need to save initial coordinates
  {
    c = new Coordinates[m.natoms()];
    m.get_coords(c);
  }

  int ntries = 0;

  // only update nproduced when we produce a new version
  for (auto nproduced = 0; nproduced < number_transformations_per_molecule; )
  {
    if (nproduced > 0)    // restore initial coordinates
      m.setxyz(c);

    ntries++;

    if (ntries > number_transformations_per_molecule && 0 == ntries % 100)
    {
      cerr << "Warning: " << m.name() << " has made " << ntries << " attempts to generate a geometry passing the bump check\n";
      if (ntries > 100 * number_transformations_per_molecule)
        return 0;
    }

    int changes = 0;

    if (arbitrary_rotations)
    {
      do_arbitrary_rotation(m);
      changes++;
    }

    if (arbitrary_molecular_translations)
    {
      if (do_arbitrary_molecular_translation(m))
        changes++;
    }


    if (invert_molecule)
    {
      do_invert_molecule(m);
      changes++;
    }

    if (arbitrary_atomic_translations)
    {
      if (do_arbitrary_atomic_translations(m))
        changes++;
    }

    if (bump_check > 0.0 && ! ok_bump_check(m, bump_check))
      continue;

    if (0 == changes && ! write_unchanged_molecules)
      continue;

    nproduced++;

    output.write(m);
  }

  if (nullptr != c)
    delete [] c;

  return 1;
}

static int
arbitrary_geometric_changes (data_source_and_type<Molecule> & input,
                             Molecule_Output_Object & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule ()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! arbitrary_geometric_changes(*m, output))
      return 0;
  }

  return output.good();
}

static int
arbitrary_geometric_changes(const char * fname,
                            FileType input_type,
                            Molecule_Output_Object & output)
{
  assert (nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return arbitrary_geometric_changes(input, output);
}
static int
arbitrary_geometric_changes(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lo:S:T:ern:t:b:ca:s:z");

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

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c'))
  {
    set_append_coordinates_after_each_atom(1);

    if (verbose)
      cerr << "Will append atomic coordinates after each atom in a smiles\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', number_transformations_per_molecule) || number_transformations_per_molecule < 1)
    {
      cerr << "The number of changes per molecule must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will perform " << number_transformations_per_molecule << " transformations per molecule\n";
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', bump_check) || bump_check < 0.0)
    {
      cerr << "The bump check option (-b) must be a valid distance\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will reject molecules with atoms closer than " << bump_check << endl;
  }

  if (cl.option_present('r'))
  {
    arbitrary_rotations = 1;
    if (verbose)
      cerr << "Will perform arbitrary rotations\n";
  }

  if (cl.option_present('e'))
  {
    invert_molecule = 1;
    if (verbose)
      cerr << "Will invert the molecule\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', max_molecular_translation) || max_molecular_translation <= 0.0)
    {
      cerr << "The maximum molecular translation value must be a valid distance\n";
      usage(4);
    }

    if (verbose)
      cerr << "Molecular translations as far as " << max_molecular_translation << " performed\n";

    arbitrary_molecular_translations = 1;
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', max_atomic_translation) || max_atomic_translation <= 0.0)
    {
      cerr << "The maximum atomic_translation value must be a valid distance\n";
      usage(4);
    }

    if (verbose)
      cerr << "Atomic translations as far as " << max_atomic_translation << " performed\n";

    arbitrary_atomic_translations = 1;
  }

  if (cl.option_present('a'))
  {
    if (! cl.value('a', max_atoms_translated_per_molecule) || max_atoms_translated_per_molecule < 1)
    {
      cerr << "The number of atoms translated at once must be a valid +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will move a maximum of " << max_atoms_translated_per_molecule << " atoms\n";

    arbitrary_atomic_translations++;
  }

  if (cl.option_present('s'))
  {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i)
    {
      Substructure_Query * q = new Substructure_Query();
      if (! q->create_from_smarts(s))
      {
        delete q;
        cerr << "Invalid smarts '" << s << "'\n";
        return 1;
      }

      only_translate_these_atoms.add(q);
    }

    if (verbose)
      cerr << "Defined " << only_translate_these_atoms.number_elements() << " queries for atoms to translate\n";

    arbitrary_atomic_translations++;
  }

  if (cl.option_present('z'))
  {
    write_unchanged_molecules = 0;

    if (verbose)
      cerr << "Unchanged molecules will not be written\n";
  }

  if (0 == arbitrary_rotations &&        // nothing specified
      0 == invert_molecule &&
      0 == arbitrary_molecular_translations &&
      0 == arbitrary_atomic_translations)
  {
    arbitrary_rotations = 1;
    arbitrary_molecular_translations = 1;
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
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SDF);
  else if (! output.determine_output_types(cl, 'o'))
  {
    cerr << "Cannot determine output type(s)\n";
    usage(4);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (! cl.option_present('S')) {
    output.new_stem("-");
  }
  else
  {
    const_IWSubstring s = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Canot overwrite input file(s) '" << s << "'\n";
      return 4;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot initialise output stream '" << s << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Output written to '" << s << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! arbitrary_geometric_changes(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = arbitrary_geometric_changes(argc, argv);

  return rc;
}
