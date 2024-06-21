// Idea from Niraj for doing a spatial displacement of
// part of a molecule by another group.
// Input must be 3d structures.

#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace spatial_replacement {

using std::cerr;

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Performs a spatial displacement, by overlaying a replacement fragment on an existing\n";
  cerr << "molecule, and forming all possible bonds.\n";
  cerr << " -R <fname>  file containing replacement fragments\n";
  cerr << "             explicit Hydrogens are removed from the replacement\n";
  cerr << " -X <rad>    remove starting molecule atoms within <rad> of replacement atoms\n";
  cerr << " -r <rad>    make bonds between all starting molecule atoms with <rad> of replacement atom\n";
  cerr << " -m <natoms> discard replacements with fewer than <natoms> heavy atoms\n";
  cerr << " -M <natoms> discard replacements with more  than <natoms> heavy atoms\n";
  cerr << " -e <natoms> the maximum number of atoms that can be lost from the starting molecule\n";
  cerr << " -W ...      control what gets written, enter '-W help'\n";
  cerr << " -t <maxdist> randomly translate the replacement by as much as <maxdist>\n";
  cerr << " -h          just before writing, make implicit Hydrogens explicit\n";
  cerr << " -j          set the bond length the first time a replacement atom is bonded\n";
  cerr << " -I          put isotopic labels on the replacements\n";
  cerr << " -u <dist>   discard products where any pair of atoms are closer than <dist>\n";
  cerr << " -v          verbose output\n";
  // clang-format on

  ::exit(rc);
}

/*
  There are two very different cases of how the fragment is
  embedded onto the starting molecule.
  It can be singly connected, or doubly connected.

  In the singly connected case, we should be able to identify
  a bond that has zero removed atoms on one side, and all the
  removed bonds on the other side.

  In the doubly connected case, we identify an atom that
  will be replaced
*/

// clang-format off
// How does the replacement fragment overlap with the starting
// molecule. Two cases:
//   1. the closest starting molecule atom is 2 connected.
//   2. the closest starting molecule atom is 3 connected.
// clang-format on
struct Overlap {
  // The case of just one bond being broken. We will create
  // a bond from this atom to `join_to`.
  int retained = -1;

  // The case where a starting molecule atom is removed
  int starting_atom_lost;
  // The two remaining starting molecules atoms that were
  // attached to starting_atom_lost. They will be bonded
  // to `join_to`.
  int retained1 = -1;
  int retained2 = -1;

  // In both cases, the atom in the replacement to which we
  // will be making bond(s).
  int join_to = -1;
};

std::ostream&
operator<<(std::ostream& output, const Overlap& overlap) {
  if (overlap.retained >= 0) {
    output << "single attachment via " << overlap.retained;
  } else {
    output << "remove " << overlap.starting_atom_lost << " keep " << overlap.retained1
           << ' ' << overlap.retained2;
  }
  if (overlap.join_to >= 0) {
    output << " join to " << overlap.join_to;
  }

  return output;
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
 private:
  int _verbose = 0;

  FileType _input_type = FILE_TYPE_INVALID;

  int _reduce_to_largest_fragment = 0;

  int _remove_chirality = 0;

  Chemical_Standardisation _chemical_standardisation;

  // Not a part of all applications, just an example...
  Element_Transformations _element_transformations;

  int _molecules_read = 0;

  int _molecules_formed;

  float _exclusion_radius;

  float _bond_formation_radius;

  // If no atoms in the replacement are close enough to
  // displace atoms in the scaffold.
  int _molecules_not_close_to_replacement;

  // When we limit the number of starting molecule atoms that
  // can be lost, we keep track of that.
  int _molecules_losing_too_many_atoms;

  // Once we remove all the atoms in the scaffold, maybe there
  // are no atoms close enough to form a bond.
  int _closest_pair_too_far_apart;

  // The replacement fragments, the -R option. For each
  // -R instance, there will be a set of molecules.
  resizable_array_p<resizable_array_p<Molecule>> _replacement;

  // We can impose atom count limits on the replacements.
  int _min_replacement_size;
  int _max_replacement_size;

  // In order to avoid chopping off chunks of the starting molecule
  // that are too large, we can impose a limit on the max number of
  // atoms lost.
  int _max_atoms_lost;

  int _make_implicit_hydrogens_explicit;

  int _write_replacement_to_output;
  int _write_scaffold_to_output;

  int _molecules_written = 0;

  // For each starting molecule, record how many replacement
  // molecules are generated.
  extending_resizable_array<int> variants_per_starting_molecule;

  // Count the number of atoms lost per starting molecule.
  extending_resizable_array<int> _atoms_lost;

  // If requested, we can set a bond length each time a bond
  // is added to a replacement atom.
  bool _set_bond_lengths;

  float _default_bond_length;

  float _bump_check;
  int _failed_bump_check;

  // We can generate more possibilities by randomly translating
  // the replacement by as far as `_random_translations` in any
  // direction.
  int _number_random_translations;
  float _random_translations;
  std::mt19937 _rng;

  // If we are doing random perturbations of the replacement
  // we need to keep track of what we have already formed.
  IW_STL_Hash_Set _seen;

  // Private functions.
  int Process(Molecule& m, const resizable_array<Molecule*>& replacement,
              Molecule_Output_Object& output);
  int ReadReplacementMoleculesFileOfFiles(IWString& fname,
                                          resizable_array_p<Molecule>& destination);
  void DisplayDashWQualifiers(std::ostream& output) const;

  std::optional<std::tuple<int, int>> IdentifyKeyBond(Molecule& m,
                                                      const float* distance_matrix,
                                                      int initial_natoms,
                                                      int* remove_atom);
  int PlausibleKeyAtom(const Molecule& m, atom_number_t zatom,
                       const float* distance_matrix, int initial_natoms);
  int BondedToTwoRetainedAtoms(const Molecule& m, atom_number_t zatom,
                               const int* remove_atom) const;
  int IdentifyDoublyConnected(const Molecule& m, int initial_natoms,
                              const int* remove_atom, int* visited,
                              Overlap& overlap) const;
  int IdentifySinglyConnectedReplacement(Molecule& m, int initial_natoms,
                                         const int* remove_atom, int* visited,
                                         Overlap& overlap) const;
  std::optional<Overlap> IdentifySite(Molecule& m, int initial_natoms,
                                      const int* remove_atom, int* visited) const;

  void DoRandomTranslation(Molecule& m);

  int OkBumpCheck(const Molecule& m, int initial_matoms);

  int MaybeWriteScaffoldAndReplacement(Molecule& starting_molecule,
                                       const resizable_array<Molecule*>& replacement,
                                       Molecule_Output_Object& output);

 public:
  Options();

  // Get user specified command line directives.
  int Initialise(Command_Line& cl);

  int verbose() const {
    return _verbose;
  }

  // After each molecule is read, but before any processing
  // is attempted, do any preprocessing transformations.
  int Preprocess(Molecule& m);

  int ReadReplacementMolecules(data_source_and_type<Molecule>& input,
                               resizable_array_p<Molecule>& destination);
  int ReadReplacementMolecules(IWString& fname, resizable_array_p<Molecule>& destination);

  // Helpful when the -i option is not given.
  int MaybeDiscernInputType(const char* fname);

  FileType input_type() const {
    return _input_type;
  }

  // The function that actually does the processing,
  // and may write to `output`.
  // You may instead want to use a Molecule_Output_Object if
  // Molecules are being written.
  // You may choose to use a std::ostream& instead of
  // IWString_and_File_Descriptor.
  int Process(Molecule& mol, Molecule_Output_Object& output);

  // After processing, report a summary of what has been done.
  int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _molecules_formed = 0;

  _exclusion_radius = 1.0;
  _bond_formation_radius = 2.0;

  _min_replacement_size = 0;
  _max_replacement_size = std::numeric_limits<int>::max();

  _max_atoms_lost = std::numeric_limits<int>::max();

  _molecules_not_close_to_replacement = 0;
  _molecules_losing_too_many_atoms = 0;
  _closest_pair_too_far_apart = 0;

  _number_random_translations = 0;
  _random_translations = 0.0f;
  _molecules_written = 0;

  _make_implicit_hydrogens_explicit = 0;

  _write_scaffold_to_output = 0;
  _write_replacement_to_output = 0;

  _set_bond_lengths = false;
  _default_bond_length = 1.43;

  _bump_check = -1.0f;
  _failed_bump_check = 0;

  std::random_device rd;
  _rng.seed(rd());
}

void
Options::DisplayDashWQualifiers(std::ostream& output) const {
  output << " -W rep       write the replacement before each product\n";
  output << " -W once      write the starting molecule once\n";
  output << " -W each      write the starting molecule before each product\n";
  ::exit(0);
}

int
Options::Initialise(Command_Line& cl) {
  set_copy_name_in_molecule_copy_constructor(1);

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T')) {
      Usage(8);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (!cl.option_present('R')) {
    cerr << "Must specify replacement molecules via the -R option\n";
    Usage(1);
  }

  if (cl.option_present('X')) {
    if (!cl.value('X', _exclusion_radius) || _exclusion_radius < 0.0f) {
      cerr << "The initial exclusion radius (-X) must be a valid distance\n";
      return 1;
    }
    if (_verbose) {
      cerr << "Will remove starting molecule atoms within " << _exclusion_radius
           << " of an atom in the replacement\n";
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', _bond_formation_radius) || _bond_formation_radius < 0.0f) {
      cerr << "The bond formation radius (-r) must be a valid distance\n";
      return 1;
    }
    if (_verbose) {
      cerr << "Will make bonds between atoms if closer than " << _bond_formation_radius
           << '\n';
    }
  }

  if (cl.option_present('W')) {
    const_IWSubstring w;
    for (int i = 0; cl.value('W', w, i); ++i) {
      if (w == "rep") {
        _write_replacement_to_output = 1;
      } else if (w == "once") {
        _write_scaffold_to_output = 1;
      } else if (w == "each") {
        _write_scaffold_to_output = 2;
      } else if (w == "help") {
        DisplayDashWQualifiers(cerr);
      } else {
        cerr << "Unrecognised -w qualifier '" << w << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('t')) {
    const_IWSubstring t;
    for (int i = 0; cl.value('t', t, i); ++i) {
      if (t.starts_with("n=")) {
        t.remove_leading_chars(2);
        if (!t.numeric_value(_number_random_translations) ||
            _number_random_translations < 1) {
          cerr << "The number of random translations must be a whole +ve number\n";
          return 0;
        }
      } else if (t.starts_with("seed=")) {
        t.remove_leading_chars(5);
        unsigned int seed;
        if (!t.numeric_value(seed)) {
          cerr << "The rng seed must be a valid unsigned\n";
          return 0;
        }
        _rng.seed(seed);
        if (_verbose) {
          cerr << "using seed " << seed << '\n';
        }
      } else if (!t.numeric_value(_random_translations) || _random_translations <= 0.0f) {
        cerr << "The magnitude of random translations must be > 0.0\n";
        return 0;
      }
    }

    if (_random_translations == 0.0f) {
      cerr << "No magnitude specified for random translations\n";
      return 0;
    }
    if (_number_random_translations == 0) {
      _number_random_translations = 10;
    }

    if (_verbose) {
      cerr << "Will randomly translate the replacement as much as "
           << _random_translations << " angstroms " << _number_random_translations
           << " times\n";
    }
  }

  if (cl.option_present('h')) {
    _make_implicit_hydrogens_explicit = 1;
    if (_verbose) {
      cerr << "Will make implicit Hydrogens explicit\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', _min_replacement_size) || _min_replacement_size < 0) {
      cerr << "The minimum replacement molecule size (-m) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Replacements must have at least " << _min_replacement_size << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', _max_replacement_size) ||
        _max_replacement_size < _min_replacement_size) {
      cerr << "The minimum replacement molecule size (-m) must be a whole number > "
           << _min_replacement_size << "\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Replacements must have at most " << _max_replacement_size << " atoms\n";
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', _max_atoms_lost) || _max_atoms_lost < 0) {
      cerr << "The max starting molecule atoms lost value (-e) must be a whole +ve "
              "number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Starting molecules can lose a max of " << _max_atoms_lost << " atoms\n";
    }
  }

  if (cl.option_present('i')) {
    if (!process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SDF;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.option_present('R')) {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      std::unique_ptr<resizable_array_p<Molecule>> r =
          std::make_unique<resizable_array_p<Molecule>>();

      if (!ReadReplacementMolecules(fname, *r)) {
        cerr << "Cannot read replacement molecules from '" << fname << "'\n";
        return 0;
      }
      _replacement << r.release();
    }

    if (_verbose) {
      cerr << "Read " << _replacement.size() << " sets of replacement fragments\n";
      extending_resizable_array<int> ratoms;
      for (const auto* mols : _replacement) {
        for (const Molecule* m : *mols) {
          ++ratoms[m->natoms()];
        }
      }
      for (int i = 0; i < ratoms.number_elements(); ++i) {
        if (ratoms[i]) {
          cerr << ratoms[i] << " replacement molecules had " << i << " atoms\n";
        }
      }
    }

    if (_replacement.empty()) {
      cerr << "No replacement fragments read\n";
      return 0;
    }
    for (const auto& rep : _replacement) {
      if (rep->empty()) {
        cerr << "Set of replacements is empty\n";
        return 0;
      }
    }
  }

  if (cl.option_present('I')) {
    for (int i = 0; i < _replacement.number_elements(); ++i) {
      for (Molecule* m : *_replacement[i]) {
        const int matoms = m->natoms();
        for (int j = 0; j < matoms; ++j) {
          m->set_isotope(j, i + 1);
        }
      }
    }
  }

  if (cl.option_present('j')) {
    _set_bond_lengths = true;
    if (_verbose) {
      cerr << "Will set a bond length each time a new replacement atom is bonded\n";
    }
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', _bump_check) || _bump_check < 0.0) {
      cerr << "The bump check option (-u) must be a valid distance\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will discard molecules with atoms closer than " << _bump_check << '\n';
    }
  }

  return 1;
}

// `fname` is a file containing the names of replacement
// molecules.
int
Options::ReadReplacementMolecules(IWString& fname,
                                  resizable_array_p<Molecule>& destination) {
  if (fname.starts_with("F:")) {
    fname.remove_leading_chars(2);
    return ReadReplacementMoleculesFileOfFiles(fname, destination);
  }

  if (_input_type == FILE_TYPE_INVALID) {
    if (!MaybeDiscernInputType(fname.null_terminated_chars())) {
      return 0;
    }
  }

  data_source_and_type<Molecule> input(_input_type, fname);
  if (!input.good()) {
    cerr << "Options::ReadReplacementMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReplacementMolecules(input, destination);
}

int
Options::ReadReplacementMoleculesFileOfFiles(IWString& fname,
                                             resizable_array_p<Molecule>& destination) {
  iwstring_data_source input;
  if (!input.open(fname.null_terminated_chars())) {
    cerr << "Options::ReadReplacementMoleculesFileOfFiles:cannot open '" << fname
         << "'\n";
    return 0;
  }

  IWString buffer;
  while (input.next_record(buffer)) {
    if (!ReadReplacementMolecules(buffer, destination)) {
      cerr << "Options::ReadReplacementMoleculesFileOfFiles:cannot read '" << buffer
           << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadReplacementMolecules(data_source_and_type<Molecule>& input,
                                  resizable_array_p<Molecule>& destination) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    m->reduce_to_largest_fragment_carefully();
    // Explicit hydrogens on the replacement are too hard to handle.
    m->remove_all(1);
    if (m->natoms() < _min_replacement_size || m->natoms() > _max_replacement_size) {
      continue;
    }

    destination << m;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules, formed " << _molecules_formed
         << '\n';
  output << _molecules_not_close_to_replacement
         << " replacements not close to starting mol\n";
  output << _closest_pair_too_far_apart << " replacements not close to starting atom\n";
  output << _molecules_losing_too_many_atoms
         << " molecules discarded for losing more than " << _max_atoms_lost << " atoms\n";
  output << _failed_bump_check << " molecules failed bump check " << _bump_check << '\n';
  output << "Wrote " << _molecules_written << " molecules\n";

  for (int i = 0; i < _atoms_lost.number_elements(); ++i) {
    if (_atoms_lost[i]) {
      output << _atoms_lost[i] << " molecules lost " << i << " atoms\n";
    }
  }

  return 1;
}

// If the input type is known, return it.
// Otherwise examine the file name's suffix to
// determine the type.
int
Options::MaybeDiscernInputType(const char* fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  // Makes things too hard otherwise.
  m.remove_all(1);

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m, Molecule_Output_Object& output) {
  ++_molecules_read;

  _seen.clear();

  if (!m.valence_ok()) {
    cerr << "Starting molecule bad valence '" << m.name() << "'\n";
  }

  // Just so we can process things with a single for loop.
  if (_number_random_translations == 0) {
    _number_random_translations = 1;
  }

  // Even if not doing translations, just easier to get them.
  std::unique_ptr<float[]> unperturbed = m.GetCoordinates();

  // The number of items in each R replacement.
  std::vector<int> count;
  for (const auto* mols : _replacement) {
    count.push_back(mols->size());
  }

  combinations::Combinations combinations(count);

  std::vector<int> state;
  state.resize(_replacement.size(), 0);

  do {
    resizable_array<Molecule*> replacements;
    replacements.reserve(_replacement.number_elements());
    for (int i = 0; i < _replacement.number_elements(); ++i) {
      replacements << _replacement[i]->item(state[i]);
    }
    for (int i = 0; i < _number_random_translations; ++i) {
      if (i > 0) {
        DoRandomTranslation(m);
      }
      Process(m, replacements, output);
      if (i > 0) {
        m.SetCoordinates(unperturbed.get());
      }
    }
  } while (combinations.Next(state));

  return 1;
}

// the only reason this is not const is that the random number
// generator changes state.
void
Options::DoRandomTranslation(Molecule& m) {
  std::uniform_real_distribution<float> u(-_random_translations, _random_translations);
  coord_t x = u(_rng);
  coord_t y = u(_rng);
  coord_t z = u(_rng);

  m.translate_atoms(x, y, z);
}

// Maybe return the atom number of the closest bonded neighbour
// of `zatom`.
// Ignore neighbours that are in the same ring.
std::optional<atom_number_t>
ClosestBondedAtom(Molecule& m, atom_number_t zatom) {
  float shortest_distance = std::numeric_limits<float>::max();
  atom_number_t closest_atom = INVALID_ATOM_NUMBER;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t o = b->other(zatom);
    if (m.in_same_ring(zatom, o)) {
      continue;
    }

    float d = m.distance_between_atoms(zatom, o);
    if (d < shortest_distance) {
      shortest_distance = d;
      closest_atom = o;
    }
  }

  if (closest_atom == INVALID_ATOM_NUMBER) {
    return std::nullopt;
  }

  return closest_atom;
}

// Atoms in `remove_atom` will be removed from `m`. In order
// to make bonds, we record the closest atom to the removed
// atom. If there is no viable bonded neighbour to get the
// bond, that will be set to INVALID_ATOM_NUMBER.
// Note that `matoms` is likely smaller than m.natoms().
std::unique_ptr<int[]>
ClosestToRemoved(Molecule& m, int matoms, int* remove_atom) {
  std::unique_ptr<int[]> result(new_int(matoms, INVALID_ATOM_NUMBER));
  for (int i = 0; i < matoms; ++i) {
    if (!remove_atom[i]) {
      continue;
    }
    if (m.ncon(i) == 0) {
      continue;
    }

    std::optional<atom_number_t> adj = ClosestBondedAtom(m, i);
    if (!adj) {
      continue;
    }
    result[i] = *adj;
  }

  return result;
}

// Set the name of `destination` to a composite of the parent names.
void
SetName(const Molecule& parent1, const resizable_array<Molecule*>& replacements,
        Molecule& destination) {
  static constexpr char kSep = '!';

  IWString new_name;
  new_name << parent1.name();
  for (const Molecule* r : replacements) {
    new_name << kSep << r->name();
  }
  destination.set_name(new_name);
}

// Given a set of atoms that are to be removed from `m`
// extend the `remove_atom` to include all atoms in ring
// systems containing those initially marked atoms.
int
ExtendToRingSystems(Molecule& m, int* remove_atom) {
  // Force sssr determination.
  m.ring_membership();

  // Collect all the fused system identifiers that the
  // marked atoms touch.
  resizable_array<int> fsid;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (remove_atom[i] == 0) {
      continue;
    }
    if (!m.is_ring_atom(i)) {
      continue;
    }

    const int f = m.fused_system_identifier(i);
    fsid.add_if_not_already_present(f);
  }

  if (fsid.empty()) {
    return 0;
  }

  for (const Ring* r : m.sssr_rings()) {
    if (!fsid.contains(r->fused_system_identifier())) {
      continue;
    }
    r->set_vector(remove_atom, 1);
  }

  return 1;
}

int
Options::PlausibleKeyAtom(const Molecule& m, atom_number_t zatom,
                          const float* distance_matrix, int initial_natoms) {
  if (m.ncon(zatom) < 2) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = initial_natoms; i < matoms; ++i) {
    float d = distance_matrix[zatom * matoms + i];
    if (d < _bond_formation_radius) {
      return 1;
    }
  }

  return 0;
}

// `values` will be 0 or 1 indicating a set membership.
// For each set, count the number of times `lost` is nonzero.
std::tuple<int, int>
CountLost(const int* values, const int n, const int* lost) {
  int s1 = 0;
  int s2 = 0;
  for (int i = 0; i < n; ++i) {
    if (values[i] == 0) {
      if (lost[i] == 1) {
        ++s1;
      }
    } else {
      if (lost[i] == 1) {
        ++s2;
      }
    }
  }

  return std::make_tuple(s1, s2);
}

int
Options::BondedToTwoRetainedAtoms(const Molecule& m, atom_number_t zatom,
                                  const int* remove_atom) const {
  int retained = 0;
  const Atom& a = m.atom(zatom);
  for (const Bond* bond : a) {
    atom_number_t j = bond->other(zatom);
    if (remove_atom[j] == 0) {
      ++retained;
    }
  }

  return retained == 2;
}

// Count the number of items in `remove_atom` that are set.
int
CountRemoved(const Molecule& m, atom_number_t zatom, const int* remove_atom,
             int* visited) {
  visited[zatom] = 1;
  int result = (remove_atom[zatom] == 1);

  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    if (visited[j]) {
      continue;
    }
    result += CountRemoved(m, j, remove_atom, visited);
  }

  return result;
}

// #define DEBUG_IDENTIFYDOUBLYCONNECTED

// Identify a 3 connected atom where two of the attachments
// have zero removed atoms, and the other has removed atoms.
// Returns true of successful and fills in `overlap`.
int
Options::IdentifyDoublyConnected(const Molecule& m, int initial_natoms,
                                 const int* remove_atom, int* visited,
                                 Overlap& overlap) const {
  overlap.starting_atom_lost = -1;

  // Look for atoms that will be removed that have two retained
  // neighbours.
  for (int i = 0; i < initial_natoms; ++i) {
    if (remove_atom[i] == 0) {
      continue;
    }

    const Atom& a = m.atom(i);
    if (a.ncon() != 3) {
      continue;
    }

    if (!BondedToTwoRetainedAtoms(m, i, remove_atom)) {
      continue;
    }

    // The atoms that will be retained. They have zero atoms
    // to be removed.
    Set_of_Atoms nbrs;
    // The attached atom that contains atoms being removed.
    atom_number_t greater_zero = INVALID_ATOM_NUMBER;

    for (const Bond* bond : a) {
      atom_number_t j = bond->other(i);
      std::fill_n(visited, m.natoms(), 0);
      visited[i] = 1;
      if (int c = CountRemoved(m, j, remove_atom, visited); c == 0) {
        nbrs << j;
      } else {
        greater_zero = j;
      }
    }

#ifdef DEBUG_IDENTIFYDOUBLYCONNECTED
    cerr << "atom " << i << " has " << nbrs.size() << " zero removed nbrs, gtzer "
         << greater_zero << '\n';
#endif

    // There must be two connections with zero atoms lost and one
    // that has one or more lost atoms.
    if (nbrs.size() != 2 || greater_zero == INVALID_ATOM_NUMBER) {
      continue;
    }

    overlap.starting_atom_lost = i;
    overlap.retained1 = nbrs[0];
    overlap.retained2 = nbrs[1];
    return 1;
  }

  return 0;
}

int
Options::IdentifySinglyConnectedReplacement(Molecule& m, int initial_natoms,
                                            const int* remove_atom, int* visited,
                                            Overlap& overlap) const {
  overlap.retained = -1;

  m.ring_membership();

  for (const Bond* b : m.bond_list()) {
    if (b->nrings()) {
      continue;
    }
    if (!b->is_single_bond()) {
      continue;
    }

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    // Looking for a split bond.
    if (remove_atom[a1] == remove_atom[a2]) {
      continue;
    }
    std::fill_n(visited, initial_natoms, 0);
    visited[a1] = 1;
    const int c2 = CountRemoved(m, a2, remove_atom, visited);
    std::fill_n(visited, initial_natoms, 0);
    visited[a2] = 1;
    const int c1 = CountRemoved(m, a1, remove_atom, visited);
    // cerr << "atoms " << a1 << " and " << a2 << " counts " << c1 << " and " << c2 <<
    // '\n';

    if (c1 == 0 && c2 > 0) {
      overlap.retained = a1;
      return 1;
    } else if (c1 > 0 && c2 == 0) {
      overlap.retained = a2;
      return 1;
    }
  }

  return 0;
}

// The atoms in the starting molecule in `overlap` have been
// determined, we need to set the `join_to` value. How that gets
// done is a function of what kind of join is being done.
int
DetermineJoinTo(Molecule& m, int initial_natoms, const float* distance_matrix,
                Overlap& overlap) {
  const int matoms = m.natoms();

  atom_number_t in_starting_molecule;
  if (overlap.retained >= 0) {
    in_starting_molecule = overlap.retained;
  } else {
    in_starting_molecule = overlap.starting_atom_lost;
  }

  float closest_distance = std::numeric_limits<float>::max();
  for (int i = initial_natoms; i < matoms; ++i) {
    const float d = distance_matrix[in_starting_molecule * matoms + i];
    if (d < closest_distance) {
      closest_distance = d;
      overlap.join_to = i;
    }
  }

  return 1;
}

// We want to alternatively check first for a single connection
// replacement, and a doubly connected replacement.
std::optional<Overlap>
Options::IdentifySite(Molecule& m, int initial_natoms, const int* remove_atom,
                      int* visited) const {
  static bool single_first = false;

  // single_first = ! single_first;

  Overlap result;

  if (single_first) {
    if (IdentifySinglyConnectedReplacement(m, initial_natoms, remove_atom, visited,
                                           result)) {
      return result;
    }
    if (IdentifyDoublyConnected(m, initial_natoms, remove_atom, visited, result)) {
      return result;
    }
  } else {
    if (IdentifyDoublyConnected(m, initial_natoms, remove_atom, visited, result)) {
      return result;
    }
    if (IdentifySinglyConnectedReplacement(m, initial_natoms, remove_atom, visited,
                                           result)) {
      return result;
    }
  }

  return std::nullopt;
}

int
Options::OkBumpCheck(const Molecule& m, int initial_matoms) {
  // If no bump check in effect, must be OK.
  if (_bump_check < 0.0f) {
    return 1;
  }

  const int matoms = m.natoms();
  for (int i = 0; i < initial_matoms; ++i) {
    for (int j = initial_matoms; j < matoms; ++j) {
      float d = m.distance_between_atoms(i, j);
      if (d <= _bump_check) {
        ++_failed_bump_check;
        return 0;
      }
    }
  }

  return 1;
}

// This is not working. This usually detects a bad valence,
// but when I run the resulting molecules through fileconv
// there are no valence problems. Do not have time to track down...
void
ComplainIfBadValence(Molecule& m) {
  return;  // Turn off until this can be debugged.

  if (m.valence_ok()) {
    return;
  }

  cerr << "Warning bad valence '" << m.name() << "'\n";
  for (int i = 0; i < m.natoms(); ++i) {
    if (m.valence_ok(i)) {
      continue;
    }
    cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " bad\n";
  }
}

// #define DEBUG_SPATIAL_REPLACEMENT

// Systematically add the molecules in `replacement` to `starting_molecule`.
int
Options::Process(Molecule& starting_molecule,
                 const resizable_array<Molecule*>& replacement,
                 Molecule_Output_Object& output) {
  // Do not destroy our input.
  Molecule m(starting_molecule);

  const int initial_natoms = m.natoms();

  // Keep track of how many atoms are in `m` as the various
  // replacements are added. `initial_natoms` will be the
  // same as `atom_count[0]`.
  resizable_array<int> atom_count;

  // Add the replacement fragments.
  for (const Molecule* rep : replacement) {
    atom_count << m.natoms();
    m.add_molecule(rep);
  }

  ++_molecules_formed;

  // cerr << "After adding replacements " << m.smiles() << '\n';
  // output.write(m);

#ifdef DEBUG_SPATIAL_REPLACEMENT
  Molecule qqq(m);
  for (int i = 0; i < qqq.natoms(); ++i) {
    qqq.set_isotope(i, i);
  }
  output.write(qqq);
#endif

  // The atom count after replacement is added.
  const int matoms = m.natoms();

  std::unique_ptr<float[]> distance_matrix = std::make_unique<float[]>(matoms * matoms);
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      distance_matrix[i * matoms + j] = distance_matrix[j * matoms + i] =
          m.distance_between_atoms(i, j);
    }
  }

#ifdef DEBUG_SPATIAL_REPLACEMENT
  cerr << "Overlapped " << write_isotopically_labelled_smiles(m, false, cerr);
  cerr << ' ' << m.name() << '\n';
  for (int i = 0; i < initial_natoms; ++i) {
    float mindist = std::numeric_limits<float>::max();
    for (int j = initial_natoms; j < matoms; ++j) {
      cerr << " dist " << i << ' ' << j << ' ' << distance_matrix[i * matoms + j] << '\n';
      if (distance_matrix[i * matoms + j] < mindist) {
        mindist = distance_matrix[i * matoms + j];
      }
    }
    cerr << i << ' ' << mindist << '\n';
  }
#endif

  // Whether or not each atom in `m` will be removed.
  std::unique_ptr<int[]> remove_atom(new_int(matoms));
  // Whether or not an atom might possibly be involved in
  // forming a new bond.
  std::unique_ptr<int[]> atom_could_be_bonded(new_int(matoms * matoms));
  // Whether or not a pair of atoms might possibly be involved in
  // forming a new bond.
  std::unique_ptr<int[]> pair_could_be_bonded(new_int(matoms * matoms));

  // Count number atoms being removed.
  int removed = 0;
  int possible_bonds = 0;

  for (int i = 0; i < initial_natoms; ++i) {
    for (int j = initial_natoms; j < matoms; ++j) {
      if (distance_matrix[i * matoms + j] <= _exclusion_radius) {
        remove_atom[i] = 1;
        ++removed;
      }
      if (distance_matrix[i * matoms + j] <= _bond_formation_radius) {
        pair_could_be_bonded[i * matoms + j] = 1;
        pair_could_be_bonded[j * matoms + i] = 1;
        atom_could_be_bonded[i] = 1;
        atom_could_be_bonded[j] = 1;
        ++possible_bonds;
      }
    }
  }

  if (possible_bonds == 0) {
    cerr << m.name() << " no atoms within " << _bond_formation_radius << " for bond formation\n";
    ++_molecules_not_close_to_replacement;
    return 1;
  }

  if (removed > _max_atoms_lost) {
    ++_molecules_losing_too_many_atoms;
    if (_verbose > 1) {
      cerr << m.name() << " too many atoms lost " << removed << '\n';
      return 0;
    }
  }

#ifdef DEBUG_SPATIAL_REPLACEMENT
  for (int i = 0; i < initial_natoms; ++i) {
    if (remove_atom[i]) {
      cerr << " remove " << i << ' ' << m.smarts_equivalent_for_atom(i) << ' ' << m.x(i)
           << ',' << m.y(i) << ',' << m.z(i) << '\n';
    }
  }
#endif

  // If an atom to be removed is in a ring, remove all atoms
  // in that ring system.
  ExtendToRingSystems(m, remove_atom.get());
  // Note that we might be breaking a double bond in the starting molecule.

  // Unconnect all the atoms that will be removed. Don't delete them
  // yet because we want to preserve our atom numbering.
  for (int i = 0; i < matoms; ++i) {
    if (remove_atom[i]) {
      m.remove_bonds_to_atom(i);
    }
  }

  // Add all plausible bonds from retained atoms in the starting
  // molecule, and atoms in the replacements.

  // If setting bond lengths, we only do it once per
  // replacement atom added.
  std::unique_ptr<int[]> set_bond_length;
  if (_set_bond_lengths) {
    set_bond_length.reset(new_int(matoms, 1));
  }

  for (int i = 0; i < initial_natoms; ++i) {
    if (remove_atom[i]) {
      continue;
    }
    if (!atom_could_be_bonded[i]) {
      continue;
    }
    if (m.hcount(i) == 0) {
      continue;
    }

    for (int j = initial_natoms; j < matoms; ++j) {
      if (!atom_could_be_bonded[j]) {
        continue;
      }
      if (!pair_could_be_bonded[i * matoms + j]) {
        continue;
      }
      if (m.hcount(i) == 0 || m.hcount(j) == 0) {
        continue;
      }

      m.add_bond(i, j, SINGLE_BOND);

      if (set_bond_length && set_bond_length[j]) {
        if (!m.in_same_ring_system(i, j)) {
          m.set_bond_length(i, j, _default_bond_length);
          set_bond_length[j] = 0;
        }
      }
    }
  }

  // Just before writing, remove the old atoms.
  if (_verbose) {
    int c = count_non_zero_occurrences_in_array(remove_atom.get(), initial_natoms);
    ++_atoms_lost[c];
  }

  // Keep track of the number of atoms so we can check how many are lost.
  int atoms_before = m.natoms();

  m.remove_atoms(remove_atom.get());

  m.reduce_to_largest_fragment();

  if (int atoms_lost = atoms_before - m.natoms(); atoms_lost > _max_atoms_lost) {
    if (_verbose > 1) {
      cerr << m.name() << " loses too many atoms " << atoms_lost << '\n';
    }
    ++_molecules_losing_too_many_atoms;
    return 0;
  }

  if (_make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  if (_number_random_translations > 0) {
    if (_seen.contains(m.unique_smiles())) {
      return 0;
    }
    _seen.insert(m.unique_smiles());
  }

  MaybeWriteScaffoldAndReplacement(starting_molecule, replacement, output);

  SetName(m, replacement, m);

  ComplainIfBadValence(m);

  if (!OkBumpCheck(m, initial_natoms)) {
    return 1;
  }

  output.write(m);
  ++_molecules_written;

  return 1;
}

int
Options::MaybeWriteScaffoldAndReplacement(Molecule& starting_molecule,
                                          const resizable_array<Molecule*>& replacement,
                                          Molecule_Output_Object& output) {
  if (_write_scaffold_to_output == 2) {
    output.write(starting_molecule);
  }
  if (_write_replacement_to_output) {
    for (Molecule* m : replacement) {
      output.write(*m);
    }
  }

  return 1;
}

int
SpatialReplacement(Options& options, Molecule& m, Molecule_Output_Object& output) {
  return options.Process(m, output);
}

int
SpatialReplacement(Options& options, data_source_and_type<Molecule>& input,
                   Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!options.Preprocess(*m)) {
      continue;
    }

    if (!SpatialReplacement(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
SpatialReplacement(Options& options, const char* fname, Molecule_Output_Object& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (!input.good()) {
    cerr << "SpatialReplacement:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return SpatialReplacement(options, input, output);
}

int
SpatialReplacement(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:R:X:r:S:o:W:m:M:ht:e:Iju:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name stem with the -S option\n";
    Usage(1);
  }

  IWString output_stem = cl.string_value('S');

  Molecule_Output_Object output;
  if (!cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SDF);
  } else if (!output.determine_output_types(cl)) {
    cerr << "Cannot determine output types (-o)\n";
    return 1;
  }

  if (output.would_overwrite_input_files(cl, output_stem)) {
    cerr << "Cannot overwrite input '" << output_stem << "'\n";
    return 1;
  }

  if (!output.new_stem(output_stem)) {
    cerr << "Cannot open output '" << output_stem << "'\n";
    return 1;
  }

  for (const char* fname : cl) {
    if (!SpatialReplacement(options, fname, output)) {
      cerr << "SpatialReplacement::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace spatial_replacement

int
main(int argc, char** argv) {
  int rc = spatial_replacement::SpatialReplacement(argc, argv);

  return rc;
}
