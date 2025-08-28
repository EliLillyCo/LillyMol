// Given a set of molecules, identify all the sidechains and remove them.
// Then re-attach them to molecules in the set.

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_set.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace sidechain_switcheroo {

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
  cerr << R"(Reads a set of molecules and identifies all sidechains - as specified by a query.
Sidechains are detatched from the molecules and stored as a set.
The scaffolds retain their attachment points.
All available sites on the scaffolds are then enumerated with all combinations of sidechains.
Note that this can generate large numbers of molecules!
 -s <smarts>    smarts to define join between scaffold (first matched atom) and sidechain (second matched atom).
 -q ...         sidechain queries as query file(s).
 -z i           ignore molecules that do not match any of the -s queries.
 -x <n>         create molecules with as many as <n> sites permuted - beware combinatorics!
 -X <n>         keep only the <n> smallest sidechains found - helps to reduce the number of molecules formed.
 -I             remove isotopes from product molecules - isotopes are used internally to label atoms. 
 -V             report and discard product molecules with invalid valences - hopefully will not happen.
 -c             remove chirality as molecules are read.
 -r <n>         report progress every <n> products generated.
 -g ...         chemical standardisation options.
 -l             strip to largest fragment.
 -i ...         input options.
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

class MoleculeAndAtoms {
  private:
    // Note that we own the atom.
    Molecule* _m;

    Set_of_Atoms _atoms;

  public:
    MoleculeAndAtoms(Molecule* m);
    ~MoleculeAndAtoms();

    void set_molecule(Molecule* m) {
      _m = m;
    }
    const Molecule& molecule() const {
      return *_m;
    }
    Molecule& molecule() {
      return *_m;
    }

    void add_atom(atom_number_t a) {
      _atoms << a;
    }

    int number_reactive_sites() const {
      return _atoms.number_elements();
    }

    void SortAtoms();

    const Set_of_Atoms atoms() const {
      return _atoms;
    }
};

MoleculeAndAtoms::MoleculeAndAtoms(Molecule* m) : _m(m) {
}

MoleculeAndAtoms::~MoleculeAndAtoms() {
  delete _m;
}

void
MoleculeAndAtoms::SortAtoms() {
  if (_atoms.size() < 2) {
    return;
  }

  if (_atoms.size() == 2) {
    if (_atoms[0] > _atoms[1]) {
      _atoms.swap_elements(0, 1);
    }
    return;
  }

  std::sort(_atoms.begin(), _atoms.end(), [](int a1, int a2) {
    return a1 < a2;
  });
}

class MoleculeAndAtom {
  private:
    // Note that we own the atom.
    Molecule* _m;

    atom_number_t _atom;

  public:
    MoleculeAndAtom(Molecule* m);
    ~MoleculeAndAtom();

    void set_atom(atom_number_t s) {
      _atom = s;
    }
    atom_number_t atom() const {
      return _atom;
    }
    const Molecule& molecule() const {
      return *_m;
    }
    Molecule& molecule() {
      return *_m;
    }
};

MoleculeAndAtom::MoleculeAndAtom(Molecule* m) : _m(m) {
  _atom = kInvalidAtomNumber;
}

MoleculeAndAtom::~MoleculeAndAtom() {
  delete _m;
}

constexpr isotope_t kScaffoldAtom = 1;
constexpr isotope_t kSidechainAtom = 2;

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    resizable_array_p<Molecule> _molecules;

    // Scaffolds can have any number of attachment points.
    resizable_array_p<MoleculeAndAtoms> _scaffolds;

    // Sidechains have just one attachment point.
    resizable_array_p<MoleculeAndAtom> _sidechains;

    // We sort the _sidechains array by atom count, and can optionally
    // retain only the _max_sidechains smallest.
    uint32_t _max_sidechains;

    // Default is [aD3]-!@{a<10}[R0]
    resizable_array_p<Substructure_Query> _sidechain_query;

    int _max_scaffold_attachments;

    int _add_hydrogen_substituent;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint32_t _molecules_read = 0;

    int _ignore_molecules_not_matching_queries;
    uint32_t _molecules_not_matching_queries;

    uint64_t _molecules_written;

    absl::flat_hash_set<IWString> _product_seen;

    int _remove_isotopes_from_products;

    uint64_t _undesirable_adjacencies;

    int _discard_bad_valence;
    uint64_t _products_with_bad_valence;

    Report_Progress _report_progress;

  // Private functions

    int ReadMolecules(data_source_and_type<Molecule>& input);
    int SplitScaffoldSidechain(Molecule& m, absl::flat_hash_set<IWString>& seen);
    int SplitScaffoldSidechain(Molecule& m, const Substructure_Results& sresults, absl::flat_hash_set<IWString>& seen);

    int AddHydrogen();

    int CreateScaffold(const Molecule& m,
                        const Set_of_Atoms& scaffold_atom,
                        const Set_of_Atoms& sidechain_atom,
                        const std::vector<int>& state,
                        absl::flat_hash_set<IWString>& seen);
    uint32_t Enumerate(MoleculeAndAtoms& m,
                   IWString_and_File_Descriptor& output);
    uint32_t MakeProduct(MoleculeAndAtoms& m,
                     const Set_of_Atoms& scaffold_state,
                     const std::vector<int>& state,
                     IWString_and_File_Descriptor& output);
    uint32_t EnumerateOneScaffoldSite(MoleculeAndAtoms& m,
                IWString_and_File_Descriptor& output);
    uint32_t MakeProduct(const MoleculeAndAtoms& scaffold,
                     const MoleculeAndAtom& sidechain,
                     IWString_and_File_Descriptor& output);
    bool ProductIsUnique(Molecule& product);
    uint32_t MaybeWrite(Molecule& product,
                    IWString_and_File_Descriptor& output);

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

    int ReadMolecules(const char* fname, FileType input_type);

    uint32_t Enumerate(IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _max_scaffold_attachments = 1;
  _max_sidechains = 0;
  _remove_isotopes_from_products = 0;
  _add_hydrogen_substituent = 0;
  _molecules_read = 0;
  _ignore_molecules_not_matching_queries = 0;
  _molecules_not_matching_queries = 0;
  _molecules_written = 0;
  _discard_bad_valence = 0;
  _undesirable_adjacencies = 0;
  _products_with_bad_valence = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
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

  if (cl.option_present('s')) {
    IWString smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _sidechain_query << qry.release();
    }
  } else if (cl.option_present('q')) {
    if (! process_queries(cl, _sidechain_query, _verbose, 'q')) {
      cerr << "Cannot process sidechain queries (-q)\n";
      return 0;
    }
  } else {
    std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
    if (! q->create_from_smarts("[aD3]-!@{a<10}[R0]")) {
      cerr << "Options::Initialise:cannot parse default smarts!\n";
      return 0;
    }
    _sidechain_query << q.release();
  }

  for (Substructure_Query* q : _sidechain_query) {
    q->set_embeddings_do_not_overlap(1);
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_scaffold_attachments) || _max_scaffold_attachments < 1) {
      cerr << "The maximum number of scaffold attachments (-x) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will allow as many as " << _max_scaffold_attachments << " scaffold replacements\n";
    }
  }

  if (cl.option_present('X')) {
    if (! cl.value('X', _max_sidechains) || _max_sidechains < 1) {
      cerr << "The maximum number of sidechains (-X) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will retain only the " << _max_sidechains <<
              " smallest sidechains\n";
    }
  }

  if (cl.option_present('z')) {
    _ignore_molecules_not_matching_queries = 1;
    if (_verbose) {
      cerr << "Will skip molecules not matching any of the breakage queries\n";
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes_from_products = 1;
    if (_verbose) {
      cerr << "Isotopes will be removed from product molecules\n";
    }
  }

  if (cl.option_present('V')) {
    _discard_bad_valence = 1;
    if (_verbose) {
      cerr << "Will discard products with bad valences\n";
    }
  }

  if (cl.option_present('h')) {
    _add_hydrogen_substituent = 1;
    if (_verbose) {
      cerr << "Will add Hydrogen as a substituent\n";
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Have " << _scaffolds.size() << " scaffolds and " << _sidechains.size() << " sidechain groups\n";
  // Other information about what has happened.
  if (_ignore_molecules_not_matching_queries) {
    output << "Ignored " <<  _molecules_not_matching_queries << " molecules not matching queries\n";
  }
  if (_discard_bad_valence) {
    output << _products_with_bad_valence << " products with bad valence discarded\n";
  }
  output << _undesirable_adjacencies << " product molecules with undesirable adjacencies\n";
  output << "Wrote " << _molecules_written << " molecules\n";

  return 1;
}

int
Options::ReadMolecules(const char* fname, FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Options::ReadMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 1) {
    input.set_verbose(1);
  }

  return ReadMolecules(input);
}

int
Options::ReadMolecules(data_source_and_type<Molecule>& input) {

  absl::flat_hash_set<IWString> seen;

  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! Preprocess(*m)) {
      continue;
    }

    _product_seen.insert(m->unique_smiles());

    if (SplitScaffoldSidechain(*m, seen)) {
      // great
    } else if (_ignore_molecules_not_matching_queries) {
      ++_molecules_not_matching_queries;
    } else {
      cerr << "Options::ReadMolecules:no query matches to\n";
      cerr << m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }

    ++_molecules_read;
  }

  if (_scaffolds.empty()) {
    cerr << "Options::ReadMolecules:no scaffolds\n";
    return 0;
  }
  if (_sidechains.empty()) {
    cerr << "Options:;ReadMolecules:no sidechains\n";
    return 0;
  }

  for (MoleculeAndAtoms* s : _scaffolds) {
    s->SortAtoms();
  }

  if (_add_hydrogen_substituent) {
    AddHydrogen();
  }

  // sort the sidechains by size.
  std::sort(_sidechains.begin(), _sidechains.end(),
    [](const MoleculeAndAtom* m1, const MoleculeAndAtom* m2) {
      return m1->molecule().natoms() < m2->molecule().natoms();
    });

  if (_max_sidechains > 0 && _sidechains.size() > _max_sidechains) {
    if (_verbose) {
      cerr << "Have " << _sidechains.size() << " sidechains, max is " << _max_sidechains << '\n';
    }
    _sidechains.resize(_max_sidechains);
  }

  if (_verbose) {
    cerr << "Read " << _scaffolds.size() << " scaffolds and " << _sidechains.size()
         << " sidechains\n";
//  for (MoleculeAndAtom* s : _sidechains) {
//    cerr << s->molecule().smiles() << ' ' << s->molecule().name() << '\n';
//  }

//  cerr << "Scaffolds\n";
//  for (MoleculeAndAtoms* s : _scaffolds) {
//    cerr << s->molecule().smiles() << ' ' << s->molecule().name() << '\n';
//    cerr << "Atoms " << s->atoms() << '\n';
//  }
  }

  return 1;
}

int
Options::AddHydrogen() {
  Molecule* hydrogen = new Molecule();
  hydrogen->add(get_element_from_atomic_number(1));
  hydrogen->set_name("H");

  MoleculeAndAtom* h = new MoleculeAndAtom(hydrogen);
  h->set_atom(0);
  // The sidechains array will be sorted and this will end up
  // at the beginning.
  _sidechains.insert_at_beginning(h);

  return 1;
}

int
Options::SplitScaffoldSidechain(Molecule& m, absl::flat_hash_set<IWString>& seen) {
  Molecule_to_Match target(&m);

  int queries_matching = 0;

  for (Substructure_Query* q : _sidechain_query) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }

    if (SplitScaffoldSidechain(m, sresults, seen)) {
      ++queries_matching;
    }
  }

  if (queries_matching == 0) {
    cerr << "No match to queries " << m.name() << '\n';
  }

  return queries_matching;
}

int
IsUnique(Molecule& m,
         absl::flat_hash_set<IWString>& seen) {
  auto iter = seen.find(m.unique_smiles());
  if (iter != seen.end()) {
    return 0;
  }

  IWString s = m.unique_smiles();
  seen.emplace(std::move(s));
  return 1;
}

int
Options::SplitScaffoldSidechain(Molecule& m, 
                const Substructure_Results& sresults,
                absl::flat_hash_set<IWString>& seen) {
  Set_of_Atoms scaffold_atom, sidechain_atom;

  for (const Set_of_Atoms* e : sresults.embeddings()) {
    // cerr << "Processing embedding " << *e << '\n';
    if (e->size() < 2) {
      cerr << "Options::SplitScaffoldSidechain:not enough matched atoms " << *e << "'\n";
      return 0;
    }

    atom_number_t a1 = e->item(0);
    atom_number_t a2 = e->item(1);
    if (! m.are_bonded(a1, a2)) {
      cerr << "Options::SplitScaffoldSidechain:matched atoms not bonded\n";
      return 0;
    }

    scaffold_atom << a1;
    sidechain_atom << a2;

    Molecule mcopy(m);
    mcopy.remove_bond_between_atoms(a1, a2);
    mcopy.set_isotope(a1, kScaffoldAtom);
    mcopy.set_isotope(a2, kSidechainAtom);

    resizable_array_p<Molecule> frags;
    mcopy.create_components(frags);
    assert(frags.size() == 2);

    if (!IsUnique(*frags[1], seen)) {
      frags.remove_item(1);
    }
    if (!IsUnique(*frags[0], seen)) {
      frags.remove_item(0);
    }

    if (frags.empty()) {
      continue;
    }

    for (Molecule* f : frags) {
      f->set_name(m.name());
    }

    // We have one or two fragments.
    // We may have one with isotope 1 and another with isotope 2.
    // Or we will have a single fragment with either isotope 1 or isotope 2.
    // This should be more straightforward...

    // First check to see if isotope 1 is on the first fragment.
    if (atom_number_t a = frags[0]->atom_with_isotope(kScaffoldAtom); a != kInvalidAtomNumber) {
      MoleculeAndAtoms* maa = new MoleculeAndAtoms(frags[0]);
      maa->add_atom(a);
      _scaffolds << maa;
      if (frags.size() == 2) {
        MoleculeAndAtom* maa = new MoleculeAndAtom(frags[1]);
        atom_number_t iso = frags[1]->atom_with_isotope(kSidechainAtom);
        maa->set_atom(iso);
        _sidechains << maa;
      }
    } else if (atom_number_t a = frags[0]->atom_with_isotope(kSidechainAtom); a != kInvalidAtomNumber) {
      MoleculeAndAtom* maa = new MoleculeAndAtom(frags[0]);
      maa->set_atom(a);
      _sidechains << maa;

      if (frags.size() == 2) {
        atom_number_t iso = frags[1]->atom_with_isotope(kScaffoldAtom);
        MoleculeAndAtoms* maa = new MoleculeAndAtoms(frags[1]);
        maa->add_atom(iso);
        _scaffolds << maa;
      }
    }

    // Very important, we have transferred ownership of the pointers.
    frags.resize_no_delete(0);
  }

  if (scaffold_atom.empty()) {
    return 1;
  }

  // No combinatorics.
  if (scaffold_atom.size() == 1) {
    return 1;
  }

  // Construct vectors consisting of some number of leading zero's and
  // then 1's, with the number of 1's being the number of sites that will 
  // be substituted.
  const int nsites = scaffold_atom.number_elements();
  // cerr << "Scaffold contains " << nsites << " sites\n";

  int istop;
  if (nsites < _max_scaffold_attachments) {
    istop = nsites;
  } else {
    istop = _max_scaffold_attachments;
  }

  for (int i = 2; i <= istop; ++i) {
    std::vector<int> state(nsites, 0);
    for (int j = 0; j < i; ++j) {
      state[nsites - j - 1] = 1;
    }
    do {
      CreateScaffold(m, scaffold_atom, sidechain_atom, state, seen);
    } while (std::next_permutation(state.begin(), state.end()));
  }

  return 1;
}

int
Options::CreateScaffold(const Molecule& m,
                        const Set_of_Atoms& scaffold_atom,
                        const Set_of_Atoms& sidechain_atom,
                        const std::vector<int>& state,
                        absl::flat_hash_set<IWString>& seen) {
  Molecule* mcopy = new Molecule(m);

  const int n = scaffold_atom.number_elements();
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }

    atom_number_t a1 = scaffold_atom[i];
    atom_number_t a2 = sidechain_atom[i];
    if (! mcopy->are_bonded(a1, a2)) {
      cerr << "Atoms not bonded " << a1 << ' ' << a2 << '\n';
      write_isotopically_labelled_smiles(*mcopy, false, cerr);
      cerr << ' ' << m.name() << '\n';
      continue;
    }
    mcopy->remove_bond_between_atoms(a1, a2);
    mcopy->set_isotope(a1, kScaffoldAtom);
    mcopy->set_isotope(a2, kSidechainAtom);
  }

  resizable_array<int> fragments_to_delete;
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }

    fragments_to_delete << mcopy->fragment_membership(sidechain_atom[i]);
  }

  mcopy->delete_fragments(fragments_to_delete);

  if (! IsUnique(*mcopy, seen)) {
    return 1;
  }

  // We continue to use the mcopy pointer even after we have 'given' it
  // to `mcc`. That is OK.
  std::unique_ptr<MoleculeAndAtoms> mcc = std::make_unique<MoleculeAndAtoms>(mcopy);
  const int matoms = mcopy->natoms();
  for (int i = 0; i < matoms; ++i) {
    if (mcopy->isotope(i) == kScaffoldAtom) {
      mcc->add_atom(i);
    }
  }

  // cerr << "Added scaffold " << mcc->molecule().smiles() << '\n';

  _scaffolds << mcc.release();

  return 1;
}

uint32_t
Options::Enumerate(IWString_and_File_Descriptor& output) {
  uint32_t rc = 0;
  for (MoleculeAndAtoms* scaffold : _scaffolds) {
    rc += Enumerate(*scaffold, output);
  }

  return rc;
}

uint32_t
Options::Enumerate(MoleculeAndAtoms& m,
                   IWString_and_File_Descriptor& output) {
  const int number_scaffold_sites = m.number_reactive_sites();

  if (number_scaffold_sites == 1) {
    return EnumerateOneScaffoldSite(m, output);
  }

  const int number_sidechains = _sidechains.number_elements();

  std::vector<int> sidechain_state(number_sidechains, 0);
  for (int i = 0; i < number_scaffold_sites; ++i) {
    sidechain_state[number_sidechains - i - 1] = 1;
  }

  // Note local copy.
  // This is cycled through multiple times.
//  Set_of_Atoms scaffold_state = m.atoms();

  uint32_t rc = 0;
  do {
    Set_of_Atoms scaffold_state = m.atoms();
    do {
      rc += MakeProduct(m, scaffold_state, sidechain_state, output);
    } while (std::next_permutation(scaffold_state.begin(), scaffold_state.end()));
  } while (std::next_permutation(sidechain_state.begin(), sidechain_state.end()));

  return rc;
}

uint32_t
Options::MakeProduct(MoleculeAndAtoms& m,
                     const Set_of_Atoms& scaffold_state,
                     const std::vector<int>& sidechain_state,
                     IWString_and_File_Descriptor& output) {
  const int nsites = scaffold_state.number_elements();

  const int nsidechains = _sidechains.number_elements();

  // Index into _sidechains array. We need to find those
  // members of that array for which sidechain_state[i] is nonzero.
  int j = 0;

  Molecule product(m.molecule());

  // cerr << "Generating combinatoric product\n";

  for (int i = 0; i < nsites; ++i) {
    for ( ; j < nsidechains; ++j) {
      if (sidechain_state[j] == 0) {
        continue;
      }

      if (product.hcount(scaffold_state[i]) == 0) {
        continue;
      }

      const int initial_natoms = product.natoms();
      product.add_molecule(&_sidechains[j]->molecule());
      atom_number_t in_sidechain = initial_natoms + _sidechains[j]->atom();
      product.add_bond(scaffold_state[i], in_sidechain, SINGLE_BOND);
      product.unset_all_implicit_hydrogen_information(scaffold_state[i]);
      product.unset_all_implicit_hydrogen_information(in_sidechain);
      // cerr << "To scaffold site " << scaffold_state[i] << " added sidechain " << j << '\n';
      ++j;
      break;
    }
  }

  return MaybeWrite(product, output);
}

// Return true if there are bonded, isotopically labelled heteroatoms.
// Also reject halogens to aliphatic atoms.
bool
UndesirableAtomAdjacencies(Molecule& m) {
  for (const Bond* b : m.bond_list()) {
    const Atom& a1 = m[b->a1()];
    const Atom& a2 = m[b->a2()];

    if (a1.isotope() == 0 || a2.isotope() == 0) {
      continue;
    }

    atomic_number_t z1 = a1.atomic_number();
    atomic_number_t z2 = a2.atomic_number();

    if (z1 == 6 && z2 == 6) {
      continue;
    }

    // Both heteroatoms, definitely bad.
    if (z1 != 6 && z2 != 6) {
      return true;
    }

    // At least one is carbon. Reject if this is a halogen to aliphatic carbon
    atom_number_t at1 = b->a1();
    atom_number_t at2 = b->a2();

    // Make sure z1 is the carbon.
    if (z2 == 6) {
      std::swap(z1, z2);
      std::swap(at1, at2);
    }

    // Some common cases.
    if (z2 == 7 || z2 == 8 || z2 == 16) {
      continue;
    }

    if (z2 == 9 || z2 == 17 || z2 == 35 || z2 == 53) {
      // If the carbon is not aromatic, reject.
      if (m.ring_bond_count(at1) == 0 || m.hcount(at1) > 1 || ! m.is_aromatic(at1)) {
        return true;
      }
    }
  }

  // Nothing bad found.
  return false;
}

uint32_t
Options::MaybeWrite(Molecule& product,
                    IWString_and_File_Descriptor& output) {
  if (_report_progress()) {
    cerr << "Written " << _molecules_written << " products generated\n";
  }

  // This check depends on isotopes being present.
  if (UndesirableAtomAdjacencies(product)) {
    if (_verbose > 1) {
      cerr << product.smiles() << " UndesirableAtomAdjacencies\n";
    }
    ++_undesirable_adjacencies;
    return 0;
  }

  if (_remove_isotopes_from_products) {
    product.unset_isotopes();
  }

  if (_add_hydrogen_substituent) {
    product.remove_all(lillymol::kHydrogen);
  }

  if (_discard_bad_valence && ! product.valence_ok()) {
    cerr << product.smiles() << ' ' << product.name() << " bad valence\n";
    ++_products_with_bad_valence;
    return 0;
  }

  if (! ProductIsUnique(product)) {
    return 0;
  }

  ++_molecules_written;

  output << product.smiles() << ' ' << product.name() << '.' << _molecules_written << '\n';

  output.write_if_buffer_holds_more_than(4196);

  return 1;
}

// The common case of a scaffold with just one attachment site is handled by just
// adding each of the sidechains.
uint32_t
Options::EnumerateOneScaffoldSite(MoleculeAndAtoms& m,
                IWString_and_File_Descriptor& output) {
  uint32_t rc = 0;
  for (const MoleculeAndAtom* sidechain : _sidechains) {
    rc += MakeProduct(m, *sidechain, output);
  }
  return rc;
}

uint32_t
Options::MakeProduct(const MoleculeAndAtoms& scaffold,
                     const MoleculeAndAtom& sidechain,
                     IWString_and_File_Descriptor& output) {
  assert(scaffold.atoms().size() == 1);

  Molecule product(scaffold.molecule());
  const int initial_natoms = product.natoms();
  product.add_molecule(&sidechain.molecule());
  atom_number_t a1 = scaffold.atoms()[0];
  atom_number_t a2 = initial_natoms + sidechain.atom();
  product.add_bond(a1, a2, SINGLE_BOND);
  product.unset_all_implicit_hydrogen_information(a1);
  product.unset_all_implicit_hydrogen_information(a2);

  if (! product.valence_ok()) {
    cerr << "Bad valence generated from \n";
    Molecule q(scaffold.molecule());
    cerr << q.smiles() << '\n';
    Molecule j(sidechain.molecule());
    cerr << j.smiles() << '\n';
    cerr << product.smiles() << " product\n";
  }

  return MaybeWrite(product, output);
}

// If we are not removing isotopes, we need to temporarily remove them before
// checking uniqueness.
bool
Options::ProductIsUnique(Molecule& product) {
  // Check current form - regardless of isotopes.
  if (! IsUnique(product, _product_seen)) {
    return 0;
  }

  // If isotopes removed the previous check was definitive. It was
  // not seen, so is unique.
  if (_remove_isotopes_from_products) {
    return 1;
  }

  // Temporarily remove isotopes.
  std::unique_ptr<isotope_t[]> iso = product.GetIsotopes();
  product.unset_isotopes();
  // No need to reset isotopes if a duplicate, product will be discarded.
  if (! IsUnique(product, _product_seen)) {
    return 0;
  }

  product.set_isotopes(iso.get());

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

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
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
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:q:s:x:X:z:IVhr:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (! options.ReadMolecules(cl[0], input_type)) {
    cerr << "Cannot read '" << cl[0] << "'\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  options.Enumerate(output);

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace sidechain_switcheroo

int
main(int argc, char ** argv) {

  int rc = sidechain_switcheroo::Main(argc, argv);

  return rc;
}
