// Identify substituents that are defined by matched query atoms.

#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace get_substituents {

using std::cerr;

using iw_tf_data_record::TFDataWriter;

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
  const char* msg = R"(Extract substituents based on a query and the adjacent unmatched atoms
 -s <smarts>       identify the matched atoms
 -q <query>        identify the matched atoms
 -O <atom>         from which query atom does the substituent sprout (default 0)
 -O all            process each query atom match
 -m <natoms>       min fragment size
 -M <natoms>       max fragment size
 -r <nrings>       min number of rings in the fragment
 -R <nrings>       max number of rings in the fragment
 -y <nsys>         min number of ring systems in the fragment
 -Y <nsys>         max number of ring systems in the fragment
 -I <isotope>      place <isotope> at the join atoms
 -P ...            atom typing specification - isotopes are for adjacent atom
 -f <smarts>       fragment must contain smarts (matches as `or`) 
 -F <query>        fragment must contain query specification *matches as `or`)
 -a                each of the fragment must contain queries must match (by default any)
 -n                suppress normal output (write summary via -S)
 -d                allow substituents bonded via double (or triple) bonds
 -S <fname>        write summary data Dicer.DicerFragment proto, to <fname>
 -z i              ignore molecules not matching any query
 -X ...            obscure options, enter '-X help' for info
 -v                verbose output

Specifying the query and the matched atom can be tricky. A common use case is to look
for substituents on a ring. A typical invocation might look like
-s '[cx2D3]-!@*' -O 1 -X anchor
This says that the substituent starts at matched atom 1 (-O 1), which is the
exocyclic atom * and that when the substituent is formed, the matched atom
itself is part of the substituent.
)";
  cerr << msg;
// clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    resizable_array_p<Substructure_Query> _queries;

    // With the matched atoms, which one defines the attachment point.

    int _matched_atom;

    // Alternatively we can process all attachment points.
    int _process_all_attachment_points;

    int _ignore_molecules_not_matching_any_queries;
    int _molecules_not_matching_queries;

    // We can impose limits on what gets saved.
    int _min_fragment_size;
    int _max_fragment_size;

    // We can impose limits on the number of rings in the substituents.
    int _min_rings;
    int _max_rings;

    // Limits on the ring systems in the substituent.
    int _min_ring_systems;
    int _max_ring_systems;

    // Keep track of some statistics about the fragments generated.
    extending_resizable_array<int> _atoms_in_fragment;
    extending_resizable_array<int> _rings_in_fragment;

    int _molecules_read = 0;

    int _write_fragments;
    int _accumulate_fragment_data;

    IW_STL_Hash_Map<IWString, dicer_data::DicerFragment> _fragments;

    // String appended to fragment names.
    IWString _suffix;

    // We can place an isotope to make the join points.
    int _isotope;

    int _fragments_generated;
    int _fragments_written;

    // We can put constraings on substructures in the fragment.
    resizable_array_p<Substructure_Query> _fragment_must_contain;

    // By default, if any of the _fragment_must_contain queries
    // match a fragment, that is good. By we can change that to
    // require that all of them match.
    int _all_fragment_must_have_queries_must_match;

    int _fragments_not_matching_substructure;

    // Some arrays that are used with the current molecule.
    // This means we are not thread safe.
    std::unique_ptr<int[]> _storage;
    // When we call create_subset, we need the atom cross refernce.
    std::unique_ptr<int[]> _xref;

    // Optionally write the parent molecule. Mostly helpful for
    // debugging.
    int _write_parent;

    // We can restrict ourselves to singly bonded attachments only.
    // Unfortunately the nature of the attachment point is not being
    // captured. TODO:ianwatson fix this.
    int _singly_bonded_attachments_only;

    // We can optionally record the atom types of the atom
    // to which our fragments are attached. These will be
    // placed as isotopes on the join atom.
    Atom_Typing_Specification _atom_typing;

    // It can be useful to include the anchor atom with the fragment.
    // Note there is a potential problem with this. Imagine the case
    // where the anchor atom has multiple unmatched atoms attached
    // to it. Right now, the program enumerates each unmatched
    // atom separately. But it would probably be more correct
    // to enumerate all those unmatched atoms as a single fragment.
    // TODO:ianwatson fix this.
    int _include_achor_atom_in_fragment;

    // By default we write textproto forms.
    int _write_serialized_protos;

  // private functions

    int QueryMatchToFragment(Molecule& m,
                        const Set_of_Atoms& embedding,
                        atom_number_t start_atom,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output);
    int QueryMatchToFragment(Molecule& m,
                        const Substructure_Results& sresults,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output);
    int QueryMatchToFragment(Molecule& m,
                        const Set_of_Atoms& atoms,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output);

    int MaybeAddToHash(const Molecule& parent, Molecule& frag);
    int MaybePlaceIsotope(const Molecule& parent,
                        const std::unique_ptr<uint32_t[]>& atype, 
                        atom_number_t zatom,
                        atom_number_t start_atom,
                        const int* xref, Molecule& frag) const;
    int UpdateStatistics(Molecule& frag);

    int PassesConstraints(Molecule& m);
    int AllMustHaveQueriesMatch(Molecule& m);
    int AnyMustHaveQueryMatches(Molecule& m);

    Set_of_Atoms GetUnmatchedAtoms(const Molecule& m,
                  atom_number_t zatom,
                  const int* matched) const;

    int IdentifyFragmentAtoms(const Molecule& m,
                      atom_number_t zatom,
                      int* storage,
                      int& atoms_in_fragment) const;
    int WriteSerialisedProtos(const char* fname) const;
    int WriteSerialisedProtos(TFDataWriter& output) const;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int WriteFragmentData(const char* fname) const;
    int WriteFragmentData(IWString_and_File_Descriptor& output) const;

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }

    void set_singly_bonded_attachments_only(int s) {
      _singly_bonded_attachments_only = s;
    }
    void set_include_achor_atom_in_fragment(int s) {
      _include_achor_atom_in_fragment = s;
    }
    void set_write_serialized_protos(int s) {
      _write_serialized_protos = s;
    }
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _write_fragments = 1;
  _accumulate_fragment_data = 0;

  _isotope = 0;

  _min_fragment_size = 0;
  _max_fragment_size = std::numeric_limits<int>::max();

  _min_rings = 0;
  _max_rings = std::numeric_limits<int>::max();

  _min_ring_systems = 0;
  _max_ring_systems = 0;

  _fragments_generated = 0;
  _fragments_written = 0;

  _include_achor_atom_in_fragment = 0;

  _matched_atom = 0;
  _process_all_attachment_points = 0;
  _ignore_molecules_not_matching_any_queries = 0;
  _molecules_not_matching_queries = 0;

  _fragments_not_matching_substructure = 0;

  _all_fragment_must_have_queries_must_match = 0;

  _suffix = ".frag";

  _write_parent = 0;

  _singly_bonded_attachments_only = 1;

  _write_serialized_protos = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
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

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Options::Initialise:cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << qry.release();
    }
  }

  if (_queries.empty()) {
    cerr << "Options::Initialise:no query atoms specified\n";
    return 0;
  }

  for (Substructure_Query* q : _queries) {
    q->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('F')) {
    if (! process_queries(cl, _fragment_must_contain, _verbose, 'F')) {
      cerr << "Options::Initialise:cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('f')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('f', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _fragment_must_contain << qry.release();
    }
  }

  if (cl.option_present('a')) {
    _all_fragment_must_have_queries_must_match = 1;
    if (_verbose) {
      cerr << "All fragment must have queries (-f, -F) must match\n";
    }
  }

  if (cl.option_present('n')) {
    _write_fragments = 0;
    if (_verbose) {
      cerr << "Normal output suppressed, fragments accumulated\n";
    }
    _accumulate_fragment_data = 1;
  }

  // Even though the caller handles the -S file, we take a peek.
  if (cl.option_present('S')) {
    _accumulate_fragment_data = 1;
  }

  if (cl.option_present('z')) {
    const IWString z = cl.string_value('z');
    if (z == 'i') {
      _ignore_molecules_not_matching_any_queries = 1;
      if (_verbose) {
        cerr << "Will ignore molecules not matching any query\n";
      }
    } else {
      cerr << "Unrecognised -z qualifier '" << z << "'\n";
      return 0;
    }
  }

  if (_accumulate_fragment_data && ! cl.option_present('S')) {
    cerr << "Must specify file name for fragment data via the -S option\n";
    return 0;
  }

  if (cl.option_present('I') && cl.option_present('P')) {
    cerr << "The -I and -P options are not compatible\n";
    return 0;
  }

  if (cl.option_present('I')) {
    if (! cl.value('I', _isotope) || _isotope < 0) {
      cerr << "Options::Initialise:the isotope flag (-I) must be a whole +ve number\n";
      return 0;
    }
  }

  if (cl.option_present('O')) {
    const_IWSubstring o = cl.string_value('O');
    if (o == "all") {
      _process_all_attachment_points = 1;
      if (_verbose) {
        cerr << "WIll process all attachment points\n";
      }
    } else if (! cl.value('O', _matched_atom) || _matched_atom < 0) {
      cerr << "Options::Initialise:the matched atom (-O) option must be a whole +ve number\n";
      return 0;
    } else if (_verbose) {
      cerr << "Matched atom " << _matched_atom << " will be the origin\n";
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_fragment_size) || _min_fragment_size < 0) {
      cerr << "Options::Initialise:the min fragment size option (-m) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Min fragment size " << _min_fragment_size << '\n';
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_fragment_size) || _max_fragment_size < _min_fragment_size) {
      cerr << "Options::Initialise:the max fragment size option (-M) must be a whole +ve number > "
           << _min_fragment_size << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Max fragment size " << _max_fragment_size << '\n';
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _min_rings) || _min_rings < 0) {
      cerr << "Options::Initialise:the min number of rings option (-r) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Min number of rings " << _min_rings << '\n';
    }
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _max_rings) || _max_rings < _min_rings) {
      cerr << "Options::Initialise:the max number of rings option (-R) must be a whole +ve number > "
           << _min_rings << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Max number of rings " << _max_rings << '\n';
    }
  }

  if (cl.option_present('y')) {
    if (! cl.value('y', _min_ring_systems) || _min_ring_systems < 0) {
      cerr << "Options::Initialise:the min number of ring systems option (-y) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "min number of ring sytems in a substituent " << _min_ring_systems << '\n';
    }
  }

  if (cl.option_present('Y')) {
    if (! cl.value('Y', _max_ring_systems) || _max_ring_systems < 0) {
      cerr << "Options::Initialise:the max number of ring systems option (-Y) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Max number of ring sytems in a substituent " << _max_ring_systems << '\n';
    }
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n:";
    }
  }

  if (cl.option_present('d')) {
    _singly_bonded_attachments_only = 0;
    if (_verbose) {
      cerr << "Will detect doubly bonded substituents\n";
    }
  }

  if (cl.option_present('P')) {

    const_IWSubstring p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Invalid atom type specification '" << p << "'\n";
      return 0;
    }
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _molecules_not_matching_queries << " molecules did not match any of "
         << _queries.size() << " queries " <<
         iwmisc::Fraction<float>(_molecules_not_matching_queries, _molecules_read)
         << '\n';
  output << _fragments_generated << " fragments generated\n";
  if (! _fragment_must_contain.empty()) {
    output << _fragments_not_matching_substructure << " fragments did not match substructure constraints\n";
  }
  if (_write_fragments) {
    output << _fragments_written << " fragments written\n";
  }
  output << _fragments.size() << " fragments in hash\n";

  for (uint32_t i = 0; i < _atoms_in_fragment.size(); ++i) {
    if (_atoms_in_fragment[i]) {
      output << _atoms_in_fragment[i] << " fragments had " << i << " atoms\n";
    }
  }

  for (uint32_t i = 0; i < _rings_in_fragment.size(); ++i) {
    if (_rings_in_fragment[i]) {
      output << _rings_in_fragment[i] << " fragments had " << i << " rings\n";
    }
  }

  return 1;
}

int
Options::MaybeDiscernInputType(const char * fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }
  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  std::unique_ptr<uint32_t[]> atype;
  if (_atom_typing.active()) {
    atype = std::make_unique<uint32_t[]>(m.natoms());
    _atom_typing.assign_atom_types(m, atype.get());
  }

  m.compute_aromaticity_if_needed();

  Molecule_to_Match target(&m);

  _storage.reset(new_int(m.natoms()));
  _xref.reset(new int[m.natoms()]);

  int queries_matched = 0;

  // Parent is written even if no fragments generated.
  if (_write_parent) {
    output << m.smiles() << ' ' << m.name() << '\n';
  }

  Substructure_Results sresults;
  for (Substructure_Query* q : _queries) {
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    ++queries_matched;
    if (QueryMatchToFragment(m, sresults, atype, output)) {
      return 1;
    }
  }

  if (queries_matched == 0) {
    ++_molecules_not_matching_queries;
    if (_ignore_molecules_not_matching_any_queries) {
      return 1;
    }
    cerr << "Options::Process:None of " << _queries.size() 
         << " queries matched " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  return 1;
}

int
Options::QueryMatchToFragment(Molecule& m,
                        const Substructure_Results& sresults,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output) {
  for (const Set_of_Atoms* atoms : sresults.embeddings()) {
    if (QueryMatchToFragment(m, *atoms, atype, output)) {
      // Not sure why this was returning, processing all embeddings is what we want.
      // return 1;
    }
  }

  return 0;
}

// Place a value of 2 in `storage` for every atom that
// is traversed.
// Increment `atoms_in_fragment` for every atom
// encountered.
int
Options::IdentifyFragmentAtoms(const Molecule& m,
                      atom_number_t zatom,
                      int* storage,
                      int& atoms_in_fragment) const {
  storage[zatom] = 2;
  ++atoms_in_fragment;
  if (atoms_in_fragment > _max_fragment_size) {
    return 0;
  }

  const Atom* a = m.atomi(zatom);
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (storage[j]) {
      continue;
    }
    if (! IdentifyFragmentAtoms(m, j, storage, atoms_in_fragment)) {
      return 0;
    }
  }

  return 1;
}

// Return the list of connections to `zatom` which are NOT set
// in `matched`. Do not include ring bonds, since we will be
// growing substituents from these atoms.
// If we are only processing singly bonded attachments, skip
// any attachments that do not match.
Set_of_Atoms
Options::GetUnmatchedAtoms(const Molecule& m,
                  atom_number_t zatom,
                  const int* matched) const {
  Set_of_Atoms result;

  const Atom* a = m.atomi(zatom);
  for (const Bond* b : *a) {
    if (_singly_bonded_attachments_only && ! b->is_single_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    //cerr << "From matched atom " << zatom << " to atom " << j << " matched " << matched[j] << " ring " << b->nrings() << '\n';
    if (matched[j]) {
      continue;
    }
    result << j;
  }

  return result;
}

// Branch depending on which matched atom specification is being used.
int
Options::QueryMatchToFragment(Molecule& m,
                        const Set_of_Atoms& atoms,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output) {
  // cerr << "Examining embedding " << atoms << ", _process_all_attachment_points " << _process_all_attachment_points << '\n';
  if (! _process_all_attachment_points) {
    return QueryMatchToFragment(m, atoms, _matched_atom, atype, output);
  }

  // Process all matched atoms.
  int rc = 0;
  const int n = atoms.size();
  for (int i = 0; i < n; ++i) {
    if (atoms[i] < 0) {
      continue;
    }
    rc += QueryMatchToFragment(m, atoms, i, atype, output);
  }

  return rc;
}

// We have formed a subset of `m` in `frag` driven by `xref`.
// If `zatom` was [sD3] in `m`, set the implicit hydrogen count 
// for the corresponding atom in `frag`.
static int
FixAromaticSulphur(Molecule& m, 
                   atom_number_t zatom, 
                   const int* xref,
                   Molecule& frag) {
  const Atom& s = m[zatom];
  if (s.atomic_number() != 16) {
    return 0;
  }
  if (s.ncon() != 3) {
    return 0;
  }
  if (! m.is_aromatic(zatom)) {
    return 0;
  }

  frag.set_implicit_hydrogens(xref[zatom], 1, 1);
  return 1;
}

//#define DEBUG_MAKE_FRAGMENT

// Identify the substituent defined by atoms[atom].
int
Options::QueryMatchToFragment(Molecule& m,
                        const Set_of_Atoms& embedding,
                        atom_number_t start_atom,
                        const std::unique_ptr<uint32_t[]>& atype,
                        IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  std::fill_n(_storage.get(), matoms, 0);
  embedding.set_vector(_storage.get(), 1);

  if (! embedding.ok_index(start_atom)) {
    cerr << "Options::QueryMatchToFragment:invalid matched atom\n";
    cerr << "Embedding " << embedding << " ndx " << start_atom << " invald\n";
    return 0;
  }

#ifdef DEBUG_MAKE_FRAGMENT
  if (atype != nullptr) {
    for (int i = 0; i < matoms; ++i) {
      cerr << i << " " << m.smarts_equivalent_for_atom(i) << " atype " << atype[i] << '\n';
    }
  }
#endif

  // Force sssr if needed.
  m.ring_membership();

  const atom_number_t zatom = embedding[start_atom];
  if (zatom < 0) {
    cerr << "QueryMatchToFragment:excluded atom is start atom " << embedding << '\n';
    return 0;
  }

  // The atoms attached to `zatom` which are not matched by
  // the query.
  Set_of_Atoms unmatched_attached = GetUnmatchedAtoms(m, zatom, _storage.get());
#ifdef DEBUG_MAKE_FRAGMENT
  cerr << "Embedding " << embedding << " ndx " << start_atom << " atom " << embedding[start_atom] << '\n';
  cerr << "unmatched_attached " << unmatched_attached << '\n';
#endif

  if (_include_achor_atom_in_fragment) {
    _storage[zatom] = 0;
  }

  for (atom_number_t atom : unmatched_attached) {
    for (int i = 0; i < matoms; ++i) {
      if (_storage[i] == 2) {
        _storage[i] = 0;
      }
    }
    int atoms_in_fragment = 0;
    if (! IdentifyFragmentAtoms(m, atom, _storage.get(), atoms_in_fragment)) {
      _storage[zatom] = 2;
      continue;
    }

    ++_atoms_in_fragment[atoms_in_fragment];
#ifdef DEBUG_MAKE_FRAGMENT
    for (int i = 0; i < matoms; ++i) {
      cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " in_fragment " << _storage[i] << '\n';
    }
#endif
    Molecule frag;
    m.create_subset(frag, _storage.get(), 2, _xref.get());
    if (!PassesConstraints(frag)) {
      continue;
    }
    UpdateStatistics(frag);
#ifdef DEBUG_MAKE_FRAGMENT
    // cerr << "Fragment is " << frag.smiles() << '\n';
#endif
    MaybePlaceIsotope(m, atype, atom, zatom, _xref.get(), frag);
    FixAromaticSulphur(m, atom, _xref.get(), frag);
    MaybeAddToHash(m, frag);
    if (_write_fragments) {
      ++_fragments_written;
      output << frag.unique_smiles() << ' ' << m.name() <<  _suffix << '\n';
    }
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::MaybeAddToHash(const Molecule& parent, 
                Molecule& frag) {
  if (! _accumulate_fragment_data) {
    return 1;
  }

  auto iter = _fragments.find(frag.unique_smiles());
  if (iter != _fragments.end()) {
    auto current = iter->second.n();
    iter->second.set_n(current + 1);
    return 1;
  }

  auto [iter2, _] = _fragments.emplace(frag.unique_smiles(), dicer_data::DicerFragment());
  if (_atom_typing.active()) {
    iter2->second.set_iso(dicer_data::ATYPE);
  } else if (_isotope) {
    iter2->second.set_iso(dicer_data::ATT);
  } else {
    // No isotope.
  }
  iter2->second.set_n(1);
  iter2->second.set_par(parent.name().AsString());
  iter2->second.set_smi(frag.unique_smiles().AsString());
  iter2->second.set_nat(frag.natoms());

  return 1;
}

// `frag` is a subset of `parent` related by `xref`.
// If we are placing isotopes on fragments, place
// one in `frag` on the atom that came from `zatom`;
int
Options::MaybePlaceIsotope(const Molecule& parent,
                const std::unique_ptr<uint32_t[]>& atype,
                atom_number_t zatom,
                atom_number_t start_atom,
                const int* xref,
                Molecule& frag) const {
#ifdef DEBUG_MAKE_FRAGMENT
  cerr << "MaybePlaceIsotope: _isotope " << _isotope << " atom " << zatom << " xref " << xref[zatom] << '\n';
#endif
  if (_isotope == 0) {
  } else if (_include_achor_atom_in_fragment) {
    frag.set_isotope(xref[start_atom], _isotope);
    return 1;
  } else {
    frag.set_isotope(xref[zatom], _isotope);
    return 1;
  }

  // cerr << "MaybePlaceIsotope zatom " << zatom << " xref " << xref[zatom] << " start_atom " << start_atom << " atype " << atype[start_atom] << '\n';

  if (atype) {
    frag.set_isotope(xref[zatom], atype[start_atom]);
    return 1;
  }

  return 0;
}

int
Options::UpdateStatistics(Molecule& frag) {
  ++_atoms_in_fragment[frag.natoms()];
  ++_rings_in_fragment[frag.nrings()];

  return 1;
}

int
Options::WriteFragmentData(const char* fname) const {
  if (_write_serialized_protos) {
    return WriteSerialisedProtos(fname);
  }

  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteFragmentData:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteFragmentData(output);
}

int
Options::WriteSerialisedProtos(const char* fname) const {
  TFDataWriter output;
  if (! output.Open(fname)) {
    cerr << "Options::WriteSerialisedProtos:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteSerialisedProtos(output);
}
int
Options::WriteSerialisedProtos(TFDataWriter& output) const {
  std::string buffer;
  for (const auto& [_, proto] : _fragments) {
    output.WriteSerializedProto<dicer_data::DicerFragment>(proto);
  }

  return 1;
}

int
Options::WriteFragmentData(IWString_and_File_Descriptor& output) const {

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  static constexpr char kSep = ' ';

  for (const auto& [key, proto] : _fragments) {
    std::string buffer;
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteFragmentData:cannot write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }
    output << key << kSep;
    output << buffer;
    output << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// `m` is a proposed fragment. Is it OK wrt any substructure constraints
int
Options::PassesConstraints(Molecule& m) {
  ++_fragments_generated;

  if (m.natoms() < _min_fragment_size) {
    return 0;
  }

  if (m.natoms() > _max_fragment_size) {
    return 0;
  }

  if (_min_rings > 0 && m.nrings() < _min_rings) {
    return 0;
  }

  if (_max_rings == std::numeric_limits<int>::max()) {
  } else if (m.nrings() > _max_rings) {
    return 0;
  }

  if (_min_ring_systems == 0) {
  } else if (m.number_ring_systems() < _min_ring_systems) {
    return 0;
  }

  if (_max_ring_systems == 0) {
  } else if (m.number_ring_systems() > _max_ring_systems) {
    return 0;
  }

  if (_fragment_must_contain.empty()) {
    return 1;
  }

  if (_all_fragment_must_have_queries_must_match) {
    return AllMustHaveQueriesMatch(m);
  } else {
    return AnyMustHaveQueryMatches(m);
  }
}

int
Options::AllMustHaveQueriesMatch(Molecule& m) {
  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  for (Substructure_Query* q : _fragment_must_contain) {
    if (! q->substructure_search(target, sresults)) {
      ++_fragments_not_matching_substructure;
      return 0;
    } 
  }

  return 1;
}

int
Options::AnyMustHaveQueryMatches(Molecule& m) {
  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  for (Substructure_Query* q : _fragment_must_contain) {
    if (q->substructure_search(target, sresults)) {
      return 1;
    } 
  }

  ++_fragments_not_matching_substructure;
  return 0;
}

int
GetSubstituents(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
GetSubstituents(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! GetSubstituents(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetSubstituents(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "GetSubstituents:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return GetSubstituents(options, input, output);
}

void
DisplayDashXOptions(std::ostream& output) {
  output << " -X 3d              include coordinates with smiles (needs 3D input)\n";
  output << " -X okdouble        allow doubly bonded substituents\n";
  output << " -X anchor          include the anchor atom with the fragment\n";
  ::exit(0);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:S:I:s:q:f:F:m:M:nz:O:aX:pdP:r:R:y:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "3d") {
        set_append_coordinates_after_each_atom(1);
        if (verbose) {
          cerr << "Will write 3D smiles\n";
        }
      } else if (x == "okdouble") {
        options.set_singly_bonded_attachments_only(0);
        if (verbose) {
          cerr << "Will allow multiply bonded substituents\n";
        }
      } else if (x == "anchor") {
        options.set_include_achor_atom_in_fragment(1);
        if (verbose) {
          cerr << "Will include the anchor atom with the fragment\n";
        }
      } else if (x == "binproto") {
        options.set_write_serialized_protos(1);
        if (verbose) {
          cerr << "Will write serialized protos\n";
        }
      } else if (x == "help") {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GetSubstituents(options, fname, output)) {
      cerr << "GetSubstituents::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (verbose) {
      cerr << "Writing fragment data to '" << fname << "'\n";
    }
    if (! options.WriteFragmentData(fname.null_terminated_chars())) {
      cerr << "Count not write fragment data '" << fname << "'\n";
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace get_substituents

int
main(int argc, char ** argv) {

  int rc = get_substituents::Main(argc, argv);

  return rc;
}
