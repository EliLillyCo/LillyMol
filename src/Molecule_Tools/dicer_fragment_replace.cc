// Consumes the output of dicer_to_topology_types to make changes
// We may still generate no products, but let's write the parent here
// if requested.
// to input molecules.

#include <iostream>
#include <memory>
#include <vector>

#include "re2/re2.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/molecule_filter_lib.h"
#include "Molecule_Tools/molecule_filter.pb.h"

namespace dicer_fragment_replace {

using std::cerr;
using iw_tf_data_record::TFDataReader;
using dicer_data::DicerFragment; 
using combinations::Combinations;

// Warning, this will need to be updated when we make the switch globally.
constexpr int kSingleBond = SINGLE_BOND;

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
cerr << R"(
Uses fragments from dicer_to_topology_types to make changes to molecules.
Fragment files are specified via the -F option(s). From the name, we can figure out what kind of
frgaments are in that file, and we handle them as either
  Sidechain
  2 connected linker
  3 connected linker
 -F <fname>             file(s) from the -S option to dicer_to_topology_types
 -L <fname              molecule_filter textproto applied to replacement fragments.

 -m <smarts>            smarts for features the replacement fragment must have
 -M <qry>               query file for features the replacement fragment must have
 -x <smarts>            smarts for features the replacement fragment must NOT have
 -X <qry>               query file for features the replacement fragment must NOT have

 -s <smarts>            smarts for one or more bonds to break.
 -q <query>             query form for the same thing.
                        If 1 query, sidechain replacement.
                        If 2-3 queries, linker replacement.
                        If sidechains, we assume the second matched atom is lost.
                        If a linker, we assume the first query atom is in the region being removed.
 -P ...                 atom typing specification. If used must be the same atom typing that was
                        used during the dicer run used by dicer_to_topology_types.
 -Y ...                 various other options, enter '-Y help' for info
 -V <fname>             discard invalid valences, write to <fname>, Use '-V none' to just discard
 -z ...                 ignore various failures, enter '-z help' for info
 -l                     reduce to largest fragment
 -c                     remove chirality
 -g ...                 chemical standardisation
 -A ...                 standard aromaticity options
 -E ...                 standard element options
 -v          verbose output
)";

  ::exit(rc);
}

void
DisplayDashYOptions(int rc) {
  // clang-format off
  cerr << R"(Miscellaneous options
 -Y minextra=<n>          minimum number of extra atoms that must be in a product
 -Y maxextra=<n>          maximum number of extra atoms that must be in a product
 -Y minfewer=<n>          minimum fewer atoms that must be in a product
 -Y maxfewer=<n>          maximum fewer atoms that must be in a product
 -Y rmiso                 remove isotopes from product molecules
 -Y wrparent              write the parent molecule before each calculation
 
)";
  // clang-format on

  ::exit(rc);
}

void
DisplayDashZOptions() {
  // clang-format off
  cerr << R"(Controlling how various non match and failures are handled
 -z smiles              sometimes aromatic smiles in protos cannot be decoded. Ignore that.
 -z i                   ignore the molecule if it does not match the -s/-q bond breaking queries.
)";
  // clang-format on

  ::exit(1);
}


// As we read frgaments, we can filter them.
class ReplacementSpecifications {
  private:
    resizable_array_p<Substructure_Query> _fragment_must_contain;
    resizable_array_p<Substructure_Query> _fragment_must_not_contain;

    // Stats on the impact of query matches.
    int _fail_must_contain;
    int _fail_must_not_contain;

    // We can choose to only replace fragments with a certain support level.
    uint32_t _min_support;

    // We can impose requirements on the reagents.
    molecule_filter_lib::MoleculeFilter _filter;
    // Set to true if the filter has been activated.
    int _filter_active;
    uint32_t _rejected_by_filter;


    // Aromatised proto smiles might be hard to interpret.
    int _ignore_smiles_interpretation_errors;

    // Statistics on the impact of our requirements.
    int _protos_examined;
    int _fragments_below_support_requirement;

  public:
    ReplacementSpecifications();

    int Initialise(Command_Line& cl);

    int PassesConstraints(Molecule& m, const dicer_data::DicerFragment& proto);

    int ignore_smiles_interpretation_errors() const {
      return _ignore_smiles_interpretation_errors;
    }

    int Report(std::ostream& output) const;
};

ReplacementSpecifications::ReplacementSpecifications() {
  _min_support = 0;
  _ignore_smiles_interpretation_errors = 0;
  _protos_examined = 0;
  _filter_active = 0;
  _rejected_by_filter = 0;
  _fragments_below_support_requirement = 0;
  _fail_must_contain = 0;
  _fail_must_not_contain = 0;
}

int
ReplacementSpecifications::Initialise(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  if (cl.option_present('m')) {
    const_IWSubstring m;
    for (int i = 0; cl.value('m', m, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(m)) {
        cerr << "ReplacementSpecifications::invalid smarts '" << m << "'\n";
        return 0;
      }
      _fragment_must_contain << q.release();
    }
  }

  if (cl.option_present('x')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('x', x, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(x)) {
        cerr << "ReplacementSpecifications::invalid smarts '" << x << "'\n";
        return 0;
      }
      _fragment_must_not_contain << q.release();
    }
  }

  static constexpr int kVerbose = 0;

  if (cl.option_present('M')) {
    if (! process_queries(cl, _fragment_must_contain, kVerbose, 'M')) {
      cerr << "ReplacementSpecifications::cannot process fragment must have queries\n";
      return 0;
    }
  }

  if (cl.option_present('X')) {
    if (! process_queries(cl, _fragment_must_not_contain, kVerbose, 'X')) {
      cerr << "ReplacementSpecifications::cannot process fragment must not have queries\n";
      return 0;
    }
  }

  if (cl.option_present('L')) {
    IWString fname = cl.string_value('L');
    if (! _filter.Build(fname)) {
      cerr << "ReplacementSpecifications::cannot build filter '" << fname << "'\n";
      return 0;
    }
    _filter_active = 1;
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _min_support) || _min_support < 1) {
      cerr << "The min support level (-p) must be a valid +ve number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will only select fragments with at least " << _min_support << " examplars\n";
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == "smiles") {
        _ignore_smiles_interpretation_errors = 1;
        if (verbose) {
          cerr << "Will ignore smiles interpretation errors\n";
        }
      } else if (z == "help") {
        DisplayDashZOptions();
      }
    }
  }

  return 1;
}

int
ReplacementSpecifications::Report(std::ostream& output) const {
  output << "Examined " << _protos_examined << " protos\n";
  if (_min_support > 0) {
    output << _fragments_below_support_requirement << " fragments below support " << _min_support << '\n';
  }

  if (_fragment_must_contain.size()) {
    output << _fail_must_contain << " fragments failed contain queries\n";
  }

  if (_fragment_must_not_contain.size()) {
    output << _fail_must_not_contain << " fragments failed NOT contain queries\n";
  }

  if (_filter_active) {
    output << _rejected_by_filter << " rejected by molecule filter\n";
  }

  return 1;
}

int
AnyQueryMatches(resizable_array_p<Substructure_Query> & query,
                Molecule_to_Match& target) {
  for (Substructure_Query* q : query) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

// Both _fragment_must_contain and _fragment_must_not_contain are interpreted
// as OR matches.
// If any _fragment_must_not_contain query matches `m`, it is rejected.
// In order to be accepted at least one of _fragment_must_contain must match.

int
ReplacementSpecifications::PassesConstraints(Molecule& m, 
                        const dicer_data::DicerFragment& proto) {
  if (_min_support > 0 && proto.n() < _min_support) {
    ++_fragments_below_support_requirement;
    return 0;
  }

  if (_fragment_must_contain.size() > 0 ||
      _fragment_must_not_contain.size() > 0) {
    Molecule_to_Match target(&m);
    if (_fragment_must_not_contain.size() > 0 && 
        AnyQueryMatches(_fragment_must_not_contain, target)) {
      ++_fail_must_contain;
      return 0;
    }
    if (_fragment_must_contain.size() > 0 &&
        ! AnyQueryMatches(_fragment_must_contain, target)) {
      ++_fail_must_not_contain;
      return 0;
    }
  }

  if (_filter_active) {
    if (! _filter.Ok(m)) {
      ++_rejected_by_filter;
      return 0;
    }
  }

  return 1;
}

class MoleculeOneAtom : public Molecule {
  private:
    atom_number_t _attach;

  public:
    MoleculeOneAtom();

    int IdentifyAttachmentPoints();

    atom_number_t attach() const {
      return _attach;
    }
};

MoleculeOneAtom::MoleculeOneAtom() {
  _attach = kInvalidAtomNumber;
}

int
MoleculeOneAtom::IdentifyAttachmentPoints() {
  const int matoms = natoms();
  for (int i = 0; i < matoms; ++i) {
    if (this->isotope(i) > 0) {
      if (_attach == kInvalidAtomNumber) {
        _attach = i;
      } else {
        cerr << "MoleculeOneAtom::more than 1 isotope " << this->smiles() << ' ' << this->name() << '\n';
        return 0;
      }
    }
  }

  return 1;
}

class MoleculeTwoAtoms : public Molecule {
  private:
    atom_number_t _attach1;
    atom_number_t _attach2;

  public:
    MoleculeTwoAtoms();

    int IdentifyAttachmentPoints();

    atom_number_t attach1() const { 
      return _attach1;
    }
    atom_number_t attach2() const { 
      return _attach2;
    }
};

MoleculeTwoAtoms::MoleculeTwoAtoms() {
  _attach1 = kInvalidAtomNumber;
  _attach2 = kInvalidAtomNumber;
}

int
MoleculeTwoAtoms::IdentifyAttachmentPoints() {
  const int matoms = this->natoms();
  for (int i = 0; i < matoms; ++i) {
    if (this->isotope(i) > 0) {
      if (_attach1 == kInvalidAtomNumber) {
        _attach1 = i;
      } else if (_attach2 == kInvalidAtomNumber) {
        _attach2 = i;
      } else {
        cerr << "MoleculeOneAtom::more than two isotopes " << this->smiles() << ' ' << this->name() << '\n';
        return 0;
      }
    }
  }

  return 1;
}

class MoleculeThreeAtoms : public Molecule {
  private:
    atom_number_t _attach1;
    atom_number_t _attach2;
    atom_number_t _attach3;

  public:
    MoleculeThreeAtoms();

    int IdentifyAttachmentPoints();

    atom_number_t attach1() const { 
      return _attach1;
    }
    atom_number_t attach2() const { 
      return _attach2;
    }
    atom_number_t attach3() const { 
      return _attach3;
    }
};

MoleculeThreeAtoms::MoleculeThreeAtoms() {
  _attach1 = kInvalidAtomNumber;
  _attach2 = kInvalidAtomNumber;
  _attach3 = kInvalidAtomNumber;
}

int
MoleculeThreeAtoms::IdentifyAttachmentPoints() {
  const int matoms = this->natoms();
  for (int i = 0; i < matoms; ++i) {
    if (this->isotope(i) > 0) {
      if (_attach1 == kInvalidAtomNumber) {
        _attach1 = i;
      } else if (_attach2 == kInvalidAtomNumber) {
        _attach2 = i;
      } else if (_attach3 == kInvalidAtomNumber) {
        _attach3 = i;
      } else {
        cerr << "MoleculeOneAtom::more than three isotopes " << this->smiles() << ' ' << this->name() << '\n';
        return 0;
      }
    }
  }

  return 1;
}

// Instantiated vwith M being one of the Molecule* classes above.
template <typename M>
class BaseFragment {
  protected:
    resizable_array_p<dicer_data::DicerFragment> _fragment;
    // For each proto, the Molecule for the fragment.
    resizable_array_p<M> _mol;

    // The file name from which this was read.
    IWString _fname;

  private:
     int Read(TFDataReader& reader, ReplacementSpecifications& replacement_specifications);
     int DetermineFragmentType();

  protected:
    int Read(IWString & fname, ReplacementSpecifications& replacement_specifications);

  public:
    BaseFragment();
    ~BaseFragment();

    int NumberFragments() const {
      return _fragment.number_elements();
    }
    const IWString& FileName() const {
      return _fname;
    }
    const M** cbegin() const {
      return _mol.cbegin();
    }
    const M** cend() const {
      return _mol.cend();
    }
    M** begin() const {
      return _mol.begin();
    }
    M** end() const {
      return _mol.end();
    }
};


template <typename M>
BaseFragment<M>::BaseFragment() {
}

template <typename M>
BaseFragment<M>::~BaseFragment() {
}

template <typename M>
int
BaseFragment<M>::Read(IWString& fname, ReplacementSpecifications& replacement_specifications) {
  TFDataReader reader(fname.null_terminated_chars());
  if (! reader.good()) {
    cerr << "BaseFragment::Read:cannot open '" << fname << "'\n";
    return 0;
  }

  _fname = fname;

  return Read(reader, replacement_specifications);
}

template <typename M>
int
BaseFragment<M>::Read(TFDataReader& reader, ReplacementSpecifications& replacement_specifications) {
  int protos_read = 0;
  int failed_constraints = 0;
  while (1) {
    std::unique_ptr<DicerFragment> maybe_proto = reader.ReadProtoPtr<DicerFragment>();
    if (! maybe_proto) {
      break;
    }
    ++protos_read;

    std::unique_ptr<M> m = std::make_unique<M>();
    if (! m->build_from_smiles(maybe_proto->smi())) {
      return replacement_specifications.ignore_smiles_interpretation_errors();
    }

    if (! replacement_specifications.PassesConstraints(*m, *maybe_proto)) {
      ++failed_constraints;
      continue;
    }

    if (! m->IdentifyAttachmentPoints()) {
      return 0;
    }

    m->set_name(maybe_proto->par());

    _fragment << maybe_proto.release();
    _mol << m.release();
  }

  if (_fragment.empty()) {
    cerr << "BaseFragment::Read:no fragments\n";
    cerr << "Read " << protos_read << " protos, " << failed_constraints << " failed constraints\n";
    return 0;
  }

  return 1;
}

// Fragments with 1 attachment points.
class SubstituentFragments : public BaseFragment<MoleculeOneAtom> {
  private:
    int _natoms;

  public:
    SubstituentFragments();

    int Read(IWString& fname, ReplacementSpecifications& replacement_specifications, int natoms);
};

SubstituentFragments::SubstituentFragments() {
  _natoms = 0;
}

int
SubstituentFragments::Read(IWString& fname, ReplacementSpecifications& replacement_specifications, int natoms) {
  _natoms = natoms;

  return BaseFragment::Read(fname, replacement_specifications);
}

// Fragments with two attachment points.
class Linker2Fragments : public BaseFragment<MoleculeTwoAtoms> {
  private:
    int _distance;

  public:
    Linker2Fragments();

    int Read(IWString& fname, ReplacementSpecifications& replacement_specifications, int dist);
};

Linker2Fragments::Linker2Fragments() {
  _distance = 0;
}

int
Linker2Fragments::Read(IWString& fname, ReplacementSpecifications& replacement_specifications, int dist) {
  _distance = dist;

  return BaseFragment::Read(fname, replacement_specifications);
}

// Fragments with 3 attachment points.
class Linker3Fragments : public BaseFragment<MoleculeThreeAtoms> {
  private:
    int _d1;
    int _d2;
    int _d3;

  public:
    Linker3Fragments();

    int Read(IWString& fname,
                ReplacementSpecifications& replacement_specifications,
                int d1, int d2, int d3);
};

Linker3Fragments::Linker3Fragments() {
  _d1 = 0;
  _d2 = 0;
  _d3 = 0;
}

int
Linker3Fragments::Read(IWString& fname, 
                ReplacementSpecifications& replacement_specifications,
                int d1, int d2, int d3) {
  _d1 = d1;
  _d2 = d2;
  _d3 = d3;

  return BaseFragment::Read(fname, replacement_specifications);
}

// A series of functions all dealing with the case of needing to adjust
// an atom number for the loss of a fragment.
void
AdjustForLossOfFragmentContainingAtom(Molecule& m, 
                const atom_number_t being_removed, 
                atom_number_t& result) {
  m.fragment_membership();
  const atom_number_t initial_value = result;

  int fragment_being_removed = m.fragment_membership(being_removed);

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {

    if (m.fragment_membership(i) != fragment_being_removed) {
      continue;
    }
    if (i < initial_value) {
      --result;
    }
  }

  assert(result >= 0);
}

void
AdjustForLossOfFragmentContainingAtom(Molecule& m, 
                const atom_number_t being_removed, 
                atom_number_t& result1,
                atom_number_t& result2) {
  m.fragment_membership();

  int fragment_being_removed = m.fragment_membership(being_removed);

  const atom_number_t initial_value1 = result1;
  const atom_number_t initial_value2 = result2;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.fragment_membership(i) != fragment_being_removed) {
      continue;
    }
    if (i < initial_value1) {
      --result1;
    }
    if (i < initial_value2) {
      --result2;
    }
  }

  assert(result1 >= 0);
  assert(result2 >= 0);
}

void
AdjustForLossOfFragmentContainingAtom(Molecule& m, 
                const atom_number_t being_removed, 
                atom_number_t& result1,
                atom_number_t& result2,
                atom_number_t& result3) {
  m.fragment_membership();

  int fragment_being_removed = m.fragment_membership(being_removed);

  const atom_number_t initial_value1 = result1;
  const atom_number_t initial_value2 = result2;
  const atom_number_t initial_value3 = result3;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.fragment_membership(i) != fragment_being_removed) {
      continue;
    }
    if (i < initial_value1) {
      --result1;
    }
    if (i < initial_value2) {
      --result2;
    }
    if (i < initial_value3) {
      --result3;
    }
  }

  assert(result1 >= 0);
  assert(result2 >= 0);
  assert(result3 >= 0);
}

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

    Atom_Typing_Specification _atom_typing;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    // The bonds that define how the molecule is split.
    resizable_array_p<Substructure_Query> _bond_break;

    // And what to do with molecules that do not match any of the bond breaking queries.
    int _molecules_with_no_breakable_bonds;
    int _ignore_molecules_not_matching;

    int _min_extra_atoms;
    int _max_extra_atoms;
    int _min_fewer_atoms;
    int _max_fewer_atoms;
    // We can speed up checks to not do anything unless at least one of the above are set.
    int _need_to_check_atom_counts;
    int _rejected_due_to_size_constraints;

    int _molecules_read = 0;

    resizable_array_p<SubstituentFragments> _substituent;
    resizable_array_p<Linker2Fragments> _linker2;
    resizable_array_p<Linker3Fragments> _linker3;

    int _remove_isotopes_before_writing;

    // We can optionally write the parent molecule to the output stream.
    int _write_starting_molecule;

    IWString _name_token_separator;

    uint32_t _max_products_per_starting_molecule;

    uint32_t _substituent_generated;
    uint32_t _linker2_generated;
    uint32_t _linker3_generated;

    int _molecules_not_matching_bond_break_queries;

    int _discard_invalid_valence;
    IWString_and_File_Descriptor _stream_for_invalid_valence;
    int _invalid_valence_generated;

  // Private functions
    int QueriesConsistentWithFragments() const;
    int ProcessInner(Molecule& m,
                 IWString_and_File_Descriptor& output);
    int ReplaceSubstituents(Molecule& m, const Substructure_Results& sresults,
                IWString_and_File_Descriptor& output);
    int ReplaceSubstituent(Molecule& m, const Set_of_Atoms& embedding,
                IWString_and_File_Descriptor& output);
    int ReplaceSubstituent(Molecule& m, atom_number_t attach,
                int atoms_removed,
                const SubstituentFragments& frags,
                IWString_and_File_Descriptor& output);
    int ReplaceLinker2(Molecule& m, const Substructure_Results* sresults,
                IWString_and_File_Descriptor& output);
    int ReplaceLinker2(Molecule& m,
                        atom_number_t a1, atom_number_t a2,
                        int atoms_removed,
                        const Linker2Fragments& frags,
                        IWString_and_File_Descriptor& output);
    int ReplaceLinker3(Molecule& m, const Substructure_Results* sresults,
                IWString_and_File_Descriptor& output);
    int ReplaceLinker3(Molecule& m,
                atom_number_t a1, atom_number_t a2, atom_number_t a3,
                int atoms_removed,
                const Linker3Fragments& frags,
                IWString_and_File_Descriptor& output);
    int OkValence(Molecule& m);
    int OkProduct(Molecule& m, atom_number_t a1, atom_number_t f1);
    int OkProduct(Molecule& m,
                  atom_number_t a1, atom_number_t f1,
                  atom_number_t a2, atom_number_t f2);
    int OkProduct(Molecule& m,
                  atom_number_t a1, atom_number_t f1,
                  atom_number_t a2, atom_number_t f2,
                  atom_number_t a3, atom_number_t f3);
    int OkSizeConstraints(const Molecule& m, int atoms_lost, const Molecule& frag);
    int OkSizeConstraintsInner(int atoms_list, int atoms_gained) const;

    int Write(Molecule& m, const IWString& frag_name,
               IWString_and_File_Descriptor& output) const;

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

    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _molecules_with_no_breakable_bonds = 0;
  _ignore_molecules_not_matching = 0;
  _name_token_separator = " %% ";

  _max_products_per_starting_molecule = 0;

  _min_extra_atoms = -1;
  _max_extra_atoms = -1;
  _min_fewer_atoms = -1;
  _max_fewer_atoms = -1;
  _need_to_check_atom_counts = 0;

  _write_starting_molecule = 0;

  _remove_isotopes_before_writing = 0;
  _rejected_due_to_size_constraints = 0;

  _substituent_generated = 0;
  _linker2_generated = 0;
  _linker3_generated = 0;
    
  _molecules_not_matching_bond_break_queries = 0;

  _discard_invalid_valence = 0;
  _invalid_valence_generated = 0;
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

  ReplacementSpecifications replacement_specifications;
  if (! replacement_specifications.Initialise(cl)) {
    cerr << "Cannot initialise replacement specifications\n";
    Usage(1);
  }

  if (! cl.option_present('F')) {
    cerr << "Must specify one or more sets of fragments via the -F option\n";
    Usage(1);
  }

  if (cl.option_present('F')) {
    IWString fname;
    std::unique_ptr<RE2> frag_rx_1 = std::make_unique<RE2>(R"(\S+_1_(\d+)\.tfdata$)");
    std::unique_ptr<RE2> frag_rx_2 = std::make_unique<RE2>(R"(\S+_2_(\d+)\.tfdata$)");
    std::unique_ptr<RE2> frag_rx_3 = std::make_unique<RE2>(R"(\S+_3_(\d+)\.(\d+)\.(\d+)\.tfdata$)");
    for (int i = 0; cl.value('F', fname, i); ++i) {
      absl::string_view tmp(fname.data(), fname.size());
      int d1;
      int d2;
      int d3;
      if (RE2::FullMatch(tmp, *frag_rx_1, &d1)) {
        std::unique_ptr<SubstituentFragments> s = std::make_unique<SubstituentFragments>();
        if (! s->Read(fname, replacement_specifications, d1)) {
          cerr << "Cannot read substituent '" << fname << "'\n";
          return 0;
        }
        _substituent << s.release();
      } else if (RE2::FullMatch(tmp, *frag_rx_2, &d1)) {
        std::unique_ptr<Linker2Fragments> s = std::make_unique<Linker2Fragments>();
        if (! s->Read(fname, replacement_specifications, d1)) {
          cerr << "Cannot read linker2 '" << fname << "'\n";
          return 0;
        }
        _linker2 << s.release();
      } else if (RE2::FullMatch(tmp, *frag_rx_3, &d1, &d2, &d3)) {
        std::unique_ptr<Linker3Fragments> s = std::make_unique<Linker3Fragments>();
        if (! s->Read(fname, replacement_specifications, d1, d2, d3)) {
          cerr << "Cannot read linker3 '" << fname << "'\n";
          return 0;
        }
        _linker3 << s.release();
      } else {
        cerr << "Ignoring unrecognised file name type '" << fname << '\n';
      } 
    }
  }


  if (_verbose) {
    replacement_specifications.Report(cerr);

    for (const SubstituentFragments* s : _substituent) {
      cerr << "Read " << s->NumberFragments() << " substituents from " << s->FileName() << '\n';
    }
    for (const Linker2Fragments* s : _linker2) {
      cerr << "Read " << s->NumberFragments() << " linker 2 from " << s->FileName() << '\n';
    }
    for (const Linker3Fragments* s : _linker3) {
      cerr << "Read " << s->NumberFragments() << " linker 3 from " << s->FileName() << '\n';
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }
      _bond_break << q.release();
    }
  }

  if (cl.option_present('q')) {
    static constexpr int kVerbose = 0;

    if (! process_queries(cl, _bond_break, kVerbose, 'q')) {
      cerr << "Cannot process bond break queries (-q)\n";
      return 0;
    }
  }

  if (_bond_break.empty()) {
    cerr << "Must specity 1-3 bond breaking queries via the -s and/or -q options\n";
    Usage(0);
  }

  if (_bond_break.size() > 3) {
    cerr << "Do not know how to break " << _bond_break.size() << " bonds at once\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Defined " << _bond_break.size() << " bond break queries\n";
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y.starts_with("minextra=")) {
        y.remove_leading_chars(9);
        if (! y.numeric_value(_min_extra_atoms) || _min_extra_atoms < 0) {
          cerr << "Invalid min extra atoms '" << y << "'\n";
          return 0;
        }
      } else if (y.starts_with("maxextra=")) {
        y.remove_leading_chars(9);
        if (! y.numeric_value(_max_extra_atoms) || _max_extra_atoms < 0) {
          cerr << "Invalid max extra atoms '" << y << "'\n";
          return 0;
        }
      } else if (y.starts_with("minfewer=")) {
        y.remove_leading_chars(10);
        if (! y.numeric_value(_min_fewer_atoms) || _min_fewer_atoms < 0) {
          cerr << "Invalid min fewer atoms '" << y << "'\n";
          return 0;
        }
      } else if (y.starts_with("maxfewer=")) {
        y.remove_leading_chars(10);
        if (! y.numeric_value(_max_fewer_atoms) || _max_fewer_atoms < 0) {
          cerr << "Invalid max fewer atoms '" << y << "'\n";
          return 0;
        }
      } else if (y == "rmiso") {
        _remove_isotopes_before_writing = 1;
        if (_verbose) {
          cerr << "Will remove isotopes before writing\n";
        }
      } else if (y == "wrparent") {
        _write_starting_molecule = 1;
        if (_verbose) {
          cerr << "Will write the starting molecule before products\n";
        }
      } else if (y.starts_with("maxgen=")) {
        y.remove_leading_chars(7);
        if (! y.numeric_value(_max_products_per_starting_molecule) ||
              _max_products_per_starting_molecule < 1) {
          cerr << "The max products per starting molecule generated must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will generate a max of " << _max_products_per_starting_molecule << 
                  " products per starting molecule\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(1);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(1);
      }
    }
  }

  if (_min_extra_atoms < 0 && _max_extra_atoms < 0 &&
      _min_fewer_atoms < 0 && _max_fewer_atoms < 0) {
    _need_to_check_atom_counts = 0;
    if (_verbose) {
      cerr << "Product atom counts will not be checked\n";
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == "i") {
        _ignore_molecules_not_matching = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not maching any of the bond breaking queries\n";
        }
      } else if (z == "help") {
        DisplayDashZOptions();
      }
    }
  }

  if (cl.option_present('V')) {
    _discard_invalid_valence = 1;
    IWString fname = cl.string_value('V');
    if (fname == "none") {
    } else {
      fname.EnsureEndsWith(".smi");
      cerr << "Not null, opening '" << fname << "'\n";
      if (! _stream_for_invalid_valence.open(fname.null_terminated_chars())) {
        cerr << "Cannot open stream for invalid valences '" << fname << "'\n";
        return 0;
      }

      if (_verbose) {
        cerr << "Invalid valences written to '" << fname << "'\n";
      }
    }
  }

  return QueriesConsistentWithFragments();
}

// Return true if the number of bond breaking queries is consistent with the
// fragments and substituents specified.
int
Options::QueriesConsistentWithFragments() const {
  if (_bond_break.size() == 1) {
    if (_linker2.size() > 0 || _linker3.size() > 0) {
      cerr << "One bond breaking query specified, but there are 2 and 3 connected linkers specified\n";
      return 0;
    }
    if (_substituent.empty()) {
      cerr << "One bond break query specified, but no substituents\n";
      return 0;
    }
  }

  if (_bond_break.size() == 2) {
    if (_substituent.size() > 0 || _linker3.size() > 0) {
      cerr << "Two bond breaking queries specified, but there are 1 and 3 connected linkers specified\n";
      return 0;
    }

    if (_linker2.empty()) {
      cerr << "Two bond break queries, but no linker2 substituents\n";
      return 0;
    }
  }

  if (_bond_break.size() == 3) {
    if (_substituent.size() > 0 || _linker2.size() > 0) {
      cerr << "Three bond breaking queries specified, but there are 1 and 2 connected linkers specified\n";
      return 0;
    }

    if (_linker3.empty()) {
      cerr << "Three bond break queries, but no linker3 substituents\n";
      return 0;
    }
  }

  return 1;
}

int
Options::OkSizeConstraints(const Molecule& m, int atoms_lost, const Molecule& frag) {
  if (! _need_to_check_atom_counts) {
    return 1;
  }

  if (! OkSizeConstraintsInner(atoms_lost, frag.natoms())) {
    ++_rejected_due_to_size_constraints;
    return 0;
  }

  return 1;
}

int
Options::OkSizeConstraintsInner(int atoms_lost, int atoms_gained) const {
  // Molecule is staying the same or growing.
  if (atoms_lost <= atoms_gained) {
    int extra_atoms = atoms_gained - atoms_lost;
    if (extra_atoms < _min_extra_atoms) {
      return 0;
    }
    if (extra_atoms > _max_extra_atoms) {
      return 0;
    }
  }

  // Size is staying the same or shrinking

  if (atoms_lost >= atoms_gained) {
    int fewer_atoms = atoms_lost - atoms_gained;
    if (fewer_atoms < _min_fewer_atoms) {
      return 0;
    }
    if (fewer_atoms > _max_fewer_atoms) {
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _molecules_not_matching_bond_break_queries << 
            " molecules did not match the bond breaking queries\n";
  if (_discard_invalid_valence) {
    output << "Discarded " << _invalid_valence_generated << " invalid valence products\n";
  }

  output << _rejected_due_to_size_constraints << " products rejected for size constraints\n";
  output << "Generated " << _substituent_generated << " substituent variations\n";
  output << "Generated " << _linker2_generated << " two connected linkers\n";
  output << "Generated " << _linker3_generated << " three connected linkers\n";

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
Options::Process(Molecule& m, 
                 IWString_and_File_Descriptor& output) {

  int rc = ProcessInner(m, output);

  output.write_if_buffer_holds_more_than(4092);

  return rc;
}

int
Options::ProcessInner(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  Molecule_to_Match target(&m);

  int nq = _bond_break.size();
  std::unique_ptr<Substructure_Results[]> sresults =
                std::make_unique<Substructure_Results[]>(nq);
  for (int i = 0; i < nq; ++i) {
    if (! _bond_break[i]->substructure_search(target, sresults[i])) {
      if (_verbose) {
        cerr << "No match to bond break query " << i << " in " << m.name() << '\n';
      }
      ++_molecules_not_matching_bond_break_queries;
      return _ignore_molecules_not_matching;
    }
  }

  // We may still generate no products, but let's write the parent here
  // if requested.
  if (_write_starting_molecule) {
    static constexpr char kSep = ' ';

    output << m.smiles() << kSep << m.name() << '\n';
  }

  // We have one or more embeddings for each query.

  if (nq == 1) {
    return ReplaceSubstituents(m, sresults[0], output);
  }

  if (nq == 2) {
    return ReplaceLinker2(m, sresults.get(), output);
  }

  if (nq == 3) {
    return ReplaceLinker3(m, sresults.get(), output);
  }

  cerr << "Cannot handle " << nq << " bond breaking queries\n";
  return 0;
}

int
Options::ReplaceSubstituents(Molecule& m, const Substructure_Results& sresults,
                        IWString_and_File_Descriptor& output) {
  for (const Set_of_Atoms * e : sresults.embeddings()) {
    if (! ReplaceSubstituent(m, *e, output)) {
      cerr << "Error processing " << m.smiles() << ' ' << m.name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Options::ReplaceSubstituent(Molecule& m, const Set_of_Atoms& embedding,
                IWString_and_File_Descriptor& output) {
  assert(embedding.size() >= 2);
  uint32_t rc = 0;

  const atom_number_t being_removed = embedding[1];

  Molecule mcopy(m);
  mcopy.remove_bond_between_atoms(embedding[0], being_removed);

  atom_number_t attach0 = embedding[0];
  AdjustForLossOfFragmentContainingAtom(mcopy, being_removed, attach0);

  const int atoms_removed = mcopy.remove_fragment_containing_atom(being_removed);
  assert(mcopy.ok_atom_number(attach0));

  for (const SubstituentFragments* substituent : _substituent) {
    if (ReplaceSubstituent(mcopy, attach0, atoms_removed, *substituent, output)) {
      ++rc;
      if (rc > _max_products_per_starting_molecule) {
        break;
      }
    }
  }

  _substituent_generated += rc;

  return 1;
}

int
Options::ReplaceSubstituent(Molecule& m, atom_number_t attach,
                int atoms_removed,
                const SubstituentFragments& frags,
                IWString_and_File_Descriptor& output) {
  int rc = 0;

  const int initial_natoms = m.natoms();

  if (m.hcount(attach) == 0) {
    return 0;
  }

  for (const MoleculeOneAtom* frag : frags) {
    if (! OkSizeConstraints(m, atoms_removed, *frag)) {
      continue;
    }

    m.add_molecule(frag);
    atom_number_t a2 = initial_natoms + frag->attach();
    if (m.hcount(a2) == 0) {
      continue;
    }
    m.add_bond(attach, a2, kSingleBond);
    if (OkProduct(m, attach, a2)) {
      Write(m, frag->name(), output);
      ++rc;
    }

    m.resize(initial_natoms);
  }

  _substituent_generated += rc;

  return 1;
}

// Return true if a1 and a2 are both heteroatoms.
int
AdjacentHeteroatoms(const Molecule& m, atom_number_t a1, atom_number_t a2) {
  const atomic_number_t z1 = m.atomic_number(a1);
  if (z1 == 6) {
    return 0;
  }

  const atomic_number_t z2 = m.atomic_number(a2);
  if (z2 == 6) {
    return 0;
  }

  return 1;
}

// REturn true if `zatom` could be the Nitrogen (or Oxygen) atom
// of an Aminal
int 
CouldBeNitrogenOfAminal(Molecule& m, atom_number_t zatom) {

  const Atom& nitrogen = m[zatom];

  atomic_number_t z = nitrogen.atomic_number();
  if (z == 6) {
    return 0;
  }
  if (nitrogen.ncon() == 1) {
    return 0;
  }
  if (nitrogen.unsaturation()) {
    return 0;
  }

  if (z == 7) {
  } else if (z == 8) {
  } else {
    return 0;
  }

  return 1;
}

// We are creating a bond between a1 and a2.
// Return true if this forms an aminal.
// We assume that `a1` is the carbon of the aminal, and that
// means that `a2` and another atom must be N (or O).
int
Aminal(Molecule& m, atom_number_t a1, atom_number_t a2) {
  const Atom& carbon = m[a1];
  if (carbon.atomic_number() != 6) {
    return 0;
  }
  if (carbon.ncon() < 3) {
    return 0;
  }
  if (carbon.unsaturation()) {
    return 0;
  }
  // Is this necessary?
  if (m.is_aromatic(a1)) {
    return 0;
  }

  if (! CouldBeNitrogenOfAminal(m, a2)) {
    return 0;
  }

  for (const Bond* b : carbon) {
    atom_number_t o = b->other(a1);
    if (o == a2) {
      continue;
    }
    if (CouldBeNitrogenOfAminal(m, o)) {
      return 1;
    }
  }

  return 0;
}

int
Options::OkValence(Molecule& m) {
  if (! _discard_invalid_valence) {
    return 1;
  }

  if (m.valence_ok()) {
    return 1;
  }

  m.unset_unnecessary_implicit_hydrogens_known_values();
  if (m.valence_ok()) {
    return 1;
  }

  const int matoms = m.natoms();
  for (int i = 0; i <  matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
    m.recompute_implicit_hydrogens(i);
  }

  if (m.valence_ok()) {
    return 1;
  }

  ++_invalid_valence_generated;

  if (_stream_for_invalid_valence.is_open()) {
    static constexpr char kSep = ' ';
    _stream_for_invalid_valence << m.smiles() << kSep << m.name() << '\n';
    _stream_for_invalid_valence.write_if_buffer_holds_more_than(4096);
  }

  return 0;
}

int
Options::OkProduct(Molecule& m, atom_number_t a1, atom_number_t f1) {
  if (! OkValence(m)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a1, f1)) {
    return 0;
  }

  if (Aminal(m, a1, f1)) {
    return 0;
  }
  if (Aminal(m, f1, a1)) {
    return 0;
  }

  return 1;
}

int
Options::Write(Molecule& m, const IWString& frag_name,
               IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';

  if (_remove_isotopes_before_writing) {
    m.unset_isotopes();
  }

  output << m.smiles() << kSep << m.name() << _name_token_separator << frag_name << '\n';

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}
    
// For simplicity, we take only the first embedding
int
Options::ReplaceLinker2(Molecule& m, const Substructure_Results* sresults,
                        IWString_and_File_Descriptor& output) {
  uint32_t rc = 0;

  std::vector<int> count(2);
  count[0] = sresults[0].number_embeddings();
  count[1] = sresults[1].number_embeddings();

  Combinations comb(count);
  std::vector<int> state(2);
  while (comb.Next(state)) {
    atom_number_t remove1 = sresults[0].embedding(state[0])->item(0);
    atom_number_t attach1 = sresults[0].embedding(state[0])->item(1);
    atom_number_t remove2 = sresults[1].embedding(state[1])->item(0);
    atom_number_t attach2 = sresults[1].embedding(state[1])->item(1);

    Molecule mcopy(m);

    mcopy.remove_bond_between_atoms(remove1, attach1);
    mcopy.remove_bond_between_atoms(remove2, attach2);
    write_isotopically_labelled_smiles(mcopy, true, cerr);
    cerr << " atoms " << remove1 << ',' << attach1 << ',' << remove2 << ',' << attach2 << '\n';

    // Some query matches may not set up a proper region.
    if (mcopy.fragment_membership(remove1) != mcopy.fragment_membership(remove2)) {
      continue;
    }

    AdjustForLossOfFragmentContainingAtom(mcopy, remove1, attach1, attach2);

    const int atoms_removed = mcopy.remove_fragment_containing_atom(remove1);
    cerr << "Base molecule ";
    write_isotopically_labelled_smiles(mcopy, true, cerr);
    cerr << " atoms " << remove1 << ',' << attach1 << ',' << remove2 << ',' << attach2 << '\n';

    for (const Linker2Fragments* frag : _linker2) {
      if (ReplaceLinker2(mcopy, attach1, attach2, atoms_removed, *frag, output)) {
        ++rc;
        if (rc > _max_products_per_starting_molecule) {
          break;
        }
      }
    }
    if (rc > _max_products_per_starting_molecule) {
      break;
    }
  }

  _linker2_generated += rc;

  return rc;
}

int
Options::ReplaceLinker2(Molecule& m,
                        atom_number_t a1, atom_number_t a2,
                        int atoms_removed,
                        const Linker2Fragments& frags,
                        IWString_and_File_Descriptor& output) {
  uint32_t rc = 0;

  const int initial_matoms = m.natoms();

  for (const MoleculeTwoAtoms* frag : frags) {
    if (! OkSizeConstraints(m, atoms_removed, *frag)) {
      continue;
    }
    m.add_molecule(frag);

    atom_number_t f1 = initial_matoms + frag->attach1();
    atom_number_t f2 = initial_matoms + frag->attach2();

    m.add_bond(a1, f1, kSingleBond);
    m.add_bond(a2, f2, kSingleBond);
    if (OkProduct(m, a1, f1, a2, f2)) {
      Write(m, frag->name(), output);
      ++rc;
    }

    m.remove_bond_between_atoms(a1, f1);
    m.remove_bond_between_atoms(a2, f2);

    m.add_bond(a1, f2, kSingleBond);
    m.add_bond(a2, f1, kSingleBond);

    if (OkProduct(m, a1, f2, a2, f1)) {
      Write(m, frag->name(), output);
      ++rc;
      if (rc > _max_products_per_starting_molecule) {
        break;
      }
    }

    m.resize(initial_matoms);
  }

  _linker2_generated += rc;

  return rc;
}

int
Options::ReplaceLinker3(Molecule& m, const Substructure_Results* sresults,
                        IWString_and_File_Descriptor& output) {
  uint32_t rc = 0;

  std::vector<int> count(3);
  for (int i = 0; i < 3; ++i) {
    count[i] = sresults[i].number_embeddings();
  }

  std::vector<int> state(3);
  Combinations comb(count);
  while (comb.Next(state)) {
    atom_number_t remove1 = sresults[0].embedding(state[0])->item(0);
    atom_number_t attach1 = sresults[0].embedding(state[0])->item(1);
    atom_number_t remove2 = sresults[1].embedding(state[1])->item(0);
    atom_number_t attach2 = sresults[1].embedding(state[1])->item(1);
    atom_number_t remove3 = sresults[2].embedding(state[2])->item(0);
    atom_number_t attach3 = sresults[2].embedding(state[2])->item(1);

    Molecule mcopy(m);

    mcopy.remove_bond_between_atoms(remove1, attach1);
    mcopy.remove_bond_between_atoms(remove2, attach2);
    mcopy.remove_bond_between_atoms(remove3, attach3);

    AdjustForLossOfFragmentContainingAtom(mcopy, remove1, attach1, attach2, attach3);

    const int atoms_removed = mcopy.remove_fragment_containing_atom(remove1);

    for (const Linker3Fragments* frag : _linker3) {
      if (ReplaceLinker3(mcopy, attach1, attach2, attach3, atoms_removed, *frag, output)) {
        ++rc;
      }
    }
    if (rc > _max_products_per_starting_molecule) {
      break;
    }
  }

  _linker3_generated = rc;

  return 1;
}

int
Options::OkProduct(Molecule& m,
                  atom_number_t a1, atom_number_t f1,
                  atom_number_t a2, atom_number_t f2) {
  if (! OkValence(m)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a1, f1)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a2, f2)) {
    return 0;
  }

  return 1;
}

// Enumerate linker with 3 attachments.
// a1, a2, a3 are atoms that are retained.
// frags is a set of 3 connected linkers.
int
Options::ReplaceLinker3(Molecule& m,
                atom_number_t a1, atom_number_t a2, atom_number_t a3,
                int atoms_removed,
                const Linker3Fragments& frags,
                IWString_and_File_Descriptor& output) {
  int rc = 0;

  int initial_natoms = m.natoms();

  for (const MoleculeThreeAtoms* frag : frags) {
    if (! OkSizeConstraints(m, atoms_removed, *frag)) {
      continue;
    }
    atom_number_t f1 = initial_natoms + frag->attach1();
    atom_number_t f2 = initial_natoms + frag->attach2();
    atom_number_t f3 = initial_natoms + frag->attach3();

    m.add_molecule(frag);

    m.add_bond(a1, f1, kSingleBond);
    m.add_bond(a2, f2, kSingleBond);
    m.add_bond(a3, f3, kSingleBond);

    if (OkProduct(m, a1, f1, a2, f2, a3, f3)) {
      Write(m, frag->name(), output);
      ++rc;
    }

    m.remove_bond_between_atoms(a1, f1);
    m.remove_bond_between_atoms(a2, f2);
    m.remove_bond_between_atoms(a3, f3);

    m.add_bond(a1, f2, kSingleBond);
    m.add_bond(a2, f3, kSingleBond);
    m.add_bond(a3, f1, kSingleBond);

    if (OkProduct(m, a1, f2, a2, f3, a3, f1)) {
      Write(m, frag->name(), output);
      ++rc;
    }

    m.remove_bond_between_atoms(a1, f2);
    m.remove_bond_between_atoms(a2, f3);
    m.remove_bond_between_atoms(a3, f1);

    m.add_bond(a1, f3, kSingleBond);
    m.add_bond(a2, f1, kSingleBond);
    m.add_bond(a3, f2, kSingleBond);

    if (OkProduct(m, a1, f3, a2, f1, a3, f2)) {
      Write(m, frag->name(), output);
      ++rc;
    }

    m.resize(initial_natoms);
  }

  _linker3_generated += rc;

  return rc;
}

int
Options::OkProduct(Molecule& m,
                  atom_number_t a1, atom_number_t f1,
                  atom_number_t a2, atom_number_t f2,
                  atom_number_t a3, atom_number_t f3) {
  if (! OkValence(m)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a1, f1)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a2, f2)) {
    return 0;
  }
  if (AdjacentHeteroatoms(m, a3, f3)) {
    return 0;
  }

  return 1;
}

int
ApplicationName(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! ApplicationName(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i::F:m:M:x:X:p:z:Y:s:q:L:V:");

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

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragment_replace

int
main(int argc, char ** argv) {

  int rc = dicer_fragment_replace::Main(argc, argv);

  return rc;
}
