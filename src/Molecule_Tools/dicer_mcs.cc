// Read output from dicer to find MCS

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace dicer_mcs {

using std::cerr;

int
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Consumes two textproto outputs from dicer performs an MCS between the two.

 -d <dbname>            name of database(s) containing dicer fragments - built by dicer2bdb
 -t                     input is textproto
 -c <natoms>            minimum fragment size to process - discards these fragments upon input
 -C <natoms>            maximum fragment size to process - discards these fragments upon input
 -p <count>             minimum number of exemplars required for a molecule to pass.
                        enter a negative value and all molecules will pass
 -S <dbname>            one or more smiles databases (buildsmidb_bdb). Fragments are looked up as complete
                        molecules in these databases, rather than looking in the fragment databases (-d).
 -X <fname>             write failed molecules (one or more fragments below -p) to <fname>
 -J ...                 look up the most common fragments, enter '-J help' for info
 -Y ...                 various other options, enter '-Y help' for details.
 -F <fname>             write fragment coverage (labelled smiles) to <fname>
 -v                     verbose output
)";
  // clang-format on

  ::exit(rc);
}

class Fragment {
  private:
    Molecule _mol;
    IWString _usmi;

    // So we know how to size various arrays.
    int _atoms_in_parent;

    // A mapping from every atom in `_mol` to an atom number in the parent molecule.
    int* _parent_atom;
    // The isotopic label placed on the parent atom. Same as the canonical ranking.
    isotope_t* _parent_isotope;

  public:
    Fragment();
    ~Fragment();

    int Build(const Molecule& parent, int* tmp, const dicer_data::DicerFragment& proto);

    int natoms() const {
      return _mol.natoms();
    }

    const IWString& usmi() const {
      return _usmi;
    }

    std::unique_ptr<int[]> InParentArray() const {
      return std::unique_ptr<int[]>(new_int(_atoms_in_parent));
    }

    // Atom `parent` must be our parent molecule.
    // Scan through _parent_atom and set the appropriate isotope in `parent`.
    void SetIsotopes(Molecule& parent) const;

    // Wrtie our atom count and the atom count as a fraction of atoms in parent
    void AtomsAndFraction(IWString_and_File_Descriptor& output) const;
};

Fragment::Fragment() {
  _parent_atom = nullptr;
  _parent_isotope = nullptr;
  _atoms_in_parent = 0;
}

Fragment::~Fragment() {
  if (_parent_atom != nullptr) {
    delete [] _parent_atom;
  }
  if (_parent_isotope != nullptr) {
    delete [] _parent_isotope;
  }
}

int
Fragment::Build(const Molecule& parent, int* tmp, const dicer_data::DicerFragment& proto) {
  if (!_mol.build_from_smiles(proto.smi())) {
    cerr << "Fragment::Build:invalid smiles '" << proto.smi() << "'\n";
    return 0;
  }

  std::unique_ptr<isotope_t[]> isosave = _mol.GetIsotopes();

  _mol.transform_to_non_isotopic_form();

  // After removal of isotopes, compute the canonical rank.
  const int* canonical_rank = _mol.canonical_ranks();

  _atoms_in_parent = parent.natoms();

  const int matoms = _mol.natoms();

  _parent_atom = new int[matoms];
  _parent_isotope = new isotope_t[matoms];

  for (int i = 0; i < matoms; ++i) {
    isotope_t iso = isosave[i];

    atom_number_t p = parent.atom_with_isotope(iso);
    if (p == kInvalidAtomNumber) {
      cerr << "Fragment::Build:no parent atom with isotope " << iso << '\n';
      return 0;
    }

    _parent_atom[i] = p;
    _parent_isotope[i] = canonical_rank[i];
  }

  _usmi = _mol.unique_smiles();

  return 1;
}

void
Fragment::SetIsotopes(Molecule& parent) const {
  const int matoms = _mol.natoms();

  for (int i = 0; i < matoms; ++i) {
    parent.set_isotope(_parent_atom[i], _parent_isotope[i]);
  }
}
    
void
Fragment::AtomsAndFraction(IWString_and_File_Descriptor& output) const {
  output << _mol.natoms() << ' ' << iwmisc::Fraction<float>(_mol.natoms(),
                _atoms_in_parent);
}

struct ComparisonConditions {
  // Each molecule consists of a number of fragments.
  // We scan each list looking for fragments that match.
  // If we stop at the first match, that will be a fully connected MCS.
  // If we allow disconnected MCS to be found then we need to
  // continue looking for non-overlapping fragments.
  int allow_disconnected_matches = 1;

  char record_separator = '\n';

  Accumulator_Int<uint32_t> _acc_natoms_overlap;
  Accumulator<double> _acc_fraction_overlap;
};

class DicedMolecule {
  private:
    // The parent molecule. When read in it must have isotopic labels.

    Molecule _mol;

    // Across all our fragments, a mapping from unique smiles to 
    // fragment number.
    absl::flat_hash_map<IWString, int> _usmi;

    resizable_array_p<Fragment> _fragment;

    std::unique_ptr<uint32_t[]> _atype;

    // Private functions
    int CompareInner(DicedMolecule& rhs, const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output);
    int ProcessMatch(int lhs_frag_number, DicedMolecule& rhs,
                int rhs_frag_number,
                const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output);
    int TryExpansion(DicedMolecule& rhs,
                const ComparisonConditions& cmp,
                isotope_t& next_isotope_to_assign);
    atom_number_t UnmatchedAtom(atom_number_t zatom, const Bond* b1,
                             uint32_t atype) const;

  public:
    DicedMolecule();
    ~DicedMolecule();

    int Build(const dicer_data::DicedMolecule& proto);

    void RemoveSmallFragments(int min_size);

    int size() const {
      return _fragment.number_elements();
    }

    const IWString& name() const {
      return _mol.name();
    }

    // If _usmi contains `s`, return _usmi[s].
    std::optional<int> Contains(const IWString& s) const;

    int Compare(DicedMolecule& rhs, const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output);

    int AssignAtomTypes(Atom_Typing_Specification& ats);
};

DicedMolecule::DicedMolecule() {
}

DicedMolecule::~DicedMolecule() {
}

int
DicedMolecule::Build(const dicer_data::DicedMolecule& proto) {
  if (! _mol.build_from_smiles(proto.smiles())) {
    cerr << "DicedMolecule::Build:bad smiles\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  _mol.set_name(proto.name());

  const int number_fragments = proto.fragment_size();

  if (number_fragments == 0) {
    cerr << "DicedMolecule::Build:no fragments\n";
    cerr << proto.ShortDebugString() << '\n';
    return 1;
  }

  _fragment.reserve(number_fragments);

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(_mol.natoms());
  cerr << "Paremt molecule contains " << _mol.natoms() << " atoms\n";

  for (int i = 0; i < number_fragments; ++i) {
    std::unique_ptr<Fragment> f = std::make_unique<Fragment>();

    if (! f->Build(_mol, tmp.get(), proto.fragment(i))) {
      cerr << "DicedMolecule::Build:invalid fragment\n";
      cerr << proto.fragment(i).ShortDebugString() << '\n';
      return 0;
    }

    _fragment << f.release();
  }

  std::sort(_fragment.rawdata(), _fragment.rawdata() + number_fragments, [](const Fragment* f1, const Fragment* f2) {
    return f1->natoms() > f2->natoms();
  });

  for (int i = 0; i < number_fragments; ++i) {
    _usmi[_fragment[i]->usmi()] = i;
  }

  _mol.transform_to_non_isotopic_form();

  return 1;
}

void
DicedMolecule::RemoveSmallFragments(int min_size) {
  for (int i = _fragment.number_elements(); i >= 0; --i) {
    if (_fragment[i]->natoms() >= min_size) {
      continue;
    }
    _fragment.resize(i);
    return;
  }
}

std::optional<int>
DicedMolecule::Contains(const IWString& s) const {
  const auto iter = _usmi.find(s);
  if (iter == _usmi.end()) {
    return std::nullopt;
  }

  return iter->second;
}

int
DicedMolecule::AssignAtomTypes(Atom_Typing_Specification& ats) {
  _atype = std::make_unique<uint32_t[]>(_mol.natoms());
  return ats.assign_atom_types(_mol, _atype.get());
}

int
DicedMolecule::Compare(DicedMolecule& rhs, const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output) {
//std::unique_ptr<int[]> lhs_ipa = InParentArray();
//std::unique_ptr<int[]> rhs_ipa = rhs.InParentArray();

  int rc = CompareInner(rhs, cmp, output);

  _mol.transform_to_non_isotopic_form();
  rhs._mol.transform_to_non_isotopic_form();

  return rc;
}

int
DicedMolecule::CompareInner(DicedMolecule& rhs, const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output) {

  int matches_found = 0;
  for (int i = 0; i < _fragment.number_elements(); ++i) {
    std::optional<int> maybe_rhs = rhs.Contains(_fragment[i]->usmi());
    if (! maybe_rhs) {
      continue;
    }

    ProcessMatch(i, rhs, *maybe_rhs, cmp, output);

    if (_atype) {
      isotope_t next_isotope_to_assign = std::max<isotope_t>(_mol.natoms(), rhs._mol.natoms());
      isotope_t isave = next_isotope_to_assign;
      TryExpansion(rhs, cmp, next_isotope_to_assign);
      if (next_isotope_to_assign > isave) {
        ProcessMatch(i, rhs, *maybe_rhs, cmp, output);
      }
    }

    if (! cmp.allow_disconnected_matches) {
      return 1;
    }
  }

  return 0;
}

int
DicedMolecule::ProcessMatch(int lhs_frag_number, DicedMolecule& rhs,
                int rhs_frag_number,
                const ComparisonConditions& cmp,
                IWString_and_File_Descriptor& output) {
  const Fragment& lhs_frag = *_fragment[lhs_frag_number];
  const Fragment& rhs_frag = *rhs._fragment[rhs_frag_number];

  lhs_frag.SetIsotopes(_mol);
  rhs_frag.SetIsotopes(rhs._mol);

  static constexpr char kSep = ' ';

  output << _mol.smiles() << kSep << _mol.name() << kSep << "LHS" << kSep;
  lhs_frag.AtomsAndFraction(output);
  output << cmp.record_separator;
  output << rhs._mol.smiles() << kSep << rhs._mol.name() << kSep << "RHS" << kSep;
  rhs_frag.AtomsAndFraction(output);
  output << cmp.record_separator;

  return 1;
}

int
DicedMolecule::TryExpansion(DicedMolecule& rhs,
                const ComparisonConditions& cmp,
                isotope_t& next_isotope_to_assign) {
  int rc = 0;

  const int matoms = _mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso1 = _mol.isotope(i);
    if (iso1 == 0) {
      continue;
    }

    atom_number_t a2 = rhs._mol.atom_with_isotope(iso1);
    assert(a2 != kInvalidAtomNumber);
    cerr << " lhs at " << i << " in system isotope " << _mol.isotope(i) << " ncon " << _mol.ncon(i) << '\n';

    for (const Bond * b : _mol[i]) {
      atom_number_t o = b->other(i);
      if (_mol.isotope(o) > 0) {
        continue;
      }

      atom_number_t rhs_atom = rhs.UnmatchedAtom(a2, b, _atype[o]);
      if (rhs_atom == kInvalidAtomNumber) {
        cerr << "No matchon other side\n";
        continue;
      }
      cerr << " rhs at " << o << " being added " << rhs._mol.smarts_equivalent_for_atom(rhs_atom) << '\n';

      _mol.set_isotope_no_perturb_canonical_ordering(o, next_isotope_to_assign);
      rhs._mol.set_isotope_no_perturb_canonical_ordering(rhs_atom, next_isotope_to_assign);
      ++next_isotope_to_assign;
      ++rc;
    }
  }

  if (rc == 0) {
    return 1;
  }

  return TryExpansion(rhs, cmp, next_isotope_to_assign);
}

// Look for an isotope 0 connection to `zatom`, connected via a
// bond of type `b` with the same atom type `atype`.
atom_number_t
DicedMolecule::UnmatchedAtom(atom_number_t zatom, const Bond* b1,
                             uint32_t atype) const {
  cerr << "Looking for unmatched atom from " << zatom << " isotope " << _mol.isotope(zatom) << '\n';

  cerr << " atom has " << _mol.ncon(zatom) << " connections\n";
  for (const Bond* b : _mol[zatom]) {
    cerr << "same_bond_type " << b->same_bond_type(*b1) << '\n';
    if (! b->same_bond_type(*b1)) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    cerr << " atom " << o << " isotpe " << _mol.isotope(o) << '\n';
    if (_mol.isotope(o) > 0) {
      continue;
    }

    cerr << "Atype " << _atype[o] << " cmp " << atype << '\n';
    if (_atype[o] != atype) {
      continue;
    }

    return o;
  }

  return kInvalidAtomNumber;
}

class DicedMolecules {
  private:
    int _number_molecules;

    DicedMolecule* _mol;

  public:
    DicedMolecules();
    ~DicedMolecules();

    int Build(const char* fname);
    int Build(iwstring_data_source& input);

    int number_molecules() const {
      return _number_molecules;
    }

    void RemoveSmallFragments(int min_size);

    DicedMolecule& operator[](int i) {
      return _mol[i];
    }
    
    void AssignAtomTypes(Atom_Typing_Specification& ats);
};

DicedMolecules::DicedMolecules() {
  _number_molecules = 0;
  _mol = nullptr;
}

DicedMolecules::~DicedMolecules() {
  if (_mol != nullptr) {
    delete [] _mol;
  }
}

int
DicedMolecules::Build(iwstring_data_source& input) {
  _number_molecules = input.records_remaining();
  if (_number_molecules == 0) {
    cerr << "DicedMolecules:Build: empty file\n";
    return 0;
  }

  _mol = new DicedMolecule[_number_molecules];

  const_IWSubstring buffer;
  for (int i = 0; i < _number_molecules; ++i) {
    input.next_record(buffer);
    dicer_data::DicedMolecule proto;
    google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());
    if (! google::protobuf::TextFormat::Parse(&input, &proto)) {
      cerr << "DicedMolecule::Build:::invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
    if (! _mol[i].Build(proto)) {
      cerr << "DicedMolecules::Build:invalid data\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
DicedMolecules::Build(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DicedMolecule::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

void
DicedMolecules::RemoveSmallFragments(int min_size) {
  for (int i = 0; i < _number_molecules; ++i) {
    _mol[i].RemoveSmallFragments(min_size);
  }
}

void
DicedMolecules::AssignAtomTypes(Atom_Typing_Specification& ats) {
  for (int i = 0; i < _number_molecules; ++i) {
    _mol[i].AssignAtomTypes(ats);
  }
}

class Options {
  private:
    int _verbose;

  public:
    Options();

    int Initialise(Command_Line& cl);
};

Options::Options() {
  _verbose = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:cm:P:w");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.size() != 2) {
    cerr << "Must specify one or two files of dicer_data::DicedMolecule textproto inputs\n";
    Usage(1);
  }

  DicedMolecules diced1, diced2;

  if (! diced1.Build(cl[0])) {
    cerr << "Cannot read fragments from '" << cl[0] << "'\n";
    return 1;
  }

  if (! diced2.Build(cl[1])) {
    cerr << "Cannot read fragments from '" << cl[1] << "'\n";
    return 1;
  }

  if (verbose) {
    cerr << "Read " << diced1.number_molecules() << " molecules from " << cl[0] << '\n';
    cerr << "Read " << diced2.number_molecules() << " molecules from " << cl[1] << '\n';
  }

  if (cl.option_present('m')) {
    int min_count;
    if (! cl.value('m', min_count) || min_count < 1) {
      cerr << "The min fragment size (-m) option must be a whole +ve number\n";
      return 1;
    }

    diced1.RemoveSmallFragments(min_count);
    diced2.RemoveSmallFragments(min_count);
  }

  if (cl.option_present('P')) {
    IWString p = cl.string_value('P');

    Atom_Typing_Specification ats;
    if (! ats.build(p)) {
      cerr << "INvalid atom typing specification '" << p << "'\n";
      return 1;
    }

    diced1.AssignAtomTypes(ats);
    diced2.AssignAtomTypes(ats);
  }

  bool ok_compare_same_name = false;
  if (cl.option_present('w')) {
    ok_compare_same_name = true;
  }

  IWString_and_File_Descriptor output(1);

  ComparisonConditions cmp;

  cmp.allow_disconnected_matches = 0;

  for (int i = 0; i < diced1.number_molecules(); ++i) {
    DicedMolecule& m1 = diced1[i];
    for (int j = 0; j < diced2.number_molecules(); ++j) {
      DicedMolecule& m2 = diced2[j];
      if (ok_compare_same_name) {
      } else if (m1.name() == m2.name()) {
        continue;
      }
      m1.Compare(m2, cmp, output);
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  output.flush();

  return 0;
}

}  // namespace dicer_mcs

int
main(int argc, char ** argv) {

  int rc = dicer_mcs::Main(argc, argv);

  return rc;
}
