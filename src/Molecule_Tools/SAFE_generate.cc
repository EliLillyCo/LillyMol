// Consume the output of mol2SAFE and generate random variants.

#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/matcher.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/molecule_filter_lib.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace safe_generate {

using std::cerr;

constexpr char kPercent = '%';

constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

using iwmatcher::Matcher;
using molecule_filter_lib::MoleculeFilter;

// We ignore any fragments with more than this many connections.
int max_ncon = 10;

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
  cerr << R"(Denovo generation of molecules from SAFE smiles
Consumes the output from mol2SAFE, and the -L file is also from mol2SAFE.
 -L <fname>             fragment library of SAFE fragments
 -Y <query>             queries for atoms that are allowed to change
 -N <query>             queries for atoms that are NOT allowed to change
 -x extra=<n>           the number of extra atoms in a fragment being substitued
 -x fewer=<n>           the number of fewer atoms in a fragment being substitued
 -n <n>                 number of molecules to generate
 -v                     verbose output
)";
// clang-format on

  ::exit(rc);
}

constexpr int kCarbon = 0;
constexpr int kArCarbon = 1;
constexpr int kNitrogen = 2;
constexpr int kArNitrogen = 3;
constexpr int kOxygen = 4;
constexpr int kArOxygen = 5;
constexpr int kFluorine = 6;
constexpr int kPhosphorus = 7;
constexpr int kSulphur = 8;
constexpr int kArSulphur = 9;
constexpr int kChlorine = 10;
constexpr int kBromine = 11;
constexpr int kIodine = 12;
constexpr int kOther = 13;

class MFormula {
  private:
    int _count[kOther];

  // Private functions
    void ZeroCountArray();

  public:
    MFormula();

    int Build(Molecule& m);
};

void
MFormula::ZeroCountArray() {
  std::fill_n(_count, kOther, 0);
}

MFormula::MFormula() {
  ZeroCountArray();
}

int
MFormula::Build(Molecule& m) {
  m.compute_aromaticity_if_needed();

  ZeroCountArray();

  for (int i = 0; i < m.natoms(); ++i) {
    atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      if (m.is_aromatic(i)) {
        ++_count[kArCarbon];
      } else {
        ++_count[kCarbon];
      }
    } else if (z == 7) {
      if (m.is_aromatic(i)) {
        ++_count[kArNitrogen];
      } else {
        ++_count[kNitrogen];
      }
    } else if (z == 8) {
      if (m.is_aromatic(i)) {
        ++_count[kArOxygen];
      } else {
        ++_count[kOxygen];
      }
    } else if (z == 9) {
      ++_count[kFluorine];
    } else if (z == 15) {
      ++_count[kPhosphorus];
    } else if (z == 16) {
      if (m.is_aromatic(i)) {
        ++_count[kArSulphur];
      } else {
        ++_count[kSulphur];
      }
    } else if (z == 17) {
      ++_count[kChlorine];
    } else if (z == 35) {
      ++_count[kBromine];
    } else if (z == 53) {
      ++_count[kIodine];
    } else {
      ++_count[kOther];
    }
  }

  return 1;
}

class SafeFragment {
  private:
    // the smiles of the fragment
    // Individual dot separated tokens from something like
    // [1Cl]%10.[1C]%11.[1C]%121=NO[1C]%13%11O1.[1C]%131CCCCN1.C1CC[1C]%10[1C]%12C1
    IWString _smiles;

    // The indices within _smiles where the first digit after the % signs are.
    resizable_array<int> _first_digit;

    // A molecule build from _smiles by dropping the %nn characters
    Molecule _m;

    MFormula _mformula;

    // Number of atoms in the fragment.
    int _natoms;

    // Number of connections (ring openings) in the fragment.
    int _ncon;

    // If two or more connections, the number of bonds between the closest two.
    int _distance;

    // Within _smiles the ring openings present.
    std::vector<int> _ring;

    // In a molecule consisting of individual SafeFragment's each fragment
    // can be marked as ok to select or not - depending on substructure
    // matches in the parent molecule.
    int _ok_to_select;

  // Private functions
    int ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination);

  public:
    SafeFragment();

    int DebugPrint(std::ostream& output) const;

    const IWString& name() const {
      return _m.name();
    }

    int natoms() const {
      return _natoms;
    }

    int ncon() const {
      return _ncon;
    }

    int Build(const const_IWSubstring& buffer);
    int Build(const dicer_data::DicerFragment& proto);

    int ok_to_select() const {
      return _ok_to_select;
    }
    void set_ok_to_select(int s) {
      _ok_to_select = s;
    }

    const IWString& smiles () const {
      return _smiles;
    }

    Molecule& mol() {
      return _m;
    }

    // Place into `new_smiles` the smiles of `f2` with the ring numbers
    // from `this`.
    int SameNumbers(const SafeFragment& f2, IWString& new_smiles) const;
};

SafeFragment::SafeFragment() {
  _natoms = 0;
  _ncon = 0;
  _distance = 0;
}

int
SafeFragment::DebugPrint(std::ostream& output) const {
  output << "SafeFragment smiles " << _smiles << " natoms " << _natoms << " ncon " << _ncon << '\n';
  output << "first";;
  for (int f : _first_digit) {
    output << ' ' << f << " '" << _smiles[f] << "'";
  }
  output << '\n';

  return output.good();
}

#ifdef NO_TUSED
int
SafeFragment::Build(const const_IWSubstring& line) {
  google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "SafeFragment:Build:cannot parse proto " << line << '\n';
    return 0;
  }

  return Build(proto);
}
#endif

int
SafeFragment::Build(const const_IWSubstring& smi) {
  const int nchars = smi.length();
  // cerr << "SafeFragment::Build from '" << smi << "'\n";
  // Must have at least an atomic symbol followed by % and a two digit ring number.
  if (nchars < 4) {
    cerr << "SafeFragment::Build:too short '" << smi << "'\n";
    return 0;
  }

  IWString smiles;
  smiles.reserve(smi.length());

  for (int i = 0; i < nchars; ++i) {
    const char c = smi[i];
    if (c == kPercent) {
      int n = (smi[i + 1] - '0') * 10 + smi[i + 2] - '0';
      _ring.push_back(n);
      _first_digit << (i + 1);
      i += 2;
    } else {
      smiles << c;
    }
  }

  if (! _m.build_from_smiles(smiles)) {
    cerr << "SafeFragment::Build:invalid smiles '" << smiles << "' from '" << smi << "'\n";
    return 0;
  }

  _natoms = _m.natoms();
  _ncon = _first_digit.number_elements();
  _smiles = smi;

  const int niso = _m.number_isotopic_atoms();
  if (niso == 0) {
    cerr << "SafeFragment::Build:no isotopes '" << smiles << "' from '" << smi << "'\n";
    return 0;
  }

  if (niso > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (_m.isotope(i) == 0) {
        continue;
      }
      for (int j = i + 1; j < _natoms; ++j) {
        if (_m.isotope(j) == 0) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  _mformula.Build(_m);

  return 1;
}

int
SafeFragment::ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination) {
  assert(smi[i] == kOpenSquareBracket);
  destination << kOpenSquareBracket;
  ++i;

  // Should not happen.
  if (i == smi.length()) {
    return 0;
  }

  int got_isotope = 0;
  if (std::isdigit(smi[i])) {
    got_isotope = 1;
  }

  for ( ;i < smi.length(); ++i) {
    destination << smi[i];
    if (smi[i] == kCloseSquareBracket) {
      break;
    }
  }

  if (got_isotope) {
    // Skip over any ring closures already present.
    for ( ; (i + 1) < smi.length() && std::isdigit(smi[i + 1]); ++i) {
      destination << smi[i + 1];
    }
    destination << '%';
    _first_digit << destination.length();
    destination << next_ring;
    ++next_ring;
  }

  return 1;
}

int
SafeFragment::Build(const dicer_data::DicerFragment& proto) {
  IWString smi = proto.smi();
  if (smi.empty()) {
    cerr << "SafeFragment::Build:empty smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (! _m.build_from_smiles(smi)) {
    cerr << "SafeFragment::Build:invalid smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _m.set_name(proto.par());

  _mformula.Build(_m);

  _natoms = _m.natoms();
  // cerr << "SafeFragment::build:smiles " << smi << " natoms " << _natoms << '\n';

  _ncon = _m.number_isotopic_atoms();
  if (_ncon == 0) {
    cerr << "SafeFragment::Build:no isotopes '" << smi << "' from '" <<
            proto.ShortDebugString() << "'\n";
    return 0;
  }

  if (_ncon > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (_m.isotope(i) == 0) {
        continue;
      }
      for (int j = i + 1; j < _natoms; ++j) {
        if (_m.isotope(j) == 0) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  // The smiles that is read does not have the %nn ring opening/closing motifs.
  // We insert those after the isotopic labels.

  IWString tmp;
  tmp.reserve(smi.length() + 9);
  int next_ring = 11;
  for (int i = 0; i < smi.length(); ++i) {
    const char c = smi[i];
    if (c == kOpenSquareBracket) {
      ProcessSquareBracket(smi, i, next_ring, tmp);
    } else {
      tmp << c;
    }
  }

  // cerr << "from " << smi << " transform to '" << tmp << " atoms " << _natoms << '\n';
  _smiles = tmp;

  return 1;
}

class SafedMolecule {
  private:
    Molecule _m;

    // the SAFE smiles from which we are built.
    IWString _smiles;

    // the indices where each fragment starts in _smiles
    resizable_array<int> _frag_start_smiles;
    // the indices where each fragment starts in _m;
    resizable_array<int> _frag_start_atom;

    resizable_array_p<SafeFragment> _frag;

    // When random items are selected, we limit the number of tries.
    int _max_attempts;

    std::mt19937 _rng;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _dist;

    // Private functions
    int GetFragment(atom_number_t zatom) const;

  public:
    SafedMolecule();

    int Build(const const_IWSubstring& buffer);

    int DebugPrint(std::ostream& output) const;

    Molecule& mol() {
      return _m;
    }

    const Molecule& mol() const {
      return _m;
    }

    const IWString& name() const {
      return _m.name();
    }

    const IWString& smiles() const {
      return _smiles;
    }

    int number_fragments() const {
      return _frag.number_elements();
    }

    const SafeFragment* fragment(int ndx) const {
      return _frag[ndx];
    }

    int SetupRng();

    // Given a set of queries, perform the searches. For every match examine
    // the matched atoms, and any fragment that contains a matched atom is
    // marked for being OK to change.
    int IdentifyChangingFragments(resizable_array_p<Substructure_Query>& queries);
    int IdentifyUnChangingFragments(resizable_array_p<Substructure_Query>& queries);

    SafeFragment* ChooseFragment(int ncon, const Matcher<int>& natoms);

    // Return the index of a random fragment.
    std::optional<int> RandomFragment();

    int NewSmiles(int f1_ndx, const SafeFragment& f2, IWString& new_smiles) const;
};

SafedMolecule::SafedMolecule() {
  std::random_device rd;
  _rng.seed(rd());
  _max_attempts = 10;  // arbitrary number.
}

int
SafedMolecule::SetupRng() {
  if (_frag.empty()) {
    cerr << "SafedMolecule::SetupRng:empty\n";
    return 0;
  }

  _dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _frag.size() - 1);

  return 1;
}

int
SafedMolecule::DebugPrint(std::ostream& output) const {
  output << "SafedMolecule with " << _frag.size() << " fragments\n";
  output << _smiles << ' ' << _m.name() << '\n';

  for (int i = 0; i < _frag.number_elements(); ++i) {
    output << i << ' ';
    _frag[i]->DebugPrint(output);
  }
  return 1;
}

int
SafedMolecule::Build(const const_IWSubstring& buffer) {
  const_IWSubstring smiles, id;
  if (! buffer.split(smiles, ' ', id) ||
        smiles.empty() || id.empty()) {
    cerr << "SafedMolecule::Build:cannot split to smiles id '" << buffer << "'\n";
    return 0;
  }

  if (! _m.build_from_smiles(smiles)) {
    cerr << "SafedMolecule::build_from_smiles:invalid smiles '" << smiles << "' from '" <<
             buffer << "'\n";
    return 0;
  }

  _smiles = smiles;

  _m.set_name(id);

  const_IWSubstring token;
  int i = 0;
  int prev_i = i;
  int atom_count = 0;
  while (buffer.nextword(token, i, '.')) {
    std::unique_ptr<SafeFragment> f = std::make_unique<SafeFragment>();
    if (! f->Build(token)) {
      cerr << "SafedMolecule::Build:invalid fragment '" << token << "' from '" <<
              buffer << "'\n";
      return 0;
    }

    atom_count += f->natoms();
    _frag << f.release();
    _frag_start_smiles << prev_i;
    _frag_start_atom << atom_count;
    prev_i = i;
  }

#ifdef DEBUG_BUILD_SAFE_MOLECULE
  cerr << "From\n" << buffer << '\n';
  for (int s : _frag_start_smiles) {
    cerr << ' ' << s;
  }
  cerr << '\n';
#endif

  return 1;
}

#ifdef NOT_USED
int
SafedMolecule::Build(const dicer_data::DicerFragment& proto) {
  IWString smi = proto.smi();

  if (smi.empty()) {
    cerr << "SafedMolecule::Build:no smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _m.build_from_smiles(smi)) {
    cerr << "SafedMolecule::Build:invalid smiles '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  return 1;
}
#endif

// Return the fragment number where `zatom` is found.
// This could be made more efficient with an atom_to_fragment
// array, but this is not expected to be a bottleneck.
int
SafedMolecule::GetFragment(atom_number_t zatom) const {
  for (int i = _frag_start_atom.number_elements() - 1; i >= 0; --i) {
    if (zatom >= _frag_start_atom[i]) {
      return i;
    }
  }

  cerr << "SafedMolecule::GetFragment:this should not happen\n";
  return -1;
}

int
SafedMolecule::IdentifyChangingFragments(resizable_array_p<Substructure_Query>& queries) {
  for (SafeFragment* f : _frag) {
    f->set_ok_to_select(0);
  }

  int rc = 0;

  Molecule_to_Match target(&_m);

  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* c : sresults.embeddings()) {
      for (atom_number_t a : *c) {
        int f = GetFragment(a);
        _frag[f]->set_ok_to_select(1);
        rc = 1;
      }
    }
  }

  return rc;
}

int
SafedMolecule::IdentifyUnChangingFragments(resizable_array_p<Substructure_Query>& queries) {
  int rc = 0;

  Molecule_to_Match target(&_m);

  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* c : sresults.embeddings()) {
      for (atom_number_t a : *c) {
        int f = GetFragment(a);
        _frag[f]->set_ok_to_select(0);
        rc = 1;
      }
    }
  }

  return rc;
}


#ifdef QWEQWEQWE
SafeFragment*
SafedMolecule::ChooseFragment(int ncon, const Matcher<int>& natoms) {
  for (int i = 0; i < _max_attempts; ++i) {
    const int ndx = (*_dist)(_rng);
    SafeFragment* f = _frag[ndx];
    if (! f->ok_to_select()) {
      continue;
    }
    if (f->ncon() != ncon) {
      continue;
    }
    if (! natoms.matches(f->natoms())) {
      continue;
    }

    return f;
  }

  return nullptr;
}
#endif

int
SafeFragment::SameNumbers(const SafeFragment& f2, IWString& new_smiles) const {
  new_smiles = f2._smiles;

  // Now copy our ring numbers to `new_smiles`.
  const int nf = _first_digit.number_elements();
  for (int i = 0; i < nf; ++i) {
    new_smiles[f2._first_digit[i]] = _smiles[_first_digit[i]];
    new_smiles[f2._first_digit[i] + 1] = _smiles[_first_digit[i] + 1];
  }

  return 1;
}

std::optional<int>
SafedMolecule::RandomFragment() {
  for (int i = 0; i < _max_attempts; ++i) {
    const int ndx = (*_dist)(_rng);
    SafeFragment* f = _frag[ndx];
    if (f->ok_to_select()) {
      return ndx;
    }
  }

  return std::nullopt;
}

int
SafedMolecule::NewSmiles(int f1_ndx, const SafeFragment& f2, IWString& new_smiles) const {
  const SafeFragment* f1 = _frag[f1_ndx];
  f1->SameNumbers(f2, new_smiles);

  return 1;
}

class Library {
  private:
    // A mapping from <ncon, natoms> to the molecules fitting that description.
    absl::flat_hash_map<std::tuple<uint32_t, uint32_t>, resizable_array_p<SafeFragment>> _ff;

    // As we read the library, we keep track of the number of fragments
    // we discard because of too many connections
    int _excessive_connections;

    std::mt19937 _rng;

    int _max_attempts;

  public:
    Library();
    ~Library();

    int Build(IWString& fname);
    int Build(iwstring_data_source& input);
    int BuildMember(const const_IWSubstring& line);

    int SetupRng();

    // Change the isotopic label of each atom that has an isotopic label to `iso`.
    int SetIsotope(isotope_t iso);

    uint32_t size() const;

    const SafeFragment* GetFragment(int ncon, int natoms);
};

Library::Library() {
  _max_attempts = 10;  // arbitrary number.
  _excessive_connections = 0;
}

Library::~Library() {
}

int
Library::SetupRng() {
  std::random_device rd;
  _rng.seed(rd());

  return 1;
}

uint32_t
Library::size() const {
  uint32_t rc = 0;

  for (const auto& [k, v] : _ff) {
    rc += v.size();
  }

  return rc;
}

const SafeFragment*
Library::GetFragment(int ncon, int natoms) {
  std::tuple<uint32_t, uint32_t> key(ncon, natoms);
  auto iter = _ff.find(key);
  if (iter == _ff.end()) {
    return nullptr;
  }

  std::uniform_int_distribution<uint32_t> u(0, iter->second.size() - 1);
  uint32_t ndx = u(_rng);

  return iter->second.item(ndx);
}

int
Library::SetIsotope(isotope_t iso) {
  for (const auto& [_, v] : _ff) {
    for (SafeFragment* f : v) {
      Molecule& m = f->mol();
      const int matoms = m.natoms();
      for (int i = 0; i < matoms; ++i) {
        if (m.isotope(i)) {
          m.set_isotope(i, iso);
        }
      }
    }
  }

  return 1;
}

int
Library::Build(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Library::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
Library::Build(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    // cerr << "Library building from '" << buffer << "'\n";
    if(! BuildMember(buffer)) {
      cerr << "Library::Build:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Library::BuildMember(const const_IWSubstring& line) {
  google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "Library:BuildMember:cannot parse proto " << line << '\n';
    return 0;
  }

  // cerr << "Proto built " << proto.ShortDebugString() << '\n';
  std::unique_ptr<SafeFragment> f = std::make_unique<SafeFragment>();
  if (! f->Build(proto)) {
    cerr << "line::BuildMember:cannot parse '" << line << "'\n";
    return 0;
  }

  // cerr << "Fragment contains " << f->ncon() << " connections, cmp " << max_ncon << '\n';
  if (f->ncon() > max_ncon) {
    ++_excessive_connections;
    return 1;
  }

  std::tuple<uint32_t, uint32_t> key(f->ncon(), f->natoms());
  auto iter = _ff.find(key);
  if (iter != _ff.end()) {
    iter->second << f.release();
    return 1;
  }

  resizable_array_p<SafeFragment> tmp;
  tmp << f.release();
  _ff.emplace(key, std::move(tmp));

  return 1;
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

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    resizable_array_p<Library> _library;

    resizable_array_p<SafedMolecule> _mols;

    // The queries instantiated from the -Y and -N options.
    resizable_array_p<Substructure_Query> _can_change;
    resizable_array_p<Substructure_Query> _cannot_change;

    int _ignore_molecules_not_matching_queries;

    uint32_t _extra_atoms;
    uint32_t _fewer_atoms;

    // The -X option.
    resizable_array_p<Substructure_Query> _discard_if_match;

    std::mt19937 _rng;

    // when asked to generate fragments, this is how many times we
    // try.
    int _max_attempts;

    // Random number generators for the _mols and _library arrays.
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _mols_dist;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _libs_dist;

    // we can impose limits on the size of fragments that are selected
    // for replacement.
    int _min_atoms_in_fragment;
    int _max_atoms_in_fragment;

    absl::flat_hash_set<IWString> _seen;

    MoleculeFilter _filter;

    uint64_t _new_molecules_formed;
    uint64_t _rejected_by_bad_valence;
    uint64_t _rejected_by_discard_queries;
    uint64_t _rejected_by_seen_before;
    uint64_t _rejected_by_filter;
    uint64_t _rejected_by_adjacent_atoms;

  // Private functions.
    int ReadLibrary(IWString& fname);
    int ReadMolecules(iwstring_data_source& input);
    int ReadMolecule(const const_IWSubstring& line);

    int IdentifyChangingFragments();
    int IdentifyChangingFragments(SafedMolecule& m);
    void IdentifyUnChangingFragments();
    void IdentifyUnChangingFragments(SafedMolecule& m);

    const SafeFragment* GetFragment(int ncon, int natoms);

    int Generate(SafedMolecule& m, Library& lib, IWString_and_File_Descriptor& output);
    int Generate(SafedMolecule& m,
                  const int f1_ndx, const SafeFragment& f2,
                  IWString_and_File_Descriptor& output);
    int AnyDiscardQueriesMatch(Molecule& m);
    int SeenBefore(Molecule& m);
    int BadBonds(Molecule& m) const;
    int OkAtomCount(const int natoms) const;
    int ProcessNewMolecule(Molecule& m, const IWString& name1,
                        const SafeFragment& f2,
                        IWString_and_File_Descriptor& output);

    int Breed(IWString_and_File_Descriptor& output);
    int Breed(SafedMolecule& m1, SafedMolecule& m2, IWString_and_File_Descriptor& output);
    int SelectFragments(const SafedMolecule& m1, const SafedMolecule& m2,
                int& f1, int& f2);

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

    int ReadMolecules(const char* fname);

    const resizable_array_p<SafedMolecule>& mols() const {
      return _mols;
    }

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;

    int SetupRng();

    int DoSubstructureSearches();

    int SetIsotope(isotope_t iso);

    int Generate(int ngenerate, IWString_and_File_Descriptor& output);
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _ignore_molecules_not_matching_queries = 0;

  _extra_atoms = 0;
  _fewer_atoms = 0;

  _max_attempts = 100;

  std::random_device rd;
  _rng.seed(rd());

  _min_atoms_in_fragment = 0;
  _max_atoms_in_fragment = std::numeric_limits<int>::max();

  _new_molecules_formed = 0;
  _rejected_by_bad_valence = 0;
  _rejected_by_discard_queries = 0;
  _rejected_by_seen_before = 0;
  _rejected_by_filter = 0;
  _rejected_by_adjacent_atoms = 0;
}

int
ReadQueries(Command_Line& cl, char flag, int verbose,
            resizable_array_p<Substructure_Query>& destination) {
  IWString q;
  for (int i = 0; cl.value(flag, q, i); ++i) {
    if (! process_cmdline_token(flag, q, destination, verbose)) {
      cerr << "ReadQueries:cannot process '" << q << "'\n";
      return 0;
    }
  }

  return destination.size();
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

  if (cl.option_present('C')) {
    if (! cl.value('C', max_ncon) || max_ncon < 1) {
      cerr << "Options::Initialise:invalid max ncon value (-C)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only consider fragments with at most " << max_ncon <<
              " connections\n";
    }
  }

  if (cl.option_present('Y')) {
    if (! ReadQueries(cl, 'Y', _verbose, _can_change)) {
      cerr << "Options::Initialise:cannot initialise -Y options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _can_change.size() << " -Y queries\n";
    }
  }

  if (cl.option_present('N')) {
    if (! ReadQueries(cl, 'N', _verbose, _cannot_change)) {
      cerr << "Options::Initialise:cannot initialise -N options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _cannot_change.size() << " -N queries\n";
    }
  }

  if (cl.option_present('X')) {
    if (! ReadQueries(cl, 'X', _verbose, _discard_if_match)) {
      cerr << "Options::Initialise:cannot initialise -X options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _discard_if_match.size() << " -X queries\n";
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_queries = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching any of the can match queries\n";
        }
      } else  {
      }
    }
  }

  if (cl.option_present('L')) {
    IWString fname;
    for (int i = 0; cl.value('L', fname, i); ++i) {
      cerr << "Reading '" << fname << "'\n";
      if (! ReadLibrary(fname)) {
        cerr << "Options::Initialise:cannot read library file '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      uint32_t nfrag = 0;
      for (const Library* lib : _library) {
        nfrag += lib->size();
      }
      cerr << "Read " << nfrag << " Library fragments from " <<
              _library.size() << " library files\n";
    }

    for (Library* lib : _library) {
      lib->SetupRng();
    }
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! _filter.Build(fname)) {
      cerr << "Options::Initialise:cannot initialise filter '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Filter initialised '" << fname << "'\n";
    }
  } 

  if (cl.option_present('x')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('x', x, i); ++i) {
      if (x.starts_with("extra=")) {
        x.remove_leading_chars(6);
        if (! x.numeric_value(_extra_atoms)) {
          cerr << "Invalid 'extra=" << x << " directive\n";
          return 0;
        }
      } else if (x .starts_with("fewer=")) {
        x.remove_leading_chars(6);
        if (! x.numeric_value(_fewer_atoms)) {
          cerr << "Invalid 'fewer=" << x << " directive\n";
          return 0;
        }
      } else {
        cerr << "Unrecognised -x directive '" << x << "'\n";
        return 0;
      }
    }
  }

  return 1;
}

int
Options::ReadLibrary(IWString& fname) {
  std::unique_ptr<Library> lib = std::make_unique<Library>();
  if (! lib->Build(fname)) {
    cerr << "Options::ReadLibrary:cannot read '" << fname << "'\n";
    return 0;
  }

  _library << lib.release();

  return _library.size();
}

int
Options::ReadMolecules(const char* fname) {
  if (_mols.empty()) {
    _mols.resize(100000);
  }

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadMolecules(input);
}

int
Options::ReadMolecules(iwstring_data_source& input) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    if (! ReadMolecule(line)) {
      cerr << "Options::ReadMolecule:cannot process '" << line << "'\n";
      return 0;
    }
  }

  return _mols.number_elements();
}

int
Options::ReadMolecule(const const_IWSubstring& line) {
  std::unique_ptr<SafedMolecule> f = std::make_unique<SafedMolecule>();
  if (f == nullptr) {
    cerr << "Options::ReadMolecule:memory failure\n";
    return 0;
  }
  if (! f->Build(line)) {
    cerr << "Options::ReadMolecule:cannot process '" << line << "'\n";
    return 0;
  }

  _mols << f.release();

  return _mols.size();
}

int
Options::SetupRng() {
  if (_mols.empty()) {
    cerr << "Options::SetupRng: no molecules\n";
    return 0;
  }

  _mols_dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _mols.size() - 1);
  _libs_dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _library.size() - 1);

  for (SafedMolecule* f : _mols) {
    f->SetupRng();
  }

  return 1;
}

int
Options::SetIsotope(isotope_t iso) {
  for (Library* lib : _library) {
    lib->SetIsotope(iso);
  }

  return 1;
}

int
Options::DoSubstructureSearches() {
  if (_can_change.size() > 0) {
    if (! IdentifyChangingFragments()) {
      return 0;
    }
  }

  if (_cannot_change.size() > 0) {
    IdentifyUnChangingFragments();
  }

  return 1;
}

int
Options::AnyDiscardQueriesMatch(Molecule& m) {
  if (_discard_if_match.empty()) {
    return 0;
  }

  Molecule_to_Match target(&m);
  for (Substructure_Query* q : _discard_if_match) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

int
Options::IdentifyChangingFragments() {
  for (SafedMolecule* m : _mols) {
    if (! IdentifyChangingFragments(*m)) {
      cerr << "Options::IdentifyChangingFragments:no matches to " << m->name() << '\n';
      return _ignore_molecules_not_matching_queries;
    }
  }

  return 1;
}

int
Options::IdentifyChangingFragments(SafedMolecule& m) {
  return m.IdentifyChangingFragments(_can_change);
}

void
Options::IdentifyUnChangingFragments() {
  for (SafedMolecule* m : _mols) {
    m->IdentifyUnChangingFragments(_cannot_change);
  }
}

#ifdef NOT_NEEDED_AASD
const SafeFragment*
Options::GetFragment(int ncon, int natoms) {
  if (_library.size() == 1) {
    return _library[0]->GetFragment(ncon, natoms);
  }

  int lib = (*_libs_dist)(_rng);

  return _library[lib]->GetFragment(ncon, natoms);
}
#endif

int
Options::Report(std::ostream& output) const {

  output << "Options: generated " << _new_molecules_formed << " molecules\n";
  output << _rejected_by_bad_valence << " _rejected_by_bad_valence\n";
  output << _rejected_by_discard_queries << " _rejected_by_discard_queries\n";
  output << _rejected_by_seen_before << " _rejected_by_seen_before\n";
  output << _rejected_by_filter << " _rejected_by_filter\n";
  output << _rejected_by_adjacent_atoms << " _rejected_by_adjacent_atoms\n";

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
Options::Generate(int ngenerate, IWString_and_File_Descriptor& output) {
  _max_attempts = 10 * ngenerate;

  int generated = 0;
  for (int i = 0; i < _max_attempts && generated < ngenerate; ++i) {
    int mindex = (*_mols_dist)(_rng);
    int libindex = (*_libs_dist)(_rng);
    // cerr << "indices " << mindex << ' ' << libindex << '\n';
    generated += Generate(*_mols[mindex], *_library[libindex], output);
  }

  return 1;
}

int
Options::Generate(SafedMolecule& m, Library& lib, IWString_and_File_Descriptor& output) {
  std::optional<int> f1_ndx = m.RandomFragment();
  if (f1_ndx == std::nullopt) {
    return 0;
  }

  const SafeFragment* f1 = m.fragment(*f1_ndx);
  // cerr << "f1_ndx " << *f1_ndx << " natoms " << f1->natoms() << " smiles " << f1->smiles() << '\n';
  // f1->DebugPrint(cerr);
  if (! OkAtomCount(f1->natoms())) {
    return 0;
  }

  int natoms;
  if (_extra_atoms == 0 && _fewer_atoms == 0) {
    natoms = f1->natoms();
  } else {
    // Make sure we have a valid minimum atom count.
    uint32_t amin;
    if (_fewer_atoms >= static_cast<uint32_t>(f1->natoms())) {
      amin = 1;
    } else {
      amin = f1->natoms() - _fewer_atoms;
    }
    std::uniform_int_distribution<uint32_t> u(amin, f1->natoms() + _extra_atoms);
    natoms = u(_rng);
    // cerr << "f1->natoms() " << f1->natoms() << " choose " << natoms << '\n';
  }

  const SafeFragment* f2 = lib.GetFragment(f1->ncon(), natoms);
  if (f2 == nullptr) {
    return 0;
  }
#ifdef DEBUG_GENERATE
  cerr << f1->natoms() << " requested " << natoms << " atoms got " << f2->natoms() << " QQ " << (natoms == f2->natoms()) << '\n';
  cerr << f1->smiles() << ' ' << m.name() << '\n';
  Molecule mcopy(const_cast<SafeFragment*>(f1)->mol());
  cerr << mcopy.smiles() << " from molecule, nat " << mcopy.natoms() << " cmp " << f1->natoms() << '\n';
  cerr << "F2 smiles " << f2->smiles() << '\n';
#endif

  return Generate(m, *f1_ndx, *f2, output);
}

int
FormNewSmiles(const IWString& starting_smiles, int f1_ndx,
              const IWString& new_smiles,
              IWString& destination) {
  destination.reserve(starting_smiles.size() + 20);   // 20 is an arbitrary choice
  int frag_number = 0;
  int new_smiles_added = 0;
  for (int i = 0; i < starting_smiles.number_elements(); ++i) {
    char c = starting_smiles[i];
    if (c == '.') {
      destination << c;
      ++frag_number;
    } else if (frag_number == f1_ndx) {
      if (! new_smiles_added) {
        destination << new_smiles;
        new_smiles_added = 1;
      }
    } else {
      destination << c;
    }
  }

  return 1;
}

// `f1` is a fragment number from within `m` and `f2` is a library fragment.
// Generate a molecule by replacing `f1` with `f2`.
int
Options::Generate(SafedMolecule& m,
                  const int f1_ndx, const SafeFragment& f2,
                  IWString_and_File_Descriptor& output) {
  const SafeFragment* f1 = m.fragment(f1_ndx);
  if (f1->ncon() != f2.ncon()) {
    cerr << "Options::Generate:ncon mismatch\n";
    return 0;
  }

  IWString new_smiles;
  if (! m.NewSmiles(f1_ndx, f2, new_smiles)) {
    return 0;
  }
  // cerr << "Fragment smiles " << f1->smiles() << '\n';
  // cerr << "from " << m.smiles() << " and " << f2.smiles() << " get " << new_smiles << '\n';

  IWString tmp;
  FormNewSmiles(m.smiles(), f1_ndx, new_smiles, tmp);

  Molecule newm;
  if (! newm.build_from_smiles(tmp)) {
    cerr << "Options::Generate:invalid smiles '" << tmp << "'\n";
    return 0;
  }

  output << m.smiles() << ' ' << m.name() << " parent\n";

  return ProcessNewMolecule(newm, m.name(), f2, output);
}

int
Options::ProcessNewMolecule(Molecule& m, const IWString& name1,
                        const SafeFragment& f2,
                        IWString_and_File_Descriptor& output) {
  output.write_if_buffer_holds_more_than(4096);

  ++_new_molecules_formed;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
  }

  if (! m.valence_ok()) {
    cerr << "Options::ProcessNewMolecule:invalid valence " << m.smiles() <<
            ' ' << name1 << '\n';
    ++_rejected_by_bad_valence;
    return 0;
  }

  if (AnyDiscardQueriesMatch(m)) {
    ++_rejected_by_discard_queries;
    return 0;
  }
  if (SeenBefore(m)) {
    ++_rejected_by_seen_before;
    return 0;
  }

  if (! _filter.Ok(m)) {
    ++_rejected_by_filter;
    return 0;
  }

  if (BadBonds(m)) {
    ++_rejected_by_adjacent_atoms;
    return 0;
  }

  // Form a new name.
  m << name1;

  m << " %% " << f2.name() << '.' << f2.ncon() << '.' << f2.natoms();

  output << m.aromatic_smiles() << ' ' << m.name() << '\n';

  return 1;
}

int
Options::OkAtomCount(const int natoms) const {
  if (natoms < _min_atoms_in_fragment) {
    return 0;
  }

  if (natoms > _max_atoms_in_fragment) {
    return 0;
  }

  return 1;
}

int
Options::BadBonds(Molecule& m) const {
  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    // Bonds are only formed with isotopes
    if (m.isotope(a1) == 0 || m.isotope(a2) == 0) {
      continue;
    }

    atomic_number_t z1 = m.atomic_number(a1);
    if (z1 == 6) {
      continue;
    }
    atomic_number_t z2 = m.atomic_number(a2);
    if (z2 == 6) {
      continue;
    }

    // Heteratom joined to heteroatom.
    return 1;
  }

  return 0;
}

int
Options::SeenBefore(Molecule& m) {
  auto f = _seen.find(m.unique_smiles());
  if (f != _seen.end()) {
    return 1;
  }

  _seen.insert(m.unique_smiles());
  return 0;
}

int
Options::Breed(IWString_and_File_Descriptor& output) {
  int m1 = (*_mols_dist)(_rng);
  int m2 = 0;
  while (1) {
    m2 = (*_mols_dist)(_rng);
    if (m1 == m2) {
      break;
    }
  }

  return Breed(*_mols[m1], *_mols[m2], output);
}

// Identify a fragment number, `f1` from `m1` that is about the
// same size as a fragment `f2` from `m2`.
int
Options::SelectFragments(const SafedMolecule& m1, const SafedMolecule& m2,
                int& f1, int& f2) {
  int nf1 = m1.number_fragments();
  int nf2 = m2.number_fragments();

  if (nf1 == 1) {
    f1 = 0;
  } else {
    std::uniform_int_distribution<uint32_t> u1(0, nf1 - 1);
    f1 = u1(_rng);
  }

  int atoms_in_f1 = m1.fragment(f1)->natoms();

  // Find the fragment in `m2` closest to atoms_in_f1.
  f2 = -1;
  int min_diff = std::numeric_limits<int>::max();
  for (int i = 0; i < nf2; ++i) {
    const SafeFragment* f = m2.fragment(i);
    int d = std::abs(f->natoms() - atoms_in_f1);
    if (d < min_diff) {
      min_diff = d;
      f2 = i;
    }
  }

  return f2 >= 0;
}

int
Options::Breed(SafedMolecule& m1, SafedMolecule& m2,
               IWString_and_File_Descriptor& output) {
  int f1_ndx, f2_ndx;
  if (! SelectFragments(m1, m2, f1_ndx, f2_ndx)) {
    return 0;
  }

  const SafeFragment* f1 = m1.fragment(f1_ndx);
  const SafeFragment* f2 = m2.fragment(f2_ndx);

  IWString replacement_smiles;
  f1->SameNumbers(*f2, replacement_smiles);

  IWString tmp;
  FormNewSmiles(m1.smiles(), f1_ndx, replacement_smiles, tmp);

  Molecule m;
  if (! m.build_from_smiles(tmp)) {
    cerr << "Options::Breed:invalid smiles '" << tmp << "'\n";
    return 0;
  }

  return ProcessNewMolecule(m, m1.name(), *f2, output);
}

int
SafeGenerate(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:L:Y:N:C:z:X:F:x:n:");

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

  for (const char * fname : cl) {
    if (verbose) {
      cerr << "Reading molecules from '" << fname << "'\n";
    }
    if (! options.ReadMolecules(fname)) {
      cerr << "SafeGenerate::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << options.mols().size() << " molecules\n";
  }

  options.SetupRng();

  options.DoSubstructureSearches();

  options.SetIsotope(2);

  IWString_and_File_Descriptor output(1);

  int ngenerate = 1;
  if (cl.option_present('n')) {
    if (! cl.value('n', ngenerate) || ngenerate < 1) {
      cerr << "Invalid number to generate (-n)\n";
      return 1;
    }
    if (verbose) {
      cerr << "Will generate " << ngenerate << " variants\n";
    }
  }
  options.Generate(ngenerate, output);

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace safe_generate

int
main(int argc, char ** argv) {

  int rc = safe_generate::SafeGenerate(argc, argv);

  return rc;
}
