// take a set of diced molecules and generate subsets that
// exemplify all the diced fragments.

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "absl/container/flat_hash_map.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace parsimonious_set {

using std::cerr;

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
Generates a parsimonious subset of a set of molecules that are greedily selected to examplify as many
of the fragments in the set of molecules as possible.
This tool was developed when there is a large set of molecules for which we would like to perform
an expensive computation - frequently 3D. But there are too many molecules to run the method on all
candidates. This tool takes a fragment representation of the set of molecules and generates a
subset that exemplifies as many of the fragments as possible, given the number of molecules
to be selected.

The input is binary dicer_data::DicedMolecule serialised protos formed by dicer. For example 

dicer -I 1 -X 64 -m 3 -M 20 -M maxnr=10 -v -B nbamide -B brcb -B serialized_proto -S dicer.data -k 3 -c file.smi

The nature of the dicing is critically important, and there is no "right" answer for what
is best. Fragments should definitely be labelled, either with a constant isotope or with an atom type.

The following options are recognised

-support <n>    discard any fragment found in fewer than <n> molecules.
-nsel <n>       the number of molecules to test.
-randf          randomise the fragments before doing the selection. Makes the tool non-deterministic.
-sortf          within groups of same atom count fragments, sort by prevalence
-nosort         do NOT sort the input molecules by size. Do this is they are already sorted by activity.
-nexample <n>   for each fragment, the number of molecules exemplifying that fragment - might be impossible.
-v              verbose output
)";
  // clang-format on
  ::exit(rc);
}

class Fragment {
 private:
  // Is this fragment under consideration or not.
  // By default yes.
  int _active = 1;

  // Both of these are copied from the dicer_data::DicerFragment proto
  // from which this is constructed.
  std::string _smiles;
  uint32_t _natoms = 0;

  // Each fragment is assigned a unique id.
  uint32_t _uid = 0;

  // the number of molecules containing this feature - before selection.
  uint32_t _number_instances = 0;

  // Whether or not this fragment has been selected.
  bool _selected = false;

 public:
  Fragment(const dicer_data::DicerFragment& proto);

  int
  active() const {
    return _active;
  }

  void
  deactivate() {
    _active = 0;
  }

  uint32_t
  uid() const {
    return _uid;
  }

  void
  set_uid(uint32_t s) {
    _uid = s;
  }

  const std::string&
  smiles() const {
    return _smiles;
  }

  uint32_t
  number_instances() const {
    return _number_instances;
  }

  void
  set_number_instances(uint32_t s) {
    _number_instances = s;
  }

  void
  increment() {
    ++_number_instances;
  }

  bool
  selected() const {
    return _selected;
  }

  void
  set_selected(bool s) {
    _selected = s;
  }

  uint32_t
  natoms() const {
    return _natoms;
  }
};

Fragment::Fragment(const dicer_data::DicerFragment& proto) {
  _smiles = proto.smi();
  _natoms = proto.nat();
  _number_instances = 1;
}

using SmilesToFragment = absl::flat_hash_map<std::string, Fragment*>;

class Candidate {
 private:
  // Attributes copied from the proto.
  std::string _smiles;
  std::string _name;
  uint32_t _natoms;

  Molecule _mol;

  bool _selected = false;

  // array of fragments that are in this molecule. Determined during
  // construction.
  std::vector<Fragment*> _fragment;

  // private functions

  int BuildFragment(const dicer_data::DicerFragment& proto,
                    SmilesToFragment& fragment_to_uid);

 public:
  // Build from `proto` while updating the global has of smiles to Fragment*
  int Build(const dicer_data::DicedMolecule& proto, SmilesToFragment& fragment_to_uid);

  int Initialise(Command_Line_v2& cl);

  bool
  selected() const {
    return _selected;
  }

  void
  set_selected(bool s) {
    _selected = s;
  }

  int
  natoms() const {
    return _natoms;
  }

  const std::string&
  smiles() const {
    return _smiles;
  }

  const std::string&
  name() const {
    return _name;
  }

  uint32_t
  number_fragments() const {
    return _fragment.size();
  }

  bool
  ContainsFragment(const Fragment* f) const {
    return absl::c_find(_fragment, f) != _fragment.end();
  }

  // Update the global array of where the first candidate with each fragment
  // is to be found. 
  // We are candidate number `ndx`. For each of our fragments, see if their
  // uid is smaller than what is in `first_candidate_with_fragment` and update
  // if smaller.
  void UpdateFirstCandidate(uint32_t ndx,
                            uint32_t* first_candidate_with_fragment) const;

  // This candidate molecule has been selected.
  // Update _selected and increment `selected`
  // with all our fragments.
  // Return the number of newly selected bits.
  int IAmSelected(uint32_t* selected);
};

int
Candidate::Build(const dicer_data::DicedMolecule& proto,
                 SmilesToFragment& fragment_to_uid) {
  _smiles = proto.smiles();
  _name = proto.name();
  _natoms = proto.natoms();

  if (!_mol.build_from_smiles(_smiles)) {
    cerr << "Candidate::Build:bad smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _fragment.reserve(proto.fragment_size());

  for (const auto& fragment : proto.fragment()) {
    if (!BuildFragment(fragment, fragment_to_uid)) {
      cerr << "Candidate::Build:error processing " << fragment.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Candidate::IAmSelected(uint32_t* selected) {
  int rc = 0;

  for (const Fragment* f : _fragment) {
    if (!f->active()) {
      // cerr << "Inactive fragment\n";
      continue;
    }

    uint32_t uid = f->uid();
    if (selected[uid] == 0) {
      selected[uid] = 1;
      ++rc;
    } else {
      ++selected[uid];
    }
    // cerr << "  incremented " << uid << " " << selected[uid] << '\n';
  }

  _selected = true;

  return rc;
}

void
Candidate::UpdateFirstCandidate(uint32_t ndx,
                            uint32_t* first_candidate_with_fragment) const {
  for (const Fragment* f : _fragment) {
    if (! f->active()) {
      continue;
    }
    uint32_t uid = f->uid();
    if (first_candidate_with_fragment[uid] > ndx) {
      first_candidate_with_fragment[uid] = ndx;
    }
  }
}

// Add fragment `proto` to `this`.
// Look up the smiles in `fragment_to_uid`. If already present, update the existing
// Fragment information.
// If this is a new fragment, create a new one and add to `fragment_to_uid`.
// In both cases, add the Fragment ptr to _fragment.
int
Candidate::BuildFragment(const dicer_data::DicerFragment& proto,
                         SmilesToFragment& fragment_to_uid) {
  auto iter = fragment_to_uid.find(proto.smi());

  if (iter != fragment_to_uid.end()) {
    Fragment* f = iter->second;
    f->increment();
    _fragment.push_back(f);
    return 1;
  }

  Fragment* f = new Fragment(proto);
  uint32_t s = fragment_to_uid.size();
  f->set_uid(s);

  fragment_to_uid[proto.smi()] = f;
  _fragment.push_back(f);

  return 1;
}

class Candidates {
 private:
  resizable_array_p<Candidate> _candidates;

  // Extract from each DicedMolecule the list of fragments.
  resizable_array_p<Fragment> _fragment;

  SmilesToFragment _fragment_to_uid;

  // Initialised from command line

  int _verbose = 0;

  // If set, we can discard fragments set in less than this many molecules.
  uint32_t _support = 0;

  // Whether or not to sort the candidate molecules by atom count.
  int _sort_candidates_by_atom_count = 1;

  // The number if items to select.
  uint32_t _nsel = 0;

  // The number of examples of each fragment we would like to get.
  // Note that this may be impossible since there will be fragments
  // that occur in fewer molecules than this number.
  uint32_t _examples_needed = 1;

  // We can make the algorithm non-deterministic by sorting the fragments
  // within bands of fragments with common atom counts.
  int _randomise_fragments = 0;

  // Within Fragments with the same atom count, sort the fragments
  // by prevalence.
  int _sort_fragments_by_prevalence = 0;

  // We can speed things up (a lot) by precomputing
  // where in the array of Fragment's the first one that exemplifies
  // a particular fragment is.
  uint32_t* _first_candidate_with_fragment = nullptr;

  Report_Progress _report_progress;

  Fraction_as_String _fraction_as_string;

  // private functions

  int ReadFragmentionData(iw_tf_data_record::TFDataReader& input);
  int FirstMoleculeWithFragment(const Fragment* f) const;
  int NotifyFragmentSelected(const Fragment* f, uint32_t* selected);
  void SortCandidates();
  void RandomiseFragments(std::vector<Fragment*>& fragments);
  void SortFragmentsByPrevalence(std::vector<Fragment*>& fragments);
  int WriteSelection(const Candidate* candidate, uint32_t features_selected, uint32_t nfrag,
             IWString_and_File_Descriptor& output) const;

 public:
  Candidates();
  ~Candidates();

  int Initialise(Command_Line_v2& cl);

  int ReadFragmentionData(const char* fname);

  int ReportSize(std::ostream& output) const;

  int ReportSelections(std::ostream& output) const;

  uint32_t
  number_candidates() const {
    return _candidates.size();
  }

  int Process(IWString_and_File_Descriptor& output);
};

Candidates::Candidates() {
  _fraction_as_string.set_include_leading_space(1);
  _fraction_as_string.initialise(0.0f, 1.0f, 3);
  _fraction_as_string.append_to_each_stored_string("\n");
}

Candidates::~Candidates() {
  for (auto& [_, value] : _fragment_to_uid) {
    delete value;
  }

  if (_first_candidate_with_fragment != nullptr) {
    delete [] _first_candidate_with_fragment;
  }
}

int
Candidates::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present("v");

  if (cl.option_present("support")) {
    cl.value("support", _support);
    if (_verbose) {
      cerr << "Will ignore fragments with fewer than " << _support << " instances\n";
    }
  }

  if (!cl.option_present("nsel")) {
    cerr << "Must specify the number of items to select with the -nsel option\n";
    return 0;
  }

  cl.value("nsel", _nsel);
  if (_verbose) {
    cerr << "Will select " << _nsel << " items\n";
  }

  if (cl.option_present("randf")) {
    _randomise_fragments = true;
    if (_verbose) {
      cerr << "Will randomise fragments within atom count bands\n";
    }
  }

  if (cl.option_present("sortf")) {
    _sort_fragments_by_prevalence = 1;
    if (_verbose) {
      cerr << "Will sort same atom count Fragment groups by prevalence\n";
    }
  }

  if (cl.option_present("nosort")) {
    _sort_candidates_by_atom_count = 0;
    if (_verbose) {
      cerr << "Will NOT sort candidate molecules by atom count\n";
    }
  }

  if (cl.option_present("nexample")) {
    cl.value("nexample", _examples_needed);
    if (_verbose) {
      cerr << "WIll attempt to get " << _examples_needed << " examples of each fragment\n";
    }
  }

  if (cl.option_present("report")) {
    uint64_t r;
    cl.value("report", r);
    _report_progress.set_report_every(r);
    if (_verbose) {
      cerr << "Will report progress every " << r << " items selected\n";
    }
  }

  return 1;
}

int
Candidates::ReportSize(std::ostream& output) const {
  output << "Contains " << _candidates.size() << " molecules and "
         << _fragment_to_uid.size() << " fragments\n";

  extending_resizable_array<uint32_t> count;

  Accumulator_Int<uint32_t> acc_natoms;

  for (const auto& [_, value] : _fragment_to_uid) {
    ++count[value->number_instances()];
    acc_natoms.extra(value->natoms());
  }

  for (uint32_t i = 0; i < count.size(); ++i) {
    if (count[i]) {
      output << count[i] << " fragments had " << i << " examples "
             << iwmisc::Fraction<float>(i, _candidates.size()) << '\n';
    }
  }

  output << "Atom counts btw " << acc_natoms.minval() << " and " << acc_natoms.maxval()
         << " ave " << acc_natoms.average() << '\n';

  return output.good();
}

int
Candidates::ReportSelections(std::ostream& output) const {
  return output.good();
}

int
Candidates::ReadFragmentionData(const char* fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (!input.good()) {
    cerr << "Candidates::ReadFragmentionData:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentionData(input);
}

int
Candidates::ReadFragmentionData(iw_tf_data_record::TFDataReader& input) {
  while (1) {
    std::optional<dicer_data::DicedMolecule> proto =
        input.ReadProto<dicer_data::DicedMolecule>();
    if (!proto) {
      return 1;
    }

    std::unique_ptr<Candidate> candidate = std::make_unique<Candidate>();
    if (!candidate->Build(*proto, _fragment_to_uid)) {
      cerr << "Candidates::ReadFragmentionData:error processint "
           << proto->ShortDebugString() << '\n';
      return 0;
    }

    _candidates << candidate.release();
  }

  return _candidates.size();
}

// Return the fraction of items in `value` that are non-zero.
float
FractionSelected(const uint32_t* value, uint32_t n) {
  uint32_t rc = 0;

  for (uint32_t i = 0; i < n; ++i) {
    if (value[i] > 0) {
      ++rc;
    }
  }

  return iwmisc::Fraction<float>(rc, n);
}

int
Candidates::Process(IWString_and_File_Descriptor& output) {
  if (_nsel >= _candidates.size()) {
    cerr << "Candidates::Process:cannot select " << _nsel << " candidates from "
         << _candidates.size() << '\n';
    return 0;
  }

  if (_sort_candidates_by_atom_count) {
    SortCandidates();
  }

  // The number of fragments being processed may be lower after we impose
  // any support requirement.
  uint32_t nfrag = _fragment_to_uid.size();

  std::vector<Fragment*> fragments;
  fragments.reserve(nfrag);

  if (_support == 0) {
    for (const auto& [_, frag] : _fragment_to_uid) {
      fragments.push_back(frag);
    }
  } else {
    // If we have a support requirement, we need to renumber the fragments.
    uint32_t new_uid = 0;

    for (const auto& [_, frag] : _fragment_to_uid) {
      if (frag->number_instances() < _support) {
        frag->deactivate();
        continue;
      } 

      frag->set_uid(new_uid);
      ++new_uid;

      fragments.push_back(frag);
    }
  }

  nfrag = fragments.size();
  if (nfrag == 0) {
    cerr << "Candidates::Process:no fragments meet support level " << _support << '\n';
    return 0;
  }

  if (_verbose) {
    cerr << nfrag << " fragments\n";
  }

  // Sort by fragment size, largest first.
  absl::c_sort(fragments, [](const Fragment* f1, const Fragment* f2) {
    return f1->natoms() > f2->natoms();
  });

  if (_randomise_fragments) {
    RandomiseFragments(fragments);
  } else if (_sort_fragments_by_prevalence) {
    SortFragmentsByPrevalence(fragments);
  }

  _first_candidate_with_fragment = new uint32_t[nfrag];
  std::fill_n(_first_candidate_with_fragment, nfrag, nfrag);
  for (uint32_t i = 0; i < _candidates.size(); ++i) {
    _candidates[i]->UpdateFirstCandidate(i, _first_candidate_with_fragment);
  }

  // How many times has a fragment been selected.
  std::unique_ptr<uint32_t[]> selected = std::make_unique<uint32_t[]>(nfrag);
  std::fill_n(selected.get(), nfrag, 0);

  // The number selected.
  uint32_t nsel = 0;

  for (uint32_t i = 0; i < nfrag; ++i) {
    int j = FirstMoleculeWithFragment(fragments[i]);
    if (j < 0) {
      cerr << i << " huh, no molecules for " << fragments[i]->uid() << '\n';
    }
  }

  uint32_t features_selected = 0;
  // For each fragment, select the first molecule that contains that fragment.
  for (uint32_t i = 0; i < nfrag && nsel < _nsel; ++i) {
    Fragment* f = fragments[i];
    // If this fragment has already been selected - as a result of a previous
    // molecule being selected, skip.
    uint32_t uid = f->uid();
    if (selected[uid] >= _examples_needed) {
      continue;
    }
    int j = FirstMoleculeWithFragment(f);
    if (j < 0) {
      continue;
    }

    features_selected += _candidates[j]->IAmSelected(selected.get());

    ++nsel;

    WriteSelection(_candidates[j], features_selected, nfrag, output);
    
    if (_report_progress()) {
      cerr << "Selected " << nsel << " of " << _nsel << " candidates\n";
    }
  }

  // The loop above can leave some bits not selected.
  // A molecule that is skipped may have bits set that never get queried.

  output.flush();

  if (_verbose) {
    cerr << "Selected " << nsel << " of " << _candidates.size() << " candidates\n";
    extending_resizable_array<uint32_t> sel;

    uint32_t fragments_not_represented = 0;
    // Accumulate the sizes of the fragments not represented.
    Accumulator_Int<uint32_t> acc_natoms;
    Accumulator_Int<uint32_t> acc_instances;

    absl::c_sort(fragments, [&selected](const Fragment* f1, const Fragment* f2) {
      return selected[f1->uid()] < selected[f2->uid()];
    });

    for (uint32_t i = 0; i < nfrag; ++i) {
      const Fragment* f = fragments[i];
      if (!f->active()) {
        continue;
      }
      uint32_t uid = f->uid();
      cerr << f->smiles() << " nat " << f->natoms() << " uid " << uid << " instances "
           << f->number_instances() << " sel " << selected[uid] << '\n';
      ++sel[selected[uid]];
      if (selected[uid] == 0) {
        ++fragments_not_represented;
        acc_natoms.extra(f->natoms());
        acc_instances.extra(f->number_instances());
      }
    }

    for (int i = 0; i < sel.number_elements(); ++i) {
      if (sel[i]) {
        cerr << sel[i] << " fragments represented " << i << " times\n";
      }
    }
    if (fragments_not_represented > 0) {
      cerr << fragments_not_represented << " of " << nfrag
           << " fragments not represented in " << nsel << " selections "
           << iwmisc::Fraction<float>(fragments_not_represented, nfrag) << '\n';
      cerr << " atoms btw " << acc_natoms.minval() << " and " << acc_natoms.maxval()
           << " ave " << acc_natoms.average() << '\n';
      cerr << " instances btw " << acc_instances.minval() << " and " <<
              acc_instances.maxval() << " ave " << acc_instances.average() << '\n';
    }
  }

  return 1;
}

int
Candidates::WriteSelection(const Candidate* candidate, uint32_t features_selected, uint32_t nfrag,
             IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';

  output << candidate->smiles() << kSep << candidate->name();

  const float f = iwmisc::Fraction<float>(features_selected, nfrag);
  _fraction_as_string.append_number(output, f);

  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

void
Candidates::SortCandidates() {
  // Sort by atom count
  _candidates.iwqsort_lambda([](const Candidate* c1, const Candidate* c2) {
    int m1 = c1->natoms();
    int m2 = c2->natoms();
    if (m1 < m2) {
      return -1;
    }
    if (m1 > m2) {
      return 1;
    }
    return 0;
  });
}

// Randomise within groups of fragments with the same atom count.
// Assumes that `fragments` has been sorted by atom count.
void
Candidates::RandomiseFragments(std::vector<Fragment*>& fragments) {
  uint32_t istart = 0;
  uint32_t natoms = fragments[0]->natoms();

  std::random_device rd;
  std::default_random_engine rng(rd());

  for (uint32_t i = 1; i < fragments.size(); ++i) {
    if (fragments[i]->natoms() == natoms) {
      continue;
    }

    std::shuffle(fragments.begin() + istart, fragments.begin() + i - 1, rng);
    istart = i;
    natoms = fragments[i]->natoms();
  }

  std::shuffle(fragments.begin() + istart, fragments.end(), rng);
}

// Within groups of fragments with the same atom count, sort the fragments
// by prevalence.
void
Candidates::SortFragmentsByPrevalence(std::vector<Fragment*>& fragments) {
  uint32_t istart = 0;
  uint32_t natoms = fragments[0]->natoms();

  std::random_device rd;
  std::default_random_engine rng(rd());

  for (uint32_t i = 1; i < fragments.size(); ++i) {
    if (fragments[i]->natoms() == natoms) {
      continue;
    }

    std::sort(fragments.begin() + istart, fragments.begin() + i - 1,
                [](const Fragment* f1, const Fragment* f2) {
                  return f1->number_instances() > f2->number_instances();
                });
    istart = i;
    natoms = fragments[i]->natoms();
  }

  std::sort(fragments.begin() + istart, fragments.end(), 
                [](const Fragment* f1, const Fragment* f2) {
                  return f1->number_instances() > f2->number_instances();
                });
}

int
Candidates::FirstMoleculeWithFragment(const Fragment* f) const {
  const uint32_t n = _candidates.size();
  for (uint32_t i = _first_candidate_with_fragment[f->uid()]; i < n; ++i) {
    const Candidate* c = _candidates[i];
    if (c->selected()) {
      continue;
    }

    if (_candidates[i]->ContainsFragment(f)) {
      return i;
    }
  }

  return -1;
}

int
ParsimoniousSet(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-S=s-support=ipos-nsel=ipos-nosort-randf-sortf-report=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Unsufficient arguments\n";
    Usage(1);
  }

  Candidates candidates;

  if (!candidates.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  for (const char* fname : cl) {
    if (!candidates.ReadFragmentionData(fname)) {
      cerr << "Cannot read '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    candidates.ReportSize(cerr);
  }

  IWString_and_File_Descriptor output(1);

  candidates.Process(output);

  if (verbose) {
    candidates.ReportSelections(cerr);
  }

  return 0;
}

}  // namespace parsimonious_set

int
main(int argc, char** argv) {
  int rc = parsimonious_set::ParsimoniousSet(argc, argv);

  return rc;
}
