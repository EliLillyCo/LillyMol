// We have a list of de-functionalised reagents.
// Given a set of molecule, for each reagent, perform a substructure seach
// on that molecule and identify matched atoms.

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <ranges>

#include "absl/container/flat_hash_set.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwrecordio.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/reagent_substructure_search.pb.h"
#else
#include "reagent_substructure_search.pb.h"
#endif

namespace reagent_substructure_search {

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
  cerr << R"(Reagent substructure search.
Consumed
 -R <fname>        recordio file of reagent_substructure_search::ReagentData serialized protos,
                      such as might come from Enamine2Reagents.
 -X <matches>      max number of substructure matches per molecule - for efficiency, but
                        may leave some atoms unmatched.
 -m                minimum number of atoms in reagents.
 -M                maximum number of atoms in reagents.
 -d <ncon>         allowed values for number of connections.
 -O <fname>        write labelled molecules to <fname>.
 -c                remove chirality.
 -v                verbose output.
)";
// clang-format on

  ::exit(rc);
}

#ifdef NOT_NEEDED_JJAA
class MoleculeAndFormula {
  private:
    Molecule* _m = nullptr;
    mformula::MFormula _formula;

    Molecule_to_Match _target;

  public:
    int Initialise(const std::string& smiles);

    const Molecule& mol() const {
      return *_m;
    }
    Molecule& mol() {
      return *_m;
    }

    const mformula::MFormula& formula() const {
      return _formula;
    }
    mformula::MFormula& formula() {
      return _formula;
    }

};

int
MoleculeAndFormula::Initialise(const std::string& smiles) {
  if (! _m->build_from_smiles(smiles)) {
    cerr << "MoleculeAndFormula::Initialise:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  _formula.Build(*_m);

  _target.initialise_molecule(_m);
}
#endif

// As we identify matched atoms in a target molecule, we need to keep
// track of which atoms are hit by the queries. For each of the _natoms
// atoms in the molecule, keep track of which Reagent matches.
class Coverage {
  private:
    const int _natoms;

    absl::flat_hash_set<uint32_t>* _by_atom;

    int _atoms_covered;

  public:
    Coverage(int n);
    ~Coverage();

    // An indication that reagent number `reagent` has matched `atom`.
    void Extra(int atom, int reagent);

    void Update(const Substructure_Results& sresults, int ndx);

    float FractionAtomsCovered() const;

    uint32_t Count(int atom) const {
      return _by_atom[atom].size();
    }

    int atoms_covered() const {
      return _atoms_covered;
    }

    bool all_atoms_covered() const {
      return _atoms_covered == _natoms;
    }
};

Coverage::Coverage(int n) : _natoms(n) {
  _by_atom = new absl::flat_hash_set<uint32_t>[_natoms];
  _atoms_covered = 0;
}

Coverage::~Coverage() {
  delete [] _by_atom;
}

void
Coverage::Extra(int atom, int reagent) {
  if (_by_atom[atom].empty()) {
    ++_atoms_covered;
  }

  _by_atom[atom].insert(reagent);
}

float
Coverage::FractionAtomsCovered() const {
  int atoms_matched = 0;
  for (int i = 0; i < _natoms; ++i) {
    if (! _by_atom[i].empty()) {
      ++atoms_matched;
    }
  }

  return iwmisc::Fraction<float>(atoms_matched, _natoms);
}

void
Coverage::Update(const Substructure_Results& sresults, int ndx) {
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    for (atom_number_t a : *e) {
      Extra(a, ndx);
    }
  }
}

class Reagent {
  private:
    ReagentData _proto;

    // Derived form the smiles, will always be present.
    mformula::MFormula _formula_smiles;

    // will come into existence if we need the aromatic formula.
    std::unique_ptr<Molecule> _m;

    std::unique_ptr<mformula::MFormula> _formula;

    // If we decide we need to do searches, this will come into existence.
    std::unique_ptr<Substructure_Query> _query;

    // Our position in some external array.
    uint32_t _ndx;

  // private functions.

  public:
    Reagent();

    // Build from `proto`.
    // Return 1 if successful.
    // Return 0 if the constraints on min/max atoms are not met.
    // Return -1 if there is an error.
    // Should do an enum, but too lazy.
    // `proto` is destroyed.
    int BuildFromProto(ReagentData& proto,
                uint32_t min_natoms,
                uint32_t max_natoms,
                const bool* _ok_ncon);

    void set_ndx(int s) {
      _ndx = s;
    }

    uint32_t ndx() const {
      return _ndx;
    }

    bool has_formula() const {
      return _formula != nullptr;
    }
    bool has_query() const {
      return _query != nullptr;
    }

    int natoms() const {
      return _proto.natoms();
    }

    int SubstructureSearch(Molecule_to_Match& target,
                           const mformula::MFormula& rhs_mf,
                           const mformula::MFormula& rhs_mf_arom,
                           Substructure_Results& sresults);
};

Reagent::Reagent() {
  _ndx = 0;
}

int
Reagent::BuildFromProto(ReagentData& proto,
                uint32_t min_natoms,
                uint32_t max_natoms,
                const bool* _ok_ncon) {
  if (proto.natoms() < min_natoms) {
    return 0;
  }
  if (proto.natoms() > max_natoms) {
    return 0;
  }

  if (! _ok_ncon[proto.ncon()]) {
    return 0;
  }

  if (! _formula_smiles.Build(proto.smiles())) {
    cerr << "Reagent::BuildFromProto:invalid smiles '" << proto.smiles() << "'\n";
    return -1;
  }

  _proto = std::move(proto);

  return 1;
}

int
Reagent::SubstructureSearch(Molecule_to_Match& target,
                const mformula::MFormula& rhs_mf,
                const mformula::MFormula& rhs_mf_arom,
                Substructure_Results& sresults) {
  // cerr << _proto.smiles() << " SSS\n";

  if (static_cast<int32_t>(_proto.natoms()) > target.natoms()) {
    // cerr << "TOo many atoms in reagent\n";
    return 0;
  }

  if (! _formula_smiles.IsSubset(rhs_mf)) {
    // cerr << "_formula_smiles no match\n";
    return 0;
  }

  if (! _formula) {
    _m = std::make_unique<Molecule>(_proto.natoms());
    if (! _m->build_from_smiles(_proto.smiles())) {
      cerr << "Reagent::SubstructureSearch:invalid smiles '" << _proto.smiles() << '\n';
      return 0;
    }

    _formula = std::make_unique<mformula::MFormula>();
    _formula->Build(*_m);
  }

  if (! _formula->IsElementCountSubset(rhs_mf_arom)) {
    // cerr << "Formula from molecule no match\n";
    return 0;
  }

  if (! _query) {
    static Molecule_to_Query_Specifications mqs;
    mqs.set_substituents_only_at_isotopic_atoms(1);
    mqs.set_must_have_substituent_at_every_isotopic_atom(1);
    mqs.set_atoms_conserve_ring_membership(1);
    mqs.set_non_ring_atoms_become_nrings_0(1);

    _query = std::make_unique<Substructure_Query>();
    if (! _query->create_from_molecule(*_m, mqs)) {
      cerr << "Reagent::SubstructureSearch:cannot create query\n";
      return 0;
    }
  }

  return _query->substructure_search(target, sresults);
}

constexpr int kMaxNcon = 3;

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint32_t _min_natoms;
    uint32_t _max_natoms;

    // We only read Reagents that have their ncon() value set in
    // this array.
    bool _ok_ncon[kMaxNcon + 1];

    // For efficiency, we can impose a maximum number of substructure matches.
    int _max_matches;

    // Read from the -R file.
    resizable_array_p<Reagent> _reagent;

    uint64_t _molecules_read = 0;

    extending_resizable_array<uint32_t> _matches;

    Accumulator<double> _acc_coverage;

    IWString_and_File_Descriptor _stream_for_labelled_molecules;

    Report_Progress _report_progress;

    // private functions

    int ReadReagents(IWString& fname);
    int ReadReagents(iwrecordio::IWRecordIoReader& reader);
    int AddReagent(reagent_substructure_search::ReagentData& proto);
    int ReportReagentSizes(std::ostream& output) const;
    void SortReagentsByAtomCount();

    int HandleNoMatches(Molecule& m);

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

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _min_natoms = 0;
  _max_natoms = std::numeric_limits<uint32_t>::max();

  _max_matches = std::numeric_limits<int>::max();

  std::fill_n(_ok_ncon, 4, true);
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

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_natoms)) {
      cerr << "Options::Initialise:invalid minimum atom count -m\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only consider reagents with at least " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_natoms)) {
      cerr << "Options::Initialise:invalid maximum atom count -M\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only consider reagents with at most " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('d')) {
    std::fill_n(_ok_ncon, kMaxNcon + 1, true);
    uint32_t r;
    for (int i = 0; cl.value('d', r, i); ++i) {
      if (r >= kMaxNcon) {
        cerr << "Invalid number of connections " << r << '\n';
        return 0;
      }
    }
    _ok_ncon[r] = true;
  }

  if (! cl.option_present('R')) {
    cerr << "Must specify one or more IWRecordIoReader binary files of\n";
    cerr << "reagent_substructure_search::ReagentData protos via the -R option\n";
    Usage(1);
  } else {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (! ReadReagents(fname)) {
        cerr << "Options::Initialise:cannot read reagents '" << fname << "'\n";
        return 0;
      }
    }

    SortReagentsByAtomCount();
  }

  if (_verbose) {
    cerr << "Read " << _reagent.size() << " reagents\n";
    ReportReagentSizes(cerr);
  }

  if (cl.option_present('X')) {
    if (! cl.value('X', _max_matches) || _max_matches < 1) {
      cerr << "Invalid maximum number of substructure matches (-X)\n";
      Usage(1);
    }
  }

  if (cl.option_present('O')) {
    IWString fname = cl.string_value('O');
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_labelled_molecules.open(fname)) {
      cerr << "Options::Initialise:cannot open '" << fname << "' -O\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Labelled molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Options::Initialise:cannot initialise progress reporting -r\n";
      return 0;
    }
  }

  return 1;
}

void
Options::SortReagentsByAtomCount() {
  // Sort so that the largest reagents are first.
  _reagent.iwqsort_lambda([](const Reagent* r1, const Reagent* r2) {
      if (r1->natoms() > r2->natoms()) {
        return -1;
      } else if (r1->natoms() < r2->natoms()) {
        return 1;
      } else {
        return 0;
      }
    });

  for (uint32_t i = 0; i < _reagent.size(); ++i) {
    _reagent[i]->set_ndx(i);
  }
}

int
Options::ReportReagentSizes(std::ostream& output) const {
  extending_resizable_array<int> reagent_size;
  for (const Reagent* r : _reagent) {
    ++reagent_size[r->natoms()];
  }

  Accumulator<uint32_t> acc_natoms;
  for (int i = 1; i < reagent_size.number_elements(); ++i) {
    if (reagent_size[i]) {
      output << reagent_size[i] << " reagents had " << i << " atoms\n";
      acc_natoms.extra(i, reagent_size[i]);
    }
  }

  output << "Ave " << static_cast<float>(acc_natoms.average()) << '\n';

  return output.good();
}

int
Options::ReadReagents(IWString& fname) {
  iwrecordio::IWRecordIoReader reader(fname);
  if (! reader.good()) {
    cerr << "Options::ReadReagents:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReagents(reader);
}

int
Options::ReadReagents(iwrecordio::IWRecordIoReader& reader) {
  while (1) {
    std::optional<reagent_substructure_search::ReagentData> maybe_proto =
      reader.ReadProto<reagent_substructure_search::ReagentData>();
    if (maybe_proto) {
    } else if (reader.eof()) {
      return 1;
    } else {
      cerr << "Options::ReadReagents:error reading proto\n";
      return 0;
    }

    if (! AddReagent(*maybe_proto)) {
      cerr << "Options::ReadReagent:invalid proto\n";
      // Will fail if it has been destroyed
      cerr << maybe_proto->ShortDebugString() << '\n';
      return 0;
    }
  }

  return _reagent.size();
}

int
Options::AddReagent(reagent_substructure_search::ReagentData& proto) {
  std::unique_ptr<Reagent> rgnt = std::make_unique<Reagent>();

  if (int rc = rgnt->BuildFromProto(proto, _min_natoms, _max_natoms, _ok_ncon); rc == 1) {
    // Good, processed below
  } else if (rc == 0) {  // constraints violated.  // constraints violated.
    return 1;
  } else {
    cerr << "Cannot build Reagent from proto\n";
    return 0;
  }

  rgnt->set_ndx(_reagent.size());  // This later gets overwritten when the _reagent array is sorted.

  _reagent << rgnt.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";

  if (_acc_coverage.n() > 1) {
    output << "Coverage btw " << _acc_coverage.minval() << " and " <<
              _acc_coverage.maxval() << " ave " << _acc_coverage.average() << '\n';
  }

  for (int i = 0; i < _matches.number_elements(); ++i) {
    if (_matches[i]) {
      output << _matches[i] << " molecules matched " << i << " reagents\n";
    }
  }

  uint32_t has_formula = 0;
  uint32_t has_query = 0;
  for (const Reagent* r : _reagent) {
    if (! r->has_formula()) {
      continue;
    }
    ++has_formula;
    if (r->has_query()) {
      ++has_query;
    }
  }

  output << "Across " << _reagent.size() << " reagents\n";
  output << has_formula << " have a Molecule and formula computed\n";
  output << has_query << " have a Substructure_Query instantiated\n";

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
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {

  ++_molecules_read;

  if (_report_progress()) {
    cerr << "Processed " << _molecules_read << " molecules\n";
  }

  const int matoms = m.natoms();

  Molecule_to_Match target(&m);

  Coverage coverage(matoms);

  mformula::MFormula formula;
  formula.set_consider_aromatic(0);
  formula.Build(m);

  mformula::MFormula formula_arom;
  formula_arom.set_consider_aromatic(1);
  formula_arom.Build(m);

  int matches_found = 0;
  for (Reagent* r : _reagent) {
    Substructure_Results sresults;
    sresults.set_save_query_atoms_matched(0);

    if (! r->SubstructureSearch(target, formula, formula_arom, sresults)) {
      continue;
    }

    ++matches_found;

    coverage.Update(sresults, r->ndx());

    if (matches_found > _max_matches) {
      break;
    }
  }

  if (_verbose) {
    ++_matches[matches_found];
  }

  if (matches_found == 0) {
    return HandleNoMatches(m);
  }

  float f = coverage.FractionAtomsCovered();

  _acc_coverage.extra(f);

  if (_stream_for_labelled_molecules.is_open()) {
    // Make a copy - do we really need to????
    Molecule tmp(m);
    for (int i = 0; i < matoms; ++i) {
      isotope_t iso = coverage.Count(i);
      if (iso > 0) {
        tmp.set_isotope(i, iso);
      }
    }

    const float f = iwmisc::Fraction<float>(tmp.number_isotopic_atoms(), matoms);

    _stream_for_labelled_molecules << tmp.aromatic_smiles() << ' ' << m.name() << ' ' << f << '\n';

    _stream_for_labelled_molecules.write_if_buffer_holds_more_than(32768);
  }

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
Options::HandleNoMatches(Molecule& m) {
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
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:m:M:r:O:R:X:f");

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

  if (cl.option_present('f')) {
//  _exit(0);
  }

  return 0;
}

}  // namespace reagent_substructure_search

int
main(int argc, char ** argv) {

  int rc = reagent_substructure_search::Main(argc, argv);

  return rc;
}
