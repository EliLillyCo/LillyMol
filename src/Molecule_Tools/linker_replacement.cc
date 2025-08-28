// Given a set of substituents and linkers, enumerate the possibilities

#include <iostream>
#include <memory>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_set.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/mformula.h"
#include "Molecule_Tools/molecule_filter_lib.h"
#include "Molecule_Tools/set_of_molecules.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace linker_replacement {

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
  cerr << R"(Given a set of linkers and substituents R1 and R2, enumerate the products.
Possible workflow
1. Identify substituents
rgroup -s <smarts> -r 0.1 -r 0.5 ... file.smi > file.r1r2.smi
2. Convert that output into separate files for substituents
rgroup_to_linker_replace.sh -v -S R file.r1r2.smi
3. Use those resulting files with this tool
linker_replacement -L l1.textproto -L l2.textproto -S new_molecules R1.smi R2.smi

 -L <fname>             One or more dicerdata::DicerFragment textproto files with linkers.
                        These must have the right number of isotopic labels.
 -m <natoms>            Discard linker molecules with fewer than <natoms> atoms.
 -M <natoms>            Discard linker molecules with more  than <natoms> atoms.
 -p <support>           Discard linker molecules with fewer than <support> examples.
 -y <smarts>            Only process linkers matching <smarts>.
 -Y <query>             Only process linkers matching <query>.
 -n <smarts>            Do NOT process linkers matching <smarts>.
 -N <query>             do Not process linkers matching <query>.
 -O <fname>             A file containing a single molecule from which a reference molecular formula is derived.
 -O maxdiff=<n>         Only use linkers whose molecular formula differs by at most <n> from the reference formula.
 -F <fname>             a molecule_filter_data::Requirements textproto for filtering the products.
 -I                     remove isotopes from product molecules.
 -X ...                 miscellaneous options, enter '-X help' for info.
 -c                     remove chirality from all input molecules.
 -v                     verbose output.
)";
// clang-format on

  ::exit(rc);
}

class Linker : public Molecule {
  private:
    atom_number_t _a1;
    atom_number_t _a2;

  public:
    Linker();

    // Succeeds if it can parse the smiles and identify two atoms.
    int Build(const dicer_data::DicerFragment& proto);

    atom_number_t a1() const {
      return _a1;
    }
    atom_number_t a2() const {
      return _a2;
    }
};

Linker::Linker() {
  _a1 = kInvalidAtomNumber;
  _a2 = kInvalidAtomNumber;
}

int
Linker::Build(const dicer_data::DicerFragment& proto) {
  if (! proto.has_smi()) {
    cerr << "Linker::Build:no smiles in proto '" << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (! this->build_from_smiles(proto.smi())) {
    cerr << "Linker::Build:cannot parse smiles '" << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_par()) {
    this->set_name(proto.par());
  }

  const int matoms = this->natoms();
  for (int i = 0; i < matoms; ++i) {
    if (this->isotope(i) == 0) {
      continue;
    }
    if (_a1 == kInvalidAtomNumber) {
      _a1 = i;
    } else {
      _a2 = i;
      break;
    }
  }

  if (_a2 == kInvalidAtomNumber) {
    cerr << "Linker::Build:not enough isotopes " << proto.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

class Sidechain : public Molecule {
  private:
    atom_number_t _a1;

  public:
    Sidechain();

    atom_number_t a1() const {
      return _a1;
    }

    // Returns true if it finds an atom with an isotope.
    int IdentifyMatchedAtom();
};

Sidechain::Sidechain() {
  _a1 = kInvalidAtomNumber;
}

int
Sidechain::IdentifyMatchedAtom() {

  const int matoms = natoms();
  for (int i = 0; i < matoms; ++i) {
    if (isotope(i) == 0) {
      continue;
    }

    _a1 = i;
    return 1;
  }

  return 0;
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

    // Min and max atom counts for linkers.
    uint32_t _min_natoms;
    uint32_t _max_natoms;
    uint32_t _rejected_for_atom_count;

    // We discard any linker with an invalid valence.
    uint32_t _linker_rejected_for_bad_valence;

    resizable_array_p<Substructure_Query> _linker_must_have;
    resizable_array_p<Substructure_Query> _linker_must_not_have;
    uint32_t _rejected_for_substructures;

    // We can impose a threshold for how many instances of a linker must be found.
    uint32_t _min_support;
    // A counter for those rejected for low support
    uint64_t _rejected_for_below_support;

    molecule_filter_lib::MoleculeFilter _filter;
    uint64_t _rejected_by_filter;

    Chemical_Standardisation _chemical_standardisation;

    resizable_array_p<Linker> _linker;

    absl::flat_hash_set<IWString> _seen;
    uint64_t _duplicates_rejected;

    resizable_array_p<mformula::MFormula> _mformula;
    uint32_t _max_formula_difference;
    uint32_t _rejected_for_formula_difference;

    uint32_t _product_rejected_for_bad_valence;
    int _report_bad_valence;

    // The -I option.
    int _remove_isotopes_from_products;

    uint64_t _molecules_written;

  // Private functions.
    int ReadLinkers(IWString& fname);
    int ReadLinkers(iwstring_data_source& input);

    int ReadFormula(IWString& fname);
    int ReadFormula(iwstring_data_source& input);

    int OkSubstructures(Molecule& m);
    int OkFormulaDifference(Molecule& m);
 
    int Process(const Sidechain& r1, const Sidechain& r2,
                 IWString_and_File_Descriptor& output);
    int MakeProduct(Molecule& product,
                atom_number_t a1, atom_number_t a2,
                atom_number_t s1, atom_number_t s2,
                IWString_and_File_Descriptor& output);

    int OkAdjacentAtoms(const Molecule& m, atom_number_t a1, atom_number_t a2) const;
    int IsUnique(Molecule& m);
    int OkMoleculeFilter(Molecule& m);
    int OkValence(Molecule& m);
    int WriteProduct(Molecule& m, const Sidechain& r1, const Sidechain& r2,
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

    int Process(set_of_molecules::SetOfMolecules<Sidechain>& r1,
                 set_of_molecules::SetOfMolecules<Sidechain>& r2,
                 IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _min_natoms = 0;
  _max_natoms = 0;
  _rejected_for_atom_count = 0;

  _linker_rejected_for_bad_valence = 0;

  _min_support = 0;
  _rejected_for_below_support = 0;

  _max_formula_difference = 0;
  _rejected_for_formula_difference = 0;

  _rejected_by_filter = 0;

  _rejected_for_substructures = 0;

  _product_rejected_for_bad_valence = 0;
  _report_bad_valence = 1;

  _duplicates_rejected = 0;

  _remove_isotopes_from_products = 0;

  _molecules_written = 0;
}

void
DisplayDashXOptions(std::ostream& output) {
  output << "The following -X options are recognised\n";
  output << " -X nwv            do NOT warn about invalid valences\n";

  ::exit(0);
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

  if (cl.option_present('X')) {
    IWString x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "nwv") {
        _report_bad_valence = 0;
        if (_verbose) {
          cerr << "Will NOT report bad valence errors\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_natoms) || _min_natoms < 2) {
      cerr << "Options::Initialise:the minumum atom count (-m) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard linkers with fewer than " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_natoms) || _max_natoms < _min_natoms) {
      cerr << "Options::Initialise:the maxumum atom count (-M) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard linkers with more than " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _min_support) || _min_support < 1) {
      cerr << "The minimum support level (-p) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard linkers unless they have at least " << _min_support << " instances\n";
    }
  }

  if (cl.option_present('O')) {
    IWString o;
    IWString fname;
    for (int i = 0; cl.value('O', o, i); ++i) {
      if (o.starts_with("maxdiff=")) {
        o.remove_leading_chars(8);
        if (! o.numeric_value(_max_formula_difference)) {
          cerr << "Options::Initialise:invalid maxdiff= directive '" << o << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will discard linkers with formula differences larger than " <<
                  _max_formula_difference << " from target\n";
        }
      } else {
        fname = o;
      }
    }

    if (fname.empty()) {
      cerr << "Options::Initialise:no formula difference file\n";
      Usage(1);
    }

    if (! ReadFormula(fname)) {
      cerr << "Cannot read formula molecule from '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Read " << _mformula.size() << " molecular formula targets\n";
    }
  }

  if (! cl.option_present('L')) {
    cerr << "Must specify one or more dicer_data::DicerFragment textproto files via the -L option\n";
    Usage(1);
  }

  if (cl.option_present('L')) {
    IWString l;
    for (int i = 0; cl.value('L', l, i); ++i) {
      if (_verbose) {
        cerr << "Reading linkers from '" << l << "'\n";
      }
      if (! ReadLinkers(l)) {
        cerr << "Cannot read linkers from '" << l << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _linker.size() << " linkers\n";
      if (_min_support > 0) {
        cerr << _rejected_for_below_support << " linkers rejected for support " << _min_support << '\n';
      }
      cerr << _rejected_for_substructures << " rejected for substructure requirements\n";
      if (_mformula.size() > 0) {
        cerr << _rejected_for_formula_difference << " rejected for formula difference\n";
      }
      cerr << _rejected_for_atom_count << " rejected for atom count\n";
      cerr << _linker_rejected_for_bad_valence << " linkers rejected for valence problems\n";
    }
  }

  if (cl.option_present('y')) {
    IWString smarts;
    for (int i = 0; cl.value('y', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _linker_must_have << q.release();
    }
  }

  if (cl.option_present('Y')) {
    if (! process_queries(cl, _linker_must_have, _verbose, 'Y')) {
      cerr << "Options::Initialise:cannot process linker must have queries (-Y)\n";
      return 0;
    }
  }

  if (cl.option_present('n')) {
    IWString smarts;
    for (int i = 0; cl.value('n', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _linker_must_not_have << q.release();
    }
  }

  if (cl.option_present('N')) {
    if (! process_queries(cl, _linker_must_not_have, _verbose, 'N')) {
      cerr << "Options::Initialise:cannot process linker must not have queries (-N)\n";
      return 0;
    }
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! _filter.Build(fname)) {
      cerr << "Cannot initialise molecule filter (-F)\n";
      return 0;
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes_from_products = 1;
    if (_verbose) {
      cerr << "Will remove isotopes from product molecules\n";
    }
  }

  return 1;
}

int
Options::ReadLinkers(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadLinkers:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadLinkers(input);
}

int
Options::ReadLinkers(iwstring_data_source& input) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
    dicer_data::DicerFragment proto;
    if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
      cerr << "Options:ReadLinkers;cannot parse proto " << line << '\n';
      return 0;
    }

    if (_min_support > 0 && proto.has_n() && proto.n() < _min_support) {
      ++_rejected_for_below_support;
      continue;
    }

    if (_min_natoms > 0 && proto.has_nat() && proto.nat() < _min_natoms) {
      ++_rejected_for_atom_count;
      continue;
    }

    if (_max_natoms > 0 && proto.has_nat() && proto.nat() > _max_natoms) {
      ++_rejected_for_atom_count;
      continue;
    }

    std::unique_ptr<Linker> lkr = std::make_unique<Linker>();
    if (! lkr->Build(proto)) {
      cerr << "Options::ReadLinkers:cannot parse proto '" << line << "'\n";
      return 0;
    }

    if (! lkr->valence_ok()) {
      ++_product_rejected_for_bad_valence;
      continue;
    }

    if (! OkSubstructures(*lkr)) {
      ++_rejected_for_substructures;
      continue;
    }

    if (! _mformula.empty() && ! OkFormulaDifference(*lkr)) {
      ++_rejected_for_formula_difference;
      continue;
    }

    _linker << lkr.release();
  }

  return _linker.number_elements();
}

int
Options::ReadFormula(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFormula:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFormula(input);
}

int
Options::ReadFormula(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    Molecule m;
    if (! m.build_from_smiles(buffer)) {
      cerr << "Options::ReadFormula:invalid smiles '" << buffer << "'\n";
      return 0;
    }

    std::unique_ptr<mformula::MFormula> f = std::make_unique<mformula::MFormula>();
    if (! f->Build(m)) {
      cerr << "Options::ReadFormula:cannot compute formula '" << buffer << "'\n";
      return 0;
    }
    _mformula << f.release();
  }

  return _mformula.number_elements();
}

int
Options::OkSubstructures(Molecule& m) {
  if (_linker_must_have.empty() && _linker_must_not_have.empty()) {
    return 1;
  }

  Molecule_to_Match target(&m);
  if (!_linker_must_have.empty()) {
    if (! lillymol::AnyQueryMatches(target, _linker_must_have)) {
      return 0;
    }
  }
  if (!_linker_must_not_have.empty()) {
    if (lillymol::AnyQueryMatches(target, _linker_must_not_have)) {
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Generated " << _molecules_written << " molecules\n";
  output << _duplicates_rejected << " duplicates rejected\n";
  output << _product_rejected_for_bad_valence << " invalid valence products discarded\n";
  if (_filter.active()) {
    output << _rejected_by_filter << " rejected by the product molecule filter\n";
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

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return 1;
}

int
Options::Process(set_of_molecules::SetOfMolecules<Sidechain>& r1,
                 set_of_molecules::SetOfMolecules<Sidechain>& r2,
                 IWString_and_File_Descriptor& output) {

  const uint32_t n = r1.size();
  assert(r2.size() == n);

  for (uint32_t i = 0; i < n; ++i) {
    if (! Process(*r1[i], *r2[i], output)) {
      cerr << "Options::Process:error on combination " << i << '\n';
      return 0;
    }
  }

  return 1;
}

int
RemoveBondIfNeeded(Molecule& m, atom_number_t a1, atom_number_t a2) {
  if (m.are_bonded(a1, a2)) {
    m.remove_bond_between_atoms(a1, a2);
    return 1;
  }

  return 0;
}

int
Options::Process(const Sidechain& r1, const Sidechain& r2,
                 IWString_and_File_Descriptor& output) {
  for (const Linker * linker : _linker) {
    Molecule product(*linker);
    const int initial_natoms = product.natoms();
    product += r1;
    product += r2;

    atom_number_t r1a = initial_natoms + r1.a1();
    atom_number_t r2a = initial_natoms + r1.natoms() + r2.a1();
    // cerr << "atoms " << linker->a1() << ' ' << r1a << ' ' << linker->a2() << ' ' << r2a << '\n';
    if (MakeProduct(product, linker->a1(), r1a, linker->a2(), r2a, output)) {
      WriteProduct(product, r1, r2, output);
    }

    RemoveBondIfNeeded(product, linker->a1(), r1a);
    RemoveBondIfNeeded(product, linker->a2(), r2a);
    // cerr << product.smiles() << " after removing bonds\n";

    if (MakeProduct(product, linker->a2(), r1a, linker->a1(), r2a, output)) {
      WriteProduct(product, r1, r2, output);
    }
  }

  return 1;
}

int
Options::MakeProduct(Molecule& product,
                atom_number_t a1, atom_number_t s1,
                atom_number_t a2, atom_number_t s2,
                IWString_and_File_Descriptor& output) {
  if (!OkAdjacentAtoms(product, a1, s1)) {
    return 0;
  }
  if (!OkAdjacentAtoms(product, a2, s2)) {
    return 0;
  }

  product.add_bond(a1, s1, SINGLE_BOND);
  product.add_bond(a2, s2, SINGLE_BOND);

  product.unset_all_implicit_hydrogen_information(a1);
  product.unset_all_implicit_hydrogen_information(s1);
  product.unset_all_implicit_hydrogen_information(a2);
  product.unset_all_implicit_hydrogen_information(s2);

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(product);
  }

  // Before checking for uniqueness.
  if (_remove_isotopes_from_products) {
    product.unset_isotopes();
  }

  if (! IsUnique(product)) {
    return 0;
  }

  if (! OkMoleculeFilter(product)) {
    return 0;
  }

  if (! OkValence(product)) {
    return 0;
  }

  return 1;
}

// Return true if it would be OK to add a bond between `a1` and `s2`.
int
Options::OkAdjacentAtoms(const Molecule& m, atom_number_t a1, atom_number_t a2) const {
  atomic_number_t z1 = m.atomic_number(a1);
  atomic_number_t z2 = m.atomic_number(a2);

  if (z1 == 6 || z2 == 6) {
    return 1;
  }

  if (z1 == 8 && z2 == 8) {
    return 0;
  }

  if (z1 == 7 && z2 == 7) {
    return 0;
  }

  if (z1 == 16 && z2 == 16) {
    return 0;
  }

  // Do not form O-N bonds
  if (z1 == 7 && z2 == 8) {
    return 0;
  }
  if (z1 == 8 && z2 == 7) {
    return 0;
  }

  return 1;
}

int
Options::IsUnique(Molecule& m) {
  if (const auto iter = _seen.find(m.unique_smiles()); iter != _seen.end()) {
    ++_duplicates_rejected;
    return 0;
  }

  _seen.insert(m.unique_smiles());

  return 1;
}

// Return true if the molecular formula of `m` is closer than _max_formula_difference
// to any of the formulae in _mformula.
int
Options::OkFormulaDifference(Molecule& m) {
  mformula::MFormula formula;
  formula.Build(m);

  for (mformula::MFormula* f : _mformula) {
    uint32_t d = formula.Diff(*f);
    if (d <= _max_formula_difference) {
      return 1;
    }
  }

  return 0;
}

int
Options::OkMoleculeFilter(Molecule& m) {
  if (! _filter.active()) {
    return 1;
  }

  if (_filter.Ok(m)) {
    return 1;
  }

  ++_rejected_by_filter;

  return 0;
}

int
Options::OkValence(Molecule& m) {
  if (m.valence_ok()) {
    return 1;
  }

  ++_product_rejected_for_bad_valence;

  if (_report_bad_valence) {
    cerr << m.smiles() << ' ' << m.name() << " bad valence\n";
  }

  return 0;
}


int
Options::WriteProduct(Molecule& m, const Sidechain& r1, const Sidechain& r2,
                IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

  // Get rid of unique smiles.
  m.invalidate_smiles();

  output << m.smiles() << kSep << m.name() << kSep << r1.name() << kSep << r2.name() << '\n';

  ++_molecules_written;

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

// Remove any bad valences from the substituents.
int
RemoveBadValence(set_of_molecules::SetOfMolecules<Sidechain>& r1,
      set_of_molecules::SetOfMolecules<Sidechain>& r2) {
  assert(r1.number_elements() == r2.number_elements());

  int rc = 0;

  for (int i = r1.number_elements() - 1; i >= 0; --i) {
    if (r1[i]->valence_ok() && r2[i]->valence_ok()) {
      continue;
    }

    r1.remove_item(i);
    r2.remove_item(i);
    ++rc;
  }

  return rc;
}

// If we have duplicate R1 R2 combinations, those will form the same product.
// Remove them.
int
DeDup(set_of_molecules::SetOfMolecules<Sidechain>& r1,
      set_of_molecules::SetOfMolecules<Sidechain>& r2) {
  assert(r1.number_elements() == r2.number_elements());

  absl::flat_hash_set<IWString> seen;

  int rc = 0;

  for (int i = r1.number_elements() - 1; i >= 0; --i) {
    IWString s = r1[i]->unique_smiles();
    s << '.' << r2[i]->unique_smiles();
    if (auto iter = seen.find(s); iter != seen.end()) {
      r1.remove_item(i);
      r2.remove_item(i);
      ++rc;
    } else {
      seen.insert(s);
    }
  }

  return rc;
}

int
LinkerReplacement(Options& options,
                set_of_molecules::SetOfMolecules<Sidechain>& r1,
                set_of_molecules::SetOfMolecules<Sidechain>& r2,
                IWString_and_File_Descriptor& output) {
  return options.Process(r1, r2, output);
}

int
IdentifyIsotopicLabel(set_of_molecules::SetOfMolecules<Sidechain>& mols) {
  for (Sidechain* m : mols) {
    if (!m->IdentifyMatchedAtom()) {
      cerr << m->smiles() << ' ' << m->name() << " no attachment point\n";
      return 0;
    }
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:g:clL:m:M:p:y:Y:n:N:F:O:I:X:");

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

  if (cl.size() != 2) {
    cerr << "Must specify two input files, R1 and R2, each with isotopically labelled smiles\n";
    Usage(1);
  }

  set_of_molecules::SetOfMolecules<Sidechain> r1, r2;
  IWString tmp(cl[0]);
  if (! r1.Build(tmp)) {
    cerr << "Cannot read '" << cl[0] << "'\n";
    return 1;
  }
  tmp = cl[1];
  if (! r2.Build(tmp)) {
    cerr << "Cannot read '" << cl[1] << "'\n";
    return 1;
  }

  if (r1.size() != r2.size()) {
    cerr << "Must have equal numbers of R1 " << r1.size() << " and R2 " << r2.size() << '\n';
    return 1;
  }

  if (! IdentifyIsotopicLabel(r1) || ! IdentifyIsotopicLabel(r2)) {
    cerr << "Cannot identify labelled atoms in sidechains\n";
    return 1;
  }

  if (int dedup = DeDup(r1, r2); dedup > 0) {
    if (verbose) {
      cerr << "Removed " << dedup << " duplicate forming combinations\n";
    }
  }

  if (int badv = RemoveBadValence(r1, r2); badv > 0) {
    if (verbose) {
      cerr << "Removed " << badv << " sidechains with possible invalid valences\n";
    }
  }

  if (verbose) {
    cerr << "For each linker will build " << r1.size() << " products\n";
  }

  if (r1.size() != r2.size()) {
    cerr << "Internal error on R1R2 size " << r1.size() << " vs " << r2.size() << '\n';
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  if (! LinkerReplacement(options, r1, r2, output)) {
    cerr << "LinkerReplacement::fatal error\n";
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace linker_replacement

int
main(int argc, char ** argv) {

  int rc = linker_replacement::Main(argc, argv);

  return rc;
}
