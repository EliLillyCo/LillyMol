// Generate pharmacophore-like queries for 2D molecules

#include <iostream>
#include <limits>
#include <memory>
#include <optional>

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/pharmacophore_2d.pb.h"

namespace pharmacophore_2d {

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
  cerr << R"(Generates query files based on 2D pharmacophore features.
Queries define pharmacophore features and query files are generated that describe the topological relationships
between the pharmacophoric features in the starting molecules. One query per starting molecule.
Parameters should most be specified via the -C option, but some are also available via command line options.
 -C <fname>       pharmacophore2d::Pharmacophore2DConfig configuration textproto
 -S <stem>        file name stem for generated query files - multiple textproto files are generated.
 -G <fname>       write names of query files generated to <fname> - use with -q PROTOFILE:<fname>
 -t               write serialized SubstructureSearch::SubstructureQuery protos to -S value. One file generated -q TFPROTO:<fname>

 -s <smarts>      define functional groups via smarts - or use the functional_group attribute in the -C file.
 -q <query>       define functional groups via query files - or use the functional_group attribute in the -C file.
 -n               ncon values become min_ncon values in the query\n";
 -d <dist>        min separation between atoms\n";
 -D <dist>        max separation between atoms\n";
 -Y ...           other options, enter '-Y help' for info\n";
 -v               verbose output\n";
)";
  // clang-format on

  ::exit(rc);
}

// The atomic properties that can be transferred to the query.
// We OR these into an uint32_t
static constexpr uint32_t kPropAtomicNumber = 1;
static constexpr uint32_t kPropRingBondCount = 2;
static constexpr uint32_t kPropAromatic = 4;
static constexpr uint32_t kPropNcon = 8;
static constexpr uint32_t kPropHasPiElectron = 16;
static constexpr uint32_t kPropPiElectronCount = 32;
static constexpr uint32_t kPropUnsaturation = 64;
static constexpr uint32_t kPropIsotope = 128;
static constexpr uint32_t kPropRingSize = 256;
static constexpr uint32_t kPropSpinach = 512;
static constexpr uint32_t kPropFusedSystemSize = 1024;
static constexpr uint32_t kPropHcount = 2048;

class PerMoleculeData {
 private:
  Molecule& _m;

  const int _matoms;

  int* _rotbond_between;

  int* _functional_group;

  quick_rotbond::QuickRotatableBonds _rotbond;

 public:
  PerMoleculeData(Molecule& m);
  ~PerMoleculeData();

  int*
  functional_group() {
    return _functional_group;
  }

  int functional_group(atom_number_t a) const {
    return _functional_group[a];
  }

  int& functional_group(atom_number_t a) {
    return _functional_group[a];
  }

  int rotbond_between(atom_number_t a1, atom_number_t a2) const {
    return _rotbond_between[a1 * _matoms + a2];
  }
};

PerMoleculeData::PerMoleculeData(Molecule& m) : _m(m), _matoms(m.natoms()) {
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);

  _rotbond_between = _rotbond.RotatableBondsBetween(m).release();

  _functional_group = new_int(_matoms, -1);
}

PerMoleculeData::~PerMoleculeData() {
  delete[] _rotbond_between;
}

class Output {
 private:
  int _verbose;

  int _nfiles;

  IWString _stem;

  std::vector<IWString> _files_generated;

  std::unique_ptr<iw_tf_data_record::TFDataWriter> _tfdata;

 public:
  Output();

  int Initialise(Command_Line& cl);

  int Write(const SubstructureSearch::SubstructureQuery& qry);

  int nfiles() const {
    return _nfiles;
  }

  int WriteFilesGenerated(IWString& fname) const;
  int WriteFilesGenerated(IWString_and_File_Descriptor& fname) const;
};

Output::Output() {
  _verbose = 0;

  _nfiles = 0;

  _files_generated.reserve(10000);
}

int
Output::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (!cl.option_present('S')) {
    cerr << "Output::Initialise:must specify output file name stem with the -S option\n";
    Usage(0);
  }

  cl.value('S', _stem);

  if (cl.option_present('t')) {
    _stem.EnsureEndsWith(".dat");
    _tfdata = std::make_unique<iw_tf_data_record::TFDataWriter>();
    if (!_tfdata->Open(_stem)) {
      cerr << "Output::Initialise:cannot open TFDataRecord file '" << _stem << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "TFDataRecord serialized protos written to '" << _stem << "'\n";
    }

    return 1;
  }

  if (_verbose) {
    cerr << "Files written with stem '" << _stem << "'\n";
  }

  return 1;
}

int
Output::Write(const SubstructureSearch::SubstructureQuery& qry) {
  if (_tfdata) {
    return _tfdata->WriteSerializedProto<SubstructureSearch::SubstructureQuery>(qry);
  }

  ++_nfiles;
  IWString fname;
  fname << _stem << _nfiles << ".textproto";

  _files_generated.push_back(fname);

  return iwmisc::WriteTextProto(qry, fname);
}

int
Output::WriteFilesGenerated(IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Output::WriteFilesGenerated:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteFilesGenerated(output);
}

int
Output::WriteFilesGenerated(IWString_and_File_Descriptor& output) const {
  if (_verbose) {
    cerr << "writing " << _files_generated.size() << " generated query file names\n";
  }

  for (const IWString& fname : _files_generated) {
    output << fname << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

class Options {
 private:
  int _verbose;

  int _molecules_read;
  int _no_functional_groups;
  int _only_one_functional_group;

  resizable_array_p<Substructure_Query> _atoms_to_ignore;

  // We can specify certain arrangements of atoms that will be placed
  // in functional groups, which might not necessarily be discovered
  // otherwise. Acids are an example.
  resizable_array_p<Substructure_Query> _external_functional_groups;

  pharmacophore2d::Pharmacophore2DConfig _proto;

  // Statistics on how the functional groups match.
  extending_resizable_array<int> _functional_group_count;
  Accumulator<float> _fraction_atoms_in_functional_groups;

  // When generating distance constraints, ignore atoms outside the
  // range specified here.
  int _min_separation = 1;
  int _max_separation = std::numeric_limits<int>::max();

  // When two atoms are separated by `d` bonds, the min_bonds_between
  // will be `d - delta_shorter` (if positive) and the max_bonds_between
  // will be `d + delta_longer`.
  int _delta_shorter = 0;
  int _delta_longer = 0;

  int _ncon_becomes_min_ncon = 0;

  int _ring_bond_count_becomes_min_ring_bond_count = 0;

  IWString_and_File_Descriptor _stream_for_labelled_smiles;

  // If writing labelled smiles, we can write as either atom map numbers, default
  // or as isotopic labels. This variable controls that behaviour.
  int _label_with_isotopes = 0;

  // The atomic properties that will be transferred to the query atom.
  // These are or'd values of the various kProp* variables.
  uint32_t _atomic_properties;

  quick_rotbond::QuickRotatableBonds _rotbond;

  int _reduce_to_largest_fragment;

  // Private functions.
  void DisplayDashYOptions(std::ostream& output);

  int Pharmacophore2d(Molecule& m, PerMoleculeData& pmd, int number_functional_groups,
                      IWString& fname);
  std::optional<SubstructureSearch::SubstructureQuery> Pharmacophore2d(
      Molecule& m, PerMoleculeData& pmd, int number_functional_groups);

  int IdentifyAtomsToIgnore(Molecule& m, int* ignore_atom);
  std::tuple<int, int> IdentifyFunctionalGroups(Molecule& m, int* fg, const int* ignore);
  int WriteLabelledSmiles(const Molecule& m, const int* functional_group);
  void IdentifyExternallySpecified(Molecule& m, int* fg, const int* ignore,
                                   int& atoms_in_functional_groups, int& group_number);
  int BuildQuery(Molecule& m, const int* functional_group, int group_number,
                 const resizable_array<atom_number_t>& atom_order_in_smiles,
                 SubstructureSearch::SingleSubstructureQuery& query);
  void CopyAttributes(Molecule& m, atom_number_t zatom,
                      SubstructureSearch::SubstructureAtom& query_atom) const;
  int AddDistanceConstraints(Molecule& m, PerMoleculeData& pmd,
                             SubstructureSearch::SingleSubstructureQuery& query);
  int AddRotbond(Molecule& m, const PerMoleculeData& pmd, atom_number_t a1,
                 atom_number_t a2, SubstructureSearch::SeparatedAtoms& separated_atoms);

 public:
  Options();

  int Initialise(Command_Line& cl);

  void Preprocess(Molecule& m);

  int Process(Molecule& m, IWString& output_fname);
  int Process(Molecule& m, IWString_and_File_Descriptor& output);
  std::optional<SubstructureSearch::SubstructureQuery> Process(Molecule& m);

  int verbose() const {
    return _verbose;
  }

  int molecules_read() const {
    return _molecules_read;
  }

  int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;

  _molecules_read = 0;
  _no_functional_groups = 0;
  _only_one_functional_group = 0;

  _min_separation = 1;
  _max_separation = std::numeric_limits<int>::max();

  _delta_shorter = 0;
  _delta_longer = 0;

  _ncon_becomes_min_ncon = 0;

  _ring_bond_count_becomes_min_ring_bond_count = 0;

  _label_with_isotopes = 0;

  _atomic_properties = (kPropAtomicNumber | kPropAromatic);

  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);

  _reduce_to_largest_fragment = 0;
}

void
Options::DisplayDashYOptions(std::ostream& output) {
  output << " -Y iso            label functional groups in the -L file as isotopes "
            "(default is atom map)\n";
  output << " -Y fname=<fname>  file name to which labelled smiles are written\n";

  ::exit(0);
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to the largest fragment\n";
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<pharmacophore2d::Pharmacophore2DConfig> p =
        iwmisc::ReadTextProto<pharmacophore2d::Pharmacophore2DConfig>(fname);

    if (!p) {
      cerr << "Cannot read textproto config file '" << fname << "'\n";
      return 0;
    }

    _proto = *p;
    if (_verbose) {
      cerr << "Read proto config from '" << fname << "'\n";
    }

    if (_proto.has_min_separation()) {
      _min_separation = _proto.min_separation();
    }
    if (_proto.has_max_separation()) {
      _max_separation = _proto.max_separation();
    }
    if (_proto.has_ncon_becomes_min_ncon()) {
      _ncon_becomes_min_ncon = _proto.ncon_becomes_min_ncon();
    }
    if (_proto.has_ring_bond_count_becomes_min_ring_bond_count()) {
      _ring_bond_count_becomes_min_ring_bond_count =
          _proto.ring_bond_count_becomes_min_ring_bond_count();
    }

    if (_proto.has_delta_shorter()) {
      _delta_shorter = _proto.delta_shorter();
    }
    if (_proto.has_delta_longer()) {
      _delta_longer = _proto.delta_longer();
    }

    for (const std::string& q : _proto.functional_group()) {
      IWString tmp(q);
      if (!process_cmdline_token('*', tmp, _external_functional_groups, _verbose)) {
        cerr << "Invalid functional group specification '" << q << "'\n";
        return 0;
      }
    }

    for (const std::string& q : _proto.atoms_to_ignore()) {
      IWString tmp(q);
      if (!process_cmdline_token('*', tmp, _atoms_to_ignore, _verbose)) {
        cerr << "Invalid atoms to ignore specification '" << q << "'\n";
        return 0;
      }
    }

    if (_proto.atomic_property_size() > 0) {
      _atomic_properties = 0;
      for (auto ap : _proto.atomic_property()) {
        switch (ap) {
          case pharmacophore2d::ATOMIC_NUMBER:
            _atomic_properties |= kPropAtomicNumber;
            break;
          case pharmacophore2d::NCON:
            _atomic_properties |= kPropNcon;
            break;
          case pharmacophore2d::AP_AROMATIC:
            _atomic_properties |= kPropAromatic;
            break;
          case pharmacophore2d::RING_BOND_COUNT:
            _atomic_properties |= kPropRingBondCount;
            break;
          case pharmacophore2d::HAS_PI_ELECTRON:
            _atomic_properties |= kPropHasPiElectron;
            break;
          case pharmacophore2d::PI_ELECTRON_COUNT:
            _atomic_properties |= kPropPiElectronCount;
            break;
          case pharmacophore2d::UNSATURATION:
            _atomic_properties |= kPropUnsaturation;
            break;
          case pharmacophore2d::ISOTOPE:
            _atomic_properties |= kPropIsotope;
            break;
          case pharmacophore2d::RING_SIZE:
            _atomic_properties |= kPropRingSize;
            break;
          case pharmacophore2d::SPINACH:
            _atomic_properties |= kPropSpinach;
            break;
          case pharmacophore2d::FUSED_SYSTEM_SIZE:
            _atomic_properties |= kPropFusedSystemSize;
            break;
          case pharmacophore2d::HCOUNT:
            _atomic_properties |= kPropHcount;
            break;
        }
      }
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smt;
    for (int i = 0; cl.value('s', smt, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (!q->create_from_smarts(smt)) {
        cerr << "Invalid functional group smarts '" << smt << "'\n";
        return 1;
      }

      _external_functional_groups << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, _external_functional_groups, _verbose, 'q')) {
      cerr << "Cannot read external queries (-q)\n";
      return 0;
    }
  }

  if (_verbose && _external_functional_groups.size() > 0) {
    cerr << "Defined " << _external_functional_groups.size()
         << " externally specified functional groups\n";
  }

  if (cl.option_present('d')) {
    if (!cl.value('d', _min_separation) || _min_separation < 1) {
      cerr << "Invalid minimum separation (-d)\n";
      Usage(1);
    }
    if (_verbose) {
      cerr << "Will ignore atoms separated by less than " << _min_separation
           << " bonds\n";
    }
  }

  if (cl.option_present('D')) {
    if (!cl.value('D', _max_separation) || _max_separation < _min_separation) {
      cerr << "Invalid maximum separation (-D)\n";
      Usage(1);
    }
    if (_verbose) {
      cerr << "Will ignore atoms separated by more than " << _max_separation
           << " bonds\n";
    }
  }

  if (cl.option_present('n')) {
    _ncon_becomes_min_ncon = 1;
    if (_verbose) {
      cerr << "Ncon values become min_ncon query attributes\n";
    }
  }

  if (cl.option_present('Y')) {
    IWString fname;
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "iso") {
        _label_with_isotopes = 1;
        if (_verbose) {
          cerr << "Will write functional group labels as isotopes\n";
        }
      } else if (y.starts_with("fname=")) {
        y.remove_leading_chars(6);
        fname = y;
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }

    if (fname.empty()) {
      cerr << "If writing labelled smiles must specify file name\n";
      return 0;
    }

    fname.EnsureEndsWith(".smi");
    if (!_stream_for_labelled_smiles.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for labelled smiles '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Labelled molecules written to '" << fname << "'\n";
    }
  }

  return 1;
}

void
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
}

int
Options::AddDistanceConstraints(Molecule& m, PerMoleculeData& pmd,
                                SubstructureSearch::SingleSubstructureQuery& query) {
  const int matoms = m.natoms();

  const int* functional_group = pmd.functional_group();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (functional_group[i] < 0) {
      continue;
    }
    for (int j = 0; j < i - 1; ++j) {
      if (functional_group[j] < 0) {
        continue;
      }
      if (functional_group[i] == functional_group[j]) {
        continue;
      }
      // Atoms in two different functional groups.
      const int d = m.bonds_between(i, j);
      if (d < _min_separation || d > _max_separation) {
        continue;
      }

      ++rc;

      SubstructureSearch::SeparatedAtoms* separated_atoms = query.add_separated_atoms();
      separated_atoms->set_a1(i);
      separated_atoms->set_a2(j);
      if (_delta_longer == 0 && _delta_shorter == 0) {
        separated_atoms->add_bonds_between(d);
      } else {
        if (_delta_shorter > 0 && d - _delta_shorter > 0) {
          separated_atoms->set_min_bonds_between(d - _delta_shorter);
        }
        if (_delta_longer > 0) {
          separated_atoms->set_max_bonds_between(d + _delta_longer);
        }
      }

      if (_proto.preserve_rotbond() || _proto.has_extra_rotbond() ||
          _proto.has_less_rotbond()) {
        AddRotbond(m, pmd, i, j, *separated_atoms);
      }
    }
  }

  return rc;
}

int
Options::AddRotbond(Molecule& m, const PerMoleculeData& pmd, atom_number_t a1,
                    atom_number_t a2,
                    SubstructureSearch::SeparatedAtoms& separated_atoms) {
  int rb = pmd.rotbond_between(a1, a2);

  // If we set a min or max, we do not need to set a value.
  int add_rotbond = 1;

  if (_proto.has_extra_rotbond()) {
    separated_atoms.set_max_rotbond(rb + _proto.extra_rotbond());
    add_rotbond = 0;
  }

  if (_proto.has_less_rotbond()) {
    int r = rb - _proto.less_rotbond();
    if (r > 0) {
      separated_atoms.set_min_rotbond(r);
    }
    add_rotbond = 0;
  }

  if (add_rotbond) {
    separated_atoms.add_rotbond(rb);
  }

  return 1;
}

#ifdef NO_LONGER_USED_JJJ
void
Options::AddDistanceConstraints(Molecule& m, const int* functional_group, int g1, int g2,
                                SubstructureSearch::SingleSubstructureQuery& query) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }
      if (functional_group[i] != g1 || functional_group[j] != g2) {
        continue;
      }
      const int d = m.bonds_between(i, j);
      if (d < _min_separation || d > _max_separation) {
        continue;
      }

      SubstructureSearch::SeparatedAtoms* separated_atoms = query.add_separated_atoms();
      separated_atoms->set_a1(i);
      separated_atoms->set_a2(j);
      if (_delta_longer == 0 && _delta_shorter == 0) {
        separated_atoms->add_bonds_between(d);
      } else {
        if (d - _delta_shorter > 0) {
          separated_atoms->set_min_bonds_between(d - _delta_shorter);
        }
        separated_atoms->set_max_bonds_between(d + _delta_longer);
      }
    }
  }
}
#endif

int
SpinachAttribute(Molecule& m, atom_number_t zatom) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);

  m.identify_spinach(spinach.get());

  if (spinach[zatom]) {
    return 1;
  }

  return 0;
}

void
Options::CopyAttributes(Molecule& m, atom_number_t zatom,
                        SubstructureSearch::SubstructureAtom& query_atom) const {
  SubstructureSearch::SubstructureAtomSpecifier* qas = query_atom.add_atom_properties();
  const Atom& atom = m.atom(zatom);
  if (_atomic_properties & kPropAtomicNumber) {
    qas->add_atomic_number(atom.atomic_number());
  }

  if ((_atomic_properties & kPropNcon) == 0) {
  } else if (_ncon_becomes_min_ncon) {
    qas->set_min_ncon(atom.ncon());
  } else {
    qas->add_ncon(atom.ncon());
  }

  if (_atomic_properties & kPropRingBondCount) {
    if (_ring_bond_count_becomes_min_ring_bond_count) {
      qas->set_min_ring_bond_count(m.ring_bond_count(zatom));
    } else {
      qas->add_ring_bond_count(m.ring_bond_count(zatom));
    }
  }

  // Should we also set this if the atom is not aromatic.
  if (_atomic_properties & kPropAromatic) {
    if (m.is_aromatic(zatom)) {
      qas->set_aromatic(true);
    }
  }

  if (_atomic_properties & kPropUnsaturation) {
    if (m.saturated(zatom)) {
      qas->add_unsaturation(0);
    } else {
      qas->set_min_unsaturation(1);  // Or should we copy it.
    }
  }

  if (_atomic_properties & kPropHasPiElectron) {
    cerr << "Warning kPropHasPiElectron not implemented\n";
  }
  if (_atomic_properties & kPropPiElectronCount) {
    cerr << "Warning kPropPiElectronCount not implemented\n";
  }

  if (_atomic_properties & kPropSpinach) {
    qas->set_match_spinach_only(SpinachAttribute(m, zatom));
  }

  if (_atomic_properties & kPropFusedSystemSize) {
    if (m.ring_bond_count(zatom)) {
      qas->add_fused_system_size(m.fused_system_size(zatom));
    }
  }

  if (_atomic_properties & kPropRingSize) {
    const Ring* r = m.ring_containing_atom(zatom);
    if (r != nullptr) {
      qas->add_ring_size(r->size());
    }
  }

  if (_atomic_properties & kPropHcount) {
    if (_proto.hcount_becomes_min()) {
      qas->set_min_hcount(m.hcount(zatom));
    } else {
      qas->add_hcount(m.hcount(zatom));
    }
  }

  // It is an open question as to whether we should set isotope 0 or not.
  if (_atomic_properties & kPropIsotope && m.isotope(zatom)) {
    qas->add_isotope(m.isotope(zatom));
  }
}

int
AddBonds(Molecule& m, const resizable_array<atom_number_t>& atom_order_in_smiles,
         const atom_number_t zatom, const int* functional_group, int group_number,
         SubstructureSearch::SubstructureAtom& query_atom) {
  const Atom& atom = m.atom(zatom);
  int rc = 0;
  for (const Bond* b : atom) {
    const atom_number_t other = b->other(zatom);
    if (functional_group[other] != group_number) {
      continue;
    }
    if (atom_order_in_smiles[other] > atom_order_in_smiles[zatom]) {
      continue;
    }
    //  cerr << "AddBonds:adding bond between " << zatom << " fg " <<
    //  functional_group[zatom] << " and " << other << " fg " << functional_group[other]
    //  << '\n';
    SubstructureSearch::SubstructureBond* query_bond = query_atom.add_query_bond();
    query_bond->set_other_end(other);
    if (b->is_aromatic()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_AROMATIC_BOND);
    } else if (b->is_single_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_SINGLE_BOND);
    } else if (b->is_double_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_DOUBLE_BOND);
    } else if (b->is_triple_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_TRIPLE_BOND);
    }
    rc++;
  }

  return rc;
}

// Return the atoms for which `functional_group[i] == group_number` ordered
// by `atom_order_in_smiles`.
// natoms is the size of the `functional_group` array.
Set_of_Atoms
AtomsInFunctionalGroup(int natoms, const int* functional_group, int group_number,
                       const resizable_array<atom_number_t>& atom_order_in_smiles) {
  Set_of_Atoms atoms_in_group;
  for (int i = 0; i < natoms; ++i) {
    if (functional_group[i] == group_number) {
      atoms_in_group << i;
    }
  }
  atoms_in_group.iwqsort_lambda([&atom_order_in_smiles](int a1, int a2) {
    if (atom_order_in_smiles[a1] < atom_order_in_smiles[a2]) {
      return -1;
    }
    if (atom_order_in_smiles[a1] > atom_order_in_smiles[a2]) {
      return 1;
    }
    return 0;  // Should never happen here.
  });

  return atoms_in_group;
}

int
Options::BuildQuery(Molecule& m, const int* functional_group, int group_number,
                    const resizable_array<atom_number_t>& atom_order_in_smiles,
                    SubstructureSearch::SingleSubstructureQuery& query) {
  const int matoms = m.natoms();
  const resizable_array<atom_number_t> atoms_in_functional_group = AtomsInFunctionalGroup(
      matoms, functional_group, group_number, atom_order_in_smiles);
  for (atom_number_t a : atoms_in_functional_group) {
    SubstructureSearch::SubstructureAtom* atom = query.add_query_atom();
    atom->set_id(a);
    CopyAttributes(m, a, *atom);
    AddBonds(m, atom_order_in_smiles, a, functional_group, group_number, *atom);
  }

  return 1;
}

#ifdef NOT_NEEDED_ASDASD
Set_of_Atoms
GetFunctionalGroup(const int* functional_group, int group_number, int n) {
  Set_of_Atoms result;
  for (int i = 0; i < n; ++i) {
    if (functional_group[i] == group_number) {
      result << i;
    }
  }
  return result;
}

int
Options::Pharmacophore2d(Molecule& m, PerMoleculeData& pmd, int number_functional_groups,
                         IWString& fname) {
  if (_verbose > 1) {
    cerr << m.name() << " has " << number_functional_groups << " functional groups\n";
  }

  // Force a smiles computation, and copy the atom order.
  m.smiles();
  resizable_array<atom_number_t> atom_order_in_smiles = m.atom_order_in_smiles();

  SubstructureSearch::SubstructureQuery composite_query;
  composite_query.set_name(m.name().data(), m.name().length());
  SubstructureSearch::SingleSubstructureQuery* query = composite_query.add_query();
  query->set_respect_initial_atom_numbering(true);
  int group_number = 0;
  for (int i = 1; i <= number_functional_groups; ++i, ++group_number) {
    BuildQuery(m, pmd.functional_group(), i, atom_order_in_smiles, *query);
  }

  if (!AddDistanceConstraints(m, pmd, *query)) {
    return 0;
  }

  if (_proto.has_max_extra_atoms()) {
    auto* mpr = query->mutable_required_molecular_properties();
    mpr->set_max_natoms(m.natoms() + _proto.max_extra_atoms());
  }

  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Options::Pharmacophore2d:cannot open '" << fname << "'\n";
    return 0;
  }

  using google::protobuf::io::FileOutputStream;
  using google::protobuf::io::ZeroCopyOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(
      new FileOutputStream(output.fd()));
  if (!google::protobuf::TextFormat::Print(composite_query, zero_copy_output.get())) {
    cerr << "Pharmacophore2d:cannot write\n";
    return 0;
  }

  return 1;
}
#endif

std::optional<SubstructureSearch::SubstructureQuery>
Options::Pharmacophore2d(Molecule& m, PerMoleculeData& pmd,
                         int number_functional_groups) {
  if (_verbose > 1) {
    cerr << m.name() << " has " << number_functional_groups << " functional groups\n";
  }

  // Force a smiles computation, and copy the atom order.
  m.smiles();
  resizable_array<atom_number_t> atom_order_in_smiles = m.atom_order_in_smiles();

  SubstructureSearch::SubstructureQuery result;
  result.set_name(m.name().data(), m.name().length());
  SubstructureSearch::SingleSubstructureQuery* query = result.add_query();
  query->set_respect_initial_atom_numbering(true);
  int group_number = 0;
  for (int i = 1; i <= number_functional_groups; ++i, ++group_number) {
    BuildQuery(m, pmd.functional_group(), i, atom_order_in_smiles, *query);
  }

  if (!AddDistanceConstraints(m, pmd, *query)) {
    return std::nullopt;
  }

  if (_proto.has_max_extra_atoms()) {
    auto* mpr = query->mutable_required_molecular_properties();
    mpr->set_max_natoms(m.natoms() + _proto.max_extra_atoms());
  }

  return result;
}

// Run the queries in _atoms_to_ignore and set `ignore_atom`.
int
Options::IdentifyAtomsToIgnore(Molecule& m, int* ignore_atom) {
  Molecule_to_Match target(&m);
  for (Substructure_Query* q : _atoms_to_ignore) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    sresults.each_embedding_set_vector(ignore_atom, 1);
  }

  return 1;  // Not sure what should be returned, ignored for now.
}

int
IdentifyFunctionalGroup(Molecule& m, const atom_number_t zatom, int* fg, int group_number,
                        const int* ignore) {
  fg[zatom] = group_number;
  int rc = 1;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    const atom_number_t other = b->other(zatom);
    if (ignore[other]) {
      continue;
    }
    if (fg[other] > 0) {  // Might be our group or another group.
      continue;
    }

    const atomic_number_t zother = m.atomic_number(other);
    if (b->is_single_bond() && zother == 6) {
      continue;
    }
    if (b->is_aromatic() && zother == 6) {
      continue;
    }
    // Always break at an aromatic ring.
    if (b->is_single_bond() && b->nrings() == 0 &&
        (m.is_aromatic(zatom) || m.is_aromatic(other))) {
      continue;
    }
    rc += IdentifyFunctionalGroup(m, other, fg, group_number, ignore);
  }

  return rc;
}

int
AnyMembersPositive(const Set_of_Atoms& embedding, const int* values) {
  for (atom_number_t a : embedding) {
    if (values[a] >= 0) {
      return 1;
    }
  }

  return 0;
}

void
Options::IdentifyExternallySpecified(Molecule& m, int* fg, const int* ignore,
                                     int& atoms_in_functional_groups, int& group_number) {
  Molecule_to_Match target(&m);

  uint32_t queries_matching = 0;

  // cerr << "Testing " << _external_functional_groups.size() << " functional groups\n";
  for (Substructure_Query* q : _external_functional_groups) {
    Substructure_Results sresults;
    if (!q->substructure_search(target, sresults)) {
      continue;
    }

    ++queries_matching;

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      if (AnyMembersPositive(*e, fg)) {
        continue;
      }
      if (e->any_members_set_in_array(ignore)) {
        continue;
      }
      ++group_number;
      e->set_vector(fg, group_number);
      atoms_in_functional_groups += e->number_elements();
    }
  }

  if (atoms_in_functional_groups == 0) {
    return;
  }

  // If we do not match all functional groups, appear as if nothing matched.
  if (_proto.all_functional_groups_must_match() &&
      queries_matching != _external_functional_groups.size()) {
    group_number = 0;
    atoms_in_functional_groups = 0;
  }
}

std::tuple<int, int>
Options::IdentifyFunctionalGroups(Molecule& m, int* fg, const int* ignore) {
  m.compute_aromaticity_if_needed();

  int atoms_in_functional_groups = 0;
  int group_number = 0;
  if (_external_functional_groups.size() > 0) {
    IdentifyExternallySpecified(m, fg, ignore, atoms_in_functional_groups, group_number);
    return {group_number, atoms_in_functional_groups};
  }

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (fg[i] > 0) {
      continue;
    }
    if (ignore[i]) {
      continue;
    }
    if (m.atomic_number(i) == 6) {
      continue;
    }
    group_number++;
    atoms_in_functional_groups += IdentifyFunctionalGroup(m, i, fg, group_number, ignore);
  }

  return {group_number, atoms_in_functional_groups};
}

int
Options::WriteLabelledSmiles(const Molecule& m, const int* functional_group) {
  Molecule mcopy(m);
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (functional_group[i] < 0) {
      continue;
    }
    if (_label_with_isotopes) {
      mcopy.set_isotope(i, functional_group[i] + 1);
    } else {
      mcopy.set_atom_map_number(i, functional_group[i] + 1);
    }
  }

  _stream_for_labelled_smiles << mcopy.smiles() << ' ' << m.name() << '\n';
  _stream_for_labelled_smiles.write_if_buffer_holds_more_than(8192);

  return 1;
}

std::optional<SubstructureSearch::SubstructureQuery>
Options::Process(Molecule& m) {
  ++_molecules_read;

  const int matoms = m.natoms();
  if (matoms == 0) {
    cerr << "Ignoring empty molecule " << m.name() << '\n';
    return std::nullopt;
  }

  // this is not implemented yet.
  std::unique_ptr<int[]> ignore_atoms(new_int(matoms));
  if (_atoms_to_ignore.number_elements() > 0) {
    IdentifyAtomsToIgnore(m, ignore_atoms.get());
  }

  PerMoleculeData pmd(m);

  const auto [number_functional_groups, atoms_in_functional_groups] =
      IdentifyFunctionalGroups(m, pmd.functional_group(), ignore_atoms.get());
  if (_verbose) {
    ++_functional_group_count[number_functional_groups];
    _fraction_atoms_in_functional_groups.extra(
        static_cast<float>(atoms_in_functional_groups) / static_cast<float>(matoms));
  }

  if (number_functional_groups == 0) {
    ++_no_functional_groups;
    return std::nullopt;
  }

  if (_stream_for_labelled_smiles.active()) {
    WriteLabelledSmiles(m, pmd.functional_group());
  }

  if (number_functional_groups == 1) {
    ++_only_one_functional_group;
    return std::nullopt;
  }

  return Pharmacophore2d(m, pmd, number_functional_groups);
}

#ifdef NO_LONGER_USED_JJJ
int
Options::Process(Molecule& m, IWString& fname) {
  ++_molecules_read;

  const int matoms = m.natoms();
  if (matoms == 0) {
    cerr << "Ignoring empty molecule " << m.name() << '\n';
    return 1;
  }

  std::unique_ptr<int[]> ignore_atoms(new_int(matoms));
  if (_atoms_to_ignore.number_elements() > 0) {
    IdentifyAtomsToIgnore(m, ignore_atoms.get());
  }

  PerMoleculeData pmd(m);

  const auto [number_functional_groups, atoms_in_functional_groups] =
      IdentifyFunctionalGroups(m, pmd.functional_group(), ignore_atoms.get());
  if (_verbose) {
    ++_functional_group_count[number_functional_groups];
    _fraction_atoms_in_functional_groups.extra(
        static_cast<float>(atoms_in_functional_groups) / static_cast<float>(matoms));
  }

  if (number_functional_groups == 0) {
    ++_no_functional_groups;
    return 0;
  }

  if (_stream_for_labelled_smiles.active()) {
    WriteLabelledSmiles(m, pmd.functional_group());
  }

  if (number_functional_groups == 1) {
    ++_only_one_functional_group;
    return 1;
  }

  return Pharmacophore2d(m, pmd, number_functional_groups, fname);
}
#endif

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _no_functional_groups << " molecules had no functional groups\n";
  output << _only_one_functional_group << " molecules had only one functional group\n";
  for (int i = 0; i < _functional_group_count.number_elements(); ++i) {
    if (_functional_group_count[i]) {
      output << _functional_group_count[i] << " molecules had " << i
             << " functional groups\n";
    }
  }

  return output.good();
}

int
Write(const SubstructureSearch::SubstructureQuery& qry, IWString& fname) {
  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Pharmacophore2d:cannot open '" << fname << "'\n";
    return 0;
  }

  using google::protobuf::io::FileOutputStream;
  using google::protobuf::io::ZeroCopyOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(
      new FileOutputStream(output.fd()));
  if (!google::protobuf::TextFormat::Print(qry, zero_copy_output.get())) {
    cerr << "Pharmacophore2d:cannot write\n";
    return 0;
  }

  return 1;
}

int
Pharmacophore2d(Options& options, data_source_and_type<Molecule>& input, Output& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    options.Preprocess(*m);

    std::optional<SubstructureSearch::SubstructureQuery> maybe_qry = options.Process(*m);
    if (!maybe_qry) {
      continue;
    }

    output.Write(*maybe_qry);
  }

  if (options.verbose()) {
    cerr << "Generated " << output.nfiles() << " valid files\n";
  }

  return 1;
}

int
Pharmacophore2d(Options& options, const char* fname, FileType input_type,
                Output& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  return Pharmacophore2d(options, input, output);
}

int
Pharmacophore2d(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:lS:td:D:nY:s:q:C:G:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  Options options;

  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  Output output;
  if (!output.Initialise(cl)) {
    cerr << "Cannot initialise output\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot process -i option\n";
      return 1;
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 4;
  }

  IWString output_stem;
  if (!cl.option_present('S')) {
    cerr << "Must specify the output stem (-S)\n";
    Usage(1);
  }

  cl.value('S', output_stem);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const auto* fname : cl) {
    if (!Pharmacophore2d(options, fname, input_type, output)) {
      cerr << "Error processing " << fname << "\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  if (cl.option_present('G')) {
    IWString fname = cl.string_value('G');
    output.WriteFilesGenerated(fname);
  }

  return 0;
}

}  // namespace pharmacophore_2d

int
main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return pharmacophore_2d::Pharmacophore2d(argc, argv);
}
