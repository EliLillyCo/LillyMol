#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <utility>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/xlogp.h"
#include "molecule_filter_lib.h"

namespace molecule_filter_lib {

using std::cerr;

namespace {

struct FeatureNameValue {
  std::string_view name;
  Feature feature;
};

constexpr std::array<FeatureNameValue, 35> kFeatureNameValues = {{
  {"natoms", Feature::kNatoms},
  {"amw", Feature::kAmw},
  {"hba_rdkit", Feature::kHbaRdkit},
  {"hbd_rdkit", Feature::kHbdRdkit},
  {"mw", Feature::kAmw},
  {"molecular_weight", Feature::kAmw},
  {"nrings", Feature::kNrings},
  {"heteroatom_count", Feature::kHeteroatomCount},
  {"heteroatoms", Feature::kHeteroatomCount},
  {"heteroatom_fraction", Feature::kHeteroatomFraction},
  {"aromatic_ring_count", Feature::kAromaticRingCount},
  {"aromatic_rings", Feature::kAromaticRingCount},
  {"aliphatic_ring_count", Feature::kAliphaticRingCount},
  {"aliphatic_rings", Feature::kAliphaticRingCount},
  {"rotatable_bonds", Feature::kRotatableBonds},
  {"rotbond", Feature::kRotatableBonds},
  {"max_ring_system_size", Feature::kMaxRingSystemSize},
  {"ring_system_size", Feature::kMaxRingSystemSize},
  {"aromatic_rings_in_system", Feature::kAromaticRingsInSystem},
  {"tpsa", Feature::kTpsa},
  {"alogp", Feature::kAlogp},
  {"xlogp", Feature::kXlogp},
  {"hba", Feature::kHba},
  {"hbd", Feature::kHbd},
  {"largest_ring_size", Feature::kLargestRingSize},
  {"halogen_count", Feature::kHalogenCount},
  {"halogens", Feature::kHalogenCount},
  {"max_distance", Feature::kMaxDistance},
  {"longest_path", Feature::kMaxDistance},
  {"sp3_carbon", Feature::kSp3Carbon},
  {"sp3_carbon_fraction", Feature::kSp3CarbonFraction},
  {"fraction_csp3", Feature::kSp3CarbonFraction},
  {"fcsp3", Feature::kSp3CarbonFraction},
  {"aromatic_density", Feature::kAromaticDensity},
  {"chiral", Feature::kChiral},
}};

}  // namespace

std::optional<Feature>
FeatureFromName(std::string_view name) {
  if (name == "number_fragments" || name == "nfrag" || name == "fragments") {
    return Feature::kNumberFragments;
  }

  if (name == "qed") {
    return Feature::kQed;
  }

  for (const FeatureNameValue& value : kFeatureNameValues) {
    if (name == value.name) {
      return value.feature;
    }
  }

  return std::nullopt;
}

std::string_view
FeatureName(Feature feature) {
  switch (feature) {
    case Feature::kNatoms:
      return "natoms";
    case Feature::kNrings:
      return "nrings";
    case Feature::kHeteroatomCount:
      return "heteroatom_count";
    case Feature::kHeteroatomFraction:
      return "heteroatom_fraction";
    case Feature::kAromaticRingCount:
      return "aromatic_ring_count";
    case Feature::kAliphaticRingCount:
      return "aliphatic_ring_count";
    case Feature::kRotatableBonds:
      return "rotatable_bonds";
    case Feature::kMaxRingSystemSize:
      return "max_ring_system_size";
    case Feature::kAromaticRingsInSystem:
      return "aromatic_rings_in_system";
    case Feature::kTpsa:
      return "tpsa";
    case Feature::kAlogp:
      return "alogp";
    case Feature::kXlogp:
      return "xlogp";
    case Feature::kHba:
      return "hba";
    case Feature::kHbd:
      return "hbd";
    case Feature::kLargestRingSize:
      return "largest_ring_size";
    case Feature::kHalogenCount:
      return "halogen_count";
    case Feature::kMaxDistance:
      return "max_distance";
    case Feature::kSp3Carbon:
      return "sp3_carbon";
    case Feature::kSp3CarbonFraction:
      return "sp3_carbon_fraction";
    case Feature::kAromaticDensity:
      return "aromatic_density";
    case Feature::kChiral:
      return "chiral";
    case Feature::kNumberFragments:
      return "number_fragments";
    case Feature::kQed:
      return "qed";
    case Feature::kAmw:
      return "amw";
    case Feature::kHbaRdkit:
      return "hba_rdkit";
    case Feature::kHbdRdkit:
      return "hbd_rdkit";
  }

  return "unknown";
}


int
Utility::BuildFromProto(const molecule_filter_data::Utility& proto) {
  _name.clear();
  _weight = 1.0f;
  _x.clear();
  _y.clear();

  if (! proto.has_name() || proto.name().empty()) {
    cerr << "Utility::BuildFromProto:missing name\n";
    return 0;
  }
  _name = proto.name();
  std::optional<Feature> feature = FeatureFromName(_name);
  if (! feature) {
    cerr << "Utility::BuildFromProto:unrecognised feature name '" << _name << "'\n";
    return 0;
  }
  _feature = *feature;

  if (proto.has_weight()) {
    if (! std::isfinite(proto.weight()) || proto.weight() <= 0.0f) {
      cerr << "Utility::BuildFromProto:invalid weight " << proto.weight() << "\n";
      return 0;
    }
    _weight = proto.weight();
  }

  if (proto.point_size() < 2) {
    cerr << "Utility::BuildFromProto:" << _name << " must have at least two points\n";
    return 0;
  }

  std::vector<std::pair<double, double>> points;
  points.reserve(proto.point_size());
  for (const molecule_filter_data::Point& point : proto.point()) {
    const double x = point.x();
    const double y = point.y();
    if (! std::isfinite(x)) {
      cerr << "Utility::BuildFromProto:" << _name << " non-finite x value\n";
      return 0;
    }
    if (! std::isfinite(y) || y < 0.0 || y > 1.0) {
      cerr << "Utility::BuildFromProto:" << _name << " invalid y value " << y << "\n";
      return 0;
    }
    points.emplace_back(x, y);
  }

  std::sort(points.begin(), points.end(),
            [](const auto& p1, const auto& p2) { return p1.first < p2.first; });

  _x.reserve(points.size());
  _y.reserve(points.size());
  for (const auto& [x, y] : points) {
    if (! _x.empty() && x == _x.back()) {
      cerr << "Utility::BuildFromProto:" << _name << " duplicate x value " << x << "\n";
      _x.clear();
      _y.clear();
      return 0;
    }
    _x.push_back(x);
    _y.push_back(y);
  }

  return 1;
}

double
Utility::Value(double value) const {
  if (_x.empty()) {
    return 0.0;
  }

  if (value <= _x.front()) {
    return _y.front();
  }
  if (value >= _x.back()) {
    return _y.back();
  }

  for (int i = 1; i < static_cast<int>(_x.size()); ++i) {
    if (value > _x[i]) {
      continue;
    }

    const double fraction = (value - _x[i - 1]) / (_x[i] - _x[i - 1]);
    return _y[i - 1] + fraction * (_y[i] - _y[i - 1]);
  }

  return _y.back();
}

FeatureValues::FeatureValues(Molecule& m, const int matoms, const int nrings,
                             FeatureCalculators& calculators) :
    _m(m),
    _matoms(matoms),
    _nrings(nrings),
    _calculators(calculators) {
}

int
FeatureValues::HeteroatomCount() {
  if (! _heteroatom_count) {
    _heteroatom_count = molecule_filter_lib::CountHeteroatoms(_m);
  }
  return *_heteroatom_count;
}

int
FeatureValues::AromaticRingCount() {
  if (! _aromatic_ring_count) {
    _aromatic_ring_count = molecule_filter_lib::AromaticRingCount(_m);
  }
  return *_aromatic_ring_count;
}

int
FeatureValues::RotatableBonds() {
  if (! _rotatable_bonds) {
    _rotatable_bonds = _calculators.rotbond.Process(_m);
  }
  return *_rotatable_bonds;
}

void
FeatureValues::ComputeRingSystem() {
  if (_max_ring_system_size && _aromatic_rings_in_system) {
    return;
  }

  const auto [max_ring_system_size, aromatic_rings_in_system] = MaxRingSystemSize(_m, _tmp);
  _max_ring_system_size = max_ring_system_size;
  _aromatic_rings_in_system = aromatic_rings_in_system;
}

std::optional<double>
FeatureValues::ALogP() {
  if (! _alogp_value) {
    std::optional<double> value = _calculators.alogp.LogP(_m);
    if (! value) {
      return std::nullopt;
    }
    _alogp_value = *value;
  }

  return _alogp_value;
}

std::optional<double>
FeatureValues::XLogP() {
  if (! _xlogp_value) {
    if (! _tmp) {
      _tmp = std::make_unique<int[]>(_matoms);
    }
    std::fill_n(_tmp.get(), _matoms, 0);
    std::optional<double> value = _calculators.xlogp.LogP(_m, _tmp.get());
    if (! value) {
      return std::nullopt;
    }
    _xlogp_value = *value;
  }

  return _xlogp_value;
}

void
FeatureValues::ComputeHbaHbd() {
  if (_hba && _hbd) {
    return;
  }

  int hba;
  int hbd;
  RuleOfFive(_m, hba, hbd);
  _hba = hba;
  _hbd = hbd;
}

int
FeatureValues::HalogenCount() {
  if (! _halogen_count) {
    _halogen_count = molecule_filter_lib::HalogenCount(_m);
  }
  return *_halogen_count;
}

int
FeatureValues::MaxDistance() {
  if (! _max_distance) {
    _max_distance = _m.longest_path();
  }
  return *_max_distance;
}

void
FeatureValues::ComputeSp3Carbon() {
  if (_sp3_carbon && _sp3_carbon_fraction) {
    return;
  }

  const int csp3 = molecule_filter_lib::Sp3Carbon(_m);
  _sp3_carbon = csp3;

  const int carbon = _m.natoms(6);
  if (carbon == 0) {
    _sp3_carbon_fraction = 0.0;
  } else {
    _sp3_carbon_fraction = iwmisc::Fraction<double>(csp3, carbon);
  }
}

int
FeatureValues::Sp3Carbon() {
  ComputeSp3Carbon();
  return *_sp3_carbon;
}

double
FeatureValues::Sp3CarbonFraction() {
  ComputeSp3Carbon();
  return *_sp3_carbon_fraction;
}

int
FeatureValues::NumberFragments() {
  if (! _number_fragments) {
    _number_fragments = _m.number_fragments();
  }
  return *_number_fragments;
}

// Average molecular weight, implicit hydrogens included.
// molecular_weight() refuses to work on molecules containing isotopes, and
// writes to stderr when it encounters one. Since a filter is typically run
// across arbitrary collections, use the variant that quietly treats an
// isotopic atom as its natural abundance form.
double
FeatureValues::Amw() {
  if (! _amw) {
    _amw = _m.molecular_weight_ignore_isotopes();
  }
  return *_amw;
}

// RDKit compatible hydrogen bond counts. Note that these are NOT the Lipinski
// counts returned by kHba and kHbd - they are a different, more restrictive
// definition that also counts sulfur. See Molecule::RDKitNumHAcceptors.
void
FeatureValues::ComputeRdkitHbaHbd() {
  if (_hba_rdkit && _hbd_rdkit) {
    return;
  }

  int hba;
  int hbd;
  _m.RDKitHbaHbd(hba, hbd);
  _hba_rdkit = hba;
  _hbd_rdkit = hbd;
}

std::optional<double>
FeatureValues::Qed() {
  if (_qed) {
    return _qed;
  }
  // No defensive copy needed; Qed::qed protects the molecule itself.
  std::optional<float> value = _calculators.qed.qed(_m);
  if (! value) {
    return std::nullopt;
  }
  _qed = *value;
  return _qed;
}

std::optional<double>
FeatureValues::Value(Feature feature) {
  switch (feature) {
    case Feature::kNatoms:
      return _matoms;
    case Feature::kNrings:
      return _nrings;
    case Feature::kHeteroatomCount:
      return HeteroatomCount();
    case Feature::kHeteroatomFraction:
      return iwmisc::Fraction<double>(HeteroatomCount(), _matoms);
    case Feature::kAromaticRingCount:
      return AromaticRingCount();
    case Feature::kAliphaticRingCount:
      return _nrings - AromaticRingCount();
    case Feature::kRotatableBonds:
      return RotatableBonds();
    case Feature::kMaxRingSystemSize:
      ComputeRingSystem();
      return *_max_ring_system_size;
    case Feature::kAromaticRingsInSystem:
      ComputeRingSystem();
      return *_aromatic_rings_in_system;
    case Feature::kTpsa:
      if (! _tpsa) {
        _tpsa = _calculators.tpsa.PolarSurfaceArea(_m);
      }
      return _tpsa;
    case Feature::kAlogp:
      return ALogP();
    case Feature::kXlogp:
      return XLogP();
    case Feature::kHba:
      ComputeHbaHbd();
      return *_hba;
    case Feature::kHbd:
      ComputeHbaHbd();
      return *_hbd;
    case Feature::kLargestRingSize:
      if (_nrings == 0) {
        return 0.0;
      }
      return _m.ringi(_nrings - 1)->number_elements();
    case Feature::kHalogenCount:
      return HalogenCount();
    case Feature::kMaxDistance:
      return MaxDistance();
    case Feature::kSp3Carbon:
      return Sp3Carbon();
    case Feature::kSp3CarbonFraction:
      return Sp3CarbonFraction();
    case Feature::kAromaticDensity:
      return iwmisc::Fraction<double>(_m.aromatic_atom_count(), _matoms);
    case Feature::kChiral:
      return _m.chiral_centres();
    case Feature::kNumberFragments:
      return NumberFragments();
    case Feature::kQed:
      return Qed();
    case Feature::kAmw:
      return Amw();
    case Feature::kHbaRdkit:
      ComputeRdkitHbaHbd();
      return *_hba_rdkit;
    case Feature::kHbdRdkit:
      ComputeRdkitHbaHbd();
      return *_hbd_rdkit;
  }

  return std::nullopt;
}

MoleculeFilter::MoleculeFilter() {
  _active = false;
  _qed_initialised = false;
}

int
MoleculeFilter::Build(IWString& fname) {
  std::optional<molecule_filter_data::Requirements> maybe_proto =
                iwmisc::ReadTextProto<molecule_filter_data::Requirements>(fname);
  if (! maybe_proto) {
    cerr << "MoleculeFilter::Build:cannot read '" << fname << "'\n";
    _active = false;
    _utilities.clear();
    return 0;
  }

  return Build(*maybe_proto);
}

int
MoleculeFilter::Build(const molecule_filter_data::Requirements& proto) {
  _active = false;
  _utilities.clear();

  _requirements = proto;

  InitialiseOptionalFeatures();

  if (! BuildUtilities()) {
    return 0;
  }

  if (! InitialiseQEDIfNeeded()) {
    return 0;
  }

  _active = true;

  return 1;
}

FeatureCalculators
MoleculeFilter::MakeCalculators() {
  return FeatureCalculators{_rotbond, _alogp, _xlogp, _tpsa, _qed};
}

int
MoleculeFilter::EvaluateUtilities(Molecule& m, const int matoms, const int nrings,
                                  std::vector<double>& per_feature_utility,
                                  double& overall_utility) {
  FeatureCalculators calculators = MakeCalculators();
  FeatureValues feature_values(m, matoms, nrings, calculators);

  return EvaluateUtilities(feature_values, per_feature_utility, overall_utility);
}

int
MoleculeFilter::EvaluateUtilities(FeatureValues& feature_values,
                                  std::vector<double>& per_feature_utility,
                                  double& overall_utility) {
  per_feature_utility.clear();
  overall_utility = 0.0;

  if (_utilities.empty()) {
    return 1;
  }

  per_feature_utility.reserve(_utilities.size());

  double weighted_sum = 0.0;
  double weight_sum = 0.0;
  double product = 1.0;
  double minval = 1.0;
  double maxval = 0.0;

  for (const Utility& utility : _utilities) {
    std::optional<double> raw_value = feature_values.Value(utility.feature());
    if (! raw_value) {
      cerr << "MoleculeFilter::EvaluateUtilities:cannot compute feature '" <<
              utility.name() << "'\n";
      per_feature_utility.clear();
      overall_utility = 0.0;
      return 0;
    }

    const double u = utility.Value(*raw_value);
    per_feature_utility.push_back(u);

    weighted_sum += utility.weight() * u;
    weight_sum += utility.weight();
    product *= u;
    if (u < minval) {
      minval = u;
    }
    if (u > maxval) {
      maxval = u;
    }
  }

  const molecule_filter_data::UtilityCombination combination =
      _requirements.has_utility_combination()
          ? _requirements.utility_combination()
          : molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_AVERAGE;

  switch (combination) {
    case molecule_filter_data::UTILITY_COMBINATION_UNSPECIFIED:
    case molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_AVERAGE:
      overall_utility = weighted_sum / weight_sum;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_SUM:
      overall_utility = weighted_sum;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_PRODUCT:
      overall_utility = product;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_MIN:
      overall_utility = minval;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_MAX:
      overall_utility = maxval;
      return 1;
    default:
      cerr << "MoleculeFilter::EvaluateUtilities:unrecognised combination " <<
              static_cast<int>(combination) << "\n";
      return 0;
  }
}

bool
MoleculeFilter::UsesFeature(Feature feature) const {
  for (const Utility& utility : _utilities) {
    if (utility.feature() == feature) {
      return true;
    }
  }

  return false;
}

int
MoleculeFilter::InitialiseQEDIfNeeded() {
  if (! UsesFeature(Feature::kQed)) {
    return 1;
  }
  if (_qed_initialised) {
    return 1;
  }

  if (! _qed.InitialiseFromEnvironment()) {
    cerr << "MoleculeFilter::InitialiseQEDIfNeeded:cannot initialise QED\n";
    _utilities.clear();
    return 0;
  }

  _qed_initialised = true;
  return 1;
}

int
MoleculeFilter::BuildUtilities() {
  const int nutility = _requirements.utility_size();
  if (nutility == 0) {
    return 1;
  }

  _utilities.reserve(nutility);
  for (const molecule_filter_data::Utility& proto : _requirements.utility()) {
    Utility utility;
    if (! utility.BuildFromProto(proto)) {
      _utilities.clear();
      return 0;
    }
    _utilities.push_back(std::move(utility));
  }

  return 1;
}

void
MoleculeFilter::InitialiseOptionalFeatures() {
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  _tpsa.set_display_psa_unclassified_atom_messages(0);
  _xlogp.SetIssueUnclassifiedAtomMessages(false);

  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_apply_zwitterion_correction(1);
}

int
CountHeteroatoms(const Molecule& m) {
  int rc = 0;
  m.each_atom_lambda([&rc](const Atom& a) {
    if (a.atomic_number() != 6) {
      ++rc;
    }
  });

  return rc;
}

int
AromaticRingCount(Molecule& m) {
  m.compute_aromaticity_if_needed();

  int rc = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      ++rc;
    }
  }

  return rc;
}

// Return true if we examine multiple fragments, which
// means that `largest_frag` will be different from `smiles`.
// Fragments are compared on heavy atom count. An explicit Hydrogen is not a
// reason for one fragment to be considered larger than another.
bool
LargestFragment(const const_IWSubstring& smiles,
                const_IWSubstring& largest_frag,
                int& natoms, int& nrings, int& explicit_hydrogens) {
  int max_atoms = 0;
  int i = 0;
  const_IWSubstring token;
  int fragments_examined = 0;
  explicit_hydrogens = 0;
  while (smiles.nextword(token, i, '.')) {
    ++fragments_examined;
    int nri;
    int eh;
    const int nat = lillymol::count_atoms_in_smiles(token, nri, eh) - eh;
    if (nat > max_atoms) {
      max_atoms = nat;
      nrings = nri;
      largest_frag = token;
      explicit_hydrogens = eh;
    }
  }

  natoms = max_atoms;

  if (fragments_examined == 1) {
    return false;
  }
  return true;
}

std::tuple<int, int>
MaxRingSystemSize(Molecule& m, std::unique_ptr<int[]>& tmp) {
  const int matoms = m.natoms();

  m.compute_aromaticity_if_needed();

  if (! tmp) {
    tmp.reset(new int[matoms]);
  }
  std::fill_n(tmp.get(), matoms, 0);

  const int nrings = m.nrings();

  std::unique_ptr<int[]> ring_already_done = std::make_unique<int[]>(nrings);

  int max_system_size = 0;
  int max_aromatic_rings_in_system = 0;
  for (int i = 0; i < nrings; ++i) {
    if (ring_already_done[i]) {
      continue;
    }
    const Ring* ri = m.ringi(i);
    if (! ri->is_fused()) {
      continue;
    }

    int system_size = 1;
    int aromatic_rings_in_system;
    if (ri->is_aromatic()) {
      aromatic_rings_in_system = 1;
    } else {
      aromatic_rings_in_system = 0;
    }


    for (int j = i + 1; j < nrings; ++j) {
      if (ring_already_done[j]) {
        continue;
      }

      ring_already_done[j] = 1;
      const Ring* rj = m.ringi(j);
      if (ri->fused_system_identifier() == rj->fused_system_identifier()) {
        ++system_size;
        if (rj->is_aromatic()) {
          ++aromatic_rings_in_system;
        }
      }
    }
    if (system_size > max_system_size) {
      max_system_size = system_size;
    }
    if (aromatic_rings_in_system > max_aromatic_rings_in_system) {
      max_aromatic_rings_in_system = aromatic_rings_in_system;
    }
  }

  return std::make_tuple(max_system_size, max_aromatic_rings_in_system);
}

// Retained for compatibility with existing callers. The rule of five counts
// live in Molecule_Lib - this used to be a private copy that had drifted from
// the other copies in iwdescr and the python bindings. Do not reintroduce a
// local implementation here.
void
RuleOfFive(Molecule & m, int& acceptor, int& donor) {
  m.LipinskiHbaHbd(acceptor, donor);
}

int
HalogenCount(const Molecule& m) {
  static std::vector<int> halogen = {
    0,  // 0
    0,  // 1
    0,  // 2
    0,  // 3
    0,  // 4
    0,  // 5
    0,  // 6
    0,  // 7
    0,  // 8
    1,  // 9
    0,  // 10
    0,  // 11
    0,  // 12
    0,  // 13
    0,  // 14
    0,  // 15
    0,  // 16
    1,  // 17
    0,  // 18
    0,  // 19
    0,  // 20
    0,  // 21
    0,  // 22
    0,  // 23
    0,  // 24
    0,  // 25
    0,  // 26
    0,  // 27
    0,  // 28
    0,  // 29
    0,  // 30
    0,  // 31
    0,  // 32
    0,  // 33
    0,  // 34
    1,  // 35
    0,  // 36
    0,  // 37
    0,  // 38
    0,  // 39
    0,  // 40
    0,  // 41
    0,  // 42
    0,  // 43
    0,  // 44
    0,  // 45
    0,  // 46
    0,  // 47
    0,  // 48
    0,  // 49
    0,  // 50
    0,  // 51
    0,  // 52
    1   // 53
  };

  int rc = 0;

  for (const Atom* a : m) {
    const uint32_t z = a->atomic_number();
    if (z < halogen.size()) {
      rc += halogen[z];
    }
  }

  return rc;
}

int 
Sp3Carbon(Molecule & m) {
  int rc = 0;

  for (const Atom* a : m) {
    if (a->atomic_number() != 6) {
      continue;
    }
    if (a->fully_saturated()) {
      ++rc;
    }
  }

  return rc;
}

int
Reject(RejectionReason& rejection_reason, RejectionReason reason) {
  rejection_reason = reason;
  return 0;
}

int
MoleculeFilter::Ok(Molecule& m) {
  RejectionReason rejection_reason;
  return Ok(m, rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, RejectionReason& rejection_reason) {
  return Ok(m, m.natoms(), m.nrings(), rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, const int matoms, const int nrings) {
  RejectionReason rejection_reason;
  return Ok(m, matoms, nrings, rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, const int matoms, const int nrings,
                   RejectionReason& rejection_reason) {
  FeatureCalculators calculators = MakeCalculators();
  FeatureValues feature_values(m, matoms, nrings, calculators);

  return Ok(feature_values, m, matoms, nrings, rejection_reason);
}

int
MoleculeFilter::OkAndUtilities(Molecule& m, const int matoms, const int nrings,
                               RejectionReason& rejection_reason,
                               std::vector<double>& per_feature_utility,
                               double& overall_utility) {
  per_feature_utility.clear();
  overall_utility = 0.0;

  // One FeatureValues spanning both phases. Everything the filters compute is
  // still in the cache when the utilities ask for it.
  FeatureCalculators calculators = MakeCalculators();
  FeatureValues feature_values(m, matoms, nrings, calculators);

  if (! Ok(feature_values, m, matoms, nrings, rejection_reason)) {
    return 0;
  }

  if (! EvaluateUtilities(feature_values, per_feature_utility, overall_utility)) {
    return -1;
  }

  return 1;
}

// The tests below run in increasing order of cost, using the relative_cost
// figures recorded in molecule_filter_features.json, so that a molecule is
// rejected as early and as cheaply as possible. A benchmark reproduced that
// scale closely - rotatable bonds measured 12.7 against the recorded 12.
//
// Reordering cannot change WHICH molecules pass, only which requirement gets
// the blame when a molecule fails more than one. That changes the verbose
// rejection statistics, so if you reorder again, expect those to move and
// check that the set of passing molecules does not.
//
// Blocks that share a descriptor stay together - heteroatom count with
// heteroatom fraction, hba with hbd, sp3 carbon with its fraction, ring system
// size with aromatic rings in system - since the second of each pair is free
// once the first has been computed.
//
// Descriptor values come from `feature_values` rather than being computed
// here, so that a caller doing both filtering and utility evaluation pays for
// each descriptor once. `m` is still needed for the handful of properties that
// are not Features - organic, isotopes, carbon count, aromatic atom count and
// ring size.
namespace {

/*
  Comparing a computed value against a floating point threshold is not as
  simple as < or >.

  The threshold arrives as a 32 bit float from the proto, while values are
  computed as doubles, so the number the config author wrote is perturbed by
  about 1e-7 relative before it is ever compared - and in a direction that
  varies with the number. Written as doubles, (float)0.2 is 0.20000000298 but
  (float)0.35 is 0.34999999404.

  On top of that the fraction features - heteroatom_fraction,
  sp3_carbon_fraction, aromatic_density - are ratios of small integers, so
  they land exactly on the round numbers people write. One sp3 carbon in five
  is 0.2, and whether that is inside or outside a 0.2 bound should not depend
  on binary representation.

  So min_ and max_ are inclusive within a tolerance. A value that is
  conceptually equal to the threshold is kept. The tolerance is relative, with
  a floor of one so that thresholds at or near zero get an absolute 1e-6. That
  is three orders above the float32 representation error and far below any
  difference that means anything chemically.

  Integer valued requirements do not go through here. For those, exact
  comparison is both correct and what people expect - max_natoms 50 keeps a 50
  atom molecule and rejects a 51 atom one.
*/

constexpr double kThresholdTolerance = 1.0e-6;

inline double
ThresholdSlack(double threshold) {
  return kThresholdTolerance * std::max(1.0, std::fabs(threshold));
}

// True if `value` is below `minval` by more than the tolerance.
inline bool
BelowMin(double value, float minval) {
  const double m = minval;
  return value < m - ThresholdSlack(m);
}

// True if `value` is above `maxval` by more than the tolerance.
inline bool
AboveMax(double value, float maxval) {
  const double m = maxval;
  return value > m + ThresholdSlack(m);
}

}  // namespace

int
MoleculeFilter::Ok(FeatureValues& feature_values, Molecule& m,
                   const int matoms, const int nrings,
                   RejectionReason& rejection_reason) {
  rejection_reason = RejectionReason::kPass;

  // Every integer valued Feature is always computable, so dereferencing is
  // safe. The three that can fail - tpsa, alogp and xlogp - are handled
  // individually below, and note that their failure policies differ.
  auto ivalue = [&feature_values](Feature feature) {
    return static_cast<int>(*feature_values.Value(feature));
  };

  if (matoms == 0) {
    return Reject(rejection_reason, RejectionReason::kZeroAtoms);
  }

  if (_requirements.has_min_natoms() && matoms < _requirements.min_natoms()) {
    return Reject(rejection_reason, RejectionReason::kTooFewAtoms);
  }

  if (_requirements.has_max_natoms() && matoms > _requirements.max_natoms()) {
    return Reject(rejection_reason, RejectionReason::kTooManyAtoms);
  }

  if (_requirements.has_min_nrings() && nrings < _requirements.min_nrings()) {
    return Reject(rejection_reason, RejectionReason::kTooFewRings);
  }

  if (_requirements.has_max_nrings() && nrings > _requirements.max_nrings()) {
    return Reject(rejection_reason, RejectionReason::kTooManyRings);
  }

  if (_requirements.has_exclude_isotopes() && m.number_isotopic_atoms() > 0) {
    return Reject(rejection_reason, RejectionReason::kIsotope);
  }

  if (_requirements.has_min_chiral() || _requirements.has_max_chiral()) {
    const int chiral_centres = ivalue(Feature::kChiral);
    if (_requirements.has_min_chiral() && chiral_centres < _requirements.min_chiral()) {
      return Reject(rejection_reason, RejectionReason::kTooFewChiral);
    }
    if (_requirements.has_max_chiral() && chiral_centres > _requirements.max_chiral()) {
      return Reject(rejection_reason, RejectionReason::kTooManyChiral);
    }
  }

  if (_requirements.has_min_heteroatom_count() ||
      _requirements.has_max_heteroatom_count() ||
      _requirements.has_min_heteroatom_fraction() ||
      _requirements.has_max_heteroatom_fraction()) {
    const int hac = ivalue(Feature::kHeteroatomCount);

    if (_requirements.has_min_heteroatom_count() && hac < _requirements.min_heteroatom_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHeteroatoms);
    }
    if (_requirements.has_max_heteroatom_count() && hac > _requirements.max_heteroatom_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHeteroatoms);
    }

    if (_requirements.has_min_heteroatom_fraction() ||
        _requirements.has_max_heteroatom_fraction()) {
      const double haf = *feature_values.Value(Feature::kHeteroatomFraction);
      if (_requirements.has_min_heteroatom_fraction() &&
          BelowMin(haf, _requirements.min_heteroatom_fraction())) {
        return Reject(rejection_reason, RejectionReason::kMinHeteroatomFraction);
      }
      if (_requirements.has_max_heteroatom_fraction() &&
          AboveMax(haf, _requirements.max_heteroatom_fraction())) {
        return Reject(rejection_reason, RejectionReason::kMaxHeteroatomFraction);
      }
    }
  }

  if (_requirements.has_exclude_non_organic() && ! m.organic_only()) {
    return Reject(rejection_reason, RejectionReason::kNonOrganic);
  }

  if (_requirements.has_min_halogen_count() || _requirements.has_max_halogen_count()) {
    const int h = ivalue(Feature::kHalogenCount);
    if (_requirements.has_min_halogen_count() && h < _requirements.min_halogen_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHalogen);
    }
    if (_requirements.has_max_halogen_count() && h > _requirements.max_halogen_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHalogen);
    }
  }

  if (_requirements.has_min_amw() || _requirements.has_max_amw()) {
    const double amw = *feature_values.Value(Feature::kAmw);
    if (_requirements.has_min_amw() && BelowMin(amw, _requirements.min_amw())) {
      return Reject(rejection_reason, RejectionReason::kLowAmw);
    }
    if (_requirements.has_max_amw() && AboveMax(amw, _requirements.max_amw())) {
      return Reject(rejection_reason, RejectionReason::kHighAmw);
    }
  }

  if (_requirements.has_min_largest_ring_size() ||
      _requirements.has_max_largest_ring_size()) {
    if (nrings == 0) {
      if (_requirements.has_min_largest_ring_size() &&
          _requirements.min_largest_ring_size() > 0) {
        return Reject(rejection_reason, RejectionReason::kRingTooSmall);
      }
    } else {
      const int rsze = ivalue(Feature::kLargestRingSize);
      if (_requirements.has_min_largest_ring_size() &&
          rsze < _requirements.min_largest_ring_size()) {
        return Reject(rejection_reason, RejectionReason::kRingTooSmall);
      }
      if (_requirements.has_max_largest_ring_size() &&
          rsze > _requirements.max_largest_ring_size()) {
        return Reject(rejection_reason, RejectionReason::kRingTooLarge);
      }
    }
  }

  if (_requirements.has_min_sp3_carbon() || _requirements.has_max_sp3_carbon() ||
      _requirements.has_min_sp3_carbon_fraction() ||
      _requirements.has_max_sp3_carbon_fraction()) {
    const int csp3 = ivalue(Feature::kSp3Carbon);
    if (_requirements.has_min_sp3_carbon() && csp3 < _requirements.min_sp3_carbon()) {
      return Reject(rejection_reason, RejectionReason::kTooFewSp3Carbon);
    }
    if (_requirements.has_max_sp3_carbon() && csp3 > _requirements.max_sp3_carbon()) {
      return Reject(rejection_reason, RejectionReason::kTooManySp3Carbon);
    }

    if (_requirements.has_min_sp3_carbon_fraction() ||
        _requirements.has_max_sp3_carbon_fraction()) {
      const double fcsp3 = *feature_values.Value(Feature::kSp3CarbonFraction);
      if (_requirements.has_min_sp3_carbon_fraction() &&
          BelowMin(fcsp3, _requirements.min_sp3_carbon_fraction())) {
        return Reject(rejection_reason, RejectionReason::kSp3CarbonFractionTooLow);
      }
      if (_requirements.has_max_sp3_carbon_fraction() &&
          AboveMax(fcsp3, _requirements.max_sp3_carbon_fraction())) {
        return Reject(rejection_reason, RejectionReason::kSp3CarbonFractionTooHigh);
      }
    }
  }

  if (_requirements.has_min_hba() || _requirements.has_max_hba() ||
      _requirements.has_min_hbd() || _requirements.has_max_hbd()) {
    const int hba = ivalue(Feature::kHba);
    const int hbd = ivalue(Feature::kHbd);
    if (_requirements.has_min_hba() && hba < _requirements.min_hba()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHba);
    }
    if (_requirements.has_max_hba() && hba > _requirements.max_hba()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHba);
    }
    if (_requirements.has_min_hbd() && hbd < _requirements.min_hbd()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHbd);
    }
    if (_requirements.has_max_hbd() && hbd > _requirements.max_hbd()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHbd);
    }
  }

  if (_requirements.has_min_aromatic_ring_count() ||
      _requirements.has_max_aromatic_ring_count() ||
      _requirements.has_min_aromatic_rings_in_system() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    // Check nrings first. If not enough rings if all were aromatic...
    if (_requirements.has_min_aromatic_ring_count() && nrings < _requirements.min_aromatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAromaticRings);
    }

    const int arc = ivalue(Feature::kAromaticRingCount);
    if (_requirements.has_min_aromatic_ring_count() && arc < _requirements.min_aromatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAromaticRings);
    }
    if (_requirements.has_max_aromatic_ring_count() && arc > _requirements.max_aromatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyAromaticRings);
    }
  }

  if (_requirements.has_min_aliphatic_ring_count() ||
      _requirements.has_max_aliphatic_ring_count()) {
    if (_requirements.has_min_aliphatic_ring_count() && nrings < _requirements.min_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAliphaticRings);
    }
    const int alring = ivalue(Feature::kAliphaticRingCount);
    if (_requirements.has_min_aliphatic_ring_count() && alring < _requirements.min_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAliphaticRings);
    }
    if (_requirements.has_max_aliphatic_ring_count() && alring > _requirements.max_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyAliphaticRings);
    }
  }

  if (_requirements.has_min_aromatic_density() || _requirements.has_max_aromatic_density()) {
    const double aromdens = *feature_values.Value(Feature::kAromaticDensity);
    if (_requirements.has_min_aromatic_density() &&
        BelowMin(aromdens, _requirements.min_aromatic_density())) {
      return Reject(rejection_reason, RejectionReason::kAromaticDensityTooLow);
    }
    if (_requirements.has_max_aromatic_density() &&
        AboveMax(aromdens, _requirements.max_aromatic_density())) {
      return Reject(rejection_reason, RejectionReason::kAromaticDensityTooHigh);
    }
  }

  if (_requirements.has_min_hba_rdkit() || _requirements.has_max_hba_rdkit() ||
      _requirements.has_min_hbd_rdkit() || _requirements.has_max_hbd_rdkit()) {
    const int hba = ivalue(Feature::kHbaRdkit);
    const int hbd = ivalue(Feature::kHbdRdkit);
    if (_requirements.has_min_hba_rdkit() && hba < _requirements.min_hba_rdkit()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHbaRdkit);
    }
    if (_requirements.has_max_hba_rdkit() && hba > _requirements.max_hba_rdkit()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHbaRdkit);
    }
    if (_requirements.has_min_hbd_rdkit() && hbd < _requirements.min_hbd_rdkit()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHbdRdkit);
    }
    if (_requirements.has_max_hbd_rdkit() && hbd > _requirements.max_hbd_rdkit()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHbdRdkit);
    }
  }

  if (_requirements.has_min_rotatable_bonds() || _requirements.has_max_rotatable_bonds()) {
    const int rotb = ivalue(Feature::kRotatableBonds);
    if (_requirements.has_min_rotatable_bonds() && rotb < _requirements.min_rotatable_bonds()) {
      return Reject(rejection_reason, RejectionReason::kTooFewRotatableBonds);
    }
    if (_requirements.has_max_rotatable_bonds() && rotb > _requirements.max_rotatable_bonds()) {
      return Reject(rejection_reason, RejectionReason::kTooManyRotatableBonds);
    }
  }

  if (_requirements.has_min_ring_system_size() ||
      _requirements.has_max_ring_system_size() ||
      _requirements.has_min_aromatic_rings_in_system() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    bool need_ring_system = true;
    if (! _requirements.has_min_ring_system_size() &&
        ! _requirements.has_min_aromatic_rings_in_system()) {
      need_ring_system = false;
      if (_requirements.has_max_ring_system_size() &&
          nrings >= _requirements.max_ring_system_size()) {
        need_ring_system = true;
      }
      if (_requirements.has_max_aromatic_rings_in_system() &&
          nrings >= _requirements.max_aromatic_rings_in_system()) {
        need_ring_system = true;
      }
    }

    if (need_ring_system) {
      const int max_ring_system_size = ivalue(Feature::kMaxRingSystemSize);
      const int max_aromatic_rings_in_system = ivalue(Feature::kAromaticRingsInSystem);
      if (_requirements.has_min_ring_system_size() &&
          max_ring_system_size < _requirements.min_ring_system_size()) {
        return Reject(rejection_reason, RejectionReason::kRingSystemTooSmall);
      }
      if (_requirements.has_max_ring_system_size() &&
          max_ring_system_size > _requirements.max_ring_system_size()) {
        return Reject(rejection_reason, RejectionReason::kRingSystemTooLarge);
      }
      if (_requirements.has_min_aromatic_rings_in_system() &&
          max_aromatic_rings_in_system < _requirements.min_aromatic_rings_in_system()) {
        return Reject(rejection_reason, RejectionReason::kTooFewAromaticRingsInSystem);
      }
      if (_requirements.has_max_aromatic_rings_in_system() &&
          max_aromatic_rings_in_system > _requirements.max_aromatic_rings_in_system()) {
        return Reject(rejection_reason, RejectionReason::kTooManyAromaticRingsInSystem);
      }
    }
  }

  if (_requirements.has_min_tpsa() || _requirements.has_max_tpsa()) {
    // Note that an uncomputable tpsa is a rejection, whereas an uncomputable
    // logp below is not. Preserving long standing behaviour.
    std::optional<double> tpsa = feature_values.Value(Feature::kTpsa);
    if (! tpsa) {
      return Reject(rejection_reason, RejectionReason::kLowTpsa);
    }
    if (_requirements.has_min_tpsa() && BelowMin(*tpsa, _requirements.min_tpsa())) {
      return Reject(rejection_reason, RejectionReason::kLowTpsa);
    }
    if (_requirements.has_max_tpsa() && AboveMax(*tpsa, _requirements.max_tpsa())) {
      return Reject(rejection_reason, RejectionReason::kHighTpsa);
    }
  }

  if (_requirements.has_min_distance() || _requirements.has_max_distance()) {
    if (! _requirements.has_min_distance() &&
        _requirements.has_max_distance() && matoms <= _requirements.max_distance()) {
      // No need to compute; longest path cannot exceed atom count.
    } else {
      const int d = ivalue(Feature::kMaxDistance);
      if (_requirements.has_min_distance() && d < _requirements.min_distance()) {
        return Reject(rejection_reason, RejectionReason::kTooShort);
      }
      if (_requirements.has_max_distance() && d > _requirements.max_distance()) {
        return Reject(rejection_reason, RejectionReason::kTooLong);
      }
    }
  }

  // A molecule whose logp cannot be computed passes, it is not rejected.
  if (_requirements.has_min_alogp() || _requirements.has_max_alogp()) {
    std::optional<double> x = feature_values.Value(Feature::kAlogp);
    if (x) {
      if (_requirements.has_min_alogp() && BelowMin(*x, _requirements.min_alogp())) {
        return Reject(rejection_reason, RejectionReason::kLowAlogp);
      }
      if (_requirements.has_max_alogp() && AboveMax(*x, _requirements.max_alogp())) {
        return Reject(rejection_reason, RejectionReason::kHighAlogp);
      }
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_xlogp()) {
    std::optional<double> x = feature_values.Value(Feature::kXlogp);
    if (x) {
      if (_requirements.has_min_xlogp() && BelowMin(*x, _requirements.min_xlogp())) {
        return Reject(rejection_reason, RejectionReason::kLowXlogp);
      }
      if (_requirements.has_max_xlogp() && AboveMax(*x, _requirements.max_xlogp())) {
        return Reject(rejection_reason, RejectionReason::kHighXlogp);
      }
    }
  }

  return 1;
}


}  // namespace molecule_filter_lib
