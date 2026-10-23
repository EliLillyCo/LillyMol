#include "Utilities/GFP_Tools/gfp_context.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/atom_pair_fingerprint.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/iwecfp_lib.h"
#include "Molecule_Tools/jwcats_lib.h"
#include "Molecule_Tools/maccskeys_fn5.h"
#include "Molecule_Tools/mformula.h"
#include "Molecule_Tools/mpr.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/ring_substitution.h"
#include "Molecule_Tools/xlogp.h"
#include "Utilities/GFP_Tools/dyfp.h"
#include "Utilities/GFP_Tools/molecular_property_ratio.h"

namespace gfp_context {
namespace {
using std::cerr;

constexpr char kSmilesTag[] = "$SMI<";
constexpr char kIdentifierTag[] = "PCN<";
constexpr char kMolecularPropertiesTag[] = "MPR<";
constexpr int kMaccskeysNbits = 3 * 64;
constexpr int kDefaultSparseReplicates = 9;
constexpr float kDefaultFormulaLogS = 200.0f;
constexpr float kDefaultFormulaLogD = 41.0f;
constexpr int kFormulaFingerprintNbits = mformula::kMFOther + 1;

uint64_t
Fnv1a(uint64_t hash, const void* v, int n) {
  const unsigned char* data = reinterpret_cast<const unsigned char*>(v);
  for (int i = 0; i < n; ++i) {
    hash ^= data[i];
    hash *= 1099511628211ULL;
  }

  return hash;
}

uint64_t
Fnv1a(uint64_t hash, const IWString& s) {
  return Fnv1a(hash, s.rawchars(), s.length());
}

bool
TagFromDataitem(const const_IWSubstring& item, IWString& tag) {
  const int open_angle = item.index('<');
  if (open_angle < 0) {
    return false;
  }

  tag.set(item.rawchars(), open_angle + 1);
  return true;
}

int
CountTdts(iwstring_data_source& input, const char* fname) {
  const uint64_t result = input.count_records_starting_with("|");
  if (result > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
    std::cerr << "GFPList::ReadFile:too many TDT records in '" << fname << "'\n";
    return -1;
  }

  return static_cast<int>(result);
}

int
CompareTagsCanonical(const IWString& lhs, const IWString& rhs) {
  const int n1 = lhs.ends_with('<') ? lhs.length() - 1 : lhs.length();
  const int n2 = rhs.ends_with('<') ? rhs.length() - 1 : rhs.length();
  const int n = std::min(n1, n2);

  const int cmp = strncmp(lhs.rawchars(), rhs.rawchars(), n);
  if (cmp != 0) {
    return cmp;
  }
  if (n1 < n2) {
    return -1;
  }
  if (n1 > n2) {
    return 1;
  }
  return 0;
}

bool
NearestNeighbourBetter(const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
  if (lhs.distance != rhs.distance) {
    return lhs.distance < rhs.distance;
  }
  return lhs.index < rhs.index;
}

struct NearestNeighbourMaxHeapCompare {
  bool
  operator()(const NearestNeighbour& lhs, const NearestNeighbour& rhs) const {
    // priority_queue places the highest priority item at the top. Return true
    // when `lhs` is better than `rhs`, so the worst retained neighbour is top.
    return NearestNeighbourBetter(lhs, rhs);
  }
};

using NearestNeighbourMaxHeap =
    std::priority_queue<NearestNeighbour, std::vector<NearestNeighbour>,
                        NearestNeighbourMaxHeapCompare>;

void
MaybeAddNearestNeighbour(NearestNeighbourMaxHeap& heap, int k,
                         NearestNeighbour candidate) {
  if (static_cast<int>(heap.size()) < k) {
    heap.push(candidate);
    return;
  }

  if (NearestNeighbourBetter(candidate, heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

std::vector<NearestNeighbour>
SortedNearestNeighbours(NearestNeighbourMaxHeap& heap) {
  std::vector<NearestNeighbour> result;
  result.reserve(heap.size());
  while (!heap.empty()) {
    result.push_back(heap.top());
    heap.pop();
  }
  std::sort(result.begin(), result.end(), NearestNeighbourBetter);

  return result;
}

float
FixedBinaryTanimoto(const FixedBitVector& lhs, const FixedBitVector& rhs) {
  const int bic = lhs.bits().BitsInCommon(rhs.bits());
  const int denom = lhs.nset() + rhs.nset() - bic;
  if (denom == 0) {
    return 1.0f;
  }

  return static_cast<float>(bic) / static_cast<float>(denom);
}

int
ALogPToPositiveInt(float value) {
  if (value <= -5.0f) {
    return 1;
  }

  if (value >= 10.0f) {
    return 20;
  }

  value += 5.0f;
  return static_cast<int>(value + value + 0.4999f);
}

IWString
ALogPTag(int replicates) {
  IWString result;
  result << "NCALOGP" << replicates << '<';
  return result;
}

int
XLogPToPositiveInt(double value) {
  const int result = static_cast<int>(value + 5.4999);
  if (result <= 0) {
    return 1;
  }

  return result;
}

IWString
XLogPTag(int replicates) {
  IWString result;
  result << "NCXLOGP" << replicates << '<';
  return result;
}

int
TPSAToPositiveInt(double value) {
  int result = static_cast<int>(value / 10.0 + 0.49999);
  if (result <= 0) {
    result = 1;
  }

  return result;
}

IWString
TPSATag(int replicates) {
  IWString result;
  result << "NCTPSA" << replicates << '<';
  return result;
}

IWString
FormulaTag() {
  return IWString("FCFML<");
}

IWString
CATSTag(int max_path_length, bool include_hydrophobic_pairs) {
  IWString result;
  if (include_hydrophobic_pairs) {
    result << "NCCATS";
  } else {
    result << "NCCATSP";
  }
  result << max_path_length << '<';
  return result;
}

int
AppendSanitisedAtomType(const IWString& atom_type, IWString& destination) {
  int appended = 0;
  for (int i = 0; i < atom_type.length(); ++i) {
    const char c = atom_type[i];
    if (c == ':') {
      continue;
    }
    if (!std::isalnum(static_cast<unsigned char>(c))) {
      return 0;
    }
    destination << c;
    ++appended;
  }

  return appended;
}

IWString
AtomPairTag(int min_separation, int max_separation, const IWString& atom_type,
            bool include_out_of_range) {
  IWString result;
  result << "NCAP";
  if (include_out_of_range) {
    result << 'T';
  }
  result << min_separation << 'M' << max_separation;
  if (AppendSanitisedAtomType(atom_type, result) == 0) {
    return IWString();
  }
  result << '<';
  return result;
}

IWString
ECTag(int radius, const IWString& atom_type) {
  IWString atom_type_component;
  if (radius < 0 || !AppendSanitisedAtomType(atom_type, atom_type_component)) {
    return IWString();
  }

  IWString result;
  result << "NCEC" << radius << atom_type_component << '<';
  return result;
}

IWString
SpinachTag(bool label_join_points) {
  return label_join_points ? IWString("FPSPINI<") : IWString("FPSPIN<");
}

IWString
ScaffoldTag(bool label_join_points) {
  return label_join_points ? IWString("FPSCAFI<") : IWString("FPSCAF<");
}

IWString
HexHash(uint64_t hash) {
  static constexpr char kHex[] = "0123456789ABCDEF";
  IWString result;
  for (int shift = 60; shift >= 0; shift -= 4) {
    result << kHex[(hash >> shift) & 0x0f];
  }
  return result;
}

IWString
SubstructureTag(const IWString& smarts, int radius, const IWString& atom_type) {
  IWString atom_type_component;
  if (smarts.empty() || radius < 0 ||
      !AppendSanitisedAtomType(atom_type, atom_type_component)) {
    return IWString();
  }

  uint64_t hash = 14695981039346656037ULL;
  hash = Fnv1a(hash, smarts);
  hash = Fnv1a(hash, atom_type);
  hash = Fnv1a(hash, &radius, sizeof(radius));

  IWString result;
  result << "FPSUB" << radius << atom_type_component << HexHash(hash) << '<';
  return result;
}

}  // namespace

GFPGeneratorSpec::GFPGeneratorSpec(GeneratorKind kind, bool maccs_level2, int replicates)
    : _kind(kind), _maccs_level2(maccs_level2), _replicates(replicates) {
}

GFPGeneratorSpec::GFPGeneratorSpec(GeneratorKind kind, int radius,
                                   const IWString& atom_type)
    : _kind(kind), _radius(radius), _atom_type(atom_type) {
}

GFPGeneratorSpec::GFPGeneratorSpec(GeneratorKind kind, int max_path_length,
                                   bool include_hydrophobic_pairs)
    : _kind(kind),
      _max_path_length(max_path_length),
      _include_hydrophobic_pairs(include_hydrophobic_pairs) {
}

GFPGeneratorSpec::GFPGeneratorSpec(GeneratorKind kind, int min_separation,
                                   int max_separation, const IWString& atom_type,
                                   bool include_out_of_range)
    : _kind(kind),
      _min_separation(min_separation),
      _max_separation(max_separation),
      _include_out_of_range(include_out_of_range),
      _atom_type(atom_type) {
}

GFPGeneratorSpec::GFPGeneratorSpec(GeneratorKind kind, const IWString& smarts,
                                   int radius, const IWString& atom_type,
                                   bool no_match_is_empty)
    : _kind(kind),
      _radius(radius),
      _no_match_is_empty(no_match_is_empty),
      _atom_type(atom_type),
      _smarts(smarts) {
}

GFPGeneratorSpec
GFPGeneratorSpec::MolecularProperties() {
  return GFPGeneratorSpec(GeneratorKind::kMolecularProperties);
}

GFPGeneratorSpec
GFPGeneratorSpec::IWMFingerprint() {
  return GFPGeneratorSpec(GeneratorKind::kIWMFingerprint);
}

GFPGeneratorSpec
GFPGeneratorSpec::MACCSKeys(bool level2) {
  return GFPGeneratorSpec(GeneratorKind::kMACCSKeys, level2);
}

GFPGeneratorSpec
GFPGeneratorSpec::ALogP(int replicates) {
  return GFPGeneratorSpec(GeneratorKind::kALogP, true, replicates);
}

GFPGeneratorSpec
GFPGeneratorSpec::XLogP(int replicates) {
  return GFPGeneratorSpec(GeneratorKind::kXLogP, true, replicates);
}

GFPGeneratorSpec
GFPGeneratorSpec::TPSA(int replicates) {
  return GFPGeneratorSpec(GeneratorKind::kTPSA, true, replicates);
}

GFPGeneratorSpec
GFPGeneratorSpec::FormulaFingerprint() {
  return GFPGeneratorSpec(GeneratorKind::kFormula);
}

GFPGeneratorSpec
GFPGeneratorSpec::CATS(int max_path_length, bool include_hydrophobic_pairs) {
  return GFPGeneratorSpec(GeneratorKind::kCATS, max_path_length,
                          include_hydrophobic_pairs);
}

GFPGeneratorSpec
GFPGeneratorSpec::AtomPair(int min_separation, int max_separation,
                           const IWString& atom_type, bool include_out_of_range) {
  return GFPGeneratorSpec(GeneratorKind::kAtomPair, min_separation, max_separation,
                          atom_type, include_out_of_range);
}

GFPGeneratorSpec
GFPGeneratorSpec::ECFingerprint(int radius, const IWString& atom_type) {
  return GFPGeneratorSpec(GeneratorKind::kECFingerprint, radius, atom_type);
}

GFPGeneratorSpec
GFPGeneratorSpec::RingSubstitution() {
  return GFPGeneratorSpec(GeneratorKind::kRingSubstitution);
}

GFPGeneratorSpec
GFPGeneratorSpec::SpinachFingerprint(bool label_join_points) {
  GFPGeneratorSpec result(GeneratorKind::kSpinachFingerprint);
  result._label_join_points = label_join_points;
  return result;
}

GFPGeneratorSpec
GFPGeneratorSpec::ScaffoldFingerprint(bool label_join_points) {
  GFPGeneratorSpec result(GeneratorKind::kScaffoldFingerprint);
  result._label_join_points = label_join_points;
  return result;
}

GFPGeneratorSpec
GFPGeneratorSpec::SubstructureFingerprint(const IWString& smarts, int radius,
                                          const IWString& atom_type,
                                          bool no_match_is_empty) {
  return GFPGeneratorSpec(GeneratorKind::kSubstructureFingerprint, smarts, radius,
                          atom_type, no_match_is_empty);
}

std::vector<Component>
GFPGeneratorSpec::Components() const {
  switch (_kind) {
    case GeneratorKind::kMolecularProperties:
      return {Component{ComponentKind::kMolecularProperties, IWString("MPR<"), 0, 1.0f}};
    case GeneratorKind::kIWMFingerprint:
      return {Component{ComponentKind::kFixedBinary, IWString("FPIW<"), 0, 1.0f}};
    case GeneratorKind::kMACCSKeys: {
      std::vector<Component> result;
      result.push_back(
          Component{ComponentKind::kFixedBinary, IWString("FPMK<"), 0, 1.0f});
      if (_maccs_level2) {
        result.push_back(
            Component{ComponentKind::kFixedBinary, IWString("FPMK2<"), 0, 1.0f});
      }
      return result;
    }
    case GeneratorKind::kALogP:
      return {Component{ComponentKind::kSparse, ALogPTag(_replicates), 0, 1.0f}};
    case GeneratorKind::kXLogP:
      return {Component{ComponentKind::kSparse, XLogPTag(_replicates), 0, 1.0f}};
    case GeneratorKind::kTPSA:
      return {Component{ComponentKind::kSparse, TPSATag(_replicates), 0, 1.0f}};
    case GeneratorKind::kFormula:
      return {Component{ComponentKind::kFixedCounted, FormulaTag(), 0, 1.0f}};
    case GeneratorKind::kCATS:
      return {Component{ComponentKind::kSparse,
                        CATSTag(_max_path_length, _include_hydrophobic_pairs), 0, 1.0f}};
    case GeneratorKind::kAtomPair:
      return {Component{ComponentKind::kSparse,
                        AtomPairTag(_min_separation, _max_separation, _atom_type,
                                    _include_out_of_range),
                        0, 1.0f}};
    case GeneratorKind::kECFingerprint:
      return {Component{ComponentKind::kSparse, ECTag(_radius, _atom_type), 0, 1.0f}};
    case GeneratorKind::kRingSubstitution:
      return {Component{ComponentKind::kSparse, IWString("NCRS<"), 0, 1.0f}};
    case GeneratorKind::kSpinachFingerprint:
      return {Component{ComponentKind::kFixedBinary, SpinachTag(_label_join_points), 0,
                        1.0f}};
    case GeneratorKind::kScaffoldFingerprint:
      return {Component{ComponentKind::kFixedBinary, ScaffoldTag(_label_join_points), 0,
                        1.0f}};
    case GeneratorKind::kSubstructureFingerprint:
      return {Component{ComponentKind::kFixedBinary,
                        SubstructureTag(_smarts, _radius, _atom_type), 0, 1.0f}};
  }

  return {};
}

std::string
GFPGeneratorSpec::Repr() const {
  switch (_kind) {
    case GeneratorKind::kMolecularProperties:
      return "GFP.mpr()";
    case GeneratorKind::kIWMFingerprint:
      return "GFP.iw()";
    case GeneratorKind::kMACCSKeys:
      return _maccs_level2 ? "GFP.maccs(level2=True)" : "GFP.maccs(level2=False)";
    case GeneratorKind::kALogP:
      return "GFP.alogp(replicates=" + std::to_string(_replicates) + ")";
    case GeneratorKind::kXLogP:
      return "GFP.xlogp(replicates=" + std::to_string(_replicates) + ")";
    case GeneratorKind::kTPSA:
      return "GFP.tpsa(replicates=" + std::to_string(_replicates) + ")";
    case GeneratorKind::kFormula:
      return "GFP.formula()";
    case GeneratorKind::kCATS:
      return "GFP.cats(max_path_length=" + std::to_string(_max_path_length) +
             ", include_hydrophobic_pairs=" +
             (_include_hydrophobic_pairs ? "True" : "False") + ")";
    case GeneratorKind::kAtomPair:
      return "GFP.atom_pair(min_separation=" + std::to_string(_min_separation) +
             ", max_separation=" + std::to_string(_max_separation) + ", atom_type='" +
             _atom_type.AsString() +
             "', include_out_of_range=" + (_include_out_of_range ? "True" : "False") +
             ")";
    case GeneratorKind::kECFingerprint:
      return "GFP.ec(radius=" + std::to_string(_radius) + ", atom_type='" +
             _atom_type.AsString() + "')";
    case GeneratorKind::kRingSubstitution:
      return "GFP.ring_substitution()";
    case GeneratorKind::kSpinachFingerprint:
      return std::string("GFP.spinach(label_join_points=") +
             (_label_join_points ? "True" : "False") + ")";
    case GeneratorKind::kScaffoldFingerprint:
      return std::string("GFP.scaffold(label_join_points=") +
             (_label_join_points ? "True" : "False") + ")";
    case GeneratorKind::kSubstructureFingerprint:
      return std::string("GFP.substructure(smarts='") + _smarts.AsString() +
             "', radius=" + std::to_string(_radius) + ", atom_type='" +
             _atom_type.AsString() + "', no_match='" +
             (_no_match_is_empty ? "empty" : "error") + "')";
  }

  return "GFP.unknown()";
}

struct GeneratorComponentAssignment {
  ComponentKind kind;
  int index = 0;
};

class FingerprintGeneratorImplementation {
 public:
  virtual ~FingerprintGeneratorImplementation() = default;

  virtual std::vector<Component> Components() const = 0;
  virtual int Generate(Molecule& m,
                       const std::vector<GeneratorComponentAssignment>& slots,
                       GFPFingerprint& destination) = 0;
};

class MolecularPropertiesFingerprintGenerator
    : public FingerprintGeneratorImplementation {
 private:
  Molecular_Properties_Generator _mpr;

 public:
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
MolecularPropertiesFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kMolecularProperties, IWString("MPR<"), 0, 1.0f}};
}

int
MolecularPropertiesFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kMolecularProperties) {
    std::cerr << "MolecularPropertiesFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  int properties[8];
  _mpr(m, properties);
  return destination.mutable_molecular_properties().BuildFromArray(properties, 8);
}

class IWMFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  IWMFingerprintOptions _iw_options;

 public:
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
IWMFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kFixedBinary, IWString("FPIW<"), 0, 1.0f}};
}

int
IWMFingerprintGenerator::Generate(Molecule& m,
                                  const std::vector<GeneratorComponentAssignment>& slots,
                                  GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kFixedBinary) {
    std::cerr << "IWMFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  IWMFingerprint iwfp(_iw_options);
  iwfp.construct_fingerprint(m);
  return destination.mutable_fixed_binary(slots[0].index)
      .BuildFromArray(iwfp.vector(), _iw_options.bits_per_fingerprint);
}

class MACCSFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  MACCSKeys _mk;
  bool _level2 = true;
  int _tmp[kMaccskeysNbits];

 public:
  explicit MACCSFingerprintGenerator(bool level2) : _level2(level2) {
  }

  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
MACCSFingerprintGenerator::Components() const {
  std::vector<Component> result;
  result.push_back(Component{ComponentKind::kFixedBinary, IWString("FPMK<"), 0, 1.0f});
  if (_level2) {
    result.push_back(Component{ComponentKind::kFixedBinary, IWString("FPMK2<"), 0, 1.0f});
  }
  return result;
}

int
MACCSFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.empty() || slots.size() > 2) {
    std::cerr << "MACCSFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }
  for (const GeneratorComponentAssignment& slot : slots) {
    if (slot.kind != ComponentKind::kFixedBinary) {
      std::cerr << "MACCSFingerprintGenerator::Generate:invalid slot kind\n";
      return 0;
    }
  }

  std::fill_n(_tmp, kMaccskeysNbits, 0);
  _mk(m, _tmp);
  assert(_mk.nbits() == kMaccskeysNbits);
  if (!destination.mutable_fixed_binary(slots[0].index)
           .BuildFromArray(_tmp, _mk.nbits())) {
    return 0;
  }

  if (slots.size() == 1) {
    return 1;
  }

  _mk.set_level_2_fingerprint(_tmp);
  return destination.mutable_fixed_binary(slots[1].index)
      .BuildFromArray(_tmp, _mk.nbits());
}

class ReplicatedCountFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  IWString _tag;
  int _replicates = kDefaultSparseReplicates;

 protected:
  virtual std::optional<int> Count(Molecule& m) = 0;
  virtual const char* Name() const = 0;

 public:
  ReplicatedCountFingerprintGenerator(const IWString& tag, int replicates)
      : _tag(tag), _replicates(replicates) {
  }

  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
ReplicatedCountFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kSparse, _tag, 0, 1.0f}};
}

int
ReplicatedCountFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kSparse) {
    std::cerr << Name() << "::Generate:invalid slots\n";
    return 0;
  }
  if (_replicates <= 0) {
    std::cerr << Name() << "::Generate:invalid replicates " << _replicates << '\n';
    return 0;
  }

  std::optional<int> count = Count(m);
  if (!count) {
    return 0;
  }

  return destination.mutable_sparse(slots[0].index)
      .build_from_replicates(_replicates, *count);
}

class ALogPFingerprintGenerator : public ReplicatedCountFingerprintGenerator {
 private:
  alogp::ALogP _alogp;

 protected:
  std::optional<int> Count(Molecule& m) override;

  const char*
  Name() const override {
    return "ALogPFingerprintGenerator";
  }

 public:
  explicit ALogPFingerprintGenerator(int replicates)
      : ReplicatedCountFingerprintGenerator(ALogPTag(replicates), replicates) {
    _alogp.set_display_error_messages(0);
  }
};

std::optional<int>
ALogPFingerprintGenerator::Count(Molecule& m) {
  std::optional<float> logp = _alogp.LogP(m);
  if (!logp) {
    return std::nullopt;
  }

  return ALogPToPositiveInt(*logp);
}

class XLogPFingerprintGenerator : public ReplicatedCountFingerprintGenerator {
 private:
  xlogp::XLogPCalc _xlogp;

 protected:
  std::optional<int> Count(Molecule& m) override;

  const char*
  Name() const override {
    return "XLogPFingerprintGenerator";
  }

 public:
  explicit XLogPFingerprintGenerator(int replicates)
      : ReplicatedCountFingerprintGenerator(XLogPTag(replicates), replicates) {
    _xlogp.SetIssueUnclassifiedAtomMessages(false);
  }
};

std::optional<int>
XLogPFingerprintGenerator::Count(Molecule& m) {
  std::optional<double> logp = _xlogp.LogP(m);
  if (!logp) {
    return std::nullopt;
  }

  return XLogPToPositiveInt(*logp);
}

class TPSAFingerprintGenerator : public ReplicatedCountFingerprintGenerator {
 private:
  nvrtspsa::NovartisPolarSurfaceArea _tpsa;

 protected:
  std::optional<int> Count(Molecule& m) override;

  const char*
  Name() const override {
    return "TPSAFingerprintGenerator";
  }

 public:
  explicit TPSAFingerprintGenerator(int replicates)
      : ReplicatedCountFingerprintGenerator(TPSATag(replicates), replicates) {
    _tpsa.set_display_psa_unclassified_atom_messages(0);
  }
};

std::optional<int>
TPSAFingerprintGenerator::Count(Molecule& m) {
  std::optional<double> tpsa = _tpsa.PolarSurfaceArea(m);
  if (!tpsa) {
    return std::nullopt;
  }

  return TPSAToPositiveInt(*tpsa);
}

class FormulaFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  int _count[kFormulaFingerprintNbits];

 public:
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
FormulaFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kFixedCounted, FormulaTag(), 0, 1.0f}};
}

int
FormulaFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kFixedCounted) {
    std::cerr << "FormulaFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  mformula::MFormula formula;
  formula.set_log_scaling_factors(kDefaultFormulaLogS, kDefaultFormulaLogD);
  if (!formula.Build(m)) {
    return 0;
  }
  if (!formula.ToFixedCountedFingerprint(_count, kFormulaFingerprintNbits)) {
    return 0;
  }

  return destination.mutable_fixed_counted(slots[0].index)
      .construct_from_array_of_ints(_count, kFormulaFingerprintNbits);
}

class CATSFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  int _max_path_length = 10;
  bool _include_hydrophobic_pairs = true;
  jwcats::JWCats _jwcats;

 public:
  CATSFingerprintGenerator(int max_path_length, bool include_hydrophobic_pairs)
      : _max_path_length(max_path_length),
        _include_hydrophobic_pairs(include_hydrophobic_pairs) {
  }

  int Initialise();
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

int
CATSFingerprintGenerator::Initialise() {
  if (_max_path_length < 1) {
    std::cerr << "CATSFingerprintGenerator::Initialise:invalid max path length "
              << _max_path_length << '\n';
    return 0;
  }

  _jwcats.SetMaxBondSeparation(_max_path_length);
  _jwcats.SetIncludeHydrophobicPairs(_include_hydrophobic_pairs);
  _jwcats.SetScalingType(0);

  if (!_jwcats.charge_assigner().BuildFromDefaultEnvs()) {
    std::cerr << "CATSFingerprintGenerator::Initialise:cannot initialise charge "
                 "assigner; ensure LILLYMOL_HOME is defined\n";
    return 0;
  }

  constexpr int kQuiet = 0;
  if (!_jwcats.donor_acceptor_assigner().BuildFromDefaultEnv(kQuiet)) {
    std::cerr << "CATSFingerprintGenerator::Initialise:cannot initialise donor/acceptor "
                 "assigner; ensure LILLYMOL_HOME is defined\n";
    return 0;
  }

  return _jwcats.Initialise();
}

std::vector<Component>
CATSFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kSparse,
                    CATSTag(_max_path_length, _include_hydrophobic_pairs), 0, 1.0f}};
}

int
CATSFingerprintGenerator::Generate(Molecule& m,
                                   const std::vector<GeneratorComponentAssignment>& slots,
                                   GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kSparse) {
    std::cerr << "CATSFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  jwcats::Result result;
  const jwcats::ComputeStatus status = _jwcats.Compute(m, result);
  if (status != jwcats::ComputeStatus::kOk) {
    std::cerr << "CATSFingerprintGenerator::Generate:cannot compute CATS fingerprint\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;
  const std::vector<int>& write_array_value = _jwcats.write_array_value();
  const int array_size = _jwcats.array_size();
  for (int i = 0; i < array_size; ++i) {
    if (!write_array_value[i]) {
      continue;
    }
    if (result.scaled_counts[i] == 0.0) {
      continue;
    }
    const int count = static_cast<int>(result.scaled_counts[i] + 0.01);
    if (count > 0) {
      sfc.hit_bit(i, count);
    }
  }

  return destination.mutable_sparse(slots[0].index)
      .build_from_sparse_fingerprint_creator(sfc);
}

class AtomPairFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  int _min_separation = 1;
  int _max_separation = 10;
  IWString _atom_type_string;
  bool _include_out_of_range = false;
  Atom_Typing_Specification _atom_typing;
  atom_pair_fingerprint::AtomPairFingerprint _atom_pair_fingerprint;

 public:
  AtomPairFingerprintGenerator(int min_separation, int max_separation,
                               const IWString& atom_type, bool include_out_of_range)
      : _min_separation(min_separation),
        _max_separation(max_separation),
        _atom_type_string(atom_type),
        _include_out_of_range(include_out_of_range) {
  }

  int Initialise();
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

int
AtomPairFingerprintGenerator::Initialise() {
  if (_min_separation < 0 || _max_separation < _min_separation) {
    return 0;
  }
  if (_atom_type_string.empty() || AtomPairTag(_min_separation, _max_separation,
                                               _atom_type_string, _include_out_of_range)
                                       .empty()) {
    return 0;
  }

  const_IWSubstring tmp(_atom_type_string);
  if (!_atom_typing.build(tmp)) {
    return 0;
  }

  _atom_pair_fingerprint.set_min_separation(_min_separation);
  _atom_pair_fingerprint.set_max_separation(_max_separation);
  _atom_pair_fingerprint.set_include_out_of_range_separations(_include_out_of_range);

  return 1;
}

std::vector<Component>
AtomPairFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kSparse,
                    AtomPairTag(_min_separation, _max_separation, _atom_type_string,
                                _include_out_of_range),
                    0, 1.0f}};
}

int
AtomPairFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kSparse) {
    std::cerr << "AtomPairFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  std::vector<atom_pair_fingerprint::atom_type_t> atom_type(m.natoms());
  if (!_atom_typing.assign_atom_types(m, atom_type.data())) {
    std::cerr << "AtomPairFingerprintGenerator::Generate:cannot assign atom types\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;
  if (!_atom_pair_fingerprint.Fingerprint(m, nullptr, atom_type.data(), sfc)) {
    return 0;
  }

  return destination.mutable_sparse(slots[0].index)
      .build_from_sparse_fingerprint_creator(sfc);
}

class ECFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  int _radius = 3;
  IWString _atom_type_string;
  Atom_Typing_Specification _atom_typing;
  iwecfp::Iwecfp _iwecfp;

 public:
  ECFingerprintGenerator(int radius, const IWString& atom_type)
      : _radius(radius), _atom_type_string(atom_type) {
    _iwecfp.set_max_radius(radius);
  }

  int Initialise();
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

int
ECFingerprintGenerator::Initialise() {
  if (_radius < 0 || _atom_type_string.empty()) {
    return 0;
  }

  const_IWSubstring tmp(_atom_type_string);
  return _atom_typing.build(tmp);
}

std::vector<Component>
ECFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kSparse, ECTag(_radius, _atom_type_string), 0, 1.0f}};
}

int
ECFingerprintGenerator::Generate(Molecule& m,
                                 const std::vector<GeneratorComponentAssignment>& slots,
                                 GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kSparse) {
    std::cerr << "ECFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  const int matoms = m.natoms();
  if (matoms == 0) {
    return 0;
  }

  std::unique_ptr<iwecfp::atype_t[]> atom_type =
      std::make_unique<iwecfp::atype_t[]>(matoms);
  if (!_atom_typing.assign_atom_types(m, atom_type.get())) {
    std::cerr << "ECFingerprintGenerator::Generate:cannot assign atom types\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;
  const iwecfp::FingerprintResult result = _iwecfp.Fingerprint(m, atom_type.get(), &sfc);
  if (result != iwecfp::FingerprintResult::kOk) {
    std::cerr << "ECFingerprintGenerator::Generate:fingerprint generation failed\n";
    return 0;
  }

  return destination.mutable_sparse(slots[0].index)
      .build_from_sparse_fingerprint_creator(sfc);
}

class RingSubstitutionFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  ring_substitution::RingSubstitutionGenerator _generator;

 public:
  RingSubstitutionFingerprintGenerator() {
    _generator.set_positional_information_only(false);
    _generator.set_full_atom_types(true);
    _generator.set_place_single_feature_bits(true);
  }

  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

std::vector<Component>
RingSubstitutionFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kSparse, IWString("NCRS<"), 0, 1.0f}};
}

int
RingSubstitutionFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kSparse) {
    std::cerr << "RingSubstitutionFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;
  if (m.nrings() > 0 && !_generator.Generate(m, sfc)) {
    std::cerr << "RingSubstitutionFingerprintGenerator::Generate:fingerprint generation "
                 "failed\n";
    return 0;
  }

  return destination.mutable_sparse(slots[0].index)
      .build_from_sparse_fingerprint_creator(sfc);
}

enum class FingerprintSubset : uint8_t {
  kSpinach,
  kScaffold,
};

class ScaffoldSpinachFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  FingerprintSubset _subset = FingerprintSubset::kSpinach;
  bool _label_join_points = false;
  IWString _tag;
  IWString _atom_type_string;
  Atom_Typing_Specification _atom_typing;
  IWMFingerprintOptions _iw_options;

 public:
  ScaffoldSpinachFingerprintGenerator(FingerprintSubset subset, bool label_join_points)
      : _subset(subset), _label_join_points(label_join_points) {
    _tag = (_subset == FingerprintSubset::kSpinach)
               ? SpinachTag(_label_join_points)
               : ScaffoldTag(_label_join_points);
    _atom_type_string = _label_join_points ? IWString("UST:AIRZ")
                                           : IWString("UST:ARZ");
  }

  int Initialise();
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

int
ScaffoldSpinachFingerprintGenerator::Initialise() {
  const_IWSubstring tmp(_atom_type_string);
  return _atom_typing.build(tmp);
}

std::vector<Component>
ScaffoldSpinachFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kFixedBinary, _tag, 0, 1.0f}};
}

int
ScaffoldSpinachFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kFixedBinary) {
    std::cerr << "ScaffoldSpinachFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  Molecule copy(m);
  const int matoms = copy.natoms();
  std::vector<int> include_atom(matoms, 0);
  if (!copy.IdentifySpinachLabel(include_atom.data(), _label_join_points ? 1 : 0)) {
    std::cerr << "ScaffoldSpinachFingerprintGenerator::Generate:cannot identify spinach\n";
    return 0;
  }

  if (_subset == FingerprintSubset::kScaffold) {
    for (int& include : include_atom) {
      include = include ? 0 : 1;
    }
  }

  std::vector<uint32_t> atom_type(matoms);
  if (!_atom_typing.assign_atom_types(copy, atom_type.data())) {
    std::cerr << "ScaffoldSpinachFingerprintGenerator::Generate:cannot assign atom "
                 "types\n";
    return 0;
  }

  IWMFingerprint fp(_iw_options);
  if (!fp.construct_fingerprint(copy, atom_type.data(), include_atom.data())) {
    std::cerr << "ScaffoldSpinachFingerprintGenerator::Generate:fingerprint generation "
                 "failed\n";
    return 0;
  }

  return destination.mutable_fixed_binary(slots[0].index)
      .BuildFromArray(fp.vector(), _iw_options.bits_per_fingerprint);
}

class SubstructureFingerprintGenerator : public FingerprintGeneratorImplementation {
 private:
  IWString _smarts;
  int _radius = 0;
  IWString _atom_type_string;
  bool _no_match_is_empty = true;
  IWString _tag;
  Substructure_Hit_Statistics _query;
  Atom_Typing_Specification _atom_typing;
  IWMFingerprintOptions _iw_options;

 public:
  SubstructureFingerprintGenerator(const IWString& smarts, int radius,
                                   const IWString& atom_type,
                                   bool no_match_is_empty)
      : _smarts(smarts),
        _radius(radius),
        _atom_type_string(atom_type),
        _no_match_is_empty(no_match_is_empty),
        _tag(SubstructureTag(smarts, radius, atom_type)) {
  }

  int Initialise();
  std::vector<Component> Components() const override;
  int Generate(Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
               GFPFingerprint& destination) override;
};

int
SubstructureFingerprintGenerator::Initialise() {
  if (_smarts.empty() || _radius < 0 || _atom_type_string.empty() || _tag.empty()) {
    return 0;
  }

  if (!_query.create_from_smarts(_smarts)) {
    std::cerr << "SubstructureFingerprintGenerator::Initialise:invalid smarts '"
              << _smarts << "'\n";
    return 0;
  }
  _query.set_find_unique_embeddings_only(1);

  const_IWSubstring atom_type(_atom_type_string);
  if (!_atom_typing.build(atom_type)) {
    std::cerr << "SubstructureFingerprintGenerator::Initialise:invalid atom type '"
              << _atom_type_string << "'\n";
    return 0;
  }

  return 1;
}

std::vector<Component>
SubstructureFingerprintGenerator::Components() const {
  return {Component{ComponentKind::kFixedBinary, _tag, 0, 1.0f}};
}

int
SubstructureFingerprintGenerator::Generate(
    Molecule& m, const std::vector<GeneratorComponentAssignment>& slots,
    GFPFingerprint& destination) {
  if (slots.size() != 1 || slots[0].kind != ComponentKind::kFixedBinary) {
    std::cerr << "SubstructureFingerprintGenerator::Generate:invalid slots\n";
    return 0;
  }

  const int matoms = m.natoms();
  std::vector<int> include_atom(matoms, 0);

  Molecule_to_Match target(&m);
  Substructure_Results results;
  const int nhits = _query.substructure_search(target, results);
  if (nhits == 0) {
    if (!_no_match_is_empty) {
      std::cerr << "SubstructureFingerprintGenerator::Generate:no query match for "
                << m.name() << "\n";
      return 0;
    }
  } else {
    for (int i = 0; i < nhits; ++i) {
      const Set_of_Atoms* embedding = results.embedding(i);
      embedding->set_vector(include_atom.data(), 1);
    }

    if (_radius > 0) {
      m.recompute_distance_matrix();
      std::vector<int> matched(include_atom);
      for (int i = 0; i < matoms; ++i) {
        if (!matched[i]) {
          continue;
        }
        for (int j = 0; j < matoms; ++j) {
          if (include_atom[j]) {
            continue;
          }
          if (m.bonds_between(i, j) <= _radius) {
            include_atom[j] = 1;
          }
        }
      }
    }
  }

  std::vector<uint32_t> atom_type(matoms);
  if (!_atom_typing.assign_atom_types(m, atom_type.data())) {
    std::cerr << "SubstructureFingerprintGenerator::Generate:cannot assign atom types\n";
    return 0;
  }

  IWMFingerprint fp(_iw_options);
  if (!fp.construct_fingerprint(m, atom_type.data(), include_atom.data())) {
    std::cerr << "SubstructureFingerprintGenerator::Generate:fingerprint generation "
                 "failed\n";
    return 0;
  }

  return destination.mutable_fixed_binary(slots[0].index)
      .BuildFromArray(fp.vector(), _iw_options.bits_per_fingerprint);
}

class StandardFingerprintGenerator {
 private:
  Chemical_Standardisation _chemical_standardisation;
  bool _preprocess = true;
  std::vector<std::unique_ptr<FingerprintGeneratorImplementation>> _generators;
  std::vector<std::vector<GeneratorComponentAssignment>> _assignments;

  void Preprocess(Molecule& m);

 public:
  explicit StandardFingerprintGenerator(bool preprocess);

  int AddGenerator(std::unique_ptr<FingerprintGeneratorImplementation> generator,
                   const std::vector<GeneratorComponentAssignment>& assignments);
  int Generate(Molecule& m, GFPFingerprint& destination);
};

StandardFingerprintGenerator::StandardFingerprintGenerator(bool preprocess)
    : _preprocess(preprocess) {
  _chemical_standardisation.activate_all();
}

void
StandardFingerprintGenerator::Preprocess(Molecule& m) {
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  _chemical_standardisation.process(m);
}

int
StandardFingerprintGenerator::AddGenerator(
    std::unique_ptr<FingerprintGeneratorImplementation> generator,
    const std::vector<GeneratorComponentAssignment>& assignments) {
  if (generator == nullptr || assignments.empty()) {
    return 0;
  }

  _generators.push_back(std::move(generator));
  _assignments.push_back(assignments);
  return 1;
}

int
StandardFingerprintGenerator::Generate(Molecule& input, GFPFingerprint& destination) {
  std::unique_ptr<Molecule> copy;
  Molecule* m = &input;
  if (_preprocess) {
    copy = std::make_unique<Molecule>(input);
    Preprocess(*copy);
    m = copy.get();
  }

  for (int i = 0; i < static_cast<int>(_generators.size()); ++i) {
    if (!_generators[i]->Generate(*m, _assignments[i], destination)) {
      return 0;
    }
  }

  return 1;
}

MolecularProperties::~MolecularProperties() {
  delete[] _property;
}

MolecularProperties::MolecularProperties(const MolecularProperties& rhs) {
  *this = rhs;
}

MolecularProperties&
MolecularProperties::operator=(const MolecularProperties& rhs) {
  if (this == &rhs) {
    return *this;
  }

  delete[] _property;
  _nproperties = rhs._nproperties;
  if (_nproperties == 0) {
    _property = nullptr;
    return *this;
  }

  _property = new uint8_t[_nproperties];
  std::copy(rhs._property, rhs._property + _nproperties, _property);

  return *this;
}

MolecularProperties::MolecularProperties(MolecularProperties&& rhs) noexcept {
  *this = std::move(rhs);
}

MolecularProperties&
MolecularProperties::operator=(MolecularProperties&& rhs) noexcept {
  if (this == &rhs) {
    return *this;
  }

  delete[] _property;

  _nproperties = rhs._nproperties;
  _property = rhs._property;

  rhs._nproperties = 0;
  rhs._property = nullptr;

  return *this;
}

int
MolecularProperties::Build(const const_IWSubstring& buffer) {
  IWDYFP fp;
  if (!fp.construct_from_tdt_record(buffer)) {
    std::cerr << "MolecularProperties::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  const int nproperties = fp.nbits() / 8;
  if (nproperties <= 0) {
    std::cerr << "MolecularProperties::Build:no properties in '" << buffer << "'\n";
    return 0;
  }

  if (_nproperties != nproperties) {
    delete[] _property;
    _nproperties = nproperties;
    _property = new uint8_t[_nproperties];
  }

  const unsigned char* bits = fp.bits();
  for (int i = 0; i < _nproperties; ++i) {
    _property[i] = bits[i];
  }

  return 1;
}

int
MolecularProperties::BuildFromArray(const int* properties, int nproperties) {
  if (nproperties <= 0) {
    return 0;
  }

  if (_nproperties != nproperties) {
    delete[] _property;
    _nproperties = nproperties;
    _property = new uint8_t[_nproperties];
  }

  for (int i = 0; i < _nproperties; ++i) {
    int value = properties[i];
    if (value < 0) {
      cerr << "MolecularProperties::BuildFromArray:negative value " << value <<
              " set to zero\n";
      value = 0;
    } else if (value > 255) {
      // No warning for saturation - large molecules..
      value = 255;
    }
    _property[i] = static_cast<uint8_t>(value);
  }

  return 1;
}

float
MolecularProperties::Similarity(const MolecularProperties& rhs) const {
  if (_nproperties != rhs._nproperties || _nproperties == 0) {
    return 0.0f;
  }

  return gfp_internal::molecular_property_ratio_table().similarity(_property, rhs._property,
                                                          _nproperties);
}

int
FixedBitVector::Build(const const_IWSubstring& buffer) {
  if (!_bits.ConstructFromTdtRecordNset(buffer, _nset)) {
    std::cerr << "FixedBitVector::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
FixedBitVector::BuildFromArray(const int* bits, int nbits) {
  if (nbits <= 0) {
    return 0;
  }

  if (!_bits.ConstructFromArrayIWBitsOrder(bits, nbits)) {
    return 0;
  }
  _nset = _bits.nset();

  return 1;
}

float
FixedBitVector::Tanimoto(const FixedBitVector& rhs) const {
  return FixedBinaryTanimoto(*this, rhs);
}

void
GFPFingerprint::FreeArrays() {
  delete[] _fixed_binary;
  delete[] _sparse;
  delete[] _fixed_counted;

  _fixed_binary = nullptr;
  _sparse = nullptr;
  _fixed_counted = nullptr;
}

GFPFingerprint::~GFPFingerprint() {
  FreeArrays();
}

GFPFingerprint::GFPFingerprint(GFPFingerprint&& rhs) noexcept {
  *this = std::move(rhs);
}

GFPFingerprint&
GFPFingerprint::operator=(GFPFingerprint&& rhs) noexcept {
  if (this == &rhs) {
    return *this;
  }

  FreeArrays();

  _context_hash = rhs._context_hash;
  _molecular_properties = std::move(rhs._molecular_properties);
  _fixed_binary = rhs._fixed_binary;
  _sparse = rhs._sparse;
  _fixed_counted = rhs._fixed_counted;

  rhs._context_hash = 0;
  rhs._fixed_binary = nullptr;
  rhs._sparse = nullptr;
  rhs._fixed_counted = nullptr;

  return *this;
}

int
GFPFingerprint::Allocate(const GFPContext& context) {
  FreeArrays();

  if (context.nfixed_binary()) {
    _fixed_binary = new FixedBitVector[context.nfixed_binary()];
  }
  if (context.nsparse()) {
    _sparse = new Sparse_Fingerprint[context.nsparse()];
  }
  if (context.nfixed_counted()) {
    _fixed_counted = new Fixed_Size_Counted_Fingerprint_uchar[context.nfixed_counted()];
  }

  _context_hash = context.context_hash();
  return 1;
}

int
GFPFingerprint::Build(const IW_TDT& tdt, const GFPContext& context) {
  if (_context_hash != context.context_hash()) {
    if (!Allocate(context)) {
      return 0;
    }
  }

  for (const Component& component : context.components()) {
    const_IWSubstring dataitem;
    if (!tdt.dataitem(component.tag, dataitem)) {
      std::cerr << "GFPFingerprint::Build:missing '" << component.tag << "'\n";
      return 0;
    }
    if (dataitem.ends_with('\n')) {
      dataitem.chop();
    }

    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        if (!_molecular_properties.Build(dataitem)) {
          return 0;
        }
        break;
      case ComponentKind::kFixedBinary:
        if (!_fixed_binary[component.index].Build(dataitem)) {
          return 0;
        }
        break;
      case ComponentKind::kSparse:
        if (!_sparse[component.index].construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPFingerprint::Build:cannot parse sparse '" << component.tag
                    << "'\n";
          return 0;
        }
        break;
      case ComponentKind::kFixedCounted:
        if (!_fixed_counted[component.index].construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPFingerprint::Build:cannot parse fixed counted '"
                    << component.tag << "'\n";
          return 0;
        }
        break;
    }
  }

  return 1;
}

void
GFPContext::ComputeHash() {
  uint64_t hash = 14695981039346656037ULL;
  for (const Component& component : _components) {
    const uint8_t kind = static_cast<uint8_t>(component.kind);
    hash = Fnv1a(hash, &kind, sizeof(kind));
    hash = Fnv1a(hash, &component.index, sizeof(component.index));
    hash = Fnv1a(hash, component.tag);
  }

  _context_hash = hash;
}

void
GFPContext::CanonicalizeComponents() {
  std::sort(_components.begin(), _components.end(),
            [](const Component& lhs, const Component& rhs) {
              return CompareTagsCanonical(lhs.tag, rhs.tag) < 0;
            });

  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;

  for (Component& component : _components) {
    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        component.index = 0;
        _has_molecular_properties = true;
        break;
      case ComponentKind::kFixedBinary:
        component.index = _nfixed_binary;
        ++_nfixed_binary;
        break;
      case ComponentKind::kSparse:
        component.index = _nsparse;
        ++_nsparse;
        break;
      case ComponentKind::kFixedCounted:
        component.index = _nfixed_counted;
        ++_nfixed_counted;
        break;
    }
  }
}

void
GFPContext::BuildDefaultActiveList() {
  _active.clear();
  if (_components.empty()) {
    return;
  }

  const float weight = 1.0f / static_cast<float>(_components.size());
  for (Component& component : _components) {
    component.weight = weight;
    _active.push_back(ActiveComponent{component.kind, component.index, weight});
  }
}

GFPContext::GFPContext() = default;

GFPContext::~GFPContext() = default;

namespace {

std::unique_ptr<FingerprintGeneratorImplementation>
MakeGenerator(const GFPGeneratorSpec& spec) {
  switch (spec.kind()) {
    case GeneratorKind::kMolecularProperties:
      return std::make_unique<MolecularPropertiesFingerprintGenerator>();
    case GeneratorKind::kIWMFingerprint:
      return std::make_unique<IWMFingerprintGenerator>();
    case GeneratorKind::kMACCSKeys:
      return std::make_unique<MACCSFingerprintGenerator>(spec.maccs_level2());
    case GeneratorKind::kALogP:
      return std::make_unique<ALogPFingerprintGenerator>(spec.replicates());
    case GeneratorKind::kXLogP:
      return std::make_unique<XLogPFingerprintGenerator>(spec.replicates());
    case GeneratorKind::kTPSA:
      return std::make_unique<TPSAFingerprintGenerator>(spec.replicates());
    case GeneratorKind::kFormula:
      return std::make_unique<FormulaFingerprintGenerator>();
    case GeneratorKind::kCATS:
      return std::make_unique<CATSFingerprintGenerator>(spec.max_path_length(),
                                                        spec.include_hydrophobic_pairs());
    case GeneratorKind::kAtomPair:
      return std::make_unique<AtomPairFingerprintGenerator>(
          spec.min_separation(), spec.max_separation(), spec.atom_type(),
          spec.include_out_of_range());
    case GeneratorKind::kECFingerprint:
      return std::make_unique<ECFingerprintGenerator>(spec.radius(), spec.atom_type());
    case GeneratorKind::kRingSubstitution:
      return std::make_unique<RingSubstitutionFingerprintGenerator>();
    case GeneratorKind::kSpinachFingerprint:
      return std::make_unique<ScaffoldSpinachFingerprintGenerator>(
          FingerprintSubset::kSpinach, spec.label_join_points());
    case GeneratorKind::kScaffoldFingerprint:
      return std::make_unique<ScaffoldSpinachFingerprintGenerator>(
          FingerprintSubset::kScaffold, spec.label_join_points());
    case GeneratorKind::kSubstructureFingerprint:
      return std::make_unique<SubstructureFingerprintGenerator>(
          spec.smarts(), spec.radius(), spec.atom_type(), spec.no_match_is_empty());
  }

  return nullptr;
}

bool
SameTag(const IWString& lhs, const IWString& rhs) {
  return lhs == rhs;
}

int
FindAssignedComponent(const std::vector<Component>& components, const Component& desired,
                      GeneratorComponentAssignment& assignment) {
  for (const Component& component : components) {
    if (component.kind == desired.kind && SameTag(component.tag, desired.tag)) {
      assignment.kind = component.kind;
      assignment.index = component.index;
      return 1;
    }
  }

  return 0;
}

}  // namespace

int
GFPContext::BuildFromTdt(const IW_TDT& tdt) {
  _components.clear();
  _active.clear();
  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;
  _context_hash = 0;
  _standard_generator.reset();

  int i = 0;
  const_IWSubstring item;
  while (tdt.next_dataitem(item, i)) {
    IWString tag;
    if (!TagFromDataitem(item, tag)) {
      continue;
    }

    if (tag == kMolecularPropertiesTag) {
      _has_molecular_properties = true;
      _components.push_back(Component{ComponentKind::kMolecularProperties, tag, 0, 1.0f});
    } else if (tag.starts_with("FP")) {
      _components.push_back(
          Component{ComponentKind::kFixedBinary, tag, _nfixed_binary, 1.0f});
      ++_nfixed_binary;
    } else if (tag.starts_with("NC")) {
      _components.push_back(Component{ComponentKind::kSparse, tag, _nsparse, 1.0f});
      ++_nsparse;
    } else if (tag.starts_with("FC")) {
      _components.push_back(
          Component{ComponentKind::kFixedCounted, tag, _nfixed_counted, 1.0f});
      ++_nfixed_counted;
    }
  }

  if (_components.empty()) {
    std::cerr << "GFPContext::BuildFromTdt:no fingerprint components\n";
    return 0;
  }

  CanonicalizeComponents();
  ComputeHash();
  BuildDefaultActiveList();

  return 1;
}

int
GFPContext::BuildFromSpecs(const std::vector<GFPGeneratorSpec>& specs, bool preprocess) {
  _components.clear();
  _active.clear();
  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;
  _context_hash = 0;
  _standard_generator.reset();

  if (specs.empty()) {
    std::cerr << "GFPContext::BuildFromSpecs:no generator specs\n";
    return 0;
  }

  std::vector<std::unique_ptr<FingerprintGeneratorImplementation>> generators;
  std::vector<std::vector<Component>> generator_components;
  generators.reserve(specs.size());
  generator_components.reserve(specs.size());

  for (const GFPGeneratorSpec& spec : specs) {
    if ((spec.kind() == GeneratorKind::kALogP || spec.kind() == GeneratorKind::kXLogP ||
         spec.kind() == GeneratorKind::kTPSA) &&
        spec.replicates() <= 0) {
      std::cerr
          << "GFPContext::BuildFromSpecs:invalid sparse fingerprint replicate count "
          << spec.replicates() << '\n';
      return 0;
    }
    if (spec.kind() == GeneratorKind::kCATS && spec.max_path_length() < 1) {
      std::cerr << "GFPContext::BuildFromSpecs:invalid CATS max path length "
                << spec.max_path_length() << '\n';
      return 0;
    }
    if (spec.kind() == GeneratorKind::kAtomPair) {
      if (spec.min_separation() < 0 || spec.max_separation() < spec.min_separation() ||
          spec.atom_type().empty() ||
          AtomPairTag(spec.min_separation(), spec.max_separation(), spec.atom_type(),
                      spec.include_out_of_range())
              .empty()) {
        std::cerr << "GFPContext::BuildFromSpecs:invalid atom pair spec " << spec.Repr()
                  << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kECFingerprint) {
      if (spec.radius() < 0 || spec.atom_type().empty() ||
          ECTag(spec.radius(), spec.atom_type()).empty()) {
        std::cerr << "GFPContext::BuildFromSpecs:invalid EC spec " << spec.Repr() << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kSubstructureFingerprint) {
      if (spec.smarts().empty() || spec.radius() < 0 || spec.atom_type().empty() ||
          SubstructureTag(spec.smarts(), spec.radius(), spec.atom_type()).empty()) {
        std::cerr << "GFPContext::BuildFromSpecs:invalid substructure spec "
                  << spec.Repr() << '\n';
        return 0;
      }
    }

    std::unique_ptr<FingerprintGeneratorImplementation> generator = MakeGenerator(spec);
    if (generator == nullptr) {
      std::cerr << "GFPContext::BuildFromSpecs:cannot build " << spec.Repr() << '\n';
      return 0;
    }
    if (spec.kind() == GeneratorKind::kCATS) {
      CATSFingerprintGenerator* cats =
          dynamic_cast<CATSFingerprintGenerator*>(generator.get());
      if (cats == nullptr || !cats->Initialise()) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot initialise " << spec.Repr()
                  << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kAtomPair) {
      AtomPairFingerprintGenerator* atom_pair =
          dynamic_cast<AtomPairFingerprintGenerator*>(generator.get());
      if (atom_pair == nullptr || !atom_pair->Initialise()) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot initialise " << spec.Repr()
                  << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kECFingerprint) {
      ECFingerprintGenerator* ec = dynamic_cast<ECFingerprintGenerator*>(generator.get());
      if (ec == nullptr || !ec->Initialise()) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot initialise " << spec.Repr()
                  << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kSpinachFingerprint ||
        spec.kind() == GeneratorKind::kScaffoldFingerprint) {
      ScaffoldSpinachFingerprintGenerator* scaffold_spinach =
          dynamic_cast<ScaffoldSpinachFingerprintGenerator*>(generator.get());
      if (scaffold_spinach == nullptr || !scaffold_spinach->Initialise()) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot initialise " << spec.Repr()
                  << '\n';
        return 0;
      }
    }
    if (spec.kind() == GeneratorKind::kSubstructureFingerprint) {
      SubstructureFingerprintGenerator* substructure =
          dynamic_cast<SubstructureFingerprintGenerator*>(generator.get());
      if (substructure == nullptr || !substructure->Initialise()) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot initialise " << spec.Repr()
                  << '\n';
        return 0;
      }
    }

    std::vector<Component> components = generator->Components();
    if (components.empty()) {
      std::cerr << "GFPContext::BuildFromSpecs:no components from " << spec.Repr()
                << '\n';
      return 0;
    }

    for (const Component& component : components) {
      for (const Component& existing : _components) {
        if (existing.tag == component.tag) {
          std::cerr << "GFPContext::BuildFromSpecs:duplicate fingerprint tag '"
                    << component.tag << "'\n";
          return 0;
        }
      }
      _components.push_back(component);
    }

    generators.push_back(std::move(generator));
    generator_components.push_back(std::move(components));
  }

  CanonicalizeComponents();
  ComputeHash();
  BuildDefaultActiveList();

  auto standard_generator = std::make_unique<StandardFingerprintGenerator>(preprocess);
  for (int i = 0; i < static_cast<int>(generators.size()); ++i) {
    std::vector<GeneratorComponentAssignment> assignments;
    assignments.reserve(generator_components[i].size());
    for (const Component& component : generator_components[i]) {
      GeneratorComponentAssignment assignment;
      if (!FindAssignedComponent(_components, component, assignment)) {
        std::cerr << "GFPContext::BuildFromSpecs:cannot assign component '"
                  << component.tag << "'\n";
        return 0;
      }
      assignments.push_back(assignment);
    }

    if (!standard_generator->AddGenerator(std::move(generators[i]), assignments)) {
      return 0;
    }
  }

  _standard_generator = std::move(standard_generator);
  return 1;
}

int
GFPContext::BuildStandard(bool preprocess) {
  std::vector<GFPGeneratorSpec> specs;
  specs.reserve(3);
  specs.push_back(GFPGeneratorSpec::IWMFingerprint());
  specs.push_back(GFPGeneratorSpec::MACCSKeys(true));
  specs.push_back(GFPGeneratorSpec::MolecularProperties());

  return BuildFromSpecs(specs, preprocess);
}

std::vector<std::string>
GFPContext::Tags() const {
  std::vector<std::string> result;
  result.reserve(_components.size());
  for (const Component& component : _components) {
    result.emplace_back(component.tag.rawchars(), component.tag.length());
  }

  return result;
}

int
GFPContext::SetWeight(const const_IWSubstring& tag, float weight) {
  if (weight < 0.0f) {
    std::cerr << "GFPContext::SetWeight:invalid negative weight " << weight << '\n';
    return 0;
  }

  int found = 0;
  for (Component& component : _components) {
    if (component.tag == tag) {
      component.weight = weight;
      ++found;
    }
  }

  if (!found) {
    std::cerr << "GFPContext::SetWeight:no component '" << tag << "'\n";
    return 0;
  }

  _active.clear();
  float sum_weight = 0.0f;
  for (const Component& component : _components) {
    if (component.weight <= 0.0f) {
      continue;
    }
    _active.push_back(ActiveComponent{component.kind, component.index, component.weight});
    sum_weight += component.weight;
  }

  if (sum_weight == 0.0f) {
    std::cerr << "GFPContext::SetWeight:all weights are zero\n";
    return 0;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }

  return 1;
}

int
GFPContext::UseOnly(const std::vector<IWString>& tags) {
  _active.clear();

  for (const IWString& tag : tags) {
    bool found = false;
    for (const Component& component : _components) {
      if (component.tag != tag) {
        continue;
      }
      _active.push_back(
          ActiveComponent{component.kind, component.index, component.weight});
      found = true;
      break;
    }
    if (!found) {
      std::cerr << "GFPContext::UseOnly:no component '" << tag << "'\n";
      _active.clear();
      return 0;
    }
  }

  float sum_weight = 0.0f;
  for (const ActiveComponent& active : _active) {
    sum_weight += active.weight;
  }

  if (sum_weight == 0.0f) {
    std::cerr << "GFPContext::UseOnly:zero total weight\n";
    return 0;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }

  return 1;
}

void
GFPContext::UseAll() {
  _active.clear();

  float sum_weight = 0.0f;
  for (const Component& component : _components) {
    if (component.weight <= 0.0f) {
      continue;
    }
    _active.push_back(ActiveComponent{component.kind, component.index, component.weight});
    sum_weight += component.weight;
  }

  if (sum_weight == 0.0f) {
    BuildDefaultActiveList();
    return;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }
}

int
GFPContext::Fingerprint(Molecule& m, GFPFingerprint& result) {
  if (_standard_generator == nullptr) {
    std::cerr << "GFPContext::Fingerprint:context cannot generate fingerprints\n";
    return 0;
  }

  if (result.context_hash() != _context_hash && !result.Allocate(*this)) {
    return 0;
  }

  return _standard_generator->Generate(m, result);
}

float
GFPContext::Distance(const GFPFingerprint& lhs, const GFPFingerprint& rhs) const {
  if (lhs.context_hash() != _context_hash || rhs.context_hash() != _context_hash) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float similarity = 0.0f;
  for (const ActiveComponent& active : _active) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity += active.weight *
                      lhs.molecular_properties().Similarity(rhs.molecular_properties());
        break;
      case ComponentKind::kFixedBinary:
        similarity +=
            active.weight *
            lhs.fixed_binary(active.index).Tanimoto(rhs.fixed_binary(active.index));
        break;
      case ComponentKind::kSparse:
        similarity +=
            active.weight * lhs.sparse(active.index).tanimoto(rhs.sparse(active.index));
        break;
      case ComponentKind::kFixedCounted:
        similarity +=
            active.weight *
            lhs.fixed_counted(active.index).tanimoto(rhs.fixed_counted(active.index));
        break;
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

GFPList::GFPList() : _context(std::make_shared<GFPContext>()) {
}

GFPList::GFPList(std::shared_ptr<GFPContext> context) : _context(std::move(context)) {
}

std::shared_ptr<GFPList>
GFPList::Standard(bool preprocess) {
  auto context = std::make_shared<GFPContext>();
  if (!context->BuildStandard(preprocess)) {
    return nullptr;
  }

  return std::make_shared<GFPList>(context);
}

std::shared_ptr<GFPList>
GFPList::StandardFromMolecules(const std::vector<Molecule*>& molecules, bool preprocess,
                               bool store_metadata) {
  auto result = GFPList::Standard(preprocess);
  if (result == nullptr) {
    return nullptr;
  }

  if (!result->AddMolecules(molecules, store_metadata)) {
    return nullptr;
  }

  return result;
}

int
GFPList::SetStoreMetadata(int store_metadata) {
  const int value = store_metadata ? 1 : 0;
  if (_store_metadata < 0) {
    _store_metadata = value;
    return 1;
  }

  if (_store_metadata == value) {
    return 1;
  }

  std::cerr << "GFPList::SetStoreMetadata:cannot mix metadata and no-metadata entries\n";
  return 0;
}

int
GFPList::Reserve(int size_hint) {
  if (size_hint <= 0) {
    return 1;
  }

  _molecular_properties.reserve(size_hint);
  _fixed_binary.reserve(size_hint * _context->nfixed_binary());
  _sparse.reserve(size_hint * _context->nsparse());
  _fixed_counted.reserve(size_hint * _context->nfixed_counted());
  if (_store_metadata != 0) {
    _smiles.reserve(size_hint);
    _ids.reserve(size_hint);
  }

  return 1;
}

int
GFPList::BuildContextIfNeeded(const IW_TDT& tdt) {
  if (_context->context_hash() != 0) {
    return 1;
  }

  return _context->BuildFromTdt(tdt);
}

int
GFPList::Add(const IW_TDT& tdt) {
  if (!BuildContextIfNeeded(tdt)) {
    return 0;
  }
  if (!SetStoreMetadata(1)) {
    return 0;
  }

  IWString smiles;
  if (!tdt.dataitem_value(kSmilesTag, smiles)) {
    std::cerr << "GFPList::Add:missing smiles\n";
    return 0;
  }
  IWString id;
  if (!tdt.dataitem_value(kIdentifierTag, id)) {
    std::cerr << "GFPList::Add:missing identifier\n";
    return 0;
  }

  for (const Component& component : _context->components()) {
    const_IWSubstring dataitem;
    if (!tdt.dataitem(component.tag, dataitem)) {
      std::cerr << "GFPList::Add:missing component '" << component.tag << "'\n";
      return 0;
    }
    if (dataitem.ends_with('\n')) {
      dataitem.chop();
    }

    switch (component.kind) {
      case ComponentKind::kMolecularProperties: {
        _molecular_properties.emplace_back();
        if (!_molecular_properties.back().Build(dataitem)) {
          return 0;
        }
        break;
      }
      case ComponentKind::kFixedBinary: {
        _fixed_binary.emplace_back();
        if (!_fixed_binary.back().Build(dataitem)) {
          return 0;
        }
        break;
      }
      case ComponentKind::kSparse: {
        _sparse.emplace_back();
        if (!_sparse.back().construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPList::Add:cannot parse sparse '" << component.tag << "'\n";
          return 0;
        }
        break;
      }
      case ComponentKind::kFixedCounted: {
        _fixed_counted.emplace_back();
        if (!_fixed_counted.back().construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPList::Add:cannot parse fixed counted '" << component.tag
                    << "'\n";
          return 0;
        }
        break;
      }
    }
  }

  _smiles.push_back(smiles);
  _ids.push_back(id);
  ++_size;

  return 1;
}

int
GFPList::AddFingerprintComponents(const GFPFingerprint& fp) {
  if (fp.context_hash() != _context->context_hash()) {
    std::cerr << "GFPList::Add:incompatible fingerprint context\n";
    return 0;
  }

  for (const Component& component : _context->components()) {
    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        _molecular_properties.push_back(fp.molecular_properties());
        break;
      case ComponentKind::kFixedBinary:
        _fixed_binary.push_back(fp.fixed_binary(component.index));
        break;
      case ComponentKind::kSparse:
        _sparse.push_back(fp.sparse(component.index));
        break;
      case ComponentKind::kFixedCounted:
        _fixed_counted.push_back(fp.fixed_counted(component.index));
        break;
    }
  }

  ++_size;
  return 1;
}

int
GFPList::Add(const GFPFingerprint& fp) {
  if (!SetStoreMetadata(0)) {
    return 0;
  }

  return AddFingerprintComponents(fp);
}

int
GFPList::Add(const GFPFingerprint& fp, const IWString& smiles, const IWString& id) {
  if (!SetStoreMetadata(1)) {
    return 0;
  }

  if (!AddFingerprintComponents(fp)) {
    return 0;
  }

  _smiles.push_back(smiles);
  _ids.push_back(id);

  return 1;
}

int
GFPList::Add(Molecule& m) {
  return Add(m, true);
}

int
GFPList::Add(Molecule& m, bool store_metadata) {
  GFPFingerprint fp;
  if (!_context->Fingerprint(m, fp)) {
    return 0;
  }

  if (!store_metadata) {
    return Add(fp);
  }

  IWString smiles(m.smiles());
  return Add(fp, smiles, m.name());
}

int
GFPList::AddMolecules(const std::vector<Molecule*>& molecules, bool store_metadata) {
  if (!SetStoreMetadata(store_metadata)) {
    return 0;
  }
  Reserve(molecules.size());

  for (Molecule* molecule : molecules) {
    if (molecule == nullptr) {
      std::cerr << "GFPList::AddMolecules:null molecule pointer\n";
      return 0;
    }
    if (!Add(*molecule, store_metadata)) {
      return 0;
    }
  }

  return 1;
}

int
GFPList::ReadFile(const char* fname, int size_hint) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    std::cerr << "GFPList::ReadFile:cannot open '" << fname << "'\n";
    return 0;
  }

  if (size_hint <= 0) {
    size_hint = CountTdts(input, fname);
    if (size_hint < 0) {
      return 0;
    }
  }

  IW_TDT tdt;
  if (!tdt.next(input)) {
    std::cerr << "GFPList::ReadFile:no TDT records in '" << fname << "'\n";
    return 0;
  }

  if (!BuildContextIfNeeded(tdt)) {
    return 0;
  }

  Reserve(size_hint);

  if (!Add(tdt)) {
    return 0;
  }

  while (tdt.next(input)) {
    if (!Add(tdt)) {
      return 0;
    }
  }

  return 1;
}

float
GFPList::Distance(int i, int j) const {
  float similarity = 0.0f;
  for (const ActiveComponent& active : _context->active_components()) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity +=
            active.weight * _molecular_properties[i].Similarity(_molecular_properties[j]);
        break;
      case ComponentKind::kFixedBinary: {
        const int nfixed = _context->nfixed_binary();
        similarity += active.weight * _fixed_binary[i * nfixed + active.index].Tanimoto(
                                          _fixed_binary[j * nfixed + active.index]);
        break;
      }
      case ComponentKind::kSparse: {
        const int nsparse = _context->nsparse();
        similarity += active.weight * _sparse[i * nsparse + active.index].tanimoto(
                                          _sparse[j * nsparse + active.index]);
        break;
      }
      case ComponentKind::kFixedCounted: {
        const int nfc = _context->nfixed_counted();
        similarity += active.weight * _fixed_counted[i * nfc + active.index].tanimoto(
                                          _fixed_counted[j * nfc + active.index]);
        break;
      }
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

float
GFPList::Distance(const GFPFingerprint& fp, int j) const {
  if (fp.context_hash() != _context->context_hash()) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float similarity = 0.0f;
  for (const ActiveComponent& active : _context->active_components()) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity += active.weight *
                      fp.molecular_properties().Similarity(_molecular_properties[j]);
        break;
      case ComponentKind::kFixedBinary: {
        const int nfixed = _context->nfixed_binary();
        similarity +=
            active.weight * fp.fixed_binary(active.index)
                                .Tanimoto(_fixed_binary[j * nfixed + active.index]);
        break;
      }
      case ComponentKind::kSparse: {
        const int nsparse = _context->nsparse();
        similarity +=
            active.weight *
            fp.sparse(active.index).tanimoto(_sparse[j * nsparse + active.index]);
        break;
      }
      case ComponentKind::kFixedCounted: {
        const int nfc = _context->nfixed_counted();
        similarity +=
            active.weight * fp.fixed_counted(active.index)
                                .tanimoto(_fixed_counted[j * nfc + active.index]);
        break;
      }
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

std::vector<NearestNeighbour>
GFPList::NearestNeighbours(int query, int k) const {
  if (k <= 0 || query < 0 || query >= size()) {
    return {};
  }

  NearestNeighbourMaxHeap heap;
  for (int i = 0; i < size(); ++i) {
    if (i == query) {
      continue;
    }

    MaybeAddNearestNeighbour(heap, k, NearestNeighbour{i, Distance(query, i)});
  }

  return SortedNearestNeighbours(heap);
}

std::vector<NearestNeighbour>
GFPList::NearestNeighbours(const GFPFingerprint& query, int k) const {
  if (k <= 0 || query.context_hash() != _context->context_hash()) {
    return {};
  }

  NearestNeighbourMaxHeap heap;
  for (int i = 0; i < size(); ++i) {
    MaybeAddNearestNeighbour(heap, k, NearestNeighbour{i, Distance(query, i)});
  }

  return SortedNearestNeighbours(heap);
}

std::vector<NearestNeighbour>
GFPList::NearestNeighboursWithinDistance(int query, float max_distance) const {
  std::vector<NearestNeighbour> result;
  if (query < 0 || query >= size() || max_distance < 0.0f) {
    return result;
  }

  for (int i = 0; i < size(); ++i) {
    if (i == query) {
      continue;
    }

    const float distance = Distance(query, i);
    if (distance <= max_distance) {
      result.push_back(NearestNeighbour{i, distance});
    }
  }

  std::sort(result.begin(), result.end(),
            [](const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
              }
              return lhs.index < rhs.index;
            });

  return result;
}

std::vector<NearestNeighbour>
GFPList::NearestNeighboursWithinDistance(const GFPFingerprint& query,
                                         float max_distance) const {
  std::vector<NearestNeighbour> result;
  if (query.context_hash() != _context->context_hash() || max_distance < 0.0f) {
    return result;
  }

  for (int i = 0; i < size(); ++i) {
    const float distance = Distance(query, i);
    if (distance <= max_distance) {
      result.push_back(NearestNeighbour{i, distance});
    }
  }

  std::sort(result.begin(), result.end(),
            [](const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
              }
              return lhs.index < rhs.index;
            });

  return result;
}

}  // namespace gfp_context
