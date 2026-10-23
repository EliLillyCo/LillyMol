#ifndef UTILITIES_GFP_TOOLS_GFP_CONTEXT_H_
#define UTILITIES_GFP_TOOLS_GFP_CONTEXT_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwstring/iwstring.h"

#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"
#include "Utilities/GFP_Tools/sparsefp.h"

class IW_TDT;
class Molecule;

namespace gfp_context {

class StandardFingerprintGenerator;

enum class ComponentKind : uint8_t {
  kMolecularProperties,
  kFixedBinary,
  kSparse,
  kFixedCounted,
};

enum class GeneratorKind : uint8_t {
  kMolecularProperties,
  kIWMFingerprint,
  kMACCSKeys,
  kALogP,
  kXLogP,
  kTPSA,
  kFormula,
  kCATS,
  kAtomPair,
  kECFingerprint,
  kRingSubstitution,
  kSpinachFingerprint,
  kScaffoldFingerprint,
  kSubstructureFingerprint,
};

struct Component {
  ComponentKind kind;
  IWString tag;
  int index = 0;
  float weight = 1.0f;
};

class GFPGeneratorSpec {
 private:
  GeneratorKind _kind = GeneratorKind::kMolecularProperties;
  bool _maccs_level2 = true;
  int _replicates = 0;
  int _radius = 0;
  int _min_separation = 0;
  int _max_separation = 0;
  int _max_path_length = 0;
  bool _include_hydrophobic_pairs = true;
  bool _include_out_of_range = false;
  bool _label_join_points = false;
  bool _no_match_is_empty = true;
  IWString _atom_type;
  IWString _smarts;

 public:
  GFPGeneratorSpec() = default;
  GFPGeneratorSpec(GeneratorKind kind, bool maccs_level2 = true, int replicates = 0);
  GFPGeneratorSpec(GeneratorKind kind, int radius, const IWString& atom_type);
  GFPGeneratorSpec(GeneratorKind kind, int max_path_length,
                   bool include_hydrophobic_pairs);
  GFPGeneratorSpec(GeneratorKind kind, int min_separation, int max_separation,
                   const IWString& atom_type, bool include_out_of_range);
  GFPGeneratorSpec(GeneratorKind kind, const IWString& smarts, int radius,
                   const IWString& atom_type, bool no_match_is_empty);

  static GFPGeneratorSpec MolecularProperties();
  static GFPGeneratorSpec IWMFingerprint();
  static GFPGeneratorSpec MACCSKeys(bool level2 = true);
  static GFPGeneratorSpec ALogP(int replicates = 9);
  static GFPGeneratorSpec XLogP(int replicates = 9);
  static GFPGeneratorSpec TPSA(int replicates = 9);
  static GFPGeneratorSpec FormulaFingerprint();
  static GFPGeneratorSpec CATS(int max_path_length = 10,
                               bool include_hydrophobic_pairs = true);
  static GFPGeneratorSpec AtomPair(int min_separation = 1, int max_separation = 10,
                                   const IWString& atom_type = IWString("UST:Y"),
                                   bool include_out_of_range = false);
  static GFPGeneratorSpec ECFingerprint(int radius = 3,
                                        const IWString& atom_type = IWString("UST:Z"));
  static GFPGeneratorSpec RingSubstitution();
  static GFPGeneratorSpec SpinachFingerprint(bool label_join_points = false);
  static GFPGeneratorSpec ScaffoldFingerprint(bool label_join_points = false);
  static GFPGeneratorSpec SubstructureFingerprint(
      const IWString& smarts, int radius = 0,
      const IWString& atom_type = IWString("UST:ARY"),
      bool no_match_is_empty = true);

  GeneratorKind
  kind() const {
    return _kind;
  }

  bool
  maccs_level2() const {
    return _maccs_level2;
  }

  int
  replicates() const {
    return _replicates;
  }

  int
  radius() const {
    return _radius;
  }

  const IWString&
  atom_type() const {
    return _atom_type;
  }

  int
  max_path_length() const {
    return _max_path_length;
  }

  int
  min_separation() const {
    return _min_separation;
  }

  int
  max_separation() const {
    return _max_separation;
  }

  bool
  include_hydrophobic_pairs() const {
    return _include_hydrophobic_pairs;
  }

  bool
  include_out_of_range() const {
    return _include_out_of_range;
  }

  bool
  label_join_points() const {
    return _label_join_points;
  }

  bool
  no_match_is_empty() const {
    return _no_match_is_empty;
  }

  const IWString&
  smarts() const {
    return _smarts;
  }

  std::vector<Component> Components() const;
  std::string Repr() const;
};

struct ActiveComponent {
  ComponentKind kind;
  int index = 0;
  float weight = 1.0f;
};

struct NearestNeighbour {
  int index = -1;
  float distance = 0.0f;
};

class MolecularProperties {
 private:
  int _nproperties = 0;
  uint8_t* _property = nullptr;

 public:
  MolecularProperties() = default;
  ~MolecularProperties();

  MolecularProperties(const MolecularProperties& rhs);
  MolecularProperties& operator=(const MolecularProperties& rhs);

  MolecularProperties(MolecularProperties&& rhs) noexcept;
  MolecularProperties& operator=(MolecularProperties&& rhs) noexcept;

  int Build(const const_IWSubstring& buffer);
  int BuildFromArray(const int* properties, int nproperties);

  int
  nproperties() const {
    return _nproperties;
  }

  const uint8_t*
  properties() const {
    return _property;
  }

  float Similarity(const MolecularProperties& rhs) const;
};

class FixedBitVector {
 private:
  fixed_bit_vector::FixedBitVector _bits;
  int _nset = 0;

 public:
  int Build(const const_IWSubstring& buffer);
  int BuildFromArray(const int* bits, int nbits);

  int
  nbits() const {
    return _bits.nbits();
  }

  int
  nset() const {
    return _nset;
  }

  const fixed_bit_vector::FixedBitVector&
  bits() const {
    return _bits;
  }

  float Tanimoto(const FixedBitVector& rhs) const;
};

class GFPContext;

// A compact standalone/query fingerprint. GFPList stores large collections in
// packed component arrays.
class GFPFingerprint {
 private:
  uint64_t _context_hash = 0;

  MolecularProperties _molecular_properties;
  FixedBitVector* _fixed_binary = nullptr;
  Sparse_Fingerprint* _sparse = nullptr;
  Fixed_Size_Counted_Fingerprint_uchar* _fixed_counted = nullptr;

  void FreeArrays();

 public:
  GFPFingerprint() = default;
  ~GFPFingerprint();

  GFPFingerprint(const GFPFingerprint&) = delete;
  GFPFingerprint& operator=(const GFPFingerprint&) = delete;

  GFPFingerprint(GFPFingerprint&& rhs) noexcept;
  GFPFingerprint& operator=(GFPFingerprint&& rhs) noexcept;

  int Allocate(const GFPContext& context);
  int Build(const IW_TDT& tdt, const GFPContext& context);

  uint64_t
  context_hash() const {
    return _context_hash;
  }

  const MolecularProperties&
  molecular_properties() const {
    return _molecular_properties;
  }

  MolecularProperties&
  mutable_molecular_properties() {
    return _molecular_properties;
  }

  const FixedBitVector&
  fixed_binary(int i) const {
    return _fixed_binary[i];
  }

  FixedBitVector&
  mutable_fixed_binary(int i) {
    return _fixed_binary[i];
  }

  const Sparse_Fingerprint&
  sparse(int i) const {
    return _sparse[i];
  }

  Sparse_Fingerprint&
  mutable_sparse(int i) {
    return _sparse[i];
  }

  const Fixed_Size_Counted_Fingerprint_uchar&
  fixed_counted(int i) const {
    return _fixed_counted[i];
  }

  Fixed_Size_Counted_Fingerprint_uchar&
  mutable_fixed_counted(int i) {
    return _fixed_counted[i];
  }
};

class GFPContext {
 private:
  std::vector<Component> _components;
  std::vector<ActiveComponent> _active;
  std::unique_ptr<StandardFingerprintGenerator> _standard_generator;

  int _nfixed_binary = 0;
  int _nsparse = 0;
  int _nfixed_counted = 0;
  bool _has_molecular_properties = false;

  uint64_t _context_hash = 0;

  void ComputeHash();
  void CanonicalizeComponents();
  void BuildDefaultActiveList();

 public:
  GFPContext();
  ~GFPContext();

  int BuildFromTdt(const IW_TDT& tdt);
  int BuildFromSpecs(const std::vector<GFPGeneratorSpec>& specs, bool preprocess = true);
  int BuildStandard(bool preprocess = true);

  bool
  can_generate_fingerprints() const {
    return _standard_generator != nullptr;
  }

  uint64_t
  context_hash() const {
    return _context_hash;
  }

  int
  nfixed_binary() const {
    return _nfixed_binary;
  }

  int
  nsparse() const {
    return _nsparse;
  }

  int
  nfixed_counted() const {
    return _nfixed_counted;
  }

  bool
  has_molecular_properties() const {
    return _has_molecular_properties;
  }

  const std::vector<Component>&
  components() const {
    return _components;
  }

  const std::vector<ActiveComponent>&
  active_components() const {
    return _active;
  }

  std::vector<std::string> Tags() const;

  int SetWeight(const const_IWSubstring& tag, float weight);
  int UseOnly(const std::vector<IWString>& tags);
  void UseAll();

  int Fingerprint(Molecule& m, GFPFingerprint& result);
  float Distance(const GFPFingerprint& lhs, const GFPFingerprint& rhs) const;
};

class GFPList {
 private:
  std::shared_ptr<GFPContext> _context;
  int _size = 0;
  // -1 unknown/empty, 0 no metadata stored, 1 metadata stored for every entry.
  int _store_metadata = -1;

  std::vector<MolecularProperties> _molecular_properties;
  std::vector<FixedBitVector> _fixed_binary;
  std::vector<Sparse_Fingerprint> _sparse;
  std::vector<Fixed_Size_Counted_Fingerprint_uchar> _fixed_counted;

  std::vector<IWString> _smiles;
  std::vector<IWString> _ids;

  int BuildContextIfNeeded(const IW_TDT& tdt);
  int Add(const IW_TDT& tdt);
  int AddFingerprintComponents(const GFPFingerprint& fp);
  int Add(const GFPFingerprint& fp);
  int Add(const GFPFingerprint& fp, const IWString& smiles, const IWString& id);
  int SetStoreMetadata(int store_metadata);

 public:
  GFPList();
  explicit GFPList(std::shared_ptr<GFPContext> context);

  static std::shared_ptr<GFPList> Standard(bool preprocess = true);
  static std::shared_ptr<GFPList> StandardFromMolecules(
      const std::vector<Molecule*>& molecules, bool preprocess = true,
      bool store_metadata = false);

  int ReadFile(const char* fname, int size_hint = 0);
  int Reserve(int size_hint);
  int Add(Molecule& m);
  int Add(Molecule& m, bool store_metadata);
  int AddMolecules(const std::vector<Molecule*>& molecules, bool store_metadata = false);

  int
  size() const {
    return _size;
  }

  bool
  metadata_stored() const {
    return _store_metadata == 1;
  }

  const GFPContext&
  context() const {
    return *_context;
  }

  GFPContext&
  mutable_context() {
    return *_context;
  }

  const IWString&
  smiles(int i) const {
    return _smiles[i];
  }

  const IWString&
  id(int i) const {
    return _ids[i];
  }

  float Distance(int i, int j) const;
  float Distance(const GFPFingerprint& fp, int j) const;
  std::vector<NearestNeighbour> NearestNeighbours(int query, int k) const;
  std::vector<NearestNeighbour> NearestNeighbours(const GFPFingerprint& query,
                                                  int k) const;
  std::vector<NearestNeighbour> NearestNeighboursWithinDistance(int query,
                                                                float max_distance) const;
  std::vector<NearestNeighbour> NearestNeighboursWithinDistance(
      const GFPFingerprint& query, float max_distance) const;
};

}  // namespace gfp_context

#endif  // UTILITIES_GFP_TOOLS_GFP_CONTEXT_H_
