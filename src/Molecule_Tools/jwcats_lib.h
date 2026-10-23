#ifndef MOLECULE_TOOLS_JWCATS_LIB_H_
#define MOLECULE_TOOLS_JWCATS_LIB_H_

#include <array>
#include <string>
#include <vector>

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/molecule.h"

namespace jwcats {

static constexpr int kNumberProperties = 5;
static constexpr int kPairsPerDistance = 15;

extern const char* const kPairName[kPairsPerDistance];

int PairNumber(int p1, int p2);

enum class ComputeStatus {
  kOk,
  kMissingChargeData,
  kNotInitialised,
  kError,
};

struct Result {
  std::vector<double> scaled_counts;
  std::array<int, kNumberProperties> property_count = {0, 0, 0, 0, 0};
  int number_heavy_atoms = 0;
};

class JWCats {
 private:
  Donor_Acceptor_Assigner _donor_acceptor_assigner;
  Charge_Assigner _charge_assigner;

  int _include_hydrophobic_pairs = 1;
  int _min_bond_separation = 0;
  int _max_bond_separation = 10;
  int _array_size = 0;
  int _initialised = 0;

  int _use_queries_to_determine_hydrophobicity = 0;
  resizable_array_p<Substructure_Hit_Statistics> _queries;

  int _scaling_type = 1;
  int _make_implicit_hydrogens_explicit = 0;

  std::vector<int> _write_array_value;

  void BuildWriteArray();

 public:
  JWCats();

  int Initialise();

  int
  initialised() const {
    return _initialised;
  }

  Donor_Acceptor_Assigner&
  donor_acceptor_assigner() {
    return _donor_acceptor_assigner;
  }

  Charge_Assigner&
  charge_assigner() {
    return _charge_assigner;
  }

  resizable_array_p<Substructure_Hit_Statistics>&
  queries() {
    return _queries;
  }

  void SetIncludeHydrophobicPairs(int s);

  void
  SetMinBondSeparation(int s) {
    _min_bond_separation = s;
    _initialised = 0;
  }

  void SetMaxBondSeparation(int s);

  void
  SetScalingType(int s) {
    _scaling_type = s;
    _initialised = 0;
  }

  void
  SetMakeImplicitHydrogensExplicit(int s) {
    _make_implicit_hydrogens_explicit = s;
    _initialised = 0;
  }

  void
  SetUseQueriesToDetermineHydrophobicity(int s) {
    _use_queries_to_determine_hydrophobicity = s;
    _initialised = 0;
  }

  int
  include_hydrophobic_pairs() const {
    return _include_hydrophobic_pairs;
  }

  int
  min_bond_separation() const {
    return _min_bond_separation;
  }

  int
  max_bond_separation() const {
    return _max_bond_separation;
  }

  int
  array_size() const {
    return _array_size;
  }

  int
  scaling_type() const {
    return _scaling_type;
  }

  const std::vector<int>&
  write_array_value() const {
    return _write_array_value;
  }

  std::vector<std::string> FeatureNames() const;

  // Mutates `m` if charge assignment, explicit hydrogen expansion, or
  // Gasteiger partial charge calculation is active.
  ComputeStatus Compute(Molecule& m, Result& result);
};

}  // namespace jwcats

#endif  // MOLECULE_TOOLS_JWCATS_LIB_H_
