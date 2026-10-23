#ifndef MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_
#define MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_

#include <memory>
#include <optional>
#include <tuple>
#include <string>
#include <string_view>
#include <vector>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rotbond_common.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/qed.h"
#include "Molecule_Tools/xlogp.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/molecule_filter.pb.h"
#else
#include "molecule_filter.pb.h"
#endif

namespace molecule_filter_lib {

enum class RejectionReason : int {
  kPass = 0,
  kZeroAtoms,
  kTooFewAtoms,
  kTooManyAtoms,
  kTooFewRings,
  kTooManyRings,
  kTooFewHeteroatoms,
  kTooManyHeteroatoms,
  kMinHeteroatomFraction,
  kMaxHeteroatomFraction,
  kNonOrganic,
  kIsotope,
  kTooManyChiral,
  kTooFewAromaticRings,
  kTooManyAromaticRings,
  kTooFewAliphaticRings,
  kTooManyAliphaticRings,
  kTooFewHba,
  kTooManyHba,
  kTooFewHbd,
  kTooManyHbd,
  kTooFewSp3Carbon,
  kTooManySp3Carbon,
  kSp3CarbonFractionTooLow,
  kSp3CarbonFractionTooHigh,
  kTooFewHalogen,
  kTooManyHalogen,
  kRingTooSmall,
  kRingTooLarge,
  kTooFewRotatableBonds,
  kTooManyRotatableBonds,
  kAromaticDensityTooLow,
  kAromaticDensityTooHigh,
  kTooShort,
  kTooLong,
  kLowTpsa,
  kHighTpsa,
  kRingSystemTooSmall,
  kRingSystemTooLarge,
  kTooFewAromaticRingsInSystem,
  kTooManyAromaticRingsInSystem,
  kLowAlogp,
  kHighAlogp,
  kLowXlogp,
  kHighXlogp,
  kTooFewChiral,
  kTooFewFragments,
  kTooManyFragments,
  kLowQed,
  kHighQed,
  kLowAmw,
  kHighAmw,
  kTooFewHbaRdkit,
  kTooManyHbaRdkit,
  kTooFewHbdRdkit,
  kTooManyHbdRdkit,
};

enum class Feature : int {
  kNatoms,
  kNrings,
  kHeteroatomCount,
  kHeteroatomFraction,
  kAromaticRingCount,
  kAliphaticRingCount,
  kRotatableBonds,
  kMaxRingSystemSize,
  kAromaticRingsInSystem,
  kTpsa,
  kAlogp,
  kXlogp,
  kHba,
  kHbd,
  kLargestRingSize,
  kHalogenCount,
  kMaxDistance,
  kSp3Carbon,
  kSp3CarbonFraction,
  kAromaticDensity,
  kChiral,
  kNumberFragments,
  kQed,
  kAmw,
  kHbaRdkit,
  kHbdRdkit,
};

std::optional<Feature> FeatureFromName(std::string_view name);
std::string_view FeatureName(Feature feature);

class Utility {
  private:
    std::string _name;
    Feature _feature = Feature::kNatoms;
    float _weight = 1.0f;
    std::vector<double> _x;
    std::vector<double> _y;

  public:
    Utility() = default;

    int BuildFromProto(const molecule_filter_data::Utility& proto);

    const std::string& name() const {
      return _name;
    }

    Feature feature() const {
      return _feature;
    }

    float weight() const {
      return _weight;
    }

    int npoints() const {
      return _x.size();
    }

    // Return the utility value associated with `value`. Values outside the
    // specified range are clamped to the nearest endpoint utility.
    double Value(double value) const;
};

// The calculator objects.
struct FeatureCalculators {
  quick_rotbond::QuickRotatableBonds& rotbond;
  alogp::ALogP& alogp;
  xlogp::XLogPCalc& xlogp;
  nvrtspsa::NovartisPolarSurfaceArea& tpsa;
  qed::Qed& qed;
};

class FeatureValues {
  private:
    Molecule& _m;
    int _matoms;
    int _nrings;

    FeatureCalculators& _calculators;

    std::unique_ptr<int[]> _tmp;

    std::optional<int> _heteroatom_count;
    std::optional<int> _aromatic_ring_count;
    std::optional<int> _rotatable_bonds;
    std::optional<int> _max_ring_system_size;
    std::optional<int> _aromatic_rings_in_system;
    std::optional<double> _tpsa;
    std::optional<double> _alogp_value;
    std::optional<double> _xlogp_value;
    std::optional<int> _hba;
    std::optional<int> _hbd;
    std::optional<int> _halogen_count;
    std::optional<int> _max_distance;
    std::optional<int> _sp3_carbon;
    std::optional<double> _sp3_carbon_fraction;
    std::optional<int> _number_fragments;
    std::optional<double> _qed;
    std::optional<double> _amw;
    std::optional<int> _hba_rdkit;
    std::optional<int> _hbd_rdkit;

    int HeteroatomCount();
    int AromaticRingCount();
    int RotatableBonds();
    void ComputeRingSystem();
    std::optional<double> ALogP();
    std::optional<double> XLogP();
    void ComputeHbaHbd();
    int HalogenCount();
    int MaxDistance();
    void ComputeSp3Carbon();
    int Sp3Carbon();
    double Sp3CarbonFraction();
    int NumberFragments();
    std::optional<double> Qed();
    double Amw();
    void ComputeRdkitHbaHbd();

  public:
    FeatureValues(Molecule& m, int matoms, int nrings, FeatureCalculators& calculators);

    std::optional<double> Value(Feature feature);
};

class MoleculeFilter {
  private:
    quick_rotbond::QuickRotatableBonds _rotbond;

    alogp::ALogP _alogp;

    xlogp::XLogPCalc _xlogp;

    nvrtspsa::NovartisPolarSurfaceArea _tpsa;

    qed::Qed _qed;
    bool _qed_initialised;

    molecule_filter_data::Requirements _requirements;

    std::vector<Utility> _utilities;

    bool _active;

  // private functions
    void InitialiseOptionalFeatures();
    int BuildUtilities();
    bool UsesFeature(Feature feature) const;
    int InitialiseQEDIfNeeded();

    // Bundle the calculators for handing to a FeatureValues.
    FeatureCalculators MakeCalculators();

    // The implementations. Everything public funnels through these so that a
    // caller wanting both filtering and utilities can share one FeatureValues
    // and have each descriptor computed once.
    int Ok(FeatureValues& feature_values, Molecule& m, int matoms, int nrings,
           RejectionReason& rejection_reason);
    int EvaluateUtilities(FeatureValues& feature_values,
                          std::vector<double>& per_feature_utility,
                          double& overall_utility);

  public:
    MoleculeFilter();

    // Read textproto configuration file
    int Build(IWString& fname);

    int active() const {
      return _active;
    }

    int number_utilities() const {
      return _utilities.size();
    }

    const Utility& utility(int i) const {
      return _utilities[i];
    }

    int has_utilities() const {
      return ! _utilities.empty();
    }

    // Compute per-utility values and the overall utility for `m`. Intended to
    // be called only after the molecule has passed hard filters.
    // Prefer OkAndUtilities if you are going to filter as well - calling Ok
    // and then this recomputes every descriptor the two have in common.
    int EvaluateUtilities(Molecule& m, int matoms, int nrings,
                          std::vector<double>& per_feature_utility,
                          double& overall_utility);

    // Apply the hard filters and, if `m` passes, evaluate the utilities,
    // sharing a single FeatureValues so that nothing is computed twice.
    // Returns
    //    0  rejected by a min_/max_ requirement, `rejection_reason` says which
    //    1  passed, and the utilities such as there are have been evaluated
    //   -1  passed the filters, but a utility feature could not be computed
    int OkAndUtilities(Molecule& m, int matoms, int nrings,
                       RejectionReason& rejection_reason,
                       std::vector<double>& per_feature_utility,
                       double& overall_utility);

    // Copy `proto` to _requirements and initialise.
    int Build(const molecule_filter_data::Requirements& proto);

    // Return true if `m` is consistent with the constraints set in `_requirements`.
    int Ok(Molecule& m);
    int Ok(Molecule& m, RejectionReason& rejection_reason);

    // Same as Ok(Molecule&), but uses caller-supplied atom and ring counts. This
    // lets command-line tools preserve cheap SMILES-level screening before
    // constructing a Molecule.
    int Ok(Molecule& m, int matoms, int nrings);
    int Ok(Molecule& m, int matoms, int nrings, RejectionReason& rejection_reason);
};


// Various supporting functions used in the calculation.
int CountHeteroatoms(const Molecule& m);
int AromaticRingCount(Molecule& m);
bool LargestFragment(const const_IWSubstring& smiles,
                const_IWSubstring& largest_frag,
                int& natoms, int& nrings, int& explicit_hydrogens);
std::tuple<int, int> MaxRingSystemSize(Molecule& m, std::unique_ptr<int[]>& tmp);
void RuleOfFive(Molecule & m, int& acceptor, int& donor);
int HalogenCount(const Molecule& m);
int Sp3Carbon(Molecule & m);

} //   namespace molecule_filter_lib

#endif // MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_
