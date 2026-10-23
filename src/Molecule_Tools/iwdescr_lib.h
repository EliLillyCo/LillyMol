#ifndef MOLECULE_TOOLS_IWDESCR_LIB_H_
#define MOLECULE_TOOLS_IWDESCR_LIB_H_

#include <memory>
#include <optional>
#include <span>
#include <string>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/minmaxspc.h"
#include "Foundational/iwmisc/set_or_unset.h"
#include "Foundational/iwstring/iwstring.h"

class Command_Line;
class Molecule;
class DescriptorsToCompute;
class Sparse_Fingerprint_Creator;

// 
class Descriptor : public Set_or_Unset<float> {
 private:
  IWString _name;
  int _active;
  uint64_t _zero_value_count;
  std::optional<float>_default_value;

  int _fingerprint_replicates;
  float _min;
  float _max;
  float _dy;
  int _best_fingerprint;

  Accumulator<double> _stats;

 public:
  Descriptor();

  // Setters and getters for the name.
  void set_name(const char* newname);
  const IWString& name() const { return _name; }
  const IWString& descriptor_name() const { return _name; }

  void set_best_fingerprint(int s) { _best_fingerprint = s; }
  int best_fingerprint() const { return _best_fingerprint; }

  void set_min_max_dy(float v1, float v2, float v3) {
    _min = v1;
    _max = v2;
    _dy = v3;
  }
  void set_min_max_resolution(float v1, float v2, int r);
  int set_range(float dmin, float dmax);

  int produce_fingerprint() const { return _fingerprint_replicates; }
  void set_produce_fingerprint(int s) { _fingerprint_replicates = s; }
  int bit_replicates() const { return _fingerprint_replicates; }
  void produce_fingerprint(int bitnum, Sparse_Fingerprint_Creator&) const;

  int active() const { return _active; }

  int report_statistics(std::ostream&) const;
  void update_statistics(int verbose);

  // Legacy form retained until all call sites pass IWDescr verbosity explicitly.
  void update_statistics();

  void set_default_value(float d) { _default_value = d; }
  void reset();
  void set(int s);
  void set(float s) { Set_or_Unset<float>::set(s); }
  void set(double s) { Set_or_Unset<float>::set(static_cast<float>(s)); }

  Descriptor operator++(int);

  Descriptor& operator=(int s) {
    set(s);
    return *this;
  }
  Descriptor& operator=(float s) {
    set(s);
    return *this;
  }
  Descriptor& operator=(double s) {
    set(s);
    return *this;
  }
};

class Descriptor_Filter {
 private:
  IWString _descriptor_name;

  int _descriptor_number = -1;

  Min_Max_Specifier<float> _cond;

  uint64_t _items_rejected = 0;

 public:
  Descriptor_Filter();

  // `prefix` will be something like "w_"
  // `s` will be something like `natoms.lt.50`
  // `descriptors` is used for establishing _descriptor_number which
  // is the index into the descriptor array for this feature.
  int build(const IWString& prefix, const const_IWSubstring& s,
                const std::span<const Descriptor>& descriptors);

  int Report(std::ostream&) const;

  uint64_t items_rejected() const {
    return _items_rejected;
  }

  int satisfied(const Descriptor*);
};


// Descriptor computation engine for iwdescr.
//
// IWDescr owns descriptor configuration and descriptor calculation state only.
// It does no molecule input, no molecule output, no filtering, no fingerprint
// output formatting, no preprocessing, and no command-line testing harness work.
class IWDescr {
 public:
  IWDescr();
  ~IWDescr();

  IWDescr(const IWDescr&) = delete;
  IWDescr& operator=(const IWDescr&) = delete;
  IWDescr(IWDescr&&) noexcept;
  IWDescr& operator=(IWDescr&&) noexcept;

  // The name of the i'th descriptor.
  const IWString& descriptor_name(int i) const;

  // Initialise descriptor-specific options from a LillyMol Command_Line.
  // Global process settings such as aromaticity and element handling remain
  // outside this class. Caller/output options remain in iwdescr_main.cc.
  int Initialise(Command_Line& cl);

  // Initialise all descriptor families using the standard charge and
  // donor/acceptor queries below LILLYMOL_HOME. Operates quietly.
  int InitialiseAll();

  // Compute descriptors for `m`. `results` must point to storage for at least
  // number_descriptors() float values. IWDescr does no allocation and no I/O.
  //
  // Note: this method may currently modify `m` as part of the legacy
  // descriptor calculation path. In particular, chirality descriptors are
  // computed and then chiral centres may be removed, matching historical
  // iwdescr behaviour. Future cleanup may make this non-mutating.
  //
  // Threading contract: thread-compatible, but not thread-safe. Concurrent
  // calls on different IWDescr instances are permitted. Concurrent calls on
  // the same IWDescr instance are not permitted because some internal
  // descriptor computation state is currently updated during Process().
  // Future cleanup may make this method const/thread-safe, but not today.
  int Process(Molecule& m, float* results);

  int number_descriptors() const;

  // Returns true of descriptor `ndx` is active.
  int descriptor_active(int ndx) const;

  // Read-only access to descriptor metadata. IWDescr retains ownership. This is
  // intentionally narrow so iwdescr_main can support legacy output modes that
  // need Descriptor metadata without moving Descriptor storage out of IWDescr.
  const Descriptor& descriptor(int i) const;
  const Descriptor* descriptor_data() const;
  int valid_descriptor_index(int i) const;
  Descriptor* GetDescriptor(const std::string& name) const;
  // All descriptors.
  std::span<Descriptor> Descriptors() const;

  // DescriptorsToCompute is IWDescr-owned state. These accessors are included
  // now to make the intended API boundary explicit. Callers that only compute
  // descriptors do not need this type; Python bindings may later expose either
  // this object directly or a narrower set of wrapper methods.
  DescriptorsToCompute& mutable_descriptors_to_compute();
  const DescriptorsToCompute& descriptors_to_compute() const;

  // Called at the end of a run. For each active descriptor report summary
  // statistics of the values generated.
  int ReportDescriptorStatistics(std::ostream& output) const;

 private:
  class IWDescrImpl;
  std::unique_ptr<IWDescrImpl> _impl;
};

// This function knows about the likely ranges of each descriptor and will
// fill in that information.
void FillDescriptorExtremeties(std::span<Descriptor>&  d, const int resolution);
// Given a span of Descriptor's, return the index of the one whose name is
// `name`. If `prefix` is not empty, and `dname` starts with `prefix`,
// `prefix` is removed before searching starts.
// Returns -1 if not found.
int DescriptorNumber(const IWString& prefix, const IWString& dname,
        const std::span<const Descriptor>& descriptors);

#endif  // MOLECULE_TOOLS_IWDESCR_LIB_H_
