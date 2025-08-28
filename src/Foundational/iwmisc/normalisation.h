#ifndef FOUNDATIONAL_IWMISC_NORMALISATION_H_
#define FOUNDATIONAL_IWMISC_NORMALISATION_H_

#include "Foundational/accumulator/accumulator.h"

#ifdef BUILD_BAZEL
#include "Foundational/iwmisc/normalisation.pb.h"
#else
#include "normalisation.pb.h"
#endif

class Command_Line;

namespace normalisation {

enum class NormalisationType : int {
    kNone,
    k01,
    k11,
    kUnitVariance,
    kZeroMean
};

class NormalisationBase {
  protected:
    // When instantiating from an Accumulator, these values are determined.
    // These are the raw value ranges.
    // The classes that inherit from this will use some or all of these values.
    uint64_t _n;
    double _minval;
    double _maxval;
    double _average;
    double _variance;
    // These will be copied from the Accumulator.
    double _xsum;
    double _x2sum;

    // Once _minval and _maxval are known, either from reading from a file or
    // when EstablishRanges is called, compute the range (_maxval - _minval).
    double _range;

    // During scaling requests what do we do when we encounter a value
    // that is outside our range?
    // By default, we truncate.
    int _truncate_out_of_range_scale_requests;

    // What to do during unscaling if values are out of range.
    // By default we do not truncate.
    int _truncate_out_of_range_unscale_requests;
    uint64_t _values_truncated;

  public:
    NormalisationBase();

    int BuildFromProto(const accumulator::AccumulatorData& proto);
    int Build(const Accumulator<double>& acc);

    void set_truncate_out_of_range_scale_requests(int s) {
      _truncate_out_of_range_scale_requests = s;
    }
    void set_truncate_out_of_range_unscale_requests(int s) {
      _truncate_out_of_range_unscale_requests = s;
    }
};

// Class for normalising data.
class Normalisation : public NormalisationBase {
  private:
    // Acknowledge that it is a bad design choice to have a class that for most
    // methods branches on what the value of this variable is. But this is a
    // trivial class and don't want a profusion of classes.
    // By having a base class the door is left open to implementing multiple classes.
    NormalisationType _normalisation_type;

    // The scaled range - typically [-1,1] or [0,1]
    double _scaled_min;
    double _scaled_max;

    // private functions

  public:
    Normalisation();

    // Methods relating to protos.
    int BuildFromProto(const normalisation::NormalisationData& proto);
    int Read(IWString& fname);
    int BuildProto(normalisation::NormalisationData& proto) const;
    int WriteTextproto(IWString& fname) const;
    int WriteBinary(IWString& fname) const;

    void SetScalingType(NormalisationType t);

    // Change `v` from the original range to the compressed range.
    int Scale(double& v);
    int Scale(double input, double& output);

    // Translate from the compressed range back to native.
    int Unscale(double& v);
    int Unscale(double input, double& output);

};

}  // namespace normalisation

// Crufty horrible code that needs to be re-implemented with a proto
// and without file scope static variables. TODO:ianwatson

#define NRML_MIN_TO_MAX 1
#define NRML_UNIT_VARIANCE 2
#define NRML_SPREAD_ZERO 3
#define NRML_MEAN_CENTER 4
class NColumn : public Accumulator<double> {
  private:

    IWString _descriptor_name;

    // If inactivated.
    int _skip;

    int _missing_values_encountered;

//  To avoid recomputing averages and variances, we store them

    double _average;
    double _variance;
    double _range;

//  Since we may read this info from a file, we need to replicate all
//  the values in the accumulator

    int    _n;
    double _minval;
    double _maxval;


  public:
    NColumn();
    void set_descriptor_name (const const_IWSubstring & d) { _descriptor_name = d;}
    const IWString & descriptor_name() const { return _descriptor_name;}

    int skip() const { return _skip;}
    void set_skip (int s) { _skip = s;}

    int n() const { return Accumulator<double>::n();}

    double minval() const { return _minval;}
    double maxval() const { return _maxval;}

    void   set (double mn, double mx, double r) { _minval = mn, _maxval = mx; _range = r;}

    void extra (double f);

    void extra_missing_value() { _missing_values_encountered++;}

    int establish_ranges();
    int establish_range_from_pre_existing_data(const IWString & buffer);

    int scale (double, double &) const;
    int unscale (double, double &) const;

    int report (int, std::ostream &) const;

    int write_scaling_information(int col, IWString_and_File_Descriptor & output) const;
};

extern int  NColumn_values_out_of_range();
extern int  NColumn_report_out_of_range_values();
extern void NColumn_set_report_out_of_range_values (int s);
extern void NColumn_set_truncate_out_of_range_unscalings (int s);
extern void NColumn_set_allow_out_of_range_unscalings (int s);
extern int  NColumn_set_scaling_type (int s);
extern int  NColumn_scaling_type();

extern void NColumn_set_range_min(double);
extern void NColumn_set_range_max(double);

extern int parse_normalisation_options (Command_Line & cl,
                             char flag,
                             int verbose);

extern int  parse_normalisation_file_record(const const_IWSubstring & buffer, int & fatal, int verbose = 0);

#endif  // FOUNDATIONAL_IWMISC_NORMALISATION_H_
