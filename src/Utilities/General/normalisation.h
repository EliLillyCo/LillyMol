#ifndef IWNORMALISATION_H
#define IWNORMALISATION_H

#include "Foundational/accumulator/accumulator.h"

#define NRML_MIN_TO_MAX 1
#define NRML_UNIT_VARIANCE 2
#define NRML_SPREAD_ZERO 3
#define NRML_MEAN_CENTER 4

class NColumn : private Accumulator<double>
{
  private:

    IWString _descriptor_name;

    int _skip;

    int _missing_values_encountered;

//  To avoid recomputing averages and variances, we store them

    double _average;
    double _variance;
    double _range;

//  Since we may read this info from a file, we need to replicate all
//  the values in the accumulator

    double _minval;
    double _maxval;
    int    _n;

  public:
    NColumn();

    void set_descriptor_name(const const_IWSubstring & d) { _descriptor_name = d;}
    const IWString & descriptor_name() const { return _descriptor_name;}

    int skip() const { return _skip;}
    void set_skip(int s) { _skip = s;}

    int n() const { return Accumulator<double>::n();}

    double minval() const { return _minval;}
    double maxval() const { return _maxval;}

    void set(double mn, double mx, double r) { _minval = mn, _maxval = mx; _range = r;}

    void extra(double f);

    void extra_missing_value() { _missing_values_encountered++;}

    int establish_ranges();
    int establish_range_from_pre_existing_data(const IWString & buffer);

    int scale(double, double &) const;
    int unscale(double, double &) const;

    int report(int, std::ostream &) const;

    int write_scaling_information(int col, IWString_and_File_Descriptor & output) const;
};

extern int  NColumn_values_out_of_range();
extern int  NColumn_report_out_of_range_values();
extern void NColumn_set_report_out_of_range_values(int s);
extern void NColumn_set_truncate_out_of_range_unscalings(int s);
extern void NColumn_set_allow_out_of_range_unscalings(int s);
extern int  NColumn_set_scaling_type(int s);
extern int  NColumn_scaling_type();

extern void NColumn_set_range_min(double);
extern void NColumn_set_range_max(double);

class Command_Line;

extern int parse_normalisation_options(Command_Line & cl,
                             char flag,
                             int verbose);

extern int  parse_normalisation_file_record(const const_IWSubstring & buffer, int & fatal, int verbose = 0);

#endif
