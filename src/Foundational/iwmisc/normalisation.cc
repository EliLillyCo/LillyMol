#include <cmath>
#include <iostream>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/proto_support.h"

#include "normalisation.h"

using std::cerr;

static int scaling_type = NRML_MIN_TO_MAX;

int
NColumn_scaling_type()
{
  return scaling_type;
}

static int values_out_of_range = 0;

int
NColumn_values_out_of_range()
{
  return values_out_of_range;
}

static int report_out_of_range_values = 1;

int
NColumn_report_out_of_range_values()
{
  return report_out_of_range_values;
}

void
NColumn_set_report_out_of_range_values (int s)
{
  report_out_of_range_values = s;
}

static int truncate_out_of_range_unscalings = 0;

void
NColumn_set_truncate_out_of_range_unscalings (int s)
{
  truncate_out_of_range_unscalings = s;
}

static int allow_out_of_range_unscalings = 0;

void 
NColumn_set_allow_out_of_range_unscalings (int s)
{
  allow_out_of_range_unscalings = s;
}

static double range_min = 0.0;

void 
NColumn_set_range_min(double s)
{
  range_min = s;
}

static double range_max = 1.0;

void 
NColumn_set_range_max(double s)
{
  range_max = s;
}

int 
NColumn_set_scaling_type (int s)
{
  scaling_type = s;
  if (NRML_MIN_TO_MAX == s)
  {
    range_min = 0.0;
    range_max = 1.0;
  }
  else if (NRML_UNIT_VARIANCE == s)
  {
    range_min =  0.0;
    range_max =  1.0;
  }
  else if (NRML_SPREAD_ZERO == s)
  {
    range_min = -1.0;
    range_max =  1.0;
  }
  else if (NRML_MEAN_CENTER == s)    // ranges don't matter
  {
  }
  else
  {
    cerr << "NColumn_set_scaling_type:invalid scaling type " << s << '\n';
    return 0;
  }

  scaling_type = s;

  return 1;
}

NColumn::NColumn()
{
  _skip = 0;
  _range = 0.0;
  _average = 0.0;
  _variance = 0.0;

  _missing_values_encountered = 0;

  return;
}

int
NColumn::report (int col,
                 std::ostream & output) const
{
  output << "Column " << col;

  if (_descriptor_name.length())
    output << ' ' << _descriptor_name;

  output << ' ';

  if (_missing_values_encountered)
    output << _missing_values_encountered << " missing values";

  output << Accumulator<double>::n() << " values";
  if (Accumulator<double>::n() > 1)
    output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << " ave " << Accumulator<double>::average() << " var " << Accumulator<double>::variance();

  output << '\n';

  return output.good();
}

int 
NColumn::write_scaling_information(int col,
                                   IWString_and_File_Descriptor & output) const
{
  if (_skip)
    return 1;

  int n = Accumulator<double>::n();

  if (0 == n)
    return 1;

  if (_descriptor_name.length())
    output << _descriptor_name;
  else
    output << "Column " << col;

  output << ' ';

  assert (n > 1);

  output << Accumulator<double>::minval() << ' ' << Accumulator<double>::maxval() << ' ' << Accumulator<double>::average() << ' ' << sqrt(Accumulator<double>::variance()) << " # min max average variance\n";

  return 1;
}

void
NColumn::extra (double f)
{
  Accumulator<double>::extra (f);

  return;
}

int
NColumn::establish_ranges()
{
  if (0 == _n)
    _n = Accumulator<double>::n();

  assert (_n > 1);

  _range = Accumulator<double>::maxval() - Accumulator<double>::minval();

  if (0.0 == _range)
  {
    cerr << "NColumn::establish_ranges: '" << _descriptor_name << "' no variability in column\n";
    _skip = 1;
    return 0;
  }

  _minval = Accumulator<double>::minval();
  _maxval = Accumulator<double>::maxval();
  _n = Accumulator<double>::n();

  if (NRML_UNIT_VARIANCE == scaling_type)
  {
    _average = Accumulator<double>::average();
    _variance = sqrt(Accumulator<double>::variance());
  }
  else if (NRML_MEAN_CENTER == scaling_type)
    _average = Accumulator<double>::average();

  _skip = 0;

  return 1;
}

int
NColumn::establish_range_from_pre_existing_data(const IWString & buffer)
{
  if (4 != buffer.nwords())
  {
    cerr << "NColumn::establish_ranges:must have exactly 4 words\n";
    cerr << buffer << "\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  if (! token.numeric_value(_minval))
  {
    cerr << "NColumn::establish_ranges:invalid _minval\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_maxval))
  {
    cerr << "NColumn::establish_ranges:invalid _maxval\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_average))
  {
    cerr << "NColumn::establish_ranges:invalid _ave\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_variance))
  {
    cerr << "NColumn::establish_ranges:invalid _variance\n";
    return 0;
  }

  if (_maxval <= _minval)
  {
    cerr << "NColumn::establish_ranges:invalid range specification " << _minval << " to " << _maxval << '\n';
    return 0;
  }

  _range = _maxval - _minval;

  return 1;
}

int
NColumn::scale (double f,
                double & rc) const
{
  if (NRML_MIN_TO_MAX == scaling_type || NRML_SPREAD_ZERO == scaling_type)
    rc = range_min + (f - _minval) / _range * (range_max - range_min);
  else if (NRML_UNIT_VARIANCE == scaling_type)
    rc = range_min + (f - _average) / _variance * (range_max - range_min);
  else if (NRML_MEAN_CENTER == scaling_type)
    rc = f - _average;
  else
  {
    cerr << "NColumn::scale: what kind of scaling " << scaling_type << '\n';
    abort();
  }

  return 1;
}

int
NColumn::unscale (double f,
                  double & rc) const
{
  // cerr << "Scaling type " << scaling_type << ", min " << range_min << " max " << range_max << " value " << f << '\n';
  // cerr << "min " << _minval << " max " << _maxval << '\n';

  if (NRML_MIN_TO_MAX == scaling_type || NRML_SPREAD_ZERO == scaling_type) {
    if (f < range_min || f > range_max) {
      values_out_of_range++;
      if (report_out_of_range_values)
        cerr << "NColumn::unscale:value out of range " << f << " must be " << range_min << " to " << range_max << '\n';

//    cerr << "truncate_out_of_range_unscalings? " << truncate_out_of_range_unscalings << " allow_out_of_range_unscalings? " << allow_out_of_range_unscalings << '\n';
      if (truncate_out_of_range_unscalings)
      {
        if (f < range_min)
        {
          rc = _minval;
          return 1;
        }
        else if (f > range_max)
        {
          rc = _maxval;
          return 1;
        }
      }

      if (! allow_out_of_range_unscalings)
        return 0;
    }

    if (f == range_min)
      rc = _minval;
    else if (f == range_max)
      rc = _maxval;
    else
      rc = _minval + (f - range_min) / (range_max - range_min) * (_maxval - _minval);
  }
  else if (NRML_UNIT_VARIANCE == scaling_type)
    rc = _average + (f - range_min) / (range_max - range_min) * _variance;
  else if (NRML_MEAN_CENTER == scaling_type)
    rc += _average;
  else
  {
    cerr << "NColumn::unscale: what kind of scaling " << scaling_type << '\n';
    abort();
  }

  return 1;
}

int
parse_normalisation_file_record (const const_IWSubstring & buffer,
                                 int & fatal,
                                 int verbose)
{
  const_IWSubstring mybuffer(buffer);

  if (buffer.starts_with("#NTYPE "))
  {
    mybuffer.remove_leading_words(1);
    if ("UV" == mybuffer)
    {
      NColumn_set_scaling_type(NRML_UNIT_VARIANCE);
      if (verbose)
        cerr << "Scaling information file scaling type UV\n";
    }
    else if ("11" == mybuffer)
    {
      NColumn_set_scaling_type(NRML_SPREAD_ZERO);
        if (verbose)
        cerr << "Scaling information file scaling type 11\n";
    }
    else if ("01" == mybuffer)
    {
      NColumn_set_scaling_type(NRML_MIN_TO_MAX);
      if (verbose)
        cerr << "Scaling information file scaling type 01\n";
      range_min = 0.0;
      range_max = 1.0;
    }
    else if ("MC" == mybuffer)
    {
      NColumn_set_scaling_type(NRML_MEAN_CENTER);
      if (verbose)
        cerr << "Scaling information file scaling type mean center\n";
    }
    else
    {
      cerr << "Invalid #NTYPE directive in scaling file '" << buffer << "'\n";
      fatal = 1;
      return 0;
    }

    return 1;
  }

  fatal = 0;

  return 0;
}

int
parse_normalisation_options(Command_Line & cl,
                            char flag,
                            int verbose)
{

  int i = 0;
  const_IWSubstring o;
  while (cl.value(flag, o, i++))
  {
    if ("allow" == o)
    {
      NColumn_set_allow_out_of_range_unscalings(1);
      if (verbose)
        cerr << "Out of range scalings will be allowed\n";
    }
    else if ("trunc" == o)
    {
      NColumn_set_truncate_out_of_range_unscalings(1);
      if (verbose)
        cerr << "Out of range scalings will be truncated\n";
    }
    else if ("quiet" == o)
    {
      NColumn_set_report_out_of_range_values(0);
      if (verbose)
        cerr << "Will not report out of range errors\n";
    }
    else if ("help" == o)
    {
      cerr << " -" << flag << " allow      allow out of range un-normalisations\n";
      cerr << " -" << flag << " trunc      truncate out of range un-normalisations\n";
      cerr << " -" << flag << " quiet      no warning messages about values out of range\n";
    }
    else
    {
      cerr << "Unrecoginsed -o directive '" << o << "'\n";
      return 0;
    }
  }

  return 1;
}

namespace normalisation {

NormalisationBase::NormalisationBase() {
  _n = 0;
  _minval = 0;
  _maxval = 0;
  _average = 0.0;
  _variance = 0.0;
}

int
NormalisationBase::BuildFromProto(const accumulator::AccumulatorData& proto) {
  if (proto.n() == 0) {
    cerr << "NormalisationBase::BuildFromProto:no data in proto " <<
            proto.ShortDebugString() << '\n';
    return 0;
  }

  _n = proto.n();
  _minval = proto.minval();
  _maxval = proto.maxval();
  _average = proto.xsum() / static_cast<double>(_n);
  _variance = (proto.x2sum() - _n * _average * _average) / static_cast<double>(_n - 1);

  return 1;
}

int
NormalisationBase::Build(const Accumulator<double>& acc) {
  if (acc.empty()) {
    cerr << "NormalisationBase::Initialise: empty accumulator\n";
    return 0;
  }

  _n = acc.n();
  _minval = acc.minval();
  _maxval = acc.maxval();
  _average = acc.average();
  _variance = acc.variance();
  _xsum = acc.sum();
  _x2sum = acc.sum_of_squares();

  _range = _maxval - _minval;

  return 1;
}

Normalisation::Normalisation() {
  _normalisation_type = NormalisationType::kNone;
  _n = 0;
  _truncate_out_of_range_scale_requests = 1;
  _truncate_out_of_range_unscale_requests = 0;
  _values_truncated = 0;

  _minval = 0.0;
  _maxval = 0.0;
  _range = 0.0;
  _average = 0.0;
  _variance = 0.0;

  _scaled_min = 0.0;
  _scaled_max = 0.0;
}

void
Normalisation::SetScalingType(NormalisationType t) {
  _normalisation_type = t;

  switch (t) {
    case NormalisationType::k11:
      _scaled_min = -1.0;
      _scaled_max = 1.0;
      break;
    case NormalisationType::k01:
      _scaled_min = 0.0;
      _scaled_max = 1.0;
      break;
    case NormalisationType::kNone:
      cerr << "Normalisation::SetScalingType:No Scaling passed, invalid\n";
      break;
    case NormalisationType::kZeroMean:
      break;
    case NormalisationType::kUnitVariance:
      break;
  }
}

int
Normalisation::BuildFromProto(const normalisation::NormalisationData& proto) {
  if (! proto.has_type()) {
    cerr << "Normalisation::Build:no type " << proto.ShortDebugString() << '\n';
    return 0;
  }
  switch (proto.type()) {
    case normalisation::NTYPE_11:
      _normalisation_type = NormalisationType::k11;
      _scaled_min = -1.0;
      _scaled_max = 1.0;
      break;
    case normalisation::NTYPE_01:
      _normalisation_type = NormalisationType::k01;
      _scaled_min = 0.0;
      _scaled_max = 1.0;
      break;
    case normalisation::NTYPE_UNIT_VARIANCE:
      _normalisation_type = NormalisationType::kUnitVariance;
      break;
    case normalisation::NTYPE_ZERO_MEAN:
      _normalisation_type = NormalisationType::kZeroMean;
      break;
    case normalisation::NTYPE_UNSET:
      cerr << "Normalisation::BuildFromProto:unset normalisation type " << 
              proto.ShortDebugString() << '\n';
      return 0;
    default:
      cerr << "Normalisation::BuildFromProto:unknown type " <<
              proto.ShortDebugString() << '\n';
      return 0;
  }

  if (! proto.has_acc()) {
    cerr << "Normalisation::BuildFromProto:no accumulator data " <<
            proto.ShortDebugString() << '\n';
    return 0;
  }

  const accumulator::AccumulatorData& acc_data = proto.acc();
  if (! acc_data.has_n()) {
    cerr << "Normalisation::BuildFromProto:no n in accumulator data " <<
            proto.ShortDebugString() << '\n';
    return 0;
  }
  if (! acc_data.has_minval() || ! acc_data.has_maxval() || acc_data.has_xsum() ||
      ! acc_data.has_x2sum()) {
    cerr << "Normalisation::BuildFromProto:incomplete accumulator data " <<
            proto.ShortDebugString() << '\n';
    return 0;
  }

  // Some of the normalisation types do not need all these values.
  _n = acc_data.n();
  _minval = acc_data.minval();
  _maxval = acc_data.maxval();
  _range = _maxval - _minval;
  _xsum = acc_data.xsum();
  _x2sum = acc_data.x2sum();
  _average = _xsum / static_cast<double>(_n);
  _variance = (_x2sum - _n * _average * _average) / static_cast<double>(_n - 1);

  return 1;
}

int
Normalisation::Read(IWString& fname) {
  std::optional<normalisation::NormalisationData> maybe_proto = 
        iwmisc::ReadTextProto<normalisation::NormalisationData>(fname);
  if (! maybe_proto) {
    cerr << "Normalisation::Read:cannot read textproto '" << fname << "'\n";
    return 0;
  }

  return BuildFromProto(*maybe_proto);
}

int
Normalisation::BuildProto(normalisation::NormalisationData& proto) const {
  switch(_normalisation_type) {
    case NormalisationType::k01:
      proto.set_type(normalisation::NTYPE_01);
      break;
    case NormalisationType::k11:
      proto.set_type(normalisation::NTYPE_11);
      break;
    case NormalisationType::kUnitVariance:
      proto.set_type(normalisation::NTYPE_UNIT_VARIANCE);
      break;
    case NormalisationType::kZeroMean:
      proto.set_type(normalisation::NTYPE_ZERO_MEAN);
      break;
    case NormalisationType::kNone:
      cerr << "NormalisationType::BuildProto:type not set\n";
      return 0;
    default:
      cerr << "Normalisation::BuildProto:invalid type " << static_cast<int>(_normalisation_type) << '\n';
      return 0;
  }

  cerr << "normalisation building proto N " << _n << '\n';
  accumulator::AccumulatorData& acc = *proto.mutable_acc();
  acc.set_n(_n);
  acc.set_minval(_minval);
  acc.set_maxval(_maxval);
  acc.set_xsum(_xsum);
  acc.set_x2sum(_x2sum);

  return 1;
}

int
Normalisation::Scale(double& v) {

  if (v < _minval || v > _maxval) [[ unlikely ]] {
    if (! _truncate_out_of_range_scale_requests) {
    } else if (v < _minval) {
      v = _minval;
    } else {
      v = _maxval;
    }
  }

  if (_normalisation_type == NormalisationType::k11) {
    v = -1.0 + 2.0 * (v - _minval) / _range;
  } else if (_normalisation_type == NormalisationType::k01) {
    v = 0.0 + (v - _minval) / _range;
  } else if (_normalisation_type == NormalisationType::kZeroMean) {
    v = v - _average;
  } else if (_normalisation_type == NormalisationType::kUnitVariance) {
    v = _minval + (v - _average) / _variance * _range;
    return 1;
  } else {
    cerr << "Normalisation::Scale: type not set\n";
    return 0;
  }

  return 1;
}

int
Normalisation::Unscale(double& v) {
  if (v <= _scaled_max) [[ likely ]] {
  } else if (_truncate_out_of_range_unscale_requests) {
    v = _maxval;
    return 1;
  }
  if (v >= _scaled_min) [[ likely ]] {
  } else if (_truncate_out_of_range_unscale_requests) {
    v = _minval;
    return 1;
  }

  if (_normalisation_type == NormalisationType::k11) {
    v = _minval + (v + 1.0) * 0.5 * _range;
    return 1;
  }
  if (_normalisation_type == NormalisationType::k01) {
    v = _minval + v * _range;
    return 1;
  }

  if (_normalisation_type == NormalisationType::kZeroMean) {
    v = v + _average;
    return 1;
  }

  if (_normalisation_type == NormalisationType::kUnitVariance) {
    v = _average + (v - _scaled_min) / (_scaled_max - _scaled_min) * _variance;
  }

  cerr << "Normalisation::Unscale:unknown type " << static_cast<int>(_normalisation_type) << '\n';
  return 0;
}

}  // namespace normalisation
