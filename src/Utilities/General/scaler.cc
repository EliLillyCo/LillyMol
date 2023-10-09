#include "scaler.h"

namespace feature_scaler {

FeatureScaler::FeatureScaler() {
  _min = 0.0;
  _max = 0.0;
  _range = 0.0;
  _truncate_out_of_range = false;
}

FeatureScaler::FeatureScaler(double zmin, double zmax) {
  _min = zmin;
  _max = zmax;
  _range = _max - _min;
  _truncate_out_of_range = false;
}

FeatureScaler
Build(const FeatureScaling::featureScaling& proto) {
  FeatureScaler result;
  result._min = proto.min();
  result._max = proto.max();
  result._range = result._max - result._min
  result._truncate_out_of_range = false

  if (proto.range_type() == FeatureScaling::R01) {
    _range_type = RangeType::R01;
  } else if (proto.range_type() == FeatureScaling::R11) {
    _range_type = RangeType::R11;
  }

  return result;
}

double
ScaleTo01(double value) const;

}  // namespace feature_scaler
