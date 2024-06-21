#ifndef UTILITIES_GENERAL_SCALER_H_
#define UTILITIES_GENERAL_SCALER_H_

// Rescaling floating point values to and from the [0,1] range.

#include "Utilities/General/feature_scaling.pb.h"

namespace feature_scaler {

// A class that can scale values to a given range, and then unscale scaled values.
// The default scaled range is [0,1], but can be set to anything.
template <typename T>
class FeatureScaler {
  private:
    // The min, max and range of the raw data.
    T _min;
    T _max;
    T _range;

    // The corresponding values in the scaled space.
    // For example these might be [0,1] or [-1,1]
    T _scaled_min;
    T _scaled_max;
    T _scaled_range;

    // What to do if a value is outside the range.
    bool _truncate_out_of_range;

  // private functions

    void _default_values();

  public:
    FeatureScaler();

    // Initialized with a min and max.
    FeatureScaler(T zmin, T zmax);

    int Initialise(const FeatureScaling::FeatureScaling& proto);

    // set or reset the range of unscaled values.
    int SetRange(T zmin, T zmax);

    void Initialise01Scaling();
    void Initialise11Scaling();

    int Active() const {
      return _range > T();
    }

    void set_truncate_out_of_range(bool s) {
      _truncate_out_of_range = s;
    }

    // Factory construction from a proto.
    static FeatureScaler Build(const FeatureScaling::FeatureScaling& proto);

    bool InRange(T value) const {
      return value >= _min && value <= _max;
    }

    // Once built, the scaling operations supported.
    // `raw` is converted to a value in [_scaled_min, _scaled_max]
    T Scale(T raw) const;

    // And returning to the original range.
    // `scaled` is converted from a value relative to
    // [_scaled_min, _scaled_max] to a value relative to
    // _min,_max. Note that if a `scaled` is out of range
    // what happens is governed by _truncate_out_of_range.
    T ScaleBackToOrignalRange(T scaled) const;
};

#ifdef FEATURE_SCALER_IMPLEMENTATION

template <typename T>
void
FeatureScaler<T>::_default_values() {
  _min = T();
  _max= T();
  _range = T();

  _truncate_out_of_range = false;

  // by default we scale to [0,1].
  _scaled_min = 0.0;
  _scaled_max = 1.0;
  _scaled_range = 1.0;
}

template <typename T>
FeatureScaler<T>::FeatureScaler() {
  _default_values();
}


template <typename T>
void
FeatureScaler<T>::Initialise01Scaling() {
  _scaled_min = 0.0;
  _scaled_max = 1.0;
  _scaled_range = 1.0;
}

template <typename T>
void
FeatureScaler<T>::Initialise11Scaling() {
  _scaled_min = -1.0;
  _scaled_max = 1.0;
  _scaled_range = 2.0;
}

template <typename T>
FeatureScaler<T>::FeatureScaler(T zmin, T zmax) {
  _default_values();

  _min = zmin;
  _max = zmax;
  _range = _max - _min;
}

template <typename T>
FeatureScaler<T>
FeatureScaler<T>::Build(const FeatureScaling::FeatureScaling& proto) {
  FeatureScaler<T> result;

  result._min = proto.min();
  result._max = proto.max();
  result._range = result._max - result._min;

  switch (proto.range_type()) {
    case FeatureScaling::FeatureScaling::UNSPECIFIED:
      break;
    case FeatureScaling::FeatureScaling::R01:
      result.Initialise01Scaling();
      break;
    case FeatureScaling::FeatureScaling::R11:
      result.Initialise11Scaling();
      break;
    default:
      break;
  }

  return result;
}

template <typename T>
int
FeatureScaler<T>::Initialise(const FeatureScaling::FeatureScaling& proto) {
  _default_values();

  _min = proto.min();
  _max = proto.max();
  _range = _max - _min;

  switch (proto.range_type()) {
    case FeatureScaling::FeatureScaling::UNSPECIFIED:
      break;
    case FeatureScaling::FeatureScaling::R11:
      Initialise11Scaling();
      break;
    case FeatureScaling::FeatureScaling::R01:
      Initialise01Scaling();
      break;
    // switches with proto enums are strange. No action.
    default:
      break;
  }

  return 1;
}

template <typename T>
T
FeatureScaler<T>::Scale(T value) const {
  if (value >= _min && value <= _max) {
    return _scaled_min + (value - _min) / _range * _scaled_range;
  }

#ifdef STATS_ON_OUT_OF_RANGE
  // Out of range, update stats and decide what to do.
  if (value < scaling.min()) {
    out_of_range_low++;
  } else {
    out_of_range_high++;
  }
#endif

  std::cerr << "Scaling " << value << " _truncate_out_of_range " << _truncate_out_of_range << ' ' << _min << ',' << _max << '\n';
  std::cerr << "Scaled " << _scaled_min << ' ' << _scaled_max << '\n';
  if (_truncate_out_of_range) {
    if (value < _min) {
      std::cerr << value << " < " << _min << " returning " << _scaled_min << '\n';
      return _scaled_min;
    } else {
      std::cerr << value << " > " << _max << " returning " << _scaled_max << '\n';
      return _scaled_max;
    }
  }

  std::cerr << "Writing extrapolated form\n";
  // Write extrapolated value.
  if (value < _min) {
    return _scaled_min + (value - _min) / _range * _scaled_range;
  } else {
    return _scaled_max + (value - _max) / _range * _scaled_range;
  }
}

template <typename T>
T
FeatureScaler<T>::ScaleBackToOrignalRange(T scaled) const {

  // std::cerr << "Unscaling " << scaled << " scaled range [" << _scaled_min << ',' << _scaled_max << "]\n";
  if (scaled >= _scaled_min && scaled <= _scaled_max) {
      return _min + (scaled - _scaled_min) / _scaled_range * _range;
  }
  // std::cerr << "Out of range _truncate_out_of_range " << _truncate_out_of_range << '\n';

#ifdef STATS_ON_OUT_OF_RANGE
  // Out of range, update stats and decide what to do.
  if (scaled < _scaled_min) {
    out_of_range_low++;
  } else {
    out_of_range_high++;
  }
#endif

  if (_truncate_out_of_range) {
    if (scaled < _scaled_min) {
      return _min;
    } else {
      return _max;
    }
  }

  // Write extrapolated value.
  return _min + (scaled - _scaled_min) / _scaled_range * _range;
}

#endif   // FEATURE_SCALER_IMPLEMENTATION

}  // namespace feature_scaler

#endif  // UTILITIES_GENERAL_SCALER_H_
