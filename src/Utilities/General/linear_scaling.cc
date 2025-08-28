#include <iostream>

#include "Foundational/iwmisc/proto_support.h"

#include "linear_scaling.h"

#ifdef BUILD_BAZEL
#include "Utilities/General/linear_scaling.pb.h"
#else
#include "linear_scaling.pb.h"
#endif


namespace linear_scaling {

using std::cerr;

LinearScaling::LinearScaling() {
  _intercept = 0.0;
  _slope = 0.0;
  _active = false;
}

int
LinearScaling::Build(const linear_scaling::LinearScalingData& proto) {
  if (proto.has_intercept()) {
    _intercept = proto.intercept();
    _active = true;
  }
  if (proto.has_slope()) {
    _slope = proto.slope();
    _active = true;
  }

  if (! _active) {
    cerr << "LinearScaling::Build:noting specified " << proto.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
LinearScaling::Build(IWString& fname) {
  std::optional<linear_scaling::LinearScalingData> maybe_proto = 
        iwmisc::ReadTextProto<linear_scaling::LinearScalingData>(fname);
  if (! maybe_proto) {
    cerr << "LinearScaling::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  return Build(*maybe_proto);
}

int
LinearScaling::MaybeTransform(double& value) const {
#ifdef DEBUG_MAYBE_TRANSFORM
  cerr << "active " << _active << " intercept " << _intercept << " slopt " << _slope << " value " << value << '\n';
#endif

  if (! _active) {
    return 0;
  }

  value = (value - _intercept) / _slope;

#ifdef DEBUG_MAYBE_TRANSFORM
  cerr << "     updated to " << value << '\n';
#endif

  return 1;
}

}  // namespace LinearScaling
