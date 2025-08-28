#ifndef UTILITIES_GENERAL_LINEAR_SCALING_H_
#define UTILITIES_GENERAL_LINEAR_SCALING_H_

#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#ifdef BUILD_BAZEL
#include "Utilities/General/linear_scaling.pb.h"
#else
#include "linear_scaling.pb.h"
#endif

namespace linear_scaling {

class LinearScaling {
  private:
    double _intercept;
    double _slope;

    bool _active;

  public:
    LinearScaling();

    // Read a Textproto.
    int Build(IWString& fname);
    int Build(const linear_scaling::LinearScalingData& proto);

    void set_slope(double s) {
      _slope = s;
      _active = true;
    }
    void set_intercept(double s) {
      _intercept = s;
      _active = true;
    }

    bool active() const {
      return _active;
    }
  
    int MaybeTransform(double& value) const;
};

}  // namespace linear_scaling

#endif  // UTILITIES_GENERAL_LINEAR_SCALING_H_
