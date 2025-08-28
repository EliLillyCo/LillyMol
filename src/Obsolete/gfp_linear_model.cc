// Implement a linear model from the bits in a fingerprint.

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/gfp_linear_model.pb.h"
#else
#include "gfp_linear_model.pb.h"
#endif

namespace gfp_linear_model {

using std::cerr;

void Usage(int rc) {
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Evaluates a gfp linear model build with gfp_linear_model_optimise.
 -C <fname>     a gfp_linear_model::GfpLinearModel textproto model build by gfp_linear_model_optimise.
 -F ...         standard gfp options, enter '-F help' for info.
 -v             verbose output.
)";
  // clang-format on
}

class GfpLinearModelImplementation {
  private:
    double* _weight;
    double _bias;

    IWString _response;

    // the proto from which we are built.
    gfp_linear_model::GfpLinearModel _model;

  public:
    GfpLinearModelImplementation();
    ~GfpLinearModelImplementation();

    int Build(IWString& fname);

    double Score(IW_General_Fingerprint& fp);
};

GfpLinearModelImplementation::GfpLinearModelImplementation() {
  _weight = nullptr;
  _bias = 0.0;
}

GfpLinearModelImplementation::~GfpLinearModelImplementation() {
  if (_weight != nullptr) {
    delete [] _weight;
  }
}

int
GfpLinearModelImplementation::Build(IWString& fname) {
  std::optional<gfp_linear_model::GfpLinearModel> maybe_proto =
        iwmisc::ReadTextProto<gfp_linear_model::GfpLinearModel>(fname);
  if (! maybe_proto) {
    cerr << "GfpLinearModelImplementation::Build:cannot read textproto '" << fname << "'\n";
    return 0;
  }

  int n = maybe_proto->weight_size();
  if (n == 0) {
    cerr << "GfpLinearModelImplementation::Build:no weights\n";
    return 0;
  }

  _weight = new float[n];

  for (int i = 0; i < n; ++i) {
    _weight[i] = maybe_proto->weight(i);
  }
  _bias = maybe_proto->bias();

  _response = maybe_proto->response();

  return 1;

}

int
GfpLinearModelScore(IW_General_Fingerprint& fp,
                    Options& options,
                    IWString_and_File_Descriptor& output) {
  double s = options.score(fp);
  output << gfp.id() << _output_separator << s << '\n';

  return 1;
}

int
GfpLinearModelScore(iwstring_data_source& input,
                    Options& options,
                    IWString_and_File_Descriptor& output) {

  IW_TDT tdt;
  while (tdt.next(input)) {
    int fatal;
    IW_General_Fingerprint fp;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      if (! GfpLinearModelScore(fp, options, output)) {
        return 0;
      }

      output.write_if_buffer_contains_more_than(4096);
    }
  }

  return 1;
}

int
GfpLinearModelScore(const char* fname,
                    Options& options,
                    IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpLinearModelScore:cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpLinearModelScore(input, options, output);
}

int
Main(int argc, char** argv) {
   
  Command_Line cl(argc, argv, "vF:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    ~sage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('C')) {
    cerr << "Must specify model file via the -C option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments, must specify fingerprint input file\n";
    Usage(2);
  }


  for (const char* fname : cl) {
    if (! GfpLinearModel(fname, options, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.report(cerr);
  }

  return 0;
}

}  // namespace  gfp_linear_model

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = gfp_linear_mode::Main(argc, argv);

  return rc;
}
