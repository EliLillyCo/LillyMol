// Generate a fingerprint that represents similarity to each of
// an external set of fingerprints.
// This is for the -SIM option to gfp_make.

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "gfp.h"

namespace gfp_to_collection {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Geneates fingerprint from similarity to a set of molecules.
The number of bits generated will be the number of items in the external set.
  -C <fname>    external file to which we compare. The fingerprints must match
  -f            function as a tDT filter
  -v            verbose output
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    // The comparator fingerprints.
    resizable_array_p<IW_General_Fingerprint> _pool;

    IWString _tag;

    int _fingerprints_processed;

    Accumulator<double> _acc;

  // private functions
    int ReadPool(IWString& fname);
    int ReadPool(iwstring_data_source& input);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(IW_General_Fingerprint& fp, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose =0;
  _tag = "NCSIM<";
  _fingerprints_processed = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('C')) {
    cerr << "Must specify comparator set of molecules via the -C option\n";
    return 0;
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    if (! ReadPool(fname)) {
      cerr << "Cannot read comparator fingerprints '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _pool.size() << " comparator fingerprints from '" << fname << "'\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    if (! _tag.starts_with("NC")) {
      cerr << "Tag must be non colliding form '" << _tag << "' invalid\n";
      return 0;
    }

    if (! _tag.ends_with('<')) {
      _tag << '<';
    }

    if (_verbose) {
      cerr << "Will write tag '" << _tag << "'\n";
    }
  }

  return 1;
}

int
Options::ReadPool(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadPool:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadPool(input);
}

int
Options::ReadPool(iwstring_data_source& input) {
  const int n = input.grep("^PCN<");
  if (n == 0) {
    cerr << "Options::ReadPool:no fingerprints in input\n";
    return 0;
  }

  _pool.resize(n);

  IW_TDT tdt;
  while (tdt.next(input)) {
    std::unique_ptr<IW_General_Fingerprint> fp = std::make_unique<IW_General_Fingerprint>();
    int fatal = 0;
    if (fp->construct_from_tdt(tdt, fatal)) {
      _pool << fp.release();
    } else if (fatal) {
      return 0;
    }
  }

  if (_pool.empty()) {
    cerr << "Options::ReadPool:pool is empty\n";
    return 0;
  }

  return n;
}

// Generate a sparse fingerprint based on the similarity of `fp`
// to everything in `_pool`, one bit per member of _pool.
// Write the resulting fingerprint to `output`.
int
Options::Process(IW_General_Fingerprint& fp,
                 IWString_and_File_Descriptor& output) {
  ++_fingerprints_processed;

  Sparse_Fingerprint_Creator sfc;

  const int n = _pool.number_elements();
  for (int i = 0; i < n; ++i) {
    float d = 1.0f - fp.tanimoto(*_pool[i]);
    if (_verbose) {
      _acc.extra(d);
    }
    int c = static_cast<int>(d * 254.0f + 0.4999f);
    sfc.hit_bit(i, c + 1);
  }

  return sfc.write_fingerprint(_tag, output);
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _fingerprints_processed << " fingerprints\n";
  if (_acc.n() == 0) {
    return 1;
  }

  output << "Distances btw " << _acc.minval() << " and " <<
            _acc.maxval() << " mean " << _acc.average() << '\n';
  return 1;
}

int
GfpToCollection(Options& options,
                    IW_General_Fingerprint& fp,
                    IWString_and_File_Descriptor& output) {
  return options.Process(fp, output);
}

int
GfpToCollection(Options& options,
                    iwstring_data_source& input,
                    IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    IW_General_Fingerprint fp;
    int fatal = 0;
    if (fp.construct_from_tdt(tdt, fatal)) {
      // good
    } else if (fatal) {
      cerr << "GfpToCollection error reading TDT " << tdt << '\n';
      return 0;
    } else {
      return 1;
    }

    tdt.write_all_except_vbar(output);

    if (! GfpToCollection(options, fp, output)) {
      cerr << "Error processing " << tdt << '\n';
      return 0;
    }

    output << "|\n";
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
GfpToCollection(Options& options,
                    const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpToCollection:Cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpToCollection(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:C:J:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  if (cl.option_present('f')) {
    GfpToCollection(options, cl[0], output);
  } else {
    for (const char* fname: cl) {
      if (! GfpToCollection(options, fname, output)) {
        cerr << "Error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_to_collection

int
main(int argc, char **argv) {
  int rc = gfp_to_collection::Main(argc, argv);

  return rc;
}
