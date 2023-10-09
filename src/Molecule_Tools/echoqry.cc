// Read a query and write it.
// Most useful for converting legacy query format to proto.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"

namespace echoqry {

using std::cerr;
using std::endl;

int verbose = 0;

struct JobOptions {
  int forget_originating_smarts = 1;
};

void
Usage(int rc) {
  exit(rc);
}

int
EchoQry(Substructure_Query& query,
        int outer_ndx,
        int inner_ndx,
        const IWString& output_stem) {
  IWString fname;
  fname << output_stem << '.' << outer_ndx << '.' << inner_ndx << "_qry.txtproto";
  return query.WriteProto(fname);
}

int
EchoQry(const resizable_array_p<Substructure_Query>& queries,
        int ndx,
        const IWString& output_stem) {
  for (int i = 0; i < queries.number_elements(); ++i) {
    if (! EchoQry(*queries[i], ndx, i, output_stem)) {
      cerr << "Fatal error processing '" << queries[i]->comment() << "'\n";
      return 0;
    }
  }

  return 1;
}

int
EchoQry(const char * token,
        int ndx,
	const JobOptions& options,
        const IWString& output_stem) {
  resizable_array_p<Substructure_Query> queries;
  if (! process_cmdline_token('*', token, queries, verbose)) {
    cerr << "Cannot instantiate query '" << token << "'\n";
    return 0;
  }

  cerr << "todo " << options.forget_originating_smarts << '\n';
  if (options.forget_originating_smarts) {
    for (Substructure_Query * q : queries) {
      q->ForgetOriginatingSmarts();
    }
  }

  return EchoQry(queries, ndx, output_stem);
}

int
EchoQry(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:S:O:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  JobOptions options;

  if (cl.option_present('O')) {
    const_IWSubstring o;
    for (int i = 0; cl.value('O', o, i); ++i) {
      if (o == "keep_smarts") {
	options.forget_originating_smarts = 0;
	if (verbose) {
          cerr << "Will discard originating smarts\n";
	}
      } else {
	cerr << "Unrecognised -O directive '" << o << "'\n";
	return 1;
      }
    }
  }
  if (cl.number_elements() == 0) {
    cerr << "Must specify name(s) of query file(s)\n";
    Usage(1);
  }

  IWString output_stem = "echoqry";
  if (cl.option_present('S')) {
    cl.value('S', output_stem);
  }

  for (int i = 0; i < cl.number_elements(); ++i) {
    if (! EchoQry(cl[i], i, options, output_stem)) {
      cerr << "Fatal error processing " << cl[i] << '\n';
      return i + 1;
    }
  }

  return 0;
}
}  // namespace echoqry

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return echoqry::EchoQry(argc, argv);
}

