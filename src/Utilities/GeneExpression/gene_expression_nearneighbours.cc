// Nearest neighbours for serialized Profile gene expression data

#include <algorithm>
#include <iostream>
#include <optional>
#include <queue>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/report_progress.h"

#include "gene_expression.h"
#include "needle.h"

namespace gene_expression {

using std::cerr;
using iw_tf_data_record::TFDataReader;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  cerr << R"(Near neighbour finder for serialized Profile gene expression data.
Implements ideas from:
  A simple and robust method for connecting small-molecule drugs using gene-expression signatures
  Shu-Dong Zhang* and Timothy W Gant*      BMC Bioinformatics 2008, 9:258 doi:10.1186/1471-2105-9-258

Input consists of two TFDataRecord files of serialized Profile protos, likely generated by 
`gene_expression_to_proto`.
One file, the -needles file, is a file of needles. The file(s) presented as arguments are the haystack.
Output is to stdout.

For every item in the needles file, the -n closest neighbours from the haystack are accumulated.
  -needles <fname>      file of TFDataRecord serialized Profile protos - needles.
  -Hmaxrank <n>         only consider the top <n> ranked genes in the Haystack (suggest 100 or more)
  -Nmaxrank <n>         only consider the top <n> ranked genes in the Needles
                        If neither -Hmaxrank not -Nmaxrank are specified, all genes in the input
                        files are used.
  -n <n>                number of neighbours to find.
  -maxgeneid <n>        gene ids above <n> are placed in a hash rather than an array.
                        this will lower memory consumption at the expense of compute time.
                        Suggest a value like 100000
  -v                    verbose output.
)";
  // clang-format on
  // clang-format off
  ::exit(rc);
}


using needle::Needle;

class Options {
  private:
    int _verbose;

    // There are two maxrank values. If neither are specified, we use
    // all the data in the input files.
    uint32_t _haystack_maxrank;
    uint32_t _needle_maxrank;

    int _number_needles;
    Needle* _needles;

    uint64_t _haystack_members_read;

    Report_Progress _report_progress;

    // The number of items retrieved with unit similarity.
    uint64_t _exact_matches_found;

    // We can optionally write the gene_id's of the matched genes.
    int _write_matched_genes;

    // Statistics across all needles and all nbrs
    Accumulator<double> _acc_score;

  public:
    Options();
    ~Options();

    int Initialise(Command_Line_v2& cl);
  
    int Compare(const gene_expression::Profile& haystack);

    uint64_t haystack_members_read() const {
      return _haystack_members_read;
    }

    // Non cost because it also gathers nearest neighbour statistics.
    int WriteNeighbours(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _haystack_maxrank = std::numeric_limits<uint32_t>::max();
  _needle_maxrank = std::numeric_limits<uint32_t>::max();

  _haystack_members_read = 0;

  _number_needles = 0;
  _needles = nullptr;

  _write_matched_genes = 0;

  _exact_matches_found = 0;
}

Options::~Options() {
  if (_needles != nullptr) {
    delete [] _needles;
  }
}

uint32_t
RecordsRemaining(TFDataReader& tfdata_reader) {
  for (uint32_t rc = 0; ; ++rc) {
    std::optional<const_IWSubstring> s = tfdata_reader.Next();
    if (! s) {
      tfdata_reader.seek_zero();
      return rc;
    }
  }

  cerr << "RecordsRemaining:should not come to here\n";

  // Should not come here.
  return 0;
}

int
Options::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present("Hmaxrank")) {
    if (! cl.value("Hmaxrank", _haystack_maxrank) || _haystack_maxrank < 1) {
      cerr << "Invalid -Hmaxrank " << _haystack_maxrank << '\n';
      return 0;
    }
    if (_verbose) {
      cerr << "Will consider the " << _haystack_maxrank << " highest scores from the haystack\n";
    }
  }

  if (cl.option_present("Nmaxrank")) {
    if (! cl.value("Nmaxrank", _needle_maxrank) || _needle_maxrank < 1) {
      cerr << "Invalid -Nmaxrank " << _needle_maxrank << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Will consider the " << _needle_maxrank << " highest scores from the needles\n";
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', Needle::_number_neighbours)) {
      cerr << "Invalid number of neighbours (-n)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will keep a max of " << Needle::_number_neighbours << " neighbours\n";
    }
  }

  if (cl.option_present("maxgeneid")) {
    uint32_t tmp;
    cl.value("maxgeneid", tmp);
    if (tmp < 10000) {
      cerr << "Unrealistic value for -maxgeneid " << tmp << '\n';
      return 0;
    }
    needle::set_large_gene_id_threshold(tmp);
    if (_verbose) {
      cerr << "Large gene id threshold " << tmp << '\n';
    }
  }

  if (! cl.option_present("needles")) {
    cerr << "Must specify file of needles via the -needles option\n";
    Usage(1);
  }

  if (cl.option_present("needles")) {
    IWString fname = cl.string_value("needles");

    TFDataReader tfdata_reader(fname);
    if (! tfdata_reader.good()) {
      cerr << "Cannot open needles file '" << fname << "'\n";
      return 0;
    }

    _number_needles = RecordsRemaining(tfdata_reader);
    cerr << "Find " << _number_needles << " needles in input\n";
    if (_number_needles == 0) {
      cerr << "Cannot determine number needles in '" << fname << "'\n";
      return 0;
    }

    _needles = new Needle[_number_needles];

    for (int i = 0; i < _number_needles; ++i) {
      std::optional<gene_expression::Profile> proto = tfdata_reader.ReadProto<gene_expression::Profile>();
      if (! proto) {
        cerr << "Reading needle failed " << i << '\n';
        return 0;
      }

      _needles[i].Build(*proto, _needle_maxrank);
    }

    Accumulator_Int<uint32_t> acc_hash;
    for (int i = 0; i < _number_needles; ++i) {
      acc_hash.extra(_needles[i].GeneIdsInHash());
    }

    if (_verbose) {
      cerr << "Read " << _number_needles << " needles from '" << fname << "'\n";
      cerr << "Have btw " << acc_hash.minval() << " and " << acc_hash.maxval() <<
              "  items in the large gene id hash\n";
    }
  }


  if (cl.option_present("genes")) {
    _write_matched_genes = 1;
    if (_verbose) {
      cerr << "Output will include a list of matched genes\n";
    }
  }

  if (cl.option_present("rpt")) {
    uint32_t rpt;
    cl.value("rpt", rpt);
    _report_progress.set_report_every(rpt);
    if (_verbose) {
      cerr << "Will report progress every " << rpt << " gene expression profiles read\n";
    }
  }

  return 1;
}

int
Options::Compare(const gene_expression::Profile& haystack) {
  ++_haystack_members_read;

  if (_report_progress()) {
    cerr << "Read " << _haystack_members_read << " haystack members\n";
  }

  static int first_call = true;
  if (first_call) {
    uint32_t genes_in_haystack = haystack.gene_size();
    cerr << "First call genes_in_haystack " << genes_in_haystack << '\n';
    if (_haystack_maxrank < genes_in_haystack) {
      _needles[0].SetMaxPossibleAssociation(_haystack_maxrank);
    } else {
      _needles[0].SetMaxPossibleAssociation(genes_in_haystack);
    }

    if (_write_matched_genes) {
      needle::set_accumulate_matching_genes(1);
    }

    first_call = false;
  }


  // cerr << "Haystack has " << haystack.gene_size() << " genes stored\n";
  for (int i = 0; i < _number_needles; ++i) {
    _needles[i].Compare(haystack, _needle_maxrank, _haystack_maxrank);
  }

  return 1;
}

int
Options::WriteNeighbours(IWString_and_File_Descriptor& output) {
  // first work out the max number of genes in common.
  Accumulator_Int<uint32_t> acc_genes_in_common;

  for (int i = 0; i < _number_needles; ++i) {
    uint32_t m = _needles[i].MaxNumberGenes();
    acc_genes_in_common.extra(m);
  }

  if (_verbose && acc_genes_in_common.minval() < acc_genes_in_common.maxval()) {
    cerr << "Genes in common btw " << acc_genes_in_common.minval() <<
            " and " << acc_genes_in_common.maxval() << '\n';
  }

  const uint32_t max_genes_in_common = acc_genes_in_common.maxval();

  for (int i = 0; i < _number_needles; ++i) {
    if (_verbose) {
      _needles[i].UpdateNnbrStatistics(_exact_matches_found, _acc_score);
    }
    _needles[i].WriteNeighbours(max_genes_in_common, output);
  }

  return output.good();
}

int
Options::Report(std::ostream& output) const {
  output << "Report on " << _number_needles << " needles that compared " <<
          _haystack_members_read << " haystack members\n";

  if (_verbose == 0) {
    return output.good();
  }

  output << "Found " << _exact_matches_found << " exact matches\n";
  output << "Scores btw " << _acc_score.minval() << " and " << _acc_score.maxval();
  if (_acc_score.n() > 1) {
     output << " ave " << static_cast<float>(_acc_score.average());
  }
  output << '\n';
  return output.good();
}

int
GeneExpressionNearNeighbours(Options& options,
                             const gene_expression::Profile& proto) {
  return options.Compare(proto);
}

int
GeneExpressionNearNeighbours(Options& options,
                             TFDataReader& input) {
  while (true) {
    std::optional<gene_expression::Profile> maybe_proto = input.ReadProto<gene_expression::Profile>();
    if (! maybe_proto) {
      return 1;
    }

    if (! GeneExpressionNearNeighbours(options, *maybe_proto)) {
      cerr << "Error processing " << maybe_proto->ShortDebugString() << '\n';
      return 0;
    }

    if (options.haystack_members_read() > std::numeric_limits<uint64_t>::max()) {
      cerr << "Early return\n";
      return 1;
    }
  }
}

int
GeneExpressionNearNeighbours(Options& options,
                             const char* fname) {
  TFDataReader input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GeneExpressionNearNeighbours(options, input);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-Hmaxrank=ipos-n=ipos-needles=sfile-rpt=ipos-Nmaxrank=ipos-genes-maxgeneid=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Must specify haystack file(s) as command line arguments\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    return 1;
  }

  for (const char* fname : cl) {
    if (! GeneExpressionNearNeighbours(options, fname)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);

  options.WriteNeighbours(output);

  output.flush();

  options.Report(cerr);

  return 0;
}

}  // namespace gene_expression

int
main(int argc, char** argv) {
  gene_expression::Main(argc, argv);
}