// Convert csv gene expression data to proto form.

#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GeneExpression/gene_expression.pb.h"

namespace gene_expression_to_proto {

using std::cerr;

using iw_tf_data_record::TFDataWriter;

namespace fs = std::filesystem;

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
  cerr << R"(Converts csv of gene expression data to proto form
-d <suffix>             input files are to be interpreted as directories. All files ending in <suffix> are processed.
-S <fname>              name of output file (a .dat suffix will be added).
-n <number>             max number of items per file if generating multiple files.
-fname                  store only the file name, not the full path name.
-maxrank <n>            only store the top <n> genes - can save a lot of space and time.
-rpt <n>                report progress every <n> items processed
-v                      verbose output
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    // By default, we put all output in one file, but if this is specified, we
    // can create a sequence of files.
    uint32_t _max_items_per_file;

    // We can create a sequence of files if _max_items_per_file is specified.
    // When creating a sequence of files we need to be able to form the name
    // of the next file to be created.
    int _file_index;

    // If creating multiple files, we need to know how many items have
    // been written to the current file.
    uint32_t _items_in_current_file;

    int _ignore_errors;

    std::unique_ptr<TFDataWriter> _output;

    // The file name stem.
    IWString _stem;

    IWString _suffix;

    // We often get full path names as input. If set, store only the file
    // name component.
    int _name_is_file_name;

    // Deliberate choice to make this int rather than uint32_t since
    // proto.gene_size() returns an int.
    int _maxrank;

    uint64_t _protos_written;

    // Keep track of the number of genes written with each proto.
    extending_resizable_array<uint32_t> _gene_size;

    Report_Progress _report_progress;

  // private functions
    int OpenNextFile();

  public:
    Options();
    int Initialise(Command_Line_v2& cl);

    int SetName(const IWString& fname, gene_expression::Profile& proto) const;

    int maxrank() const {
      return _maxrank;
    }

    void GeneSize(int s) {
      ++_gene_size[s];
    }

    int Write(const gene_expression::Profile& proto);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _file_index = 0;
  _items_in_current_file = 0;
  _protos_written = 0;
  _maxrank = std::numeric_limits<int>::max();
  _max_items_per_file = std::numeric_limits<uint32_t>::max();
  _name_is_file_name = 0;
  _ignore_errors = 0;

  _suffix = ".dat";
}

int
Options::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specfy output file via the -S option\n";
    Usage(1);
  }

  if (cl.option_present("ignore_errors")) {
    _ignore_errors = 1;
    if (_verbose) {
      cerr << "Will ignore otherwise fatal errors\n";
    }
  }

  cl.value('S', _stem);
  if (_verbose) {
    cerr << "Output to '" << _stem << "'\n";
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _max_items_per_file)) {
      cerr << "The max number of items per file must be a whole +ve number\n";
      Usage(1);
    }
    if (_verbose) {
      cerr << "Will write a max of " << _max_items_per_file << " items to each output file\n";
    }
  }

  if (cl.option_present("maxrank")) {
    cl.value("maxrank", _maxrank);
    if (_maxrank <= 0) {
      cerr << "Invalid maxrank " << _maxrank << '\n';
      return 0;
    }
    if (_verbose) {
      cerr << "Will write a max of " << _maxrank << " expression values\n";
    }
  }

  if (cl.option_present("fname")) {
    _name_is_file_name = 1;
    if (_verbose) {
      cerr << "Will only store the file name component\n";
    }
  }

  if (cl.option_present("rpt")) {
    uint32_t rpt;
    cl.value("rpt", rpt);
    _report_progress.set_report_every(rpt);
    if (_verbose) {
      cerr << "Will report progress every " << rpt << " protos stored\n";
    }
  }

  return 1;
}

int
Options::OpenNextFile() {
  IWString fname;
  if (_max_items_per_file == std::numeric_limits<uint32_t>::max()) {
    fname << _stem;
    fname.EnsureEndsWith(_suffix);
  } else {
    fname << _stem << _file_index << _suffix;
  }

  if (_output) {
    _output.reset();
  }

  _output = std::make_unique<TFDataWriter>();
  if (! _output->Open(fname)) {
    cerr << "Open::OpenNextFile:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Options::Write(const gene_expression::Profile& proto) {
  if (_verbose > 1) {
    cerr << "Writing " << proto.name() << '\n';
  }

  if (_report_progress()) {
    Report(cerr);
  }

  if (! _output) {
    if (! OpenNextFile()) {
      return 0;
    }
  } else if (_items_in_current_file >= _max_items_per_file) {
    if (! OpenNextFile()) {
      cerr << "Options::Write:cannot open next file " << _file_index << '\n';
      return 0;
    }
    _items_in_current_file = 0;
  }

  if (!_output->WriteSerializedProto<gene_expression::Profile>(proto)) {
    cerr << "Cannot write\n";
    return 0;
  }

  ++_items_in_current_file;
  ++_protos_written;

  return 1;
}

int
Options::SetName(const IWString& fname, gene_expression::Profile& proto) const {
  if (! _name_is_file_name) {
    proto.set_name(fname.data(), fname.size());
    return 1;
  }

  std::string tmp(fname.data(), fname.size());
  // cerr << "Setting name '" << fs::path(fname.data()).filename() << "'\n";
  proto.set_name(fs::path(tmp).filename());
  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Wrote " << _protos_written << " protos\n";
  for (int i = 0; i < _gene_size.number_elements(); ++i) {
    if (_gene_size[i]) {
      output << _gene_size[i] << " protos had " << i << " genes\n";
    }
  }
  return 1;
}

int
AddExpressionData(const const_IWSubstring& line, gene_expression::Profile& proto) {
  int i = 0;
  const_IWSubstring s_gene, s_expression;

  constexpr char kComma = ',';

  if (! line.nextword(s_gene, i, kComma) || ! line.nextword(s_expression, i, kComma) ||
      s_gene.empty() || s_expression.empty()) {
    cerr << "Options::Process:invalid input '" << line << "'\n";
    return 0;
  }

  uint32_t gene;
  if (! s_gene.numeric_value(gene)) {
    cerr << "Options::Process:invalid gene number '" << s_gene << "'\n";
    return 0;
  }

  double expression;
  if (! s_expression.numeric_value(expression)) {
    cerr << "Options::Process:invalid expression '" << s_expression << "'\n";
    return 0;
  }

  gene_expression::Gene* g = proto.add_gene();
  g->set_gene_id(gene);
  g->set_score(static_cast<float>(expression));

  return 1;
}

int
GeneExpressionToProtoInner(IWString& fname, Options& options) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "GeneExpressionToProtoInner:cannot open '" << fname << "'\n";
    return 0;
  }

  if (fname.ends_with(".csv")) {
    fname.chop(4);
  }

  gene_expression::Profile proto;
  options.SetName(fname, proto);

  const_IWSubstring line;

  if (! input.next_record(line)) {
    cerr << "Cannot read header\n";
    return 0;
  }

  while (input.next_record(line)) {
    if (! AddExpressionData(line, proto)) {
      cerr << "Error processing '" << line << "'\n";
      return 0;
    }
    if (proto.gene_size() >= options.maxrank()) {
      break;
    }
  }

  options.GeneSize(proto.gene_size());
  // cerr << "Writing proto with " << proto.gene_size() << " items\n";

  return options.Write(proto);
}

int
GeneExpressionToProtoFile(IWString& fname, Options& options) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "GeneExpressionToProtoFile:cannot open '" << fname << "'\n";
    return 0;
  }

  IWString line;
  while (input.next_record(line)) {
    if (! GeneExpressionToProtoInner(line, options)) {
      cerr << "Fatal error processing '" << line << "'\n";
      return 0;
    }
  }

  return 1;
}

int
GeneExpressionToProto(IWString& fname, Options& options) {
  if (fname.starts_with("F:")) {
    fname.remove_leading_chars(2);
    return GeneExpressionToProtoFile(fname, options);
  } else {
    return GeneExpressionToProtoInner(fname, options);
  }
}

int
GeneExpressionToProtoDir(const std::string& dirname, const IWString& suffix,
                         Options& options) {
  for (const auto& entry : fs::directory_iterator(dirname)) {
    std::filesystem::path fname = entry.path();
    IWString tmp(fname.string());
    if (! tmp.ends_with(suffix)) {
      continue;
    }
    if (! GeneExpressionToProtoInner(tmp, options)) {
      cerr << "GeneExpressionToProtoDir:error processing '" << tmp << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-S=s-n=ipos-rpt=ipos-d=s-fname-maxrank=i-ignore_errors");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  IWString directory_suffix;
  if (cl.option_present('d')) {
    cl.value('d', directory_suffix);
    if (verbose) {
      cerr << "Input files interpreted as directories\n";
    }
  }

  Options options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise output data\n";
    return 1;
  }

  for (const char* tmp : cl) {
    if (directory_suffix.empty()) {
      IWString fname(tmp);
      if (! GeneExpressionToProto(fname, options)) {
        cerr << "Error processing '" << fname << "'\n";
        return 1;
      }
    } else {
      std::string dirname(tmp);
      if (! GeneExpressionToProtoDir(dirname, directory_suffix, options)) {
        cerr << "Error processing directory '" << dirname << "'\n";
        return 1;
      }
    }

  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gene_expression_to_proto

int
main(int argc, char **argv) {
  int rc = gene_expression_to_proto::Main(argc, argv);

  return rc;
}
