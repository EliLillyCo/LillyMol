// Randomise the column order of a file

#include <iostream>
#include <numeric>
#include <random>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"

namespace random_column_order {
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
  cerr << R"(Randomise the column order of a file.
 -j             file is a descriptor file
 -N <n>         generate <n> random variants
 -S <stem>      when generating multiple variants (-N) use <stem> as the output file name stem
 -U <str>       suffix for generated files
 -i <char>      input token separator
 -o <char>      output token separator
 -v             verbose output
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;
    char _input_separator;
    char _output_separator;
    int _is_descriptor_file;

    int _number_variants;
    IWString _stem;
    IWString _suffix;

    // Reflecting current processing
    int _lines_read;

  // Private functions
    int RandomColumnOrderOne(const IWString& header,
                              int columns_in_input,
                              iwstring_data_source& input,
                              IWString_and_File_Descriptor& output);
    int RandomColumnOrderMany(const IWString& header,
                               int columns_in_input,
                               iwstring_data_source& input);
    std::unique_ptr<uint32_t[]> GenerateCrossReference(int columns_in_input);
    template <typename T>
    int WriteRecord(const T& buffer,
            int columns_in_input,
            const uint32_t* xref,
            IWString* column,
            IWString_and_File_Descriptor& output);

  public:
    Options();

    int Initialise(Command_Line& cl);

    // Process all records in `input`.
    int RandomColumnOrder(iwstring_data_source& input, 
                  std::unique_ptr<IWString_and_File_Descriptor>& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _input_separator = ' ';
  _output_separator = ' ';
  _is_descriptor_file = 0;

  _number_variants = 0;

  _lines_read = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('j')) {
    _is_descriptor_file = 1;
    if (_verbose) {
      cerr << "Will treat input as a descriptor file\n";
    }
  }

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i)) {
      cerr << "Invalid input token separator (-i) '" << i << "'\n";
      return 0;
    }
  }

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o)) {
      cerr << "Invalid output token separator (-o) '" << o << "'\n";
      return 0;
    }
  }

  if (! cl.option_present('N') && ! cl.option_present('S')) {
  } else if (cl.option_present('N') && cl.option_present('S')) {
  } else {
    cerr << "When generating multiple variants must specify both -S and -N options\n";
    Usage(1);
  }

  if (cl.option_present('N')) {
    if (! cl.value('N', _number_variants) || _number_variants < 1) {
      cerr << "The number of variants (-N) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will generate " << _number_variants << " shuffled variants\n";
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Will generate variants with file name stem '" << _stem << "'\n";
    }
  }

  if (cl.option_present('U')) {
    cl.value('U', _suffix);
    if (_verbose) {
      cerr << "Files generated with suffix '" << _suffix << "'\n";
    }
  }

  return 1;
}

int
Options::RandomColumnOrder(iwstring_data_source& input,
                  std::unique_ptr<IWString_and_File_Descriptor>& output) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Cannot read header record\n";
    return 0;
  }

  _lines_read = 1;

  int columns_in_input = buffer.nwords_single_delimiter(_input_separator);
  if (columns_in_input < 2) {
    cerr << "Only " << columns_in_input << " columns in input, cannot process\n";
    return 0;
  }

  if (_number_variants == 0) {
    return RandomColumnOrderOne(buffer, columns_in_input, input, *output.get());
  }

  return RandomColumnOrderMany(buffer, columns_in_input, input);
}

std::unique_ptr<uint32_t[]>
Options::GenerateCrossReference(int columns_in_input) {
  std::random_device rd;
  std::mt19937 gen(rd());

  // The result.
  std::unique_ptr<uint32_t[]> xref = std::make_unique<uint32_t[]>(columns_in_input);

  if (_verbose) {
    cerr << "Input contains " << columns_in_input << " columns\n";
  }

  std::iota(xref.get(), xref.get() + columns_in_input, 0);
  if (_is_descriptor_file) {
    std::shuffle(xref.get() + 1, xref.get() + columns_in_input, gen);
  } else {
    std::shuffle(xref.get(), xref.get() + columns_in_input, gen);
  }

  return xref;
}

template <typename T>
int
Options::WriteRecord(const T& buffer,
            int columns_in_input,
            const uint32_t* xref,
            IWString* column,
            IWString_and_File_Descriptor& output) {
  for (int i = 0; i < columns_in_input; ++i) {
    column[i].resize_keep_storage(0);
  }

  const_IWSubstring token;
  int i = 0;
  for (int col = 0; buffer.nextword_single_delimiter(token, i, _input_separator); ++col) {
    if (col > columns_in_input) {
      cerr << "Non tabular file detected, expect " << columns_in_input << " columns\n";
      cerr << buffer << '\n';
      return 0;
    }
    column[xref[col]] = token;
  }

  output << column[0];
  for (int i = 1; i < columns_in_input; ++i) {
    output << _output_separator << column[i];
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
Options::RandomColumnOrderOne(const IWString& header,
                              int columns_in_input,
                              iwstring_data_source& input,
                              IWString_and_File_Descriptor& output) {

  std::unique_ptr<uint32_t[]> xref = GenerateCrossReference(columns_in_input);
  std::unique_ptr<IWString[]> tmp = std::make_unique<IWString[]>(columns_in_input);

  WriteRecord(header, columns_in_input, xref.get(), tmp.get(), output);
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    ++_lines_read;
    WriteRecord(buffer, columns_in_input, xref.get(), tmp.get(), output);
  }

  return 1;
}

int
Options::RandomColumnOrderMany(const IWString& header,
                               int columns_in_input,
                               iwstring_data_source& input) {
  std::unique_ptr<IWString[]> tmp = std::make_unique<IWString[]>(columns_in_input);

  const off_t begin_file = input.tellg();

  for (int i = 0; i < _number_variants; ++i) {
    std::unique_ptr<uint32_t[]> xref = GenerateCrossReference(columns_in_input);

    IWString fname(_stem);
    fname << (i + 1) << _suffix;

    IWString_and_File_Descriptor output;
    if (! output.open(fname)) {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }

    if (! input.seekg(begin_file)) {
      cerr << "Cannot seek back to " << begin_file << '\n';
      return 0;
    }

    WriteRecord(header, columns_in_input, xref.get(), tmp.get(), output);
  
    const_IWSubstring buffer;
    while (input.next_record(buffer)) {
      ++_lines_read;
      WriteRecord(buffer, columns_in_input, xref.get(), tmp.get(), output);
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _lines_read << " lines\n";
  return 1;
}

int
RandomColumnOrder(const char* fname, Options& options,
                  std::unique_ptr<IWString_and_File_Descriptor>& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "RandomColumnOrder:cannot open '" << fname << "'\n";
    return 0;
  }

  return options.RandomColumnOrder(input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vN:S:ji:o:U:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments, must specify input file\n";
    Usage(1);
  }

  if (cl.size() > 1) {
    cerr << "Sorry, do not know how to process multiple input files\n";
    return 1;
  }

  std::unique_ptr<IWString_and_File_Descriptor> output;
  if (! cl.option_present('S')) {
    output.reset(new IWString_and_File_Descriptor(1));
  }

  // a loop despite the fact we can only process one file.
  for (const char* fname : cl) {
    if (! RandomColumnOrder(fname, options, output)) {
      cerr << "RandomColumnOrder:error processsing '" << fname << "'\n";
      return 1;
    }
  }

  if (output) {
    output->flush();
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace random_column_order

int
main(int argc, char** argv) {
  return random_column_order::Main(argc, argv);
}
