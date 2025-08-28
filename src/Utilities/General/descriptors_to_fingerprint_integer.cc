// Insert a precomputed file of descriptors into a fingerprint file.
// this is basically the same as descriptors_to_fingerprint.cc but it is
// restricted to integer values only and the columns cannot be selected.
// And it is more modern code.

#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/absl_hash.h"

namespace descriptors_to_fingerprint {

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
  cerr << R"(Inserts integer descriptors into a fingerprint file.
Takes an existing integer descriptor file and either generates a fingerprint file from that or merges it
as a fingerprint into an existing fingerprint file.

In the first case, the -S option allows specifying the smiles
  descriptors_to_fingerprint_integer -S file.smi file.dat > file.gfp
  Generates
    $SMI<C>
    PCN<methane>
    NCDSCI<...>
    |

Or when adding to an existing fingerprint file
  descriptors_to_fingerprint_integer -f -D file.dat -
  Inserts
    NCDSCI<....>
    into the existing gfp file.

 -f             Working as a TDT filter - inserting an extra record into an existing fingerprint file.
 -J <tag>       Tag for fingerprint. If starts with 'FP' will be a fixed width fingerprint, 'NC' sparse, non colliding.
 -D <fname>     If working as a pipeline, must be specified. Values come from here.
 -S <fname>     Smiles file. Identifiers in the descriptor file are retrieved from here. Not needed with -f option.
 -s <dummy>     If no smiles specified, use 'dummy' as the smiles for each molecule.
 -r <replicate> When writing sparse fingerprints, number of bit replicates.
 -v             verbose output.
)";
  // clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    IWString _tag;

    IWString _smiles_tag;
    IWString _identifier_tag;

    // The -f option.
    int _function_as_tdt_filter;

    absl::flat_hash_map<IWString, std::vector<uint32_t>> _id_to_values;

    // The -S option.
    absl::flat_hash_map<IWString, IWString> _id_to_smiles;

    // Or we can use a constant dummy smiles
    IWString _dummy_smiles;

    uint32_t _columns_in_input;

    // If writing fixed size fingerprints.
    bool _write_fixed_size_fingerprints;
    // The number of bits in fixed size fingerprints - must be a multiple of 8
    int _nbits;

    int _bit_replicates;

    char _input_separator;

  // Private functions.
    int ReadDescriptors(const char* fname);
    int ReadDescriptors(iwstring_data_source& input);
    int ReadDescriptorsRecord(const const_IWSubstring& buffer);

    int ReadSmiles(const char* fname);
    int ReadSmiles(iwstring_data_source& input);
    int ReadSmilesRecord(const const_IWSubstring& buffer);

    void SetNBits();

    int WriteFixedSizeFingerprint(const IWString& id,
                                   const uint32_t* values,
                                   IWString_and_File_Descriptor& output);
    int WriteSparseFingerprint(const IWString& id,
                                const uint32_t* values,
                                IWString_and_File_Descriptor& output);
  public:
    Options();

    int Initialise(Command_Line& cl);

    const IWString& identifier_tag() const {
      return _identifier_tag;
    }

    int function_as_tdt_filter() const {
      return _function_as_tdt_filter;
    }

    int columns_in_input() const {
      return _columns_in_input;
    }

    char input_separator() const {
      return _input_separator;
    }

    // When reading a descriptor file input, we need to examine the first record
    // in order to set _columns_in_input;
    int SetHeader(const const_IWSubstring& header);

    int WriteSmilesId(const IWString& id,
                       IWString_and_File_Descriptor& output) const;

    int Process(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output);

    int WriteFingerprint(const IWString& id,
                         const uint32_t* values,
                         IWString_and_File_Descriptor& output);

    int Report(std::ostream& output);
};

Options::Options() {
  _verbose = 0;

  _tag = "NCDSCI<";

  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";

  _function_as_tdt_filter = 0;

  _columns_in_input = 0;

  _bit_replicates = 1;

  _write_fixed_size_fingerprints = false;
  _nbits = 0;

  _input_separator = ' ';
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('f') && cl.option_present('S')) {
    cerr << "The -f (function as filter) and -S (smiles) options do not make sense together\n";
    Usage(1);
  }

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (_verbose) {
      cerr << "Will function as a TDT filter\n";
    }

    if (! cl.option_present('D')) {
      cerr << "Must specify the descriptors via the -D option\n";
      return 0;
    }

    if (cl.option_present('D')) {
      IWString fname = cl.string_value('D');
      if (! ReadDescriptors(fname)) {
        cerr << "Options::Initialise:cannot read descriptors '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _id_to_values.size() << " descriptor values\n";
    }
  }

  if (cl.option_present('S') && cl.option_present('s')) {
    cerr << "Cannot specify both a file of smiles (-S) and a dummy smiles (-s)\n";
    Usage(1);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! ReadSmiles(fname)) {
      cerr << "Options::Initialise:cannot read smiles '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _id_to_smiles.size() << " smiles\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    _tag.EnsureEndsWith('<');
    if (_tag.starts_with("FP")) {
      _write_fixed_size_fingerprints = 1;
      if (_verbose) {
        cerr << "Will write fixed width fingerprints\n";
      }
    } else if (_tag.starts_with("NC")) {
    } else {
      cerr << "Options::Initialise:invalid tag '" << _tag << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Tag " << _tag << "\n";
    }
  }

  if (cl.option_present('s')) {
    cl.value('s', _dummy_smiles);
    if (_verbose) {
      cerr << "Will write " << _dummy_smiles << " as a dummy smiles\n";
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _bit_replicates) || _bit_replicates < 1) {
      cerr << "Options::Initialise:the number of bit replicates must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will create " << _bit_replicates << " replicates of sparse fingerprints\n";
    }
  }

  if (cl.option_present('i')) {
    IWString tmp;
    cl.value('i', tmp);
    char_name_to_char(tmp, false);
    _input_separator = tmp[0];
  }

  return 1;
}

int
Options::ReadDescriptors(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadDescriptors:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDescriptors(input);
}

// _columns_in_input has been determined. From that, if writing fixed width
// fingerprints, set _nbits;
void
Options::SetNBits() {
  assert(_write_fixed_size_fingerprints);  // does not really matter.

  if ((_columns_in_input / 8) * 8 == _columns_in_input) {
    _nbits = _columns_in_input;
  } else {
    _nbits = ((_columns_in_input / 8) + 1) * 8;
  }

  if (_verbose) {
    cerr << "Options::ReadDescriptors:given " << _columns_in_input << 
            " columns in input, will write fixed width fingerprint with " <<
            _nbits << " bits\n";
  }
}

int
Options::ReadDescriptors(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadDescriptors:cannot read header\n";
    return 0;
  }

  _columns_in_input = buffer.nwords(_input_separator) - 1;

  if (_write_fixed_size_fingerprints) {
    SetNBits();
  }

  while (input.next_record(buffer)) {
    if (! ReadDescriptorsRecord(buffer)) {
      cerr << "Options:ReadDescriptors:error processing line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return _identifier_tag.size();
}

int
Options::ReadDescriptorsRecord(const const_IWSubstring& buffer) {
  IWString id;
  int i = 0;
  if (! buffer.nextword(id, i, _input_separator)) {
    cerr << "Options:;ReadDescriptorsRecord:cannot get identifier\n";
    return 0;
  }

  if (auto iter = _id_to_values.find(id); iter != _id_to_values.end()) {
    cerr << "Options::ReadDescriptorsRecord:duplicate identifier '" << id << "'\n";
    return 0;
  }

  std::vector<uint32_t> values;
  values.reserve(_columns_in_input);

  const_IWSubstring token;
  while (buffer.nextword(token, i, _input_separator)) {
    uint32_t v;
    if (! token.numeric_value(v) || v > std::numeric_limits<uint8_t>::max())  {
      cerr << "Options::ReadDescriptorsRecord:invalid numeric " << token << '\n';
      return 0;
    }
    values.push_back(v);
  }

  if (values.size() != _columns_in_input) {
    cerr << "Options::ReadDescriptorsRecord:column count mismatch, got " << values.size() <<
            " expected " << _columns_in_input << '\n';
    return 0;
  }

  _id_to_values[id] = std::move(values);

  return 1;
}

int
Options::ReadSmiles(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

int
Options::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!ReadSmilesRecord(buffer)) {
      cerr << "Options::ReadSmiles:error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return _id_to_smiles.size();
}

int
Options::ReadSmilesRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString smiles, id;
  if (! buffer.nextword(smiles, i, _input_separator) || 
      ! buffer.nextword(id, i, _input_separator) ||
      id.empty() || smiles.empty()) {
    return 0;
  }

  _id_to_smiles[id] = std::move(smiles);

  return 1;
}

int
Options::Process(const const_IWSubstring& buffer,
                 IWString_and_File_Descriptor& output) {
  IWString id(buffer);
  id.remove_leading_chars(_identifier_tag.length());
  id.pop();

  auto iter = _id_to_values.find(id);
  if (iter == _id_to_values.end()) {
    cerr << "Options::Process:no data for '" << id << "'\n";
    return 0;
  }

  const std::vector<uint32_t>& values = iter->second;

  if (_write_fixed_size_fingerprints) {
    return WriteFixedSizeFingerprint(id, values.data(), output);
  }

  return WriteSparseFingerprint(id, values.data(), output);
}

int
Options::WriteSparseFingerprint(const IWString& id,
                                const uint32_t* values,
                                IWString_and_File_Descriptor& output) {
  Sparse_Fingerprint_Creator sfc;

  if (_bit_replicates == 1) {
    for (uint32_t i = 0; i < _columns_in_input; ++i) {
      if (values[i] == 0) {
        continue;
      }
      sfc.hit_bit(i, values[i]);
    }
  } else {
    for (uint32_t i = 0; i < _columns_in_input; ++i) {
      if (values[i] == 0) {
        continue;
      }

      uint32_t offset = i * _bit_replicates;

      for (int j = 0; j < _bit_replicates; ++j) {
        sfc.hit_bit(offset + j, values[i]);
      }
    }
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tmp);
  output << _tag << tmp << ">\n";

  return 1;
}

int
Options::WriteFixedSizeFingerprint(const IWString& id,
                                   const uint32_t* values,
                                   IWString_and_File_Descriptor& output) {
  IW_Bits_Base fp(_nbits);
  for (uint32_t i = 0; i < _columns_in_input; ++i) {
    if (values[i]) {
      fp.set_bit(i);
    }
  }

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);

  output << _tag << tmp << ">\n";

  return 1;
}

int
Options::SetHeader(const const_IWSubstring& header) {
  _columns_in_input = header.nwords(_input_separator);
  if (_columns_in_input < 2) {
    cerr << "Options::SetHeader:invalid header '" << header << "'\n";
    return 0;
  }

  if (_write_fixed_size_fingerprints) {
    SetNBits();
  }

  return _columns_in_input;
}

int
DescriptorsToFingerprintsFilter(iwstring_data_source& input, Options& options,
                                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if (! buffer.starts_with(options.identifier_tag())) {
      continue;
    }

    if (! options.Process(buffer, output)) {
      return 0;
    }
  }

  return 1;
}

int
AsIntegers(const const_IWSubstring& buffer,
           int i,
           const Options& options,
           uint32_t* destination) {
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i, options.input_separator()); ++col) {
    if (col >= options.columns_in_input()) {
      cerr << "AsIntegers:column count mismatch " << options.columns_in_input() <<  '\n';
      return 0;
    }

    if (! token.numeric_value(destination[col]) || 
        destination[col] > std::numeric_limits<uint8_t>::max()) {
      cerr << "AsIntegers:invalid numeric '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::WriteFingerprint(const IWString& id,
                         const uint32_t* values,
                         IWString_and_File_Descriptor& output) {
  if (_write_fixed_size_fingerprints) {
    return WriteFixedSizeFingerprint(id, values, output);
  }

  return WriteSparseFingerprint(id, values, output);
}

int
Options::WriteSmilesId(const IWString& id,
                       IWString_and_File_Descriptor& output) const {
  if (_dummy_smiles.size() > 0) {
    output << _smiles_tag << _dummy_smiles << ">\n";
  } else if (_id_to_smiles.size() == 0) {
  } else {
    auto iter = _id_to_smiles.find(id);
    if (iter != _id_to_smiles.end()) {
      output << _smiles_tag << iter->second << ">\n";
    } else if (_id_to_smiles.size() > 0) {
      cerr << "Options::WriteSmilesId:no smiles for '" << id << "'\n";
      return 0;
    }
  }

  output << _identifier_tag << id << ">\n";

  return 1;
}

// Reading a descriptor file.
int
DescriptorsToFingerprintsRecord(const const_IWSubstring& buffer,
                          Options& options,
                          IWString_and_File_Descriptor& output) {
  IWString id;
  int i = 0;
  if (! buffer.nextword(id, i, options.input_separator()) || id.empty()) {
    return 0;
  }

  static std::unique_ptr<uint32_t[]> data = std::make_unique<uint32_t[]>(options.columns_in_input());

  if (! AsIntegers(buffer, i, options, data.get())) {
    cerr << "DescriptorsToFingerprintsRecord:invalid numeric data\n";
    return 0;
  }

  if (! options.WriteSmilesId(id, output)) {
    return 0;
  }

  int rc = options.WriteFingerprint(id, data.get(), output);

  output << "|\n";
  return rc;
}
    
// Reading a descriptor file.
int
DescriptorsToFingerprints(iwstring_data_source& input, Options& options,
                          IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "DescriptorsToFingerprints:cannot read header\n";
    return 0;
  }

  if (! options.SetHeader(buffer)) {
    cerr << "DescriptorsToFingerprints:cannot initialise header\n";
    cerr << buffer << '\n';
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! DescriptorsToFingerprintsRecord(buffer, options, output)) {
      cerr << "DescriptorsToFingerprints:error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
DescriptorsToFingerprints(const char* fname, Options& options,
                          IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DescriptorsToFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.function_as_tdt_filter()) {
    return DescriptorsToFingerprintsFilter(input, options, output);
  } else {
    return DescriptorsToFingerprints(input, options, output);
  }
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vD:fS:J:s:r:i:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl)  {
    if (! DescriptorsToFingerprints(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
  }

  return 0;
}

}  // namespace descriptors_to_fingerprint

int
main(int argc, char** argv) {

  int rc = descriptors_to_fingerprint::Main(argc, argv);

  return rc;
}


