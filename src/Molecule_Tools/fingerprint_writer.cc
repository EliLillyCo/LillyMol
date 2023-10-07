#include <algorithm>
#include <iostream>
#include <utility>

#include "fingerprint_writer.h"

namespace fingerprint_writer {

using std::cerr;

int
OkTag(const IWString& tag) {
  if (tag.starts_with("FP")) {
    return 1;
  }
  if (tag.starts_with("NC")) {
    return 1;
  }

  return 0;
}

FingerprintWriter::FingerprintWriter()  {
 _output_vector = nullptr;
 _output_type = OutputType::kSparse;
}

FingerprintWriter::~FingerprintWriter() {
  if (_output_vector != nullptr) {
    delete [] _output_vector;
  }
}

void
FingerprintWriter::DisplayUsage(char flag,
                          std::ostream& output) const {
  output << " -" << flag << " <tag>          tag for fingerprint output, same as TAG=<tag>\n";
  output << " -" << flag << " sep=<sep>      output column separator when writing a descriptor file\n";
  output << " -" << flag << " name=<name>    prefix for feature names when writing a descriptor file (default 'fp')\n";
  output << " -" << flag << " nbits=<nbits>  number of bits in fixed width fingerprint, or columns in descriptor file\n";
  output << " -" << flag << " array          output is descriptor file form\n";
  output << " -" << flag << " fixed          output is fixed width fingerprint form\n";
  output << " -" << flag << " sparse         output is sparse fingerprint form\n";
}

int
FingerprintWriter::Initialise(Command_Line& cl,
                        char flag,
                        int verbose) {
  const_IWSubstring opt;
  for (int i = 0; cl.value(flag, opt, i); ++i) {
    if (opt.starts_with("sep=")) {
        opt.remove_leading_chars(4);
        _output_column_separator = opt;
      _output_type = OutputType::kDescriptor;
    } else if (opt.starts_with("name=")) {
      opt.remove_leading_chars(5);
      _feature_prefix = opt;
      _output_type = OutputType::kDescriptor;
    } else if (opt == "help") {
      DisplayUsage(flag, cerr);
      return 0;
    } else if (opt.starts_with("nbits=")) {
      opt.remove_leading_chars(6);
      if (! opt.numeric_value(_nbits) || _nbits < 64) {
        cerr << "Invalid array output width '" << opt << "'\n";
        return 0;
      }
    } else if (opt == "array") {
      _output_type = OutputType::kDescriptor;
      if (_nbits == 0) {
        _nbits = 2048;
      }
    } else if (opt == "fixed") {
      _output_type = OutputType::kFixed;
    } else if (opt == "sparse") {
      _output_type = OutputType::kSparse;
    } else if (opt == "svml") {
      _output_type = OutputType::kSvml;
    } else if (opt.starts_with("tag=")) {
      opt.remove_leading_chars(4);
      if (! OkTag(opt)) {
        cerr << "FingerprintWriter::Initialise:invalid fingerprint '" << opt << "'\n";
        return 0;
      }
      _tag = opt;
    } else {
      if (! OkTag(opt)) {
        cerr << "FingerprintWriter::Initialise:invalid fingerprint '" << opt << "'\n";
        return 0;
      }
      _tag = opt;
      ChangeOutputTypeIfNeeded(opt);
    }
  }

  // If tag has been specified, it must be of valid form.

  if (_tag.empty()) {
    if (_output_type == OutputType::kDescriptor) {
    } else if (_output_type == OutputType::kSvml) {
    } else {
      cerr << "FingerprintWriter:Initialise:no tag, but not descriptor output\n";
      return 0;
    }
  } else if (! OkTag(_tag)) {
    cerr << "FingerprintWriter::Initialise:unrecognised fingerprint tag form '" << opt << "', must begin NC or FP\n";
    return 0;
  } else if (! _tag.ends_with("<")) {
    _tag << '<';
  }

  if (_output_type == OutputType::kDescriptor) {
    _output_vector = new int[_nbits];
    _iwdigits.initialise(100);
    _iwdigits.set_leading_string(_output_column_separator);
  } else if (_output_type == OutputType::kFixed) {
    _fixed.resize(_nbits);
  } else if (_output_type == OutputType::kSvml) {
    _iwdigits.initialise(100);
    _iwdigits.set_leading_string(':');
  }

  if (_tag.starts_with("FP") && _output_type == OutputType::kSparse) {
    cerr << "FingerprintWriter::Initialise:sparse fingerprints must start with NC\n";
    return 0;
  }

  if (_tag.starts_with("NC") && _output_type == OutputType::kFixed) {
    cerr << "FingerprintWriter::Initialise:fixed fingerprints must start with FP '" << _tag << "'\n";
    return 0;
  }

  if (_output_type == OutputType::kFixed || _output_type == OutputType::kDescriptor) {
    if (_nbits <= 0) {
      cerr << "FingerprintWriter::Initialise:request fixed or descriptor output, but nbits zero\n";
      return 0;
    }
  }

  if (_nbits > 0 && _output_type == OutputType::kSparse) {
    cerr << "FingerprintWriter::Initialise:sparse fingerprint output incompabitle with nbits " << _nbits << '\n';
    return 0;
  }

  return 1;
}

// A new tag has been specified, update the output type to reflect the
// kind of tag specified.
int
FingerprintWriter::ChangeOutputTypeIfNeeded(const IWString& new_tag) {
  if (new_tag.starts_with("FP")) {
    _output_type = OutputType::kFixed;
    if (_nbits == 0) {
      _nbits = 2048;
    }
  } else if (new_tag.starts_with("NC")) {
    _output_type = OutputType::kSparse;
  }

  return 1;
}

int
FingerprintWriter::SetSparseOutput(const char* tag) {
  _output_type = OutputType::kSparse;
  _tag = tag;
  if (! _tag.ends_with('<')) {
    _tag << '<';
  }
  
  if (! OkTag(_tag)) {
    cerr << "FingerprintWriter:SetSparseOutput:invalid tag '" << tag << "'\n";
    return 0;
  }

  return 1;
}

// NOte that this does not check for multiple initialisation
// nor does it do any checking for any other state already set.
int
FingerprintWriter::SetArrayOutput(int nfeatures) {
  _output_type = OutputType::kDescriptor;
  _nbits = nfeatures;
  _output_vector = new int[_nbits];
  _iwdigits.initialise(100);
  _iwdigits.set_leading_string(_output_column_separator);

  return 1;
}

int
FingerprintWriter::WriteHeaderIfNeeded(IWString_and_File_Descriptor& output) const {
  if (_output_type != OutputType::kDescriptor) {
    return 0;
  }

  output << "ID";
  for (int i = 0; i < _nbits; ++i) {
    output << _output_column_separator << _feature_prefix << i;
  }
  output << '\n';

  return 1;
}

int
FingerprintWriter::WriteFingerprint(const IWString& mname,
                        const Sparse_Fingerprint_Creator& sfc,
                        IWString_and_File_Descriptor& output) {
  switch (_output_type) {
    case OutputType::kSparse:
      return WriteSparseFingerprint(sfc, output);
    case OutputType::kFixed:
      return WriteFixedFingerprint(sfc, output);
    case OutputType::kDescriptor:
      return WriteDescriptors(mname, sfc, output);
    case OutputType::kSvml:
      return WriteSvml(sfc, output);
    default:
      cerr << "FingerprintWriter::WriteFixedFingerprint:what to write?\n";
      return 0;
  }
}

int
FingerprintWriter::WriteDescriptors(const IWString& mname,
                  const Sparse_Fingerprint_Creator& sfc,
                  IWString_and_File_Descriptor& output) {
  std::fill_n(_output_vector, _nbits, 0);

  for (const auto& [bit, count] : sfc.bits_found()) {
    _output_vector[bit % _nbits] += count;
  }

  append_first_token_of_name(mname, output);
  for (int i = 0; i < _nbits; ++i) {
    _iwdigits.append_number(output, _output_vector[i]);
  }

  output << '\n';

  return 1;
}

int
FingerprintWriter::WriteSparseFingerprint(
                const Sparse_Fingerprint_Creator& sfc,
                IWString_and_File_Descriptor& output) {

  _fp.resize_keep_storage(0);
  sfc.daylight_ascii_form_with_counts_encoded(_tag, _fp);
  output << _fp << '\n';

  return 1;
}

int
FingerprintWriter::WriteFixedFingerprint(
                const Sparse_Fingerprint_Creator& sfc,
                IWString_and_File_Descriptor& output) {
  _fixed.clear();

  for (const auto& [bit, _] : sfc.bits_found()) {
    _fixed.set_bit(bit % _nbits);
  }
  output << _tag << _fixed.DaylightAsciiRepresentationIncludingNsetInfo() << ">\n";

  return 1;
}

int
FingerprintWriter::WriteSvml(const Sparse_Fingerprint_Creator& sfc, IWString_and_File_Descriptor& output) {
  const auto n = sfc.nbits();

  std::unique_ptr<std::pair<uint32_t, int>[]> for_sort = std::make_unique<std::pair<uint32_t, int>[]>(n);

  int ndx = 0;
  for (const auto& [bit, count] : sfc.bits_found()) {
    for_sort[ndx].first = bit;
    for_sort[ndx].second = count;
    ++ndx;
  }

  std::sort(for_sort.get(), for_sort.get() + n, [](const std::pair<uint32_t, int>& fp1,
                                                   const std::pair<uint32_t, int>& fp2) {
    return fp1.first < fp2.first;
  });

  // By convention the first column in an svml file is the response.
  output << ".";
  for (uint32_t i = 0; i < n; ++i) {
    output << _output_column_separator << for_sort[i].first;
    _iwdigits.append_number(output, for_sort[i].second);
  }

  output << '\n';

  return 1;
}

}  // namespace fingerprint_writer
