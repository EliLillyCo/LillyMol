// Generate a composite fingerprint for a collection of molecules.

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_map.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"

namespace gfp_collection_fingerprint {

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
  cerr << R"(Computes a composite fingerprint for a collection of fingerprints.
By default all fixed and sparse fingerprints are processed. In addition, with the -p option,
pairs of fixed and/or sparse fingerprints can also be included.
 -p fixed        profile pairs of fixed  fingerprints.
 -p sparse       profile pairs of sparse fingerprints.
 -u <support>    specify support level for a feature or feature combination to be included
 -T <title>      write <title> to the PCN<> tag of the output
 -v              verbose output
)";

  // clang-format on

  ::exit(rc);
}

class CollectionData {
  private:
    enum ScalingType {
      kLinear = 0,
      kLog10 = 1
    };

  private:
    int _verbose;

    ScalingType _scaling_type;

    // For keeping track of the bits encountered.

    // Fixed width fingerprints
    std::unique_ptr<std::vector<uint32_t>[]> _fixed;
    std::unique_ptr<std::vector<uint32_t>[]> _pairwise_fixed;

    // Sparse fingerprints
    std::unique_ptr<absl::flat_hash_map<uint32_t, uint32_t>[]> _sparse;
    std::unique_ptr<absl::flat_hash_map<uint32_t, uint32_t>[]> _pairwise_sparse;

    int _first_call;
    // We cache these from number_fingerprints() and number_sparse_fingerprints()
    int _nfixed;
    int _nsparse;

    uint32_t _support;

    int _collect_pairwise_fixed_bits;
    int _collect_pairwise_sparse_bits;

    Sparse_Fingerprint_Creator _sfc;

    IWString _fixed_tag, _fixed_pairwise_tag;
    IWString _sparse_tag, _sparse_pairwise_tag;

  // private functions
    int AllocateFingerprintArrays(IW_General_Fingerprint& fp);
    int ProfileFixed(const IWDYFP& fp, std::vector<uint32_t>& collect);
    int ProfilePairwiseFixed(const IWDYFP& fp, std::vector<uint32_t>& collect);
    int ProfileSparse(const Sparse_Fingerprint& fp,
                       absl::flat_hash_map<uint32_t, uint32_t>& acc);
    int ProfilePairwiseSparse(const Sparse_Fingerprint& fp,
                       absl::flat_hash_map<uint32_t, uint32_t>& acc);
    int WriteFixed(std::vector<uint32_t>& acc, const IWString& tag, IWString_and_File_Descriptor& output);
    int WriteSparse(const absl::flat_hash_map<uint32_t, uint32_t>& acc, const IWString& tag, IWString_and_File_Descriptor& output);

    int LinearScalingFixed(std::vector<uint32_t>& data,
                              const IWString& tag,
                              IWString_and_File_Descriptor& output);
    int LogarithmicScalingFixed(std::vector<uint32_t>& data,
                        const IWString& tag,
                        IWString_and_File_Descriptor& output);
    int LinearScalingSparse(const absl::flat_hash_map<uint32_t, uint32_t>& data,
                              const IWString& tag,
                        IWString_and_File_Descriptor& output);
    int Log10ScalingSparse(const absl::flat_hash_map<uint32_t, uint32_t>& data,
                              const IWString& tag,
                        IWString_and_File_Descriptor& output);

  public:
    CollectionData();
    ~CollectionData();

    int Initialise(Command_Line& cl);

    int Profile(IW_General_Fingerprint& fp);

    int Write(IWString_and_File_Descriptor& output);
};

CollectionData::CollectionData() {
  _verbose = 0;
  _scaling_type = kLinear;
  _first_call = 1;
  _nfixed = 0;
  _nsparse = 0;
  _collect_pairwise_fixed_bits = 0;
  _collect_pairwise_sparse_bits = 0;

  _support = 1;

  _fixed_tag = "NCFIX";
  _sparse_tag = "NCSPA";
}

CollectionData::~CollectionData() {
}

int
CollectionData::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('p')) {
    const_IWSubstring p;
    for (int i = 0; cl.value('p', p, i); ++i) {
      if (p == "fixed") {
        _collect_pairwise_fixed_bits = 1;
      } else if (p == "sparse") {
        _collect_pairwise_sparse_bits = 1;
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      if (s == "linear") {
        _scaling_type = kLinear;
      } else if (s == "log10") {
        _scaling_type = kLog10;
      } else {
        cerr << "Unrecognised -s qualifier '" << s << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('u')) {
    if (! cl.value('u', _support)) {
      cerr << "Invalid support value (-u)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will ignore bits and bit combinations occurring less than " << _support << " times\n";
    }
  }

  return 1;
}

int
CollectionData::AllocateFingerprintArrays(IW_General_Fingerprint& fp) {
  _nfixed = number_fingerprints();
  _nsparse = number_sparse_fingerprints();
  if (_nfixed == 0 && _nsparse == 0) {
    cerr << "CollectionData::AllocateFingerprintArrays:no fingerprints\n";
    return 0;
  }

  // Note that we significantly over-allocate the arrays in the case of
  // pairwise fixed fingerprints, but it would be much more complex to
  // make it smaller.
  if (_nfixed) {
    _fixed = std::make_unique<std::vector<uint32_t>[]>(_nfixed);
    if (_collect_pairwise_fixed_bits) {
      _pairwise_fixed = std::make_unique<std::vector<uint32_t>[]>(_nfixed);
    }
    for (int i = 0; i < _nfixed; ++i) {
      const int nbits = fp[i].nbits();
      _fixed[i].resize(nbits);
      if (_collect_pairwise_fixed_bits) {
        _pairwise_fixed[i].resize(nbits * nbits);
      }
    }
  }

  if (_nsparse) {
    _sparse = std::make_unique<absl::flat_hash_map<uint32_t, uint32_t>[]>(_nsparse);
    _pairwise_sparse = std::make_unique<absl::flat_hash_map<uint32_t, uint32_t>[]>(_nsparse);
  }

  _first_call = 0;

  return 1;
}

int
CollectionData::Profile(IW_General_Fingerprint& fp) {
  if (_first_call) {
    if (! AllocateFingerprintArrays(fp)) {
      return 0;
    }
    _first_call = 0;
  }

  for (int i = 0; i < _nfixed; ++i) {
    ProfileFixed(fp[i], _fixed[i]);
    if (_collect_pairwise_fixed_bits) {
      ProfilePairwiseFixed(fp[i], _pairwise_fixed[i]);
    }
  }

  for (int i = 0; i < _nsparse; ++i) {
    ProfileSparse(fp.sparse_fingerprint(i), _sparse[i]);
    if (_collect_pairwise_sparse_bits) {
      ProfilePairwiseSparse(fp.sparse_fingerprint(i), _pairwise_sparse[i]);
    }
  }

  return 1;
}

int
CollectionData::ProfileFixed(const IWDYFP& fp, std::vector<uint32_t>& collect) {
  cerr << "Fingerprint with " << fp.nbits() << " bits being profiled\n";
  int ndx = 0;
  while (true) {
    int bit = fp.next_on_bit(ndx);
    // cerr << " Bit " << bit << '\n';
    if (bit < 0) {
      break;
    }
    ++collect[bit];
  }

  return 1;
}

int
CollectionData::ProfilePairwiseFixed(const IWDYFP& fp, std::vector<uint32_t>& collect) {
  int ndx1 = 0;
  while (true) {
    int bit1 = fp.next_on_bit(ndx1);
    if (bit1 < 0) {
      break;
    }

    ++collect[bit1];

    int ndx2 = ndx1 + 1;
    while (true) {
      int bit2 = fp.next_on_bit(ndx2);
      if (bit2 < 0) {
        break;
      }

      uint32_t ndx = fp.nbits() + bit1 * fp.nbits() + bit2;
      //cerr << "Index of " << bit1 << " and " << bit2 << " is " << ndx << '\n';
      ++collect[ndx];
    }
  }

  return 1;
}
int
CollectionData::ProfileSparse(const Sparse_Fingerprint& fp,
                       absl::flat_hash_map<uint32_t, uint32_t>& acc) {
  uint32_t bit;
  int count;
  int i = 0;
  while (fp.next_bit_set(i, bit, count)) {
    ++acc[bit];
  }

  return 1;
}

int
CollectionData::ProfilePairwiseSparse(const Sparse_Fingerprint& fp,
                       absl::flat_hash_map<uint32_t, uint32_t>& acc) {
  const uint32_t nbits = fp.nbits();

  std::unique_ptr<uint32_t[]> bits = std::make_unique<uint32_t[]>(nbits);

  uint32_t bit;
  int count;
  int i = 0;
  for (int ndx = 0; fp.next_bit_set(i, bit, count); ++ndx) {
    bits[ndx] = bit;
    ++acc[bit];
  }

  std::sort(bits.get(), bits.get() + nbits, [](uint32_t b1, uint32_t b2) {
      return b1 < b2;
    }
  );

  for (uint32_t i = 0; i < nbits; ++i) {
    for (uint32_t j = i + 1; j < nbits; ++j) {
      uint32_t bit = nbits * bits[i]  + bits[j];
      ++acc[bit];
    }
  }

  return 1;
}

IWString
MakeTag(const IWString& stem, const char* extra) {
  IWString result;
  result << "NC" << stem << extra << '<';
  return result;
}

int
CollectionData::Write(IWString_and_File_Descriptor& output) {
  for (int i = 0; i < _nfixed; ++i) {
    WriteFixed(_fixed[i], MakeTag(fixed_fingerprint_tag(i), ""), output);
    output.write_if_buffer_holds_more_than(32768);
    if (_collect_pairwise_fixed_bits) {
      WriteFixed(_pairwise_fixed[i], MakeTag(fixed_fingerprint_tag(i), "P"), output);
      output.write_if_buffer_holds_more_than(32768);
    }
  }
  for (int i = 0; i < _nsparse; ++i) {
    WriteSparse(_sparse[i], MakeTag(sparse_fingerprint_tag(i), ""), output);
    output.write_if_buffer_holds_more_than(32768);
    if (_collect_pairwise_sparse_bits) {
      WriteSparse(_pairwise_sparse[i], MakeTag(sparse_fingerprint_tag(i), "p"), output);
      output.write_if_buffer_holds_more_than(32768);
    }
  }

  return 1;
}

int
CollectionData::WriteFixed(std::vector<uint32_t>& data,
                const IWString& tag,
                IWString_and_File_Descriptor& output) {
  switch (_scaling_type) {
    case kLinear:
      return LinearScalingFixed(data, tag, output);
    case kLog10:
      return LogarithmicScalingFixed(data, tag, output);
    default:
      cerr << "CollectionData::WriteFixed:unknown scaling type\n";
      return 0;
  }
}

int
CollectionData::LinearScalingFixed(std::vector<uint32_t>& data,
                              const IWString& tag,
                              IWString_and_File_Descriptor& output) {
  const uint32_t maxval = *std::max_element(data.cbegin(), data.cend());

  Sparse_Fingerprint_Creator sfc;

  const uint32_t n = data.size();
  for (uint32_t i = 0; i < n; ++i) {
    uint32_t count = data[i];
    if (count < _support) {
      continue;
    }
    float c = static_cast<float>(count) / static_cast<float>(maxval) * 255.0f;
    sfc.hit_bit(i, 1 + static_cast<int>(c + 0.4999f));
  }

  return sfc.write_fingerprint(tag, output);
}

int
CollectionData::LogarithmicScalingFixed(std::vector<uint32_t>& data,
                        const IWString& tag,
                        IWString_and_File_Descriptor& output) {
  const uint32_t maxval = *std::max_element(data.cbegin(), data.cend());
  if (maxval == 0) {
    return 0;
  }

  const float fmax = std::log10(static_cast<float>(maxval));

  const uint32_t n = data.size();

  Sparse_Fingerprint_Creator sfc;

  for (uint32_t i = 0; i < n; ++i) {
    uint32_t count = data[i];
    if (count < _support) {
      continue;
    }
    // Add 1 so we do not get zero count values
    float fcount = log10(count + 1);

    float c = fcount / fmax * 256.0f;
    sfc.hit_bit(i, 1 + static_cast<int>(c + 0.4999f));
  }

  return sfc.write_fingerprint(tag, output);
}

int
CollectionData::WriteSparse(const absl::flat_hash_map<uint32_t, uint32_t>& data,
                        const IWString& tag,
                        IWString_and_File_Descriptor& output) {
  switch(_scaling_type) {
    case kLinear:
      return LinearScalingSparse(data, tag, output);
    case kLog10:
      return Log10ScalingSparse(data, tag, output);
    default:
      cerr << "CollectionData::WriteSparse:unknown scaling type\n";
      return 0;
  }

  return 1;  // never happens
}

// Return the max value in `data`.
uint32_t
MaxValue(const absl::flat_hash_map<uint32_t, uint32_t>& data) {
  auto iter = std::max_element(data.begin(), data.end(), [](const auto &x, const auto &y) {
      return x.second < y.second;
    });

  return iter->second;
}

int
CollectionData::LinearScalingSparse(const absl::flat_hash_map<uint32_t, uint32_t>& data,
                        const IWString& tag,
                        IWString_and_File_Descriptor& output) {
  const uint32_t maxval = MaxValue(data);

  Sparse_Fingerprint_Creator sfc;

  for (const auto& [bit, count] : data) {
    if (count < _support) {
      continue;
    }

    float c = static_cast<float>(count) / static_cast<float>(maxval) * 255.0;
    sfc.hit_bit(bit, 1 + static_cast<int>(c + 0.4999f));
  }

  return sfc.write_fingerprint(tag, output);
}

int
CollectionData::Log10ScalingSparse(const absl::flat_hash_map<uint32_t, uint32_t>& data,
                        const IWString& tag,
                        IWString_and_File_Descriptor& output) {
  const uint32_t maxval = MaxValue(data);
  float fmaxval = std::log10(static_cast<float>(maxval));

  Sparse_Fingerprint_Creator sfc;

  for (const auto& [bit, count] : data) {
    if (count < _support) {
      continue;
    }

    float c = std::log10(static_cast<float>(count)) / fmaxval * 255.0f;
    sfc.hit_bit(bit, 1 + static_cast<int>(c + 0.4999f));
  }


  return sfc.write_fingerprint(tag, output);
}

int
ProfileCollection(IW_General_Fingerprint& fp,
                  CollectionData& data) {
  return data.Profile(fp);
}

int
ProfileCollection(iwstring_data_source& input, CollectionData& data) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    IW_General_Fingerprint fp;
    int fatal;
    if (! fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        return 0;
      }
      return 1;
    }

    if (! ProfileCollection(fp, data)) {
      return 0;
    }
  }

  return 1;
}

int
ProfileCollection(const char* fname, CollectionData& data) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ProfileCollection::cannot open '" << fname << "'\n";
    return 0;
  }

  return ProfileCollection(input, data);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:p:s:u:T:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  CollectionData data;

  if (! data.Initialise(cl)) {
    cerr << "Cannot initialise data collection options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const char* fname : cl) {
    if (! ProfileCollection(fname, data)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('T')) {
    IWString title = cl.string_value('T');
    output << "PCN<" << title << ">\n";
  }

  data.Write(output);
  output << "|\n";

  return 0;
}

}  // namespace gfp_collection_fingerprint

int
main(int argc, char **argv) {
  int rc = gfp_collection_fingerprint::Main(argc, argv);

  return rc;
}
