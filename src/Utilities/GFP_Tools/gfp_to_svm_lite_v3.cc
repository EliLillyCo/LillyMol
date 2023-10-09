// Converts GFP files to svm lite.output_sparse_fp_tag, tmp);
// Primary function is to assign a unique feature number to each
// bit in the gfp file(s).
// Output can be either to an svm_lite input file format, or to
// a .gfp file with a single, sparse, fingerprint.

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#define ACTIVITY_DATA_IMPLEMENATION_H
#include "Foundational/iwmisc/activity_data_from_file.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Utilities/GFP_Tools/bit_subset.h"
#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"

namespace gfp_to_svm_lite
{

using std::cerr;

// The type of features that are used as feature numbers in the svm_lite input.
using feature_type_t = uint32_t;
// The type of bit numbers that arise from gfp files.
using gfp_bit_type_t = uint32_t;

int verbose = 0;

IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";
IWString activity_tag = "ACT<";

// The activity can be part of the name.
int activity_column = -1;

IWDigits bit_number;
IWDigits count;

// Output can be svm_lite form, or a gfp subset.
// by default, only svm_lite produced.
bool output_as_svm_lite = true;
bool output_as_gfp_subset = false;
bool output_as_tsv = false;

int flatten_sparse_fingerprints = 0;

char sep = ' ';

// If we are using an existing proto map file, it can only
// be processed after we have read the first gfp item.
bool use_bit_map = false;
IWString input_proto_fname;

int fingerprints_read = 0;

int write_counts_as_float = 0;

// If the resulting svml file is destined for LightGBM, then
// the identifier is not allowed.
int append_identifier_at_end_of_each_record = 1;

void
Usage(int rc)
{
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Convert a gfp file to form(s) useful for svm_lite\n";
  cerr << " -A <fname>        activity data in <fname>\n";
  cerr << " -r                no activity data, write a random value\n";
  cerr << " -c <col>          activity is column <col> in the name in the input .gfp file\n";
  cerr << " -C <fname>        create GfpBitXref and GfpBitSubset protos and use for the current input\n";
  cerr << " -X <fname>        use a previously generated GfpBitXref proto for the current input\n";
  cerr << " -U <fname>        use a previously generated GfpBitSubset proto for the current input\n";
  cerr << " -O <type>         output type, 'svml, gfp'\n";
  cerr << " -S <stem>         output file name stem, <stem>.gfp and/or <stem>.svml\n";
  cerr << " -p <n>            support level, discard bits found in <n> or fewer fingerprints\n";
  cerr << " -f                flatten counted bits to 1\n";
  cerr << " -d                write counts as floating point values - ordinal rather than categorical values\n";
  cerr << " -V                output is for LIBSVM, suppress addition of identifier\n";
  cerr << " -v                verbose output\n";
  // clang-format on

  exit(rc);
}

struct FeatureCount {
  feature_type_t feature;
  uint32_t count;
};

// A class to facilitate bit subsetting and cross referencing.
// During the formation phase, the ProfileBits* methods are
// called and ProfileBits() is called for eahc fingerprint.
// This populates the _dense_count and _sparse_count hashes.
//
// Once the input fingerprints are read,feature subset and feature
// cross reference protos are written.
//
// Then one or more output files can be created.
//   *  Subset of the GFP
//   *  svm-lite output
//   *  a single, sparse, gfp fingerprint.

// During output, the subset of bits to be written is then either
// converted from the _dense_count and _sparse_count hashes, which
// are re-used, or read from protos.

class GfpBitsRetained
{
 private:
  // A mapping from bit number in each fingerprint, to the number of molecules
  // that have that bit set. For each fingerprint.
  // By having a count, we can enable the support functionality.
  resizable_array<std::unordered_map<gfp_bit_type_t, uint32_t>> _dense_count;
  resizable_array<std::unordered_map<gfp_bit_type_t, uint32_t>> _sparse_count;

  // The names of the fingerprints - used when writing the proto.
  resizable_array<IWString> _dense_gfp_name;
  resizable_array<IWString> _sparse_gfp_name;

  // If a support requirement is imposed, we store it and write with
  // the proto.
  int _support;

  // private functions

  GfpBitSubset::BitSubset
  MakeProto(const std::unordered_map<gfp_bit_type_t, uint32_t>& count) const;

  // Remove all features that occur in fewer than `_support` molecules.
  // Returns the number of features discarded.
  int
  ImposeSupport();

  // Both of the protos produced need to have their `params` attributes filled.
  template <typename Proto>
  void
  SetParams(Proto& proto, uint32_t highest_feature_number) const;

 public:
  GfpBitsRetained();

  // Uses the globally available counts of dense and sparse fingerprints
  // to initialize.
  int
  Initialize();

  // Building from protos.
  int
  FromProto(IWString& fname);

  int
  DebugPrint(std::ostream& output) const;

  void
  set_support(int s)
  {
    _support = s;
  }

  // Iterate through fingerprints, collecting the bits set into
  // _dense_count and _sparse_count.
  int
  ProfileBits(const IWDYFP& fp, int fpnum);
  int
  ProfileBits(const Sparse_Fingerprint& fp, int fpnum);

  // In the case of having called ProfileBits to identify the bits
  // present in the inputs, convert the _dense_count and _sparse_count
  // hashes to cross references from gfp bit number to globally
  // unique feature number. During this phase any support requirement
  // is imposed.
  // Returns the maximum feature number.
  int
  ConvertToCrossReference();

  // During the output phase, convert the retained bit numbers in `fp` to
  // feature/count pairs in `features`, starting at position `ndx`, which is
  // incremented.
  int
  BitsToFeatures(const IWDYFP& fp, int fpnum, FeatureCount* features, int& ndx);
  int
  BitsToFeatures(const Sparse_Fingerprint& fp, int fpnum, FeatureCount* features,
                 int& ndx);

  // Convert to proto form.
  GfpBitSubset::GfpBitSubset
  ToBitSubsetProto() const;
  GfpBitSubset::GfpBitToFeature
  ToBitXrefProto() const;

  // Write the bit cross reference and bit subset data in proto form.
  int
  WriteBitSubsetProto(IWString& fname) const;
  int
  writeBitXrefProto(IWString& fname) const;
};

GfpBitsRetained::GfpBitsRetained()
{
  _support = 0;
}

std::tuple<int, int>
GetNumberFingerprints()
{
  return {number_fingerprints(), number_sparse_fingerprints()};
}

int
GfpBitsRetained::Initialize()
{
  const auto [nfixed, nsparse] = GetNumberFingerprints();

  _dense_count.extend(nfixed, std::unordered_map<gfp_bit_type_t, uint32_t>());
  _sparse_count.extend(nsparse, std::unordered_map<gfp_bit_type_t, uint32_t>());
  _dense_gfp_name.extend(nfixed, IWString());
  _sparse_gfp_name.extend(nsparse, IWString());

  for (int i = 0; i < nfixed; ++i) {
    _dense_gfp_name[i] = fixed_fingerprint_tag(i);
  }
  for (int i = 0; i < nsparse; ++i) {
    _sparse_gfp_name[i] = sparse_fingerprint_tag(i);
  }

  if (verbose) {
    cerr << "Initialised for " << nfixed << " fixed and " << nsparse
         << " sparse fingerprints\n";
  }

  return 1;
}

int
GfpBitsRetained::DebugPrint(std::ostream& output) const
{
  output << "GfpBitsRetained: with " << _dense_count.number_elements() << " fixed and "
         << _sparse_count.number_elements() << " fingerprints\n";
  return output.good();
}

// Given a list of bits to retain in `proto.bits`, populate `retain`.
void
ExtractBitsToRetain(const GfpBitSubset::BitSubset& proto,
                    std::unordered_map<gfp_bit_type_t, uint32_t>& retain)
{
  for (uint32_t bit : proto.bits()) {
    retain.emplace(bit, 1);  // 1 is an arbitrary choice, could be any number.
  }
}

int
GfpBitsRetained::FromProto(IWString& fname)
{
  std::optional<GfpBitSubset::GfpBitSubset> proto =
      iwmisc::ReadBinaryProto<GfpBitSubset::GfpBitSubset>(fname);
  if (!proto) {
    cerr << "GfpBitSubset::FromProto:cannot read '" << fname << "'\n";
    return 0;
  }

  Initialize();

  const auto [nfixed, nsparse] = GetNumberFingerprints();

  for (int i = 0; i < nfixed; ++i) {
    const std::string tag(_dense_gfp_name[i].data(), _dense_gfp_name[i].length());
    auto iter = proto->xref().find(tag);
    if (iter == proto->xref().end()) {
      cerr << "GfpBitsRetained::FromProto:no data for " << _dense_gfp_name[i] << "'\n";
      return 0;
    }
    ExtractBitsToRetain(iter->second, _dense_count[i]);
  }

  for (int i = 0; i < nsparse; ++i) {
    const std::string tag(_sparse_gfp_name[i].data(), _sparse_gfp_name[i].length());
    auto iter = proto->xref().find(tag);
    if (iter == proto->xref().end()) {
      cerr << "GfpBitsRetained::FromProto:no data for " << _sparse_gfp_name[i] << "'\n";
      return 0;
    }
    ExtractBitsToRetain(iter->second, _sparse_count[i]);
  }

  return 1;
}

void
UpdateHash(std::unordered_map<uint32_t, uint32_t>& bit_count, uint32_t bit)
{
  auto iter = bit_count.find(bit);
  if (iter == bit_count.end()) {
    bit_count.emplace(bit, 1);
  } else {
    iter->second++;
  }
}

// Until C++ 20
// https://en.cppreference.com/w/cpp/container/unordered_map/erase_if
//  rc += std::erase_if(mapping, [support] (const auto& kv) {
//    const auto [bit, count] = kv;
//    return count < support;
//  })
int
RemoveLowSupport(std::unordered_map<gfp_bit_type_t, uint32_t>& mapping, uint32_t support)
{
  int rc = 0;
  for (auto iter = mapping.begin(), last = mapping.end(); iter != last;) {
    if (iter->second < support) {
      mapping.erase(iter++);
      rc++;
    } else {
      ++iter;
    }
  }

  return rc;
}

int
GfpBitsRetained::ImposeSupport()
{
  int rc = 0;
  for (auto& dense : _dense_count) {
    rc += RemoveLowSupport(dense, _support);
  }

  for (auto& sparse : _sparse_count) {
    rc += RemoveLowSupport(sparse, _support);
  }

  return rc;
}

// Convert an unordered_map to a BitSubset proto.
GfpBitSubset::BitSubset
MakeSubsetProto(const std::unordered_map<gfp_bit_type_t, uint32_t>& count)
{
  GfpBitSubset::BitSubset result;
  for (const auto& [key, value] : count) {
    result.add_bits(key);
  }
  return result;
}

#ifdef NOT_NEEDED
void
UpdateHighestFeature(GfpBitSubset::BitSubset& subset, uint32_t& highest_feature_number)
{
  for (const auto feature : subset.bits()) {
    if (feature > highest_feature_number) {
      highest_feature_number = feature;
    }
  }
}
#endif

GfpBitSubset::GfpBitSubset
GfpBitsRetained::ToBitSubsetProto() const
{
  GfpBitSubset::GfpBitSubset result;

  const auto [nfixed, nsparse] = GetNumberFingerprints();

  uint32_t highest_feature_number = 0;
  for (int i = 0; i < nfixed; ++i) {
    GfpBitSubset::BitSubset counts = MakeSubsetProto(_dense_count[i]);
    //  UpdateHighestFeature(counts, highest_feature_number);

    const std::string tag_name(_dense_gfp_name[i].data(), _dense_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitSubset::BitSubset> to_insert(tag_name,
                                                                              counts);

    result.mutable_xref()->insert(to_insert);
  }

  for (int i = 0; i < nsparse; ++i) {
    GfpBitSubset::BitSubset counts = MakeSubsetProto(_sparse_count[i]);
    //  UpdateHighestFeature(counts, highest_feature_number);

    const std::string tag_name(_sparse_gfp_name[i].data(), _sparse_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitSubset::BitSubset> to_insert(tag_name,
                                                                              counts);

    result.mutable_xref()->insert(to_insert);
  }

  SetParams(result, highest_feature_number);

  return result;
}

GfpBitSubset::BitXref
MakeXrefProto(const std::unordered_map<gfp_bit_type_t, uint32_t>& count)
{
  GfpBitSubset::BitXref result;

  for (const auto& [key, value] : count) {
    google::protobuf::MapPair<uint32_t, uint32_t> to_insert(key, value);
    result.mutable_bit_to_feature()->insert(to_insert);
  }

  return result;
}

template <typename Proto>
void
GfpBitsRetained::SetParams(Proto& proto, uint32_t highest_feature_number) const
{
  if (_support > 0) {
    proto.mutable_params()->set_support(_support);
  }

  if (flatten_sparse_fingerprints) {
    proto.mutable_params()->set_flatten_sparse_fingerprints(true);
  }

  proto.mutable_params()->set_highest_feature_number(highest_feature_number);
}

void
UpdateHighestFeature(const GfpBitSubset::BitXref& xref, uint32_t& highest_feature_number)
{
  for (const auto& [_, feature] : xref.bit_to_feature()) {
    if (feature > highest_feature_number) {
      highest_feature_number = feature;
    }
  }
}

GfpBitSubset::GfpBitToFeature
GfpBitsRetained::ToBitXrefProto() const
{
  GfpBitSubset::GfpBitToFeature result;

  const auto [nfixed, nsparse] = GetNumberFingerprints();

  uint32_t highest_feature_number = 0;

  for (int i = 0; i < nfixed; ++i) {
    GfpBitSubset::BitXref xref = MakeXrefProto(_dense_count[i]);
    UpdateHighestFeature(xref, highest_feature_number);

    const std::string tag_name(_dense_gfp_name[i].data(), _dense_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitSubset::BitXref> to_insert(tag_name,
                                                                            xref);

    result.mutable_xref()->insert(to_insert);
  }

  for (int i = 0; i < nsparse; ++i) {
    GfpBitSubset::BitXref xref = MakeXrefProto(_sparse_count[i]);
    UpdateHighestFeature(xref, highest_feature_number);

    const std::string tag_name(_sparse_gfp_name[i].data(), _sparse_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitSubset::BitXref> to_insert(tag_name,
                                                                            xref);

    result.mutable_xref()->insert(to_insert);
  }

  SetParams(result, highest_feature_number);

  return result;
}

int
GfpBitsRetained::WriteBitSubsetProto(IWString& fname) const
{
  const GfpBitSubset::GfpBitSubset proto = ToBitSubsetProto();
  return iwmisc::WriteBinaryProto(proto, fname);
}

int
GfpBitsRetained::writeBitXrefProto(IWString& fname) const
{
  const GfpBitSubset::GfpBitToFeature proto = ToBitXrefProto();
  return iwmisc::WriteBinaryProto(proto, fname);
}

int
GfpBitsRetained::ProfileBits(const IWDYFP& fp, int fpnum)
{
  int ndx = 0;
  int bit;
  while ((bit = fp.next_on_bit(ndx)) >= 0) {
    UpdateHash(_dense_count[fpnum], bit);
  }

  return 1;
}

int
GfpBitsRetained::ProfileBits(const Sparse_Fingerprint& fp, int fpnum)
{
  int i = 0;
  gfp_bit_type_t bit = 0;
  int c = 0;
  while (fp.next_bit_set(i, bit, c)) {
    UpdateHash(_sparse_count[fpnum], bit);
  }
  return 1;
}

void
ToCrossReference(std::unordered_map<gfp_bit_type_t, uint32_t>& xref, int& feature_number)
{
  for (auto iter = xref.begin(); iter != xref.end(); ++iter) {
    iter->second = feature_number;
    feature_number++;
  }
}

int
GfpBitsRetained::ConvertToCrossReference()
{
  if (_support > 0) {
    const int removed = ImposeSupport();
    if (verbose) {
      cerr << "GfpBitsRetained::ImposeSupport:support " << _support << " removed "
           << removed << " bits\n";
    }
  }

  // Globally unique feature number.
  int feature_number = 1;

  const auto [nfixed, nsparse] = GetNumberFingerprints();
  for (int i = 0; i < nfixed; ++i) {
    ToCrossReference(_dense_count[i], feature_number);
  }
  for (int i = 0; i < nsparse; ++i) {
    ToCrossReference(_sparse_count[i], feature_number);
  }

  return feature_number;
}

int
GfpBitsRetained::BitsToFeatures(const IWDYFP& fp, int fpnum, FeatureCount* features,
                                int& ndx)
{
  int rc = 0;
  int bit = 0;
  while (fp.next_on_bit(bit)) {
    const auto iter = _dense_count[fpnum].find(bit);
    if (iter == _dense_count[fpnum].end()) {
      continue;
    }
    features[ndx].feature = iter->second;
    features[ndx].count = 1;
    ndx++;
    rc++;
  }

  return rc;
}

int
GfpBitsRetained::BitsToFeatures(const Sparse_Fingerprint& fp, int fpnum,
                                FeatureCount* features, int& ndx)
{
  int rc = 0;
  int i = 0;
  uint32_t bit = 0;
  int count = 0;
  while (fp.next_bit_set(i, bit, count)) {
    const auto iter = _sparse_count[fpnum].find(bit);
    if (iter == _sparse_count[fpnum].cend()) {
      continue;
    }
    features[ndx].feature = iter->second;
    if (flatten_sparse_fingerprints) {
      features[ndx].count = 1;
    } else {
      features[ndx].count = count;
    }
    ndx++;
    rc++;
  }

  return rc;
}

// Write the `nfeatures` feature/count data in `features` in svm_lite
// form to `output`.
#ifdef NOP_LONGER_USED
template <typename Activity>
int
WriteSvmLite(FeatureCount* features, int nfeatures, const IWString& name,
             Activity activity, IWString_and_File_Descriptor& output)
{
  std::sort(features, features + nfeatures,
            [](const FeatureCount& fc1, const FeatureCount& fc2) {
              return fc1.feature < fc2.feature;
            });

  output << activity;
  for (int i = 0; i < nfeatures; ++i) {
    bit_number.append_number(output, features[i].feature);
    count.append_number(output, features[i].count);
  }
  if (append_identifier_at_end_of_each_record) {
    output << " # " << name;
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}
#endif

// In order to reduce the number of arguments passed, put some
// arguments to GfpToSvmLite that are commonly passed.
struct Args {
  // The objects that enable conversion to the output forms.
  // EIther a subset or a mapping from gfp bits to svm lite feature numbers.
  bit_subset::BitSubset bit_subset;
  bit_subset::BitXref bit_xref;

  // The output stream for svml output.
  IWString_and_File_Descriptor svml_output;
  // The output stream for gfp output.
  IWString_and_File_Descriptor gfp_output;
  // Catboost needs a tab separated descriptor file.
  IWString_and_File_Descriptor tsv_output;

  // Catboost needs a feature description file.
  IWString_and_File_Descriptor feature_description_file;

  // If tsv output is needed, compute once the number of
  // columns in the output.
  uint32_t highest_feature_number;

  // Used for tsv formation.
  char output_separator = '\t';

  int write_random_experimental_values = 0;
};

template <typename Activity>
int
SvmLiteOutput(IW_General_Fingerprint& gfp, const IWString& id, Activity activity,
              Args& args)
{
  args.svml_output << activity;
  args.bit_xref.WriteSvmlFeatures(gfp, args.svml_output);
  if (append_identifier_at_end_of_each_record) {
    args.svml_output << " # " << id;
  }
  args.svml_output << '\n';

  args.svml_output.write_if_buffer_holds_more_than(8192);

  return 1;
}

template <typename Activity>
int
TsvOutput(IW_General_Fingerprint& gfp, const IWString& id, Activity activity, Args& args)
{
  args.tsv_output << activity << args.output_separator;
  args.bit_xref.WriteDsv(gfp, args.output_separator, args.tsv_output);
  args.tsv_output << '\n';

  args.tsv_output.write_if_buffer_holds_more_than(8192);
  return 1;
}

template <typename Activity>
int
GfpSubsetOutput(IW_General_Fingerprint& gfp, const IWString& id, Activity activity,
                Args& args)
{
  args.bit_subset.MakeSubset(gfp);

  cerr << "Gfp subset output not implemented\n";

  return 1;
}

// Initialise the two bit_subset objects in `args`.
int
InitialiseGfpKnown(IW_General_Fingerprint& gfp, Args& args)
{
  if (args.bit_xref.Active()) {
    if (!args.bit_xref.InitialiseGfpKnown(gfp)) {
      return 0;
    }
  }
  if (args.bit_subset.Active()) {
    if (!args.bit_subset.InitialiseGfpKnown(gfp)) {
      return 0;
    }
  }

  return 1;
}

int
NFeatures(const IW_General_Fingerprint& gfp)
{
  const auto [nfixed, nsparse] = GetNumberFingerprints();
  int nfeatures = 0;

  for (int i = 0; i < nfixed; ++i) {
    nfeatures += gfp[i].nset();
  }
  for (int i = 0; i < nsparse; ++i) {
    nfeatures += gfp.sparse_fingerprint(i).nset();
  }

  return nfeatures;
}

template <typename Activity>
int
GfpToSvmLite(IW_General_Fingerprint& gfp, const IWString& id, Activity activity,
             Args& args)
{
  static bool first_call = true;
  if (first_call) {
    first_call = false;
    if (!InitialiseGfpKnown(gfp, args)) {
      return 0;
    }
  }

  int rc = 1;
  if (output_as_svm_lite) {
    if (!SvmLiteOutput(gfp, id, activity, args)) {
      rc = 0;
    }
  }
  if (output_as_gfp_subset) {
    if (!GfpSubsetOutput(gfp, id, activity, args)) {
      rc = 0;
    }
  }
  if (output_as_tsv) {
    if (!TsvOutput(gfp, id, activity, args)) {
      rc = 0;
    }
  }

  return rc;
}

int
ProfileBits(const IW_General_Fingerprint& gfp, GfpBitsRetained& bit_xref)
{
  const auto [nfixed, nsparse] = GetNumberFingerprints();

  for (int i = 0; i < nfixed; ++i) {
    bit_xref.ProfileBits(gfp[i], i);
  }
  for (int i = 0; i < nsparse; ++i) {
    bit_xref.ProfileBits(gfp.sparse_fingerprint(i), i);
  }

  return 1;
}

int
EchoItem(const IW_TDT& tdt, const IWString& tag, IWString_and_File_Descriptor& output)
{
  IWString s;
  if (!tdt.dataitem_value(tag, s, 0)) {
    cerr << "Cannot extract " << tag << " from " << tdt << '\n';
    return 0;
  }
  output << tag << s << ">\n";
  return 1;
}

int
EchoSmilesID(const IW_TDT& tdt, IWString_and_File_Descriptor& output)
{
  if (!EchoItem(tdt, smiles_tag, output) || !EchoItem(tdt, identifier_tag, output)) {
    cerr << "Cannot echo smiles/id\n";
    return 0;
  }

  return 1;
}

template <typename Activity>
std::optional<Activity>
ActivityFromColumn(const IWString& id, int activity_column)
{
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; id.nextword(token, i); ++col) {
    if (col != activity_column) {
      continue;
    }
    Activity result;
    if (!token.numeric_value(result)) {
      cerr << "ActivityFromColumn:invalid numeric " << id << "'\n";
      return std::nullopt;
    }
    return result;
  }
  cerr << "ActivityFromColumn:did not encounter column " << activity_column << '\n';
  return std::nullopt;
}

// Get the activity for `id` from either `activity` or from
// column `activity_column` of `id`;
template <typename Activity>
std::optional<Activity>
GetActivity(const Activity_Data_From_File<Activity>& activity, const IWString& id,
            int activity_column)
{
  if (activity_column < 0 && activity.empty()) {
    return 0;
  }

  if (activity_column > 0) {
    return ActivityFromColumn<Activity>(id, activity_column);
  }

  const auto iter = activity.find(id);
  if (iter != activity.end()) {
    return iter->second;
  }

  if (id.nwords() == 1) {
    cerr << "No activity for '" << id << "'\n";
    return std::nullopt;
  }

  IWString tmp(id);
  tmp.truncate_at_first(' ');
  return GetActivity(activity, tmp, activity_column);
}

template <typename Activity>
int
GfpToSvmLite(IW_TDT& tdt, const Activity_Data_From_File<Activity>& activity, Args& args)
{
  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    if (fatal) {
      return 0;
    }
    return 1;
  }

  // Or we could/should get the id from `gfp`.
  IWString id;
  if (!tdt.dataitem_value(identifier_tag, id, 0) || id.empty()) {
    cerr << "Cannot extract" << identifier_tag << " from " << tdt << '\n';
    return 0;
  }

  if (args.write_random_experimental_values) {
    static std::random_device rd;
    std::uniform_real_distribution<float> u(0.0, 1.0);
    Activity arbitrary = u(rd);
    return GfpToSvmLite(gfp, id, arbitrary, args);
  }

  std::optional<Activity> act = GetActivity(activity, id, activity_column);
  if (!act) {
    return 0;
  }

  return GfpToSvmLite(gfp, id, *act, args);
}

int
ProfileBits(IW_TDT& tdt, GfpBitsRetained& bit_xref)
{
  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    if (fatal) {
      return 0;
    }
    return 1;
  }

  static bool first_call = true;
  if (first_call) {
    first_call = false;
    bit_xref.Initialize();
  }

  fingerprints_read++;

  return ProfileBits(gfp, bit_xref);
}

template <typename Activity>
int
GfpToSvmLite(iwstring_data_source& input,
             const Activity_Data_From_File<Activity>& activity, Args& args)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!GfpToSvmLite(tdt, activity, args)) {
      cerr << "Cannot process " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

int
ProfileBits(iwstring_data_source& input, GfpBitsRetained& bit_xref)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!ProfileBits(tdt, bit_xref)) {
      cerr << "Cannot process " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

template <typename Activity>
int
GfpToSvmLite(const char* fname, const Activity_Data_From_File<Activity>& activity,
             Args& args)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return GfpToSvmLite(input, activity, args);
}

int
ProfileBits(const char* fname, GfpBitsRetained& bit_xref)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return ProfileBits(input, bit_xref);
}

// Write both GfpBitSubset and GfpBitToFeature protos formed from
// `bit_xref` to files with stem `output_proto_stem`.
int
WriteProtos(GfpBitsRetained& bit_xref, const IWString& output_proto_stem)
{
  IWString fname;
  fname << output_proto_stem << "_subset.dat";
  bit_xref.WriteBitSubsetProto(fname);
  fname = output_proto_stem;
  fname << "_xref.dat";
  bit_xref.writeBitXrefProto(fname);

  return 1;
}

IWString
FileNameGivenSuffix(const IWString& output_file_name_stem, const char* suffix)
{
  if (output_file_name_stem == "-") {
    return IWString("-");
  }
  if (output_file_name_stem.ends_with(suffix)) {
    return IWString(output_file_name_stem);
  }

  IWString fname;
  fname << output_file_name_stem << suffix;
  return fname;
}

int
OpenWithSuffix(const IWString& output_file_name_stem, const char* suffix,
               IWString_and_File_Descriptor& output, int verbose)
{
  IWString fname = FileNameGivenSuffix(output_file_name_stem, suffix);
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "OpenWithSuffix:cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << suffix << " output to '" << fname << "'\n";
  }
  return 1;
}

int
GfpToSvmLite(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vF:P:C:X:U:A:p:fc:S:dlO:D:Vr");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise fingerprints\n";
      return 1;
    }
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  Activity_Data_From_File<float> activity;
  if (cl.option_present('c') && cl.option_present('A')) {
    cerr << "Can use just one of the -c or -A options\n";
    Usage(1);
  }

  Args args;

  if (cl.option_present('r')) {
    args.write_random_experimental_values = 1;
    set_default_iwstring_float_concatenation_precision(4);
  } else if (cl.option_present('c')) {
    if (!cl.value('c', activity_column) || activity_column < 1) {
      cerr << "Invalid activity column value (-c)\n";
      return 1;
    }
    if (verbose) {
      cerr << "Activity values extracted from column " << activity_column << '\n';
    }
    activity_column--;
  } else if (cl.option_present('A')) {
    if (!activity.construct_from_command_line(cl, 'A', verbose)) {
      cerr << "Cannot read activity data (-A)\n";
      return 1;
    }
  } else if (cl.option_present('X') || cl.option_present('U')) {
    // activity not available when doing a test set.
  } else {
    cerr << "MUst specify a source of activity via either -A or -c options\n";
    Usage(1);
  }

  if (cl.option_present('f')) {
    flatten_sparse_fingerprints = 1;
    if (verbose) {
      cerr << "Will flatten sparse fingerprints\n";
    }
  }

  const IWString output_file_name_stem = cl.string_value('S');

  if (cl.empty()) {
    cerr << "Must specify input file(s)\n";
    Usage(1);
  }

  // Specifying a support level only makes sense if creating a profile.
  if (cl.option_present('p') && !cl.option_present('C')) {
    cerr << "Specifying a support level only makes sense with the -C option\n";
    Usage(1);
  }

  if (cl.option_present('U') && cl.option_present('X')) {
    cerr << "The subset (-U) and cross reference (-X) options are mutually exclusive\n";
    Usage(1);
  }

  if (cl.option_present('O')) {
    output_as_gfp_subset = false;
    output_as_svm_lite = false;
    output_as_tsv = false;

    const_IWSubstring o;
    for (int i = 0; cl.value('O', o, i); ++i) {
      int j = 0;
      const_IWSubstring token;
      while (o.nextword(token, j, ',')) {
        if (token == "svml") {
          output_as_svm_lite = true;
        } else if (token == "gfp") {
          output_as_gfp_subset = true;
        } else if (token == "tsv") {
          output_as_tsv = true;
        } else {
          cerr << "Unrecognised -O directive '" << o << "'\n";
          return 1;
        }
      }
    }
  }

  // If we are creating the subset/cross reference here.
  GfpBitsRetained bit_xref;

  if (cl.option_present('p')) {
    int support;
    if (!cl.value('p', support) || support < 1) {
      cerr << "Invalid support value (-p)\n";
      Usage(1);
    }

    bit_xref.set_support(support);
    if (verbose) {
      cerr << "WIll only output features found in " << support << " molecules\n";
    }
  }

  // Fill in one or both of the protos in `args`.

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    for (const char* fname : cl) {
      cerr << "Profiling " << fname << '\n';
      if (!ProfileBits(fname, bit_xref)) {
        cerr << "Cannot profile bits in " << fname << '\n';
        return 1;
      }
    }
    bit_xref.ConvertToCrossReference();
    WriteProtos(bit_xref, fname);
    // Generate both protos, even if only one needed. Should be cheap.
    args.bit_subset.Build(bit_xref.ToBitSubsetProto());
    args.bit_xref.Build(bit_xref.ToBitXrefProto());
    // args.bit_xref.DebugPrint(cerr);
  } else if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    if (!args.bit_subset.Build(fname)) {
      cerr << "GfpToSvmLite::cannot initialise bit subset proto '" << fname << "'\n";
      return 1;
    }
  } else if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    if (!args.bit_xref.Build(fname)) {
      cerr << "GfpToSvmLite::cannot initialise bit xref proto '" << fname << "'\n";
      return 1;
    }
  } else {
    cerr << "Must specify whether creating a new (-C) or using (-U, -X) an existing "
            "proto\n";
    Usage(1);
  }

  // Note that there is no checking for multiple output streams trying to use stdout.
  if (output_as_svm_lite) {
    IWString fname =
        OpenWithSuffix(output_file_name_stem, ".svml", args.svml_output, verbose);
  }
  if (output_as_gfp_subset) {
    OpenWithSuffix(output_file_name_stem, ".gfp", args.gfp_output, verbose);
  }
  if (output_as_tsv) {
    OpenWithSuffix(output_file_name_stem, ".tsv", args.tsv_output, verbose);

    args.highest_feature_number = args.bit_xref.HighestFeatureNumber();

    if (cl.option_present('D')) {
      IWString fname = cl.string_value('D');
      if (!args.feature_description_file.open(fname.null_terminated_chars())) {
        cerr << "Cannot open feature description file '" << fname << "'\n";
        return 1;
      }
      if (verbose) {
        cerr << "Feature description written to '" << fname << "'\n";
      }
    }
  }

  if (cl.option_present('V')) {
    append_identifier_at_end_of_each_record = 0;
    if (verbose) {
      cerr << "Output intended for LightGBM - no identifier info\n";
    }
  }

  // Memoized output helpers for svml lite outputs.
  bit_number.set_include_leading_space(1);
  bit_number.initialise(10000);
  count.set_leading_string(':');
  count.initialise(256);
  if (cl.option_present('d')) {
    count.append_to_each_stored_string(".");
  }

  for (const char* fname : cl) {
    if (verbose) {
      cerr << "Begin processing " << fname << '\n';
    }
    if (!GfpToSvmLite(fname, activity, args)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }

  if (verbose) {
    cerr << "Processed " << fingerprints_read << " fingerprints\n";
  }

  return 0;
}

}  // namespace gfp_to_svm_lite

int
main(int argc, char** argv)
{
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return gfp_to_svm_lite::GfpToSvmLite(argc, argv);
}
