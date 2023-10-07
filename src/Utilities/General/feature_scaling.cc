// Scale and un-scale a numeric feature in a descriptor file.
// First phase is to use the -C option to profile an existing data
// file to establish the range. That is written to a proto.
// Then there are two possible modes of operation.
//  1. Apply that same scaling to another file - compress to [0,1].
//  2. Unscale some data that has been previously generated here,
//     convert from [0,1] to [min,max]

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Utilities/General/feature_scaling.pb.h"
#define FEATURE_SCALER_IMPLEMENTATION
#include "Utilities/General/scaler.h"

namespace feature_scaling {

using std::cerr;

int verbose = 0;

int header_records = 1;

int activity_column = 1;

char input_separator = ' ';

char output_separator = ' ';

bool truncate_out_of_range = true;

// We can either scale to the target range, or unscale scaled values.
enum class Action {kScaleToRange, kScaleBackToOriginal };

// If the -subset option is given.
IW_STL_Hash_Set identifiers_to_process;

Fraction_as_String fraction_as_string;

void
Usage(int rc) {
  cerr << "Scale an unscale a numeric data column\n";
  cerr << " -C <fname>        profile the input, create a scaling and apply\n";
  cerr << " -bin              create both .txt and .dat proto files with scaling information\n";
  cerr << " -U <fname>        use a previously generated scaling to scale or unscale data\n";
  cerr << " -action ...       either 'scale' or 'unscale' data using previously generated scaling data\n";
  cerr << " -c <col>          colum contining the data\n";
  cerr << " -subset <fname>   only process the identifiers in <fname>\n";
  cerr << " -scol <col>       the identifiers in the -subset file are in <col> (use 2 for a smiles file)\n";
  cerr << " -11               scale and unscale to the [-1,1] range (default is [0,1]\n";
  cerr << " -v                verbose output\n";
  exit(rc);
}

template <typename T>
int
DoWrite(T value,
        IWString_and_File_Descriptor& output) {
  if (fraction_as_string.active()) {
    fraction_as_string.append_number(output, value);
  } else {
    output << output_separator << value;
  }
  return 1;
}

int
UseScalingLine(const const_IWSubstring& buffer,
               const feature_scaler::FeatureScaler<float>& scaling,
               const Action action,
               IWString_and_File_Descriptor& output) {
  int i = 0;
  const_IWSubstring token;
  for(int col = 0; buffer.nextword(token, i, input_separator); ++col) {
    if (col > 0) {
      output << output_separator;
    }
    if (col != activity_column) {
      output << token;
      continue;
    }

    double value;
    if (! token.numeric_value(value)) {
      cerr << "UseScalingLine:invalid numeric '" << token << "\n";
      return 0;
    }


    if (action == Action::kScaleToRange) {
      value = scaling.Scale(value);
    } else {
      value = scaling.ScaleBackToOrignalRange(value);
    }
    DoWrite(value, output);
  }

  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
UseScaling(iwstring_data_source& input,
           const feature_scaler::FeatureScaler<float>& scaling,
           Action action,
           IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  for (int i = 0; i < header_records; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Cannot head header\n";
      return 0;
    }
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  while (input.next_record(buffer)) {
    if (! UseScalingLine(buffer, scaling, action, output)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
UseScaling(const char * fname,
           const feature_scaler::FeatureScaler<float>& scaling,
           Action action,
           IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return UseScaling(input, scaling, action, output);
}

int
UseScaling(const Command_Line_v2& cl,
           const feature_scaler::FeatureScaler<float>& scaling,
           Action action,
           IWString_and_File_Descriptor& output) {
  for (const char * fname : cl) {
    if (! UseScaling(fname, scaling, action, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
UseScaling(const Command_Line_v2& cl,
           IWString& proto_file_name,
           Action action,
           IWString_and_File_Descriptor& output) {
  iwstring_data_source input(proto_file_name.null_terminated_chars());
  if (! input.good()) {
    cerr << "Cannot open proto file '" << proto_file_name << "'\n";
    return 0;
  }

  google::protobuf::io::FileInputStream file_input_stream(input.fd());

  FeatureScaling::FeatureScaling scaling;

  if (! google::protobuf::TextFormat::Parse(&file_input_stream, &scaling)) {
    cerr << "UseScaling:cannot parse '" << proto_file_name << "'\n";
    return 0;
  }

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(scaling);
  scaler.set_truncate_out_of_range(truncate_out_of_range);

  if (cl.option_present("11")) {
    scaler.Initialise11Scaling();
  } else {
    scaler.Initialise01Scaling();
  }

  return UseScaling(cl, scaler, action, output);
}

// Return true if `buffer` should be processed for determining
// the range.
// If `identifiers_to_process` is empty, no subsetting is done
// so return ok.
// Otherwise return ok if the identifier in `buffer` is in
// identifiers_to_process.
int
OkSubset(const const_IWSubstring& buffer,
         const IW_STL_Hash_Set& identifiers_to_process) {

  static constexpr int kIdentifierColumn = 0;
  static constexpr char kSep = ' ';

  // No subsetting, everything is OK.
  if (identifiers_to_process.empty()) {
    return 1;
  }

  const_IWSubstring token;
  int i = 0;
  for (int col = 0; buffer.nextword(token, i, kSep); ++col) {
    if (col != kIdentifierColumn) {
      continue;
    }

    IWString tmp(token);
    return identifiers_to_process.contains(tmp);
  }

  // Kind of strange, did not find the identifier column, 
  // silently OK????
  return 1;
}

int
GatherRangeLine(const const_IWSubstring& buffer,
                 FeatureScaling::FeatureScaling& scaling) {

  if (! OkSubset(buffer, identifiers_to_process)) {
    return 1;
  }

  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i, input_separator); ++col) {
    if (col != activity_column) {
      continue;
    }
    double value;
    if (! token.numeric_value(value)) {
      cerr << "GatherRangeLine:invalid numeric '" << token << "'\n";
      return 0;
    }
    if (value < scaling.min()) {
      scaling.set_min(value);
    }
    if (value > scaling.max()) {
      scaling.set_max(value);
    }
    return 1;
  }

  cerr << "Did not find column " << activity_column << "  in '" << buffer << "'\n";
  return 0;
}

int
GatherRange(iwstring_data_source& input,
            FeatureScaling::FeatureScaling& scaling) {
  const_IWSubstring buffer;
  for (int i = 0; i < header_records; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Cannot read header record\n";
      return 0;
    }
  }

  scaling.set_min(std::numeric_limits<double>::max());
  scaling.set_max(-std::numeric_limits<double>::max());

  while (input.next_record(buffer)) {
    if (! GatherRangeLine(buffer, scaling)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
GatherRange(const char * fname,
            FeatureScaling::FeatureScaling& scaling) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GatherRange(input, scaling);
}

int
GetIdentiferSubsetLine(const_IWSubstring& line,
                   const int scol,
                   IW_STL_Hash_Set& identifiers_to_process) {
  static constexpr char kSep = ' ';
  IWString token;
  int i = 0;
  for (int col = 0; line.nextword(token, i, kSep); ++col) {
    if (col != scol) {
      continue;
    }

    identifiers_to_process.emplace(token);
    return 1;
  }

  return 0;
}

int
GetIdentiferSubset(iwstring_data_source& input,
                   int scol,
                   IW_STL_Hash_Set& identifiers_to_process) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! GetIdentiferSubsetLine(buffer, scol, identifiers_to_process)) {
      cerr << "GetIdentiferSubset:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
GetIdentiferSubset(IWString& fname,
                   int scol,
                   IW_STL_Hash_Set& identifiers_to_process) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GetIdentiferSubset:cannot open '" << fname << "'\n";
    return 0;
  }

  return GetIdentiferSubset(input, scol, identifiers_to_process);
}

std::optional<FeatureScaling::FeatureScaling>
GatherRange(const Command_Line_v2& cl,
            IWString& output_fname) {

  FeatureScaling::FeatureScaling scaling;
  for (const char * fname : cl) {
    if (! GatherRange(fname, scaling)) {
      cerr << "Fatal error profiling '" << fname << "'\n";
      return std::nullopt;
    }
  }

  if (! cl.option_present("bin")) {
    if (! iwmisc::WriteProtoAsText(scaling, output_fname)) {
      return std::nullopt;
    }
  }

  IWString fname;
  fname << output_fname << ".txt";
  if (! iwmisc::WriteProtoAsText(scaling, fname)) {
    return std::nullopt;
  }

  fname = output_fname << ".dat";
  if (! iwmisc::WriteBinaryProto(scaling, fname)) {
    return std::nullopt;
  }

  return scaling;
}

int
FeatureScaling(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-C=s-U=s-c=ipos-action=s-hdr=ipos-prec=ipos-bin-subset-11=sfile-scol=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    cl.value('c', activity_column);
    if (verbose) {
      cerr << "Activity in column " << activity_column << '\n';
    }
    activity_column--;
  }

  if (cl.option_present("hdr")) {
    cl.value("hdr", header_records);
    if (verbose)
      cerr << "Inpus has " << header_records << " header records\n";
  }

  if (cl.option_present('C') && cl.option_present("action")) {
    cerr << "The create proto (-C) and -action options are incompatible\n";
    Usage(1);
  }

  if (cl.option_present('C') && cl.option_present('U')) {
    cerr << "Cannot specify both create (-C) and use (-U) options\n";
    Usage(1);
  }

  Action action = Action::kScaleToRange;
  if (cl.option_present("action")) {
    const IWString s = cl.string_value("action");
    if (s == "scale") {
      action = Action::kScaleToRange;
    } else if (s == "unscale" ) {
      action = Action::kScaleBackToOriginal;
    } else {
      cerr << "Unrecognised action '" << s << "'\n";
      Usage(1);
    }
  }

  // Right now, this is just set up for the [0,1] range, but it
  // could be changed to be [scaling.min(),scaling.max()].
  if (cl.option_present("prec")) {
    int precision;
    cl.value("prec", precision);
    fraction_as_string.initialise(0.0, 1.0, precision);
  }

  if (cl.option_present("subset")) {
    IWString fname = cl.string_value("subset");
    int scol = 0;
    if (cl.option_present("scol")) {
      cl.value("scol", scol);
      if (scol < 1) {
        cerr << "Invalid subset column specifier " << scol << "\n";
        return 1;
      }
      if (verbose) {
        cerr << "Subset identifiers in column " << scol << " in '" << fname << "'\n";
      }
      --scol;
    }
    if (! GetIdentiferSubset(fname, scol, identifiers_to_process)) {
      cerr << "Cannot read identifier subset '" << fname << "'\n";
      return 1;
    }
    if (verbose) {
      cerr << "Read " << identifiers_to_process.size() << " identifiers to process from '" << fname << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);
  int rc = 0;
  if (cl.option_present('C')) {
    IWString cfile = cl.string_value('C');
    std::optional<FeatureScaling::FeatureScaling> scaling = GatherRange(cl, cfile);
    if (! scaling) {
      cerr << "Cannot create scaling file '" << cfile << "'\n";
      return 1;
    }

    feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(*scaling);

    scaler.set_truncate_out_of_range(truncate_out_of_range);

    if (cl.option_present("11")) {
      scaler.Initialise11Scaling();
    } else {
      scaler.Initialise01Scaling();
    }

    rc = UseScaling(cl, scaler, Action::kScaleToRange, output);
  } else if (cl.option_present('U')) {
    IWString ufile = cl.string_value('U');
    rc = UseScaling(cl, ufile, action, output);
  } else {
    cerr << "Must specify one of -C (create) or -U (use) options\n";
    Usage(1);
  }

  if (! rc) {
    return 1;
  }

  return 0;
}

}  // namespace feature_scaling

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return feature_scaling::FeatureScaling(argc, argv);
}
