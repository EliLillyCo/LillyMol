// Tests for TFDataRecord
#include <cstdio>

#include "Foundational/cmdline/cmdline.h"

#include "tfdatarecord.h"
#include "Foundational/data_source/proto_for_testing.pb.h"

namespace test_tfdata {

const char * prog_name = nullptr;

using std::cerr;

void
Usage(int rc) {
  exit(rc);
}

int
TestWriteRead(const std::string& fname) {
  iw_tf_data_record::TFDataWriter writer;
  if (! writer.Open(fname)) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }
  constexpr int nitmes = 10;
  for (int i = 0; i < nitmes; ++i) {
    JustForTesting::TestMessage message;
    message.set_i(-i);
    message.set_u(i);

    std::string serialized;
    message.SerializeToString(&serialized);
    writer.Write(serialized.data(), serialized.size());
  }
  writer.Close();
  cerr << "Wrote " << nitmes << " to " << fname << '\n';

  int rc = 1;
  int items_read = 0;
  iw_tf_data_record::TFDataReader reader(fname);
  for (; ! reader.eof(); ++items_read) {
    std::optional<const_IWSubstring> maybe_data = reader.Next();
    if (! maybe_data) {
      break;
    }
    JustForTesting::TestMessage message;
    std::string as_string(maybe_data->data(), maybe_data->length());
    if (! message.ParseFromString(as_string)) {
      cerr << "Invalid data for i=" << items_read << '\n';
      return 0;
    }
    if (message.i() != -items_read) {
      cerr << "Mismatch on signed value, got " << message.i() << " expected -" << items_read << '\n';
      rc = 0;
    }
    if (message.u() != static_cast<uint32_t>(items_read)) {
      cerr << "Mismatch on unsigned value, got " << message.u() << " expected " << items_read << '\n';
      rc = 0;
    }
  }
  cerr << "Read " << items_read << " items\n";
  return rc;
} 

int
TestTfDataRecord(const char * fname,
                 IWString_and_File_Descriptor& output) {
  iw_tf_data_record::TFDataReader reader(fname);
  const_IWSubstring data;
  int records_read = 0;

  for ( ; ! reader.eof() ; ++records_read) {
    std::optional<const_IWSubstring> maybe_data = reader.Next();
    if (! maybe_data) {
      break;
    }
//  cerr << "Read item with " << maybe_data->length() << " bytes\n";
  }
  cerr << "Read " << records_read << " records from " << fname << '\n';

  return 1;
}

int
TestTfDataRecord(int argc, char** argv) {
  Command_Line cl(argc, argv, "vw");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options\n";
    Usage(1);
  }

  if (cl.option_present('w')) {
    const std::string fname = std::tmpnam(nullptr);
    TestWriteRead(fname);
  }

  IWString_and_File_Descriptor output(1);
  for (const char * fname : cl) {
    if (! TestTfDataRecord(fname, output)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }
  return 0;
}
}  // namespace test_tfdata

int
main(int argc, char ** argv) {
  test_tfdata::prog_name = argv[0];

  int rc = test_tfdata::TestTfDataRecord(argc, argv);

  return rc;
}
