#ifndef FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
#define FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
// Functions to support operating with protos.

#include <fcntl.h>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>

#include "google/protobuf/text_format.h"
#include "google/protobuf/util/json_util.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iwstring.h"

namespace iwmisc {

using std::cerr;

// A lightweight class to open a file descrptor and make sure
// it gets closed.
class AFile {
  private:
    int _fd;

  public:
    AFile(IWString& fname, int mode);  // O_RDONLY or O_WRONLY
    ~AFile();

    int good() const {
      return _fd >= 0;
    }

    int fd() const {
      return _fd;
    }
};

// Write a proto to `fname` using Text_Format.
template <typename Proto>
int
WriteProtoAsText(const Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);
  if (! output.good()) {
    std::cerr << "WriteProtoAsText:cannot open " << fname << '\n';
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get())) {
    std::cerr << "WriteProtoAsText:cannot write " << fname << "\n";
    return 0;
  }

  return 1;
}

template <typename Proto>
int
WriteBinaryProto(const Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);
  if (! output.good()) {
    cerr << "WriteBinaryProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return proto.SerializeToFileDescriptor(output.fd());
}

template <typename Proto>
std::optional<Proto>
ReadBinaryProto(IWString& fname) {
  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return std::nullopt;
  }

  Proto result;
  if (! result.ParseFromFileDescriptor(input.fd())) {
    cerr << "ReadBinaryProto:cannot parse '" << fname << "'\n";
    return std::nullopt;
  }

  return result;
}

// Read existing proto.
template <typename Proto>
bool
ReadBinaryProto(IWString& fname, Proto& proto) {
  proto.Clear();

  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return false;
  }

  if (! proto.ParseFromFileDescriptor(input.fd())) {
    cerr << "ReadBinaryProto:cannot parse '" << fname << "'\n";
    return false;
  }

  return true;
}

template <typename Proto>
std::optional<Proto>
ReadTextProtoJson(const IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ReadTextProtoJson:cannot open '" << fname << "'\n";
    return std::nullopt;
  }
  IWString file_contents;
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    file_contents.append_with_spacer(buffer, ' ');
  }

  std::string_view tmp(file_contents.data(), file_contents.size());
  Proto result;
  google::protobuf::util::JsonParseOptions options;
  // This makes me very nervous, we should not be using an internal namespace.
  if (auto status = google::protobuf::util::JsonStringToMessage(tmp, &result, options);
      status != absl::OkStatus()) {
    return std::nullopt;
  }

  return result;
}

template <typename Proto>
std::optional<Proto>
ReadTextProto(IWString& fname) {
  if (fname.ends_with(".json")) {
    return ReadTextProtoJson<Proto>(fname);
  }

  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadTextProto:cannot open '" << fname << "'\n";
    return std::nullopt;
  }

  using google::protobuf::io::ZeroCopyInputStream;
  using google::protobuf::io::FileInputStream;
  std::unique_ptr<FileInputStream> zero_copy_input(new FileInputStream(input.fd()));

  Proto result;
  if (! google::protobuf::TextFormat::Parse(zero_copy_input.get(), &result)) {
    cerr << "ReadTextProto:cannot read '" << fname << "'\n";
    return std::nullopt;
  }

  return result;
}

template <typename Proto>
std::unique_ptr<Proto>
ReadTextProtoPtr(IWString& fname) {
  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadTextProto:cannot open '" << fname << "'\n";
    return nullptr;
  }

  using google::protobuf::io::ZeroCopyInputStream;
  using google::protobuf::io::FileInputStream;
  std::unique_ptr<FileInputStream> zero_copy_input(new FileInputStream(input.fd()));

  std::unique_ptr<Proto> result = std::make_unique<Proto>();
  if (! google::protobuf::TextFormat::Parse(zero_copy_input.get(), result.get())) {
    cerr << "ReadTextProto:cannot read '" << fname << "'\n";
    return nullptr;
  }

  return result;
}

// TextFormat does not allow comments. We can enable them.
// Note that it is possible to imagine a failure where a text
// value gets wrapped onto a separate line??? Probably not.
template <typename Proto>
std::optional<Proto>
ReadTextProtoCommentsOK(IWString& fname) {
  if (fname.ends_with(".json")) {
    return ReadTextProtoJson<Proto>(fname);
  }

  iwstring_data_source input(fname);
  if (! input.good()) {
    std::cerr <<"ReadTextProtoCommentsOK:cannot open '" << fname << "'\n";
    return std::nullopt;
  }

  input.set_ignore_pattern("^# ");

  std::string file_contents;
  const_IWSubstring line;
  while (input.next_record(line)) {
    file_contents.append(line.data(), line.length());
    file_contents.append(" ");
  }

  Proto result;
  if (! google::protobuf::TextFormat::ParseFromString(file_contents, &result)) {
    cerr << "ReadTextProtoCommentsOK:cannot read '" << fname << "'\n";
    return std::nullopt;
  }

  return result;
}


template <typename Proto>
int
WriteTextProto(Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);

  if (! output.good()) {
    cerr << "WriteTextProto:cannot open '" << fname << "'\n";
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get())) {
    cerr << "WriteTextProto:cannot write '" << fname << "'\n";
    return 0;
  }

  return 1;
}

}  // namespace iwmisc

#endif // FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
