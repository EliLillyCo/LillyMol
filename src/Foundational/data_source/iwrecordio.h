#ifndef FOUNDATIONAL_DATA_SOURCE_IWRECORDIO_H_
#define FOUNDATIONAL_DATA_SOURCE_IWRECORDIO_H_

#include <memory>
#include <optional>
#include <string_view>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"

namespace iwrecordio {

// These have the same public methods as the corresponding TFdata classes.
// Could make them inherit from a common base class, but seems un-necessary.

class IWRecordIoReader {
  private:
    // The file descriptor from which data is retrieved.
    int _fd;
    // Status flags.

    bool _good;
    bool _eof;

    int _items_read;

    // An index into _read_buffer where the next item starts.
    uint64_t _next;

    // Data is read from _fd into _read_buffer.
    // Note that the resizable_array class is not 64 bit capable.
    resizable_array<char> _read_buffer;

    // private functions.

    void DefaultValues();

    bool OpenFile(const char * fname);

    std::optional<uint64_t> GetLength();

    bool FillReadBuffer(uint64_t bytes_needed = 8192);

  public:
    IWRecordIoReader();
    IWRecordIoReader(const char * fname);
    IWRecordIoReader(IWString & fname);
    IWRecordIoReader(const const_IWSubstring & fname);
    ~IWRecordIoReader();

    int Open(const char * fname);
    int Open(IWString & fname);
    int Open(const const_IWSubstring & fname);

    bool IsOpen() const { return _fd >= 0;}

    bool good() const { return _good;}
    bool eof() const { return _eof;}

    int seek_zero();

    int Close();

    // Will not reduce the size of _read_buffer, but may
    // increase it.
    int SetBufferSize(int buf_size);

    int items_read() const { return _items_read;}

    // The primary method for this class. If possible, return
    // the next item.
    std::optional<const_IWSubstring> Next(); 

    // Read serialized proto of type P and return decoded form.
    template <typename P>
    std::optional<P>
    ReadProto();

    // Read serialized proto type P and return a unique ptr to
    // the newly read proto.
    template <typename T>
    std::unique_ptr<T> ReadProtoPtr();
};

class IWRecordIoWriter {
  private:
    // Handy abstraction that already has a file descriptor and
    // methods for writing.
    IWString_and_File_Descriptor _output;

  // private functions

  ssize_t WriteData(const void * data, uint64_t nbytes);
  ssize_t WriteLength(const uint64_t nbytes);
  int CommonWrite(const void * data, const uint64_t nbytes);

  public:
    IWRecordIoWriter();
    IWRecordIoWriter(int fd);

    int Open(const char * fname);
    int Open(IWString& fname);
    int Open(const const_IWSubstring& fname);

    int Close();

    ssize_t Write(const void * data, uint64_t nbytes);
    ssize_t Write(const const_IWSubstring& data);
    ssize_t Write(const IWString& data);

    // Write a serialized proto.
    template <typename P>
    int WriteSerializedProto(const P& proto);
};

template <typename P>
std::optional<P>
IWRecordIoReader::ReadProto() {
  std::optional<const_IWSubstring> data = Next();
  if (! data) {
    std::cerr << "No data returned from Next\n";
    return std::nullopt;
  }

  const std::string_view as_string(data->data(), data->length());
  P proto;
  if (! proto.ParseFromString(as_string)) {
    std::cerr << "IWRecordIoReader::ReadProto:cannot parse serialized form\n";
    return std::nullopt;
  }

  return proto;
}

template <typename T>
std::unique_ptr<T>
IWRecordIoReader::ReadProtoPtr() {
  std::optional<const_IWSubstring> data = Next();
  if (! data) {
    return nullptr;
  }
  const std::string_view as_string(data->data(), data->length());
  std::unique_ptr<T> result = std::make_unique<T>();
  if (! result->ParseFromString(as_string)) {
    std::cerr << "IWRecordIoReader::ReadProto:cannot parse serialized form\n";
    return nullptr;
  }

  return result;
}

template <typename P>
int
IWRecordIoWriter::WriteSerializedProto(const P& proto) {
  const std::string as_string = proto.SerializeAsString();
  return Write(as_string.data(), as_string.size());
}

}  // namespace iwrecordio

#endif // FOUNDATIONAL_DATA_SOURCE_IWRECORDIO_H_
