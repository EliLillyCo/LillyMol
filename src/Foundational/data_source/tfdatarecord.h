#ifndef FOUNDATIONAL_DATA_SOURCE_TFDATARECORD_H
#define FOUNDATIONAL_DATA_SOURCE_TFDATARECORD_H
// Read and write TFDataRecord

#include <optional>

#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

namespace iw_tf_data_record {

class TFDataReader {
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
    TFDataReader();
    TFDataReader(const char * fname);
    TFDataReader(IWString & fname);
    TFDataReader(const const_IWSubstring & fname);
    ~TFDataReader();

    int Open(const char * fname);
    int Open(IWString & fname);
    int Open(const const_IWSubstring & fname);

    bool IsOpen() const { return _fd >= 0;}

    bool good() const { return _good;}
    bool eof() const { return _eof;}

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
};

// Writes data files that can subsequently be read by TFDataReader.
class TFDataWriter {
  private:
    // Handy abstraction that already has a file descriptor and
    // methods for writing.
    IWString_and_File_Descriptor _output;

  // private functions

  ssize_t WriteData(const void * data, uint64_t nbytes);
  ssize_t WriteLength(const uint64_t nbytes);
  int CommonWrite(const void * data, const uint64_t nbytes);

  public:
    TFDataWriter();
    TFDataWriter(int fd);

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
TFDataReader::ReadProto() {
  std::optional<const_IWSubstring> data = Next();
  if (! data) {
    return std::nullopt;
  }
  const std::string as_string(data->data(), data->length());
  P proto;
  if (! proto.ParseFromString(as_string)) {
    std::cerr << "TFDataReader::ReadProto:cannot parse serialized form\n";
    return std::nullopt;
  }

  return proto;
}

template <typename P>
int
TFDataWriter::WriteSerializedProto(const P& proto) {
  const std::string as_string = proto.SerializeAsString();
  return Write(as_string.data(), as_string.size());
}

}  // namespace iw_tf_data_record

#endif // FOUNDATIONAL_DATA_SOURCE_TFDATARECORD_H
