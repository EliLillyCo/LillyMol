#include <algorithm>
#include <cstdint>

#include <endian.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// Building crc32c with Bazel is difficult.
// For expediency, we can disable the crc based error checking.

#define LILLYMOL_HAS_CRC32C

#ifdef LILLYMOL_HAS_CRC32C
// if using crc32c directly     #include "crc32c/crc32c.h"
#include "absl/crc/crc32c.h"
#endif

#include "Foundational/iwmisc/misc.h"

#include "tfdatarecord.h"

namespace iw_tf_data_record {

using std::cerr;

unsigned int default_read_buffer_size = 4096; 

constexpr uint64_t sizeof_length = sizeof(uint64_t);
constexpr uint64_t sizeof_crc = sizeof(uint32_t);

#ifdef LILLYMOL_HAS_CRC32C
uint32_t
MaskedCrc(uint32_t crc) {
  return ((crc >> 15) | (crc << 17)) + 0xa282ead8ul;
}
#endif

void
TFDataReader::DefaultValues() {
  _fd = -1;
  _good = true;
  _eof = false;
  _next = 0;

  _read_buffer.resize(default_read_buffer_size);

//_compression_type = kUncompressed;

  _items_read = 0;
}

TFDataReader::TFDataReader() {
  DefaultValues();
}

TFDataReader::TFDataReader(const char * fname) {
  DefaultValues();
  OpenFile(fname);
}

TFDataReader::TFDataReader(IWString& fname) {
  DefaultValues();
  OpenFile(fname.null_terminated_chars());
}

TFDataReader::TFDataReader(const const_IWSubstring& fname) {
  DefaultValues();
  IWString tmp(fname);
  OpenFile(tmp.null_terminated_chars());
}

TFDataReader::~TFDataReader() {
}

int
TFDataReader::Open(const char * fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(IWString& fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(const const_IWSubstring& fname) {
  IWString tmp(fname);
  return OpenFile(tmp);
}

bool
TFDataReader::OpenFile(const char * fname) {
  if (_fd >= 0) {
    cerr << "TFDataReader::OpenFile:already open " << _fd << ", no action\n";
    return 0;
  }

  _fd = IW_FD_OPEN(fname, O_RDONLY);
#ifdef DEBUG_IW_TF_DATA
  cerr << "Openfile " << fname << " returned " << _fd << '\n';
#endif
  if (_fd < 0) {
    _good = false;
    return false;
  }

  _good = true;
  _eof = false;
  return true;
}

// _next is pointing at length data. Retrieve that length.
std::optional<uint64_t>
TFDataReader::GetLength() {
  const uint64_t* lptr = reinterpret_cast<const uint64_t*>(_read_buffer.rawdata() + _next);
  const uint32_t* crc = reinterpret_cast<const uint32_t*>(_read_buffer.rawdata() + _next + sizeof_length);

#ifdef DEBUG_IW_TF_DATA
  cerr << "Length " << *lptr  << " crc " << *crc << '\n';
#else
  (void) crc;  // keep the compiler quiet.
#endif
#ifdef LILLYMOL_HAS_CRC32C
  // If using crc32c directly
  // uint32_t result = crc32c::Crc32c(reinterpret_cast<const char*>(lptr), sizeof_length);
  // result = MaskedCrc(result);

  const absl::string_view tmp(reinterpret_cast<const char*>(lptr), sizeof_length);
  const absl::crc32c_t hash = absl::ComputeCrc32c(tmp);
  const uint32_t result = MaskedCrc(static_cast<uint32_t>(hash));
  if (result != *crc) {
    cerr << "TFDataReader::GetLength:crc fails, length " << *lptr << '\n';
    cerr << result << " vs " << *crc << '\n';
    _good = 0;
    return std::nullopt;
  }
#endif
  _next += sizeof_length + sizeof_crc;
  return *lptr;
}

std::optional<const_IWSubstring>
TFDataReader::Next() {
  if (_eof || ! _good) {
    cerr << "TFDataReader::Next:eof or not good\n";
    return std::nullopt;
  }

  // First task is to read the size of the next item.
  // At a minimum, we need 8+4 bytes for size of next and the crc
  if(_read_buffer.size() - _next < (sizeof_length + sizeof_crc)) {
    if (! FillReadBuffer()) {
      return std::nullopt;
    }
    if (_read_buffer.size() - _next < (sizeof_length + sizeof_crc)) {
      _eof = 1;
      return std::nullopt;
    }
  }

  std::optional<uint64_t> length = GetLength();
  if (! length) {
    return std::nullopt;
  }
#ifdef DEBUG_IW_TF_DATA
  cerr << "TFDataReader::Next: to read " << *length << " bytes\n";
#endif

  if (length == 0) {
    ++_items_read;
    return const_IWSubstring("");
  }

  if (_next + *length + sizeof_crc > static_cast<uint64_t>(_read_buffer.number_elements())) {
    if (! FillReadBuffer(*length + sizeof_crc)) {
      return std::nullopt;
    }
  }

  const_IWSubstring result(_read_buffer.rawdata() + _next, *length);

#ifdef LILLYMOL_HAS_CRC32C
  // If using crc32c directly
  // const uint32_t crc_data = crc32c::Crc32c(result.data(), *length);
  // const uint32_t masked_crc = MaskedCrc(crc_data);

  const absl::crc32c_t hash = absl::ComputeCrc32c(absl::string_view(result.data(), *length));
  const uint32_t masked_crc = MaskedCrc(static_cast<uint32_t>(hash));

  const uint32_t * crc = reinterpret_cast<const uint32_t*>(_read_buffer.rawdata() + _next + *length);
  if (*crc != masked_crc) {
    cerr << "TFDataReader::Next:Invalid data crc " << *length << " bytes\n";
    _good = 0;
    return std::nullopt;
  }
#endif

  _next += *length + sizeof_crc;

  ++_items_read;
  return result;
}

// The class needs to be able to read `bytes_needed` into an item.
// Upon exit, the next bytes_needed of data will be in the buffer.
bool
TFDataReader::FillReadBuffer(uint64_t bytes_needed) {
  // If the next bytes_needed bytes are already in _read_buffer we are done.
  if (_next + bytes_needed <= static_cast<uint64_t>(_read_buffer.number_elements())) {
    return true;
  }
  
  // We need to shift the data and maybe resize.
  if (_next > 0) {
    _read_buffer.remove_from_to(0, _next);
    _next = 0;
  }

  // At this stage, the next item will start at zero.

  if (static_cast<uint64_t>(_read_buffer.elements_allocated()) < bytes_needed) {
    _read_buffer.resize(bytes_needed);
  }
  const int bytes_already_present = _read_buffer.number_elements();
  const int to_read = _read_buffer.elements_allocated() - bytes_already_present;

  int bytes_read = IW_FD_READ(_fd, _read_buffer.rawdata() + bytes_already_present, to_read);
#ifdef DEBUG_IW_TF_DATA
  cerr << "Reading " << to_read << " got " << bytes_read << " bytes from file " << _fd << '\n';
#endif
  if (bytes_read < 0) {
    _good = false;
    return false;
  }
  if (bytes_read == 0) {
    _eof = true;
    return false;
  }

  _read_buffer.set_number_elements(bytes_already_present + bytes_read);
  return true;
}

TFDataWriter::TFDataWriter() {
}

TFDataWriter::TFDataWriter(int fd) : _output(fd) {
}

int
TFDataWriter::Open(IWString& fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname.null_terminated_chars())) {
    return 0;
  }

  return 1;
}

int
TFDataWriter::Open(const char* fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname)) {
    return 0;
  }

  return 1;
}

int
TFDataWriter::Open(const const_IWSubstring& fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }
  IWString tmp(fname);
  if (! _output.open(tmp.null_terminated_chars())) {
    return 0;
  }
  return 1;
}

int
TFDataWriter::Close() {
  return _output.close();
}

ssize_t
TFDataWriter::Write(const void * data, uint64_t nbytes) {
  if (! WriteLength(nbytes)) {
    cerr << "TFDataWriter::Write:cannot write " << nbytes << " length\n";
    return 0;
  }
  return WriteData(data, nbytes);
}

// Write `nbytes` from `value` as well as the crc of `data`.
int
TFDataWriter::CommonWrite(const void * data, const uint64_t nbytes) {
  const char * p = reinterpret_cast<const char *>(data);
  if (! _output.write(p, nbytes)) {
    cerr << "TFDataWriter::CommonWrite:cannot write " << nbytes << " bytes\n";
    return 0;
  }

#ifdef LILLYMOL_HAS_CRC32C
  // If using crc32c directly.
  // uint32_t crc = crc32c::Crc32c(p, nbytes);
  // crc = MaskedCrc(crc);

  const absl::crc32c_t hash = absl::ComputeCrc32c(absl::string_view(p, nbytes));
  const uint32_t crc = MaskedCrc(static_cast<uint32_t>(hash));
#else
  uint32_t crc = 0;
#endif
  p = reinterpret_cast<const char *>(&crc);
  if (! _output.write(p, sizeof_crc)) {
    cerr << "TFDataWriter::CommonWrite:cannot write crc, nbytes " << nbytes << '\n';
    return 0;
  }

  return _output.write_if_buffer_holds_more_than(8192);
}

ssize_t
TFDataWriter::WriteLength(const uint64_t nbytes) {
  const char * p = reinterpret_cast<const char *>(&nbytes);
  return CommonWrite(p, sizeof_length);
}

ssize_t
TFDataWriter::WriteData(const void * data, uint64_t nbytes) {
  if (nbytes == 0) {
    return 1;
  }
  return CommonWrite(data, nbytes);
}

}  // namespace iw_tf_data_record
