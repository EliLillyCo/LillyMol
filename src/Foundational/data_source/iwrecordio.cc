// Implementation of recordio
#include <algorithm>
#include <cstdint>

#ifdef __APPLE__
#include <machine/endian.h>
#else
#include <endian.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "Foundational/iwmisc/misc.h"

#include "iwrecordio.h"

namespace iwrecordio {

using std::cerr;

constexpr char kNewLine = '\n';

unsigned int default_read_buffer_size = 4096; 

constexpr uint64_t sizeof_crc = sizeof(uint32_t);

// The length of each item is stored as a plain text integer followed
// by a newline. We impose a maximum length of how much data we can
// read. Typically this will be used for serialized protos, so a much
// smaller number would be OK.
constexpr int64_t max_length_length = 7 + 1;

void
IWRecordIoReader::DefaultValues() {
  _fd = -1;
  _good = true;
  _eof = false;
  _next = 0;

  _read_buffer.resize(default_read_buffer_size);

//_compression_type = kUncompressed;

  _items_read = 0;
}

IWRecordIoReader::IWRecordIoReader() {
  DefaultValues();
}

IWRecordIoReader::IWRecordIoReader(const char * fname) {
  DefaultValues();
  OpenFile(fname);
}

IWRecordIoReader::IWRecordIoReader(IWString& fname) {
  DefaultValues();
  OpenFile(fname.null_terminated_chars());
}

IWRecordIoReader::IWRecordIoReader(const const_IWSubstring& fname) {
  DefaultValues();
  IWString tmp(fname);
  OpenFile(tmp.null_terminated_chars());
}

IWRecordIoReader::~IWRecordIoReader() {
}

int
IWRecordIoReader::Open(const char * fname) {
  return OpenFile(fname);
}

int
IWRecordIoReader::Open(IWString& fname) {
  return OpenFile(fname);
}

int
IWRecordIoReader::Open(const const_IWSubstring& fname) {
  IWString tmp(fname);
  return OpenFile(tmp);
}

bool
IWRecordIoReader::OpenFile(const char * fname) {
  if (_fd >= 0) {
    cerr << "IWRecordIoReader::OpenFile:already open " << _fd << ", no action\n";
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
IWRecordIoReader::GetLength() {
  const unsigned char* lptr = reinterpret_cast<const unsigned char*>(_read_buffer.rawdata() + _next);

  uint64_t len = 0;
  for (int i = 0; i < max_length_length; ++i) {
    char c = lptr[i];
    // cerr << "Character '" << c << "'\n";
    if (c == kNewLine) {
      _next += (i + 1);
      return len;
    }
    if (c >= '0' && c <= '9') {
      len = 10 * len + c - '0';
    } else {
      cerr << "IWRecordIoReader::GetLength:unrecognised character '" << c << "'\n";
      return std::nullopt;
    }
  }

  cerr << "IWRecordIoReader:GetLength:no newline found after " << max_length_length << " characters\n";
  return std::nullopt;
}

std::optional<const_IWSubstring>
IWRecordIoReader::Next() {
  if (_eof || ! _good) {
    cerr << "IWRecordIoReader::Next:eof or not good\n";
    return std::nullopt;
  }

  // First task is to read the size of the next item.
  // Make sure we have enough data to get the length.
  if(_read_buffer.size() - _next < max_length_length) {
    if (! FillReadBuffer()) {
      return std::nullopt;
    }
    if (_read_buffer.size() - _next < max_length_length) {
      _eof = 1;
      return std::nullopt;
    }
  }

  std::optional<uint64_t> length = GetLength();
  if (! length) {
    return std::nullopt;
  }
#ifdef DEBUG_IW_TF_DATA
  cerr << "IWRecordIoReader::Next: to read " << *length << " bytes\n";
#endif

  if (length == 0) {
    ++_items_read;
    return const_IWSubstring("");
  }
  // cerr << "_next " << _next << " size " << *length << " buffer " << _read_buffer.size()  << '\n';

  if (_next + *length > static_cast<uint64_t>(_read_buffer.number_elements())) {
    if (! FillReadBuffer(*length)) {
      return std::nullopt;
    }
  }

  const_IWSubstring result(_read_buffer.rawdata() + _next, *length);

  _next += *length;

  ++_items_read;
  return result;
}

// The class needs to be able to read `bytes_needed` into an item.
// Upon exit, the next bytes_needed of data will be in the buffer.
bool
IWRecordIoReader::FillReadBuffer(uint64_t bytes_needed) {
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

int
IWRecordIoReader::seek_zero() {
  if (_fd < 0) {
    return 0;
  }

  off_t rc = IW_FD_LSEEK(_fd, 0, SEEK_SET);

  if (rc < 0) {
    cerr << "IWRecordIoReader::seek_zero:cannot seek back to start of file\n";
    _good = false;
    return 0;
  }

  _good = true;
  _eof = false;
  _next = 0;
  _read_buffer.resize_keep_storage(0);
  return 1;
}


IWRecordIoWriter::IWRecordIoWriter() {
}

IWRecordIoWriter::IWRecordIoWriter(int fd) : _output(fd) {
}

int
IWRecordIoWriter::Open(IWString& fname) {
  if (_output.is_open()) {
    std::cerr << "IWRecordIoWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname.null_terminated_chars())) {
    return 0;
  }

  return 1;
}

int
IWRecordIoWriter::Open(const char* fname) {
  if (_output.is_open()) {
    std::cerr << "IWRecordIoWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname)) {
    return 0;
  }

  return 1;
}

int
IWRecordIoWriter::Open(const const_IWSubstring& fname) {
  if (_output.is_open()) {
    std::cerr << "IWRecordIoWriter::Open:already open\n";
    return 0;
  }
  IWString tmp(fname);
  if (! _output.open(tmp.null_terminated_chars())) {
    return 0;
  }
  return 1;
}

int
IWRecordIoWriter::Close() {
  return _output.close();
}

ssize_t
IWRecordIoWriter::Write(const void * data, uint64_t nbytes) {
  if (! WriteLength(nbytes)) {
    cerr << "IWRecordIoWriter::Write:cannot write " << nbytes << " length\n";
    return 0;
  }
  return WriteData(data, nbytes);
}

// Write `nbytes` from `value` as well as the crc of `data`.
int
IWRecordIoWriter::CommonWrite(const void * data, const uint64_t nbytes) {
  const char * p = reinterpret_cast<const char *>(data);
  // cerr << "Writing " << nbytes << " bytes\n";
  if (! _output.write(p, nbytes)) {
    cerr << "IWRecordIoWriter::CommonWrite:cannot write " << nbytes << " bytes\n";
    return 0;
  }

  return _output.write_if_buffer_holds_more_than(8192);
}

ssize_t
IWRecordIoWriter::WriteLength(const uint64_t nbytes) {
  IWString tmp;
  tmp << nbytes << '\n';
  return CommonWrite(tmp.data(), tmp.length());
}

ssize_t
IWRecordIoWriter::WriteData(const void * data, uint64_t nbytes) {
  if (nbytes == 0) {
    return 1;
  }
  return CommonWrite(data, nbytes);
}

}  // namespace iwrecordio
