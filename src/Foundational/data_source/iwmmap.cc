#include<sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <iostream>
using std::cerr;
using std::endl;

#define STORAGE_READER_IMPLEMENTATION

#include "iwmmap.h"

IW_MMapd_File::IW_MMapd_File()
{
  _fd = -1;
  _good = 1;

  return;
}

IW_MMapd_File::IW_MMapd_File(const char * fname)
{
  _fd = -1;
  _good = 1;

  _open_file(fname);
}

IW_MMapd_File::~IW_MMapd_File()
{
  _close_file();
}

int
IW_MMapd_File::_close_file()
{
  int rc = 1;

  if (nullptr != _p)
  {
    const auto tmp = munmap((void *)_p, _len);    // horrible old style cast

    if (0 != tmp)
    {
      cerr << "IWString_Data_Source_MMAP::destructor:munmap failed " << tmp << endl;
      rc = 0;
    }
    _p = nullptr;
    _len = 0;
  }

  if (_fd > 0)
  {
    const auto tmp = close(_fd);
    if (0 != tmp)
    {
      cerr << "IW_MMapd_File::destructor:cannot close " << _fd << " rc " << tmp << endl;
      rc = 0;
    }
    _fd = -1;
  }

  return rc;
}

int
IW_MMapd_File::open (const char * fname)
{
  return _open_file(fname);
}

int
IW_MMapd_File::_open_file(const char * fname)
{
  if (_fd > 0)
    _close_file();

  _fd = ::open(fname, O_RDONLY);

#ifdef IWCYGWIN
  if (_fd >= 0)
#else
  if (_fd > 0)
#endif
  {
    _good = 1;
  }
  else
  {
    cerr << "IW_MMapd_File::_open_file:cannot open '" << fname << "' rc " << _fd << '\n';
    return 0;
  }

  struct stat sb;

  const auto tmp = fstat(_fd, &sb);
  if (tmp < 0)
  {
    cerr << "IW_MMapd_File::_cannot stat '" << fname << "' rc " << tmp << endl;
    return 0;
  }

  if (!S_ISREG (sb.st_mode))
  {
    cerr << "IW_MMapd_File::_open_file:not a regulat file '" << fname << "'\n";
    return 0;
  }

  _len = sb.st_size;

  void * v = mmap (0, _len, PROT_READ, MAP_SHARED, _fd, 0);
  if (MAP_FAILED == v)
  {
    cerr << "IW_MMapd_File:_open_file:mmap failed " << fname << " ";
    perror("mmap");
    return 0;
  }

  _p = static_cast<const unsigned char *>(v);

  _good = 1;

  _fname = fname;

  return 1;
}

IWString_Data_Source_MMAP::IWString_Data_Source_MMAP()
{
  _default_values();

  return;
}
int
IW_MMapd_File::do_madvise (const int flag)
{
  if (! _good || nullptr == _p)
    return 0;

  const auto tmp = ::madvise((void *) _p, _len, flag);
  if (0 == tmp)
    return 1;

  cerr << "IWString_and_File_Descriptor::do_madvise:failure\n";
  perror("madvise");
  return 0;
}

#ifdef IMPLEMENT_SOMETIMEQWEQWE
int
IWString_Data_Source_MMAP::grep (char)
{
}

int
IWString_Data_Source_MMAP::grep (const const_IWSubstring & rx)
{
}
#endif

IWString_Data_Source_MMAP::IWString_Data_Source_MMAP(char const * fname) : IW_Storage_Reader<IW_MMapd_File>(fname)
{
  return;
}
//IWString_Data_Source_MMAP::IWString_Data_Source_MMAP(const char * fname) : T(fname)
//{
//}

template int IW_Storage_Reader<IW_MMapd_File>::records_remaining(int) const;
template int IW_Storage_Reader<IW_MMapd_File>::next_record(const_IWSubstring&);
template size_t IW_Storage_Reader<IW_MMapd_File>::tellg() const;
template int IW_Storage_Reader<IW_MMapd_File>::seekg(unsigned long, int);
template IW_Storage_Reader<IW_MMapd_File>::~IW_Storage_Reader();
template int IW_Storage_Reader<IW_MMapd_File>::next_record(IWString&);
template int IW_Storage_Reader<IW_MMapd_File>::echo_records(IWString_and_File_Descriptor&, int);
template size_t IW_Storage_Reader<IW_MMapd_File>::copy_raw_bytes(void*, unsigned long) const;
template int IW_Storage_Reader<IW_MMapd_File>::count_records_starting_with(const_IWSubstring const&);
