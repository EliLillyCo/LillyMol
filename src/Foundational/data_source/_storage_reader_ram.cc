#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#define STORAGE_READER_IMPLEMENTATION

#include "iwmmap.h"

template int IW_Storage_Reader<IW_File_Contents>::records_remaining(int) const;
template int IW_Storage_Reader<IW_File_Contents>::next_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_File_Contents>::seekg(unsigned long, int);
template size_t IW_Storage_Reader<IW_File_Contents>::tellg() const;
template IW_Storage_Reader<IW_File_Contents>::~IW_Storage_Reader();

IWString_Data_Source_RAM::IWString_Data_Source_RAM(char const* fname) : IW_Storage_Reader<IW_File_Contents>(fname)
{
  return;
}

IW_File_Contents::IW_File_Contents(const char * f)
{
  _p = nullptr;
  _len = 0;

  _good = 0;

  int fd;
  int is_stdin = 0;

  if (0 == strcmp(f, "-"))
  {
    fd = 0;
    is_stdin = 1;
  }
  else
  {
    fd = IW_FD_OPEN(f, O_RDONLY);

#ifdef IWCYGWIN
    if (fd >= 0)
#else
    if (fd > 0)
#endif
      ;
    else
    {
      cerr << "IW_File_Contents::constructor:cannot open '" << f << "'\n";
  
      return;
    }
  }

  const int rc =  _read_the_data(fd);

  if (! is_stdin)
    IW_FD_CLOSE(fd);

  if (! rc)
  {
    cerr << "IW_File_Contents::constructor:cannot process '" << f << "\'n";
    _good = 0;

    return;
  }

  _good = 1;
  _fname = f;

  return;
}

int
IW_File_Contents::_read_the_data(const int fd)
{
  _acc.resize(8192);

  const size_t lrecl = 4096;
  unsigned char * buffer = new unsigned char[lrecl]; std::unique_ptr<unsigned char[]> free_buffer(buffer);

//cerr << "IW_File_Contents::_read_the_data:reading from " << fd << endl;

  size_t bytes_read;
  while ((bytes_read = IW_FD_READ(fd, buffer, lrecl)) > 0)
  {
//  cerr << "Read " << bytes_read << " bytes now have " << _acc.number_elements() << endl;
    _acc.add(buffer, bytes_read);
    _acc.add('\n');
  }

  _p = _acc.rawdata();
  _len = _acc.number_elements();

  _good = 1;

  return 1;
}

IW_File_Contents::~IW_File_Contents()
{
  _good = 0;

  return;
}

template IW_Storage_Reader<IW_File_Contents>::IW_Storage_Reader(char const*);
template IW_Storage_Reader<IW_MMapd_File>::IW_Storage_Reader(char const*);
template IW_Storage_Reader<IW_MMapd_File>::IW_Storage_Reader();
template void IW_Storage_Reader<IW_MMapd_File>::_default_values();
template void IW_Storage_Reader<IW_File_Contents>::_default_values();
template int IW_Storage_Reader<IW_File_Contents>::_next_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_File_Contents>::_fetch_previous_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_File_Contents>::_find_next_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_MMapd_File>::eof() const;
template int IW_Storage_Reader<IW_MMapd_File>::skip_records(IW_Regular_Expression_Template<IW_grep_25_regex>&, int);
template int IW_Storage_Reader<IW_MMapd_File>::echo(IWString_and_File_Descriptor&, unsigned long);
template int IW_Storage_Reader<IW_MMapd_File>::most_recent_record(IWString&);
template int IW_Storage_Reader<IW_MMapd_File>::set_ignore_pattern(const_IWSubstring const&);
template int IW_Storage_Reader<IW_MMapd_File>::skip_past(char const*);
template int IW_Storage_Reader<IW_MMapd_File>::skip_to(char const*);
template int IW_Storage_Reader<IW_MMapd_File>::skip_records(int);
template int IW_Storage_Reader<IW_MMapd_File>::_fetch_previous_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_MMapd_File>::_next_record_inner(const_IWSubstring&);
template int IW_Storage_Reader<IW_MMapd_File>::_find_next_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_MMapd_File>::_next_record(const_IWSubstring&);
template int IW_Storage_Reader<IW_File_Contents>::_next_record_inner(const_IWSubstring&);
