#ifndef ACTIVITY_DATA_FF_H
#define ACTIVITY_DATA_FF_H

/*
  I have several programmes that read activity data from a file
*/

#include "cmdline.h"
#include "iw_stl_hash_map.h"
#include "iwstring_data_source.h"

template <typename T>
class Activity_Data_From_File : public IW_STL_Hash_Map<IWString, T>
{
  private:
    int _activity_column;
    int _strip_leading_zeros;
    IWString _activity_column_header;

// private functions

    int _read_activity_data_record (const const_IWSubstring & buffer);
    int _read_activity_data (iwstring_data_source & input);

  public:
    Activity_Data_From_File();

    void set_activity_column(int c) { _activity_column = c;}
    void set_strip_leading_zeros(int s) {_strip_leading_zeros = s;}

    const IWString & activity_column_header() const { return _activity_column_header;}

    int read_activity_data (const const_IWSubstring & fname);

    int get_activity(const const_IWSubstring &, T &) const;

    template <typename C> int construct_from_command_line (const C & cl, char, int verbose);
};

#ifdef ACTIVITY_DATA_IMPLEMENATION_H

template <typename T>
Activity_Data_From_File<T>::Activity_Data_From_File()
{
  _activity_column = 1;

  _strip_leading_zeros = 0;

  return;
}

template <typename T> template <typename C>
int
Activity_Data_From_File<T>::construct_from_command_line(const C & cl,
                                                        char flag,
                                                        int verbose)
{
  if (! cl.option_present(flag))
  {
    cerr << "Activity_Data_From_File::construct_from_command_line: no option '" << flag << "' specified\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring fname;
  while (cl.value (flag, fname, i++))
  {
    if (! _read_activity_data(fname))
    {
      cerr << "Activity_Data_From_File::construct_from_command_line:could not read activity data from '" << fname << "'\n";
      return 0;
    }
  }

  return static_cast<int>((*this).size());
}

template <typename T>
int
Activity_Data_From_File<T>::read_activity_data (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Activity_Data_From_File::_read_activity_data:cannot open '" << fname << "'\n";
    return 0;
  }

  return _read_activity_data (input);
}


template <typename T>
int
Activity_Data_From_File<T>::_read_activity_data (iwstring_data_source & input)
{
  input.set_dos(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! _read_activity_data_record(buffer))
    {
      cerr << "Activity_Data_From_File::_read_activity_data_record:cannot read activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return static_cast<int>((*this).size());
}


template <typename T>
int
Activity_Data_From_File<T>::_read_activity_data_record (const const_IWSubstring & buffer)
{
  int i = 0;
  IWString id;

  if (! buffer.nextword(id, i))
  {
    cerr << "Cannot extract identifier\n";
    return 0;
  }

  IWString token;

  if (! buffer.nextword(token, i))
  {
    cerr << "Not enough tokens on experimental data record\n";
    return 0;
  }

  if (1 != _activity_column)
  {
    if (! buffer.word(_activity_column, token))
    {
      cerr << "Activity_Data_From_File::_read_activity_data_record:cannot extract column '" << (_activity_column + 1) << " from record\n";
      return 0;
    }
  }

  float a;
  if (token.numeric_value(a))
    ;
  else if (0 == (*this).size())
    _activity_column_header = token;
  else
  {
    cerr << "Activity_Data_From_File::_read_activity_data_record:non numeric activity value '" << token << "', id '" << id << "'\n";
    return 0;
  }

  if (_strip_leading_zeros)
    id.remove_leading_chars('0');
  
  (*this)[id] = a;

//cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token << "'\n";

  return 1;
}

template <typename T>
int
Activity_Data_From_File<T>::get_activity(const const_IWSubstring & id,
                                         T & a) const
{
  typename Activity_Data_From_File<T>::const_iterator f = (*this).find(id);

  if (f != (*this).end())
  {
    a = (*f).second;
    return 1;
  }

  if (! _strip_leading_zeros)
    return 0;

  IWString tmp(id);
  tmp.remove_leading_chars('0');

  f = (*this).find(tmp);

  if (f == (*this).end())
    return 0;

  a = (*f).second;

  return 1;
}

#endif    /* implementation */

#endif
