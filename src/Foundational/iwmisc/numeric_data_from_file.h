#ifndef NUMBERS_FROM_FILE_H
#define NUMBERS_FROM_FILE_H

/*
  Just read a set of numbers from a file
*/

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

template <typename T>
class Numeric_Data_From_File : public IW_STL_Hash_Map<IWString, T>
{
  private:
    int _column;
    int _strip_leading_zeros;
    IWString _header;

// private functions

    int _read_data_record (const const_IWSubstring & buffer);

  public:
    Numeric_Data_From_File();

    void set_column(int c) { _column = c;}
    void set_strip_leading_zeros(int s) {_strip_leading_zeros = s;}

    int read_data (iwstring_data_source & input);
    int read_data (const const_IWSubstring & fname);

    const IWString & header() const { return _header;}

    int get_data(const const_IWSubstring &, T &) const;
};

#ifdef NUMERIC_DATA_FILE_IMPLEMENTATION

template <typename T>
Numeric_Data_From_File<T>::Numeric_Data_From_File()
{
  _column = 1;

  _strip_leading_zeros = 0;

  return;
}

template <typename T>
int
Numeric_Data_From_File<T>::read_data (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    std::cerr << "Numeric_Data_From_File::_read_data:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_data (input);
}

template <typename T>
int
Numeric_Data_From_File<T>::read_data (iwstring_data_source & input)
{
  input.set_dos(1);
  input.set_translate_tabs(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! _read_data_record(buffer))
    {
      std::cerr << "Numeric_Data_From_File::_read_data:cannot read data, line " << input.lines_read() << '\n';
      std::cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return static_cast<int>((*this).size());
}

template <typename T>
int
Numeric_Data_From_File<T>::_read_data_record (const const_IWSubstring & buffer)
{
  int i = 0;
  IWString id;

  if (! buffer.nextword(id, i))
  {
    std::cerr << "Numeric_Data_From_File::_read_data_record:cannot extract identifier\n";
    return 0;
  }

  IWString token;

  if (! buffer.nextword(token, i))
  {
    std::cerr << "Numeric_Data_From_File::_read_data_record:not enough tokens\n";
    return 0;
  }

  if (1 != _column)
  {
    if (! buffer.word(_column, token))
    {
      std::cerr << "Numeric_Data_From_File::_read_data_record:cannot extract column '" << (_column + 1) << " from record\n";
      return 0;
    }
  }

  float a;
  if (token.numeric_value(a))
    ;
  else if (0 == (*this).size())
    _header = token;
  else
  {
    std::cerr << "Numeric_Data_From_File::_read_data_record:non numeric value '" << token << "', id '" << id << "'\n";
    return 0;
  }

  if (_strip_leading_zeros)
    id.remove_leading_chars('0');
  
  if (this->contains(id))
    std::cerr << "Duplicate data for '" << id << "', ignored\n";
  else
    (*this)[id] = a;

//std::cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token << "'\n";

  return 1;
}

template <typename T>
int
Numeric_Data_From_File<T>::get_data(const const_IWSubstring & id,
                                    T & a) const
{
  typename Numeric_Data_From_File<T>::const_iterator f = (*this).find(id);

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
