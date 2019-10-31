#ifndef IWDM_BASE_H
#define IWDM_BASE_H

#include <iostream>
#include <stdint.h>

#if (__GNUC__ == 3) && (__GNUC_MINOR__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

#if (__GNUC__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

/*
  Objects for storing pair-wise relationships between objects
*/

#include "iw_stl_hash_map.h"
#include "iwstring_data_source.h"
#include "accumulator.h"

/*
  We can convert a distance matrix into a .nn file suitable for nplotnn.
  To do that, we need various conditions about what is to be written.
  Rather than creating lots of arguments, create an object that holds
  all the conditions
*/

template <typename T>
class DM_to_NN_Conditions
{
  private:
    int _min_neighbours;
    int _max_neighbours;
    T   _min_distance;
    T   _max_distance;

//  We can emulate the output from gfp_nearneighbours -h -o

    int _dash_ho;

  public:
    DM_to_NN_Conditions ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int recognised (const IWString & n, int error_occurred);   // parsing command line directives

    int max_neighbours () const { return _max_neighbours;}
    void set_max_neighbours (int m) { _max_neighbours = m;}

    int min_neighbours () const { return _min_neighbours;}
    void set_min_neighbours (int m) { _min_neighbours = m;}

    T max_distance () const { return _max_distance;}
    void set_max_distance (T m) { _max_distance = m;}

    int dash_ho () const { return _dash_ho;}
};

/*
  When forming a neighbour list, we want an object that holds the index
  of the neighbour and the distance
*/

template <typename T>
class ID_Distance
{
  protected:
    unsigned int _id;
    T _distance;

  public:
    ID_Distance ();

    unsigned int id () const { return _id;}
    void set_id (unsigned int i) { _id = i;}

    void set_id_and_distance (int i, T d) { _id = i; _distance = d;}

    T distance () const { return _distance;}
    void set_distance (T d) { _distance = d;}

    int iwqsortcompare (const ID_Distance<T> &);
};

template <typename T>
class IWDistanceMatrixBase : public IW_STL_Hash_Map_int
{
  public:
    typedef T & reference;
    typedef const T & const_reference;

    typedef T * iterator;
    typedef const T * const_iterator;

  protected:
    unsigned int _number_molecules;

    T * _zdata;

    const T * _end;

//  When creating a new matrix, we need an initialiser

    T _initialiser;

//  We save time if we figure out in advance where each row starts

    uint64_t * _row_offset;

//  If we are reading a file, we may know the ID's of our components

    IWString * _id;

// private functions

    uint64_t _offset (int, int) const;

    void _free_all_dynamically_allocated_arrays ();

    void _build_row_offset (uint64_t un);

    int _fetch_neighbours (int target, ID_Distance<T> * tmp) const;

    void _process_size_directive (const const_IWSubstring & svalue, int & error);
    void _process_id_directive (const_IWSubstring directive, const const_IWSubstring & svalue, int & error);
    int _read_the_data (iwstring_data_source & input);
    int _do_read (iwstring_data_source &, double &, double &);

    int _neighbours_sorted_by_distance (int, resizable_array<int> &, ID_Distance<T> *) const;

    int _create_nn_file (int, const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<T> &,
                        ID_Distance<T> *,
                        std::ostream &) const;
    int _create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<T> &,
                        ID_Distance<T> *,
                        std::ostream &) const;

  public:
    IWDistanceMatrixBase ();
    ~IWDistanceMatrixBase ();

    int debug_print (std::ostream & os) const;
    int echo        (std::ostream & os) const;

    void set_initialiser (T i) { _initialiser = i;}

    int resize (int);
    int resize (int, T);    // 2nd arg is initialiser

    int active () const { return NULL != _zdata;}

    int number_molecules () const { return _number_molecules;}

    const IWString & id (int i) const { return _id[i];}
    int  set_id (int, const IWString &);

    int which_item_has_id (const const_IWSubstring &) const;

//  Normally do_read fails if it encounters unrecognised header records. If the optional
//  int argument is set, it ignores unrecognised records

    int do_read (const const_IWSubstring &, int = 0);
    int do_read (iwstring_data_source &, int = 0);

    int recognised (const const_IWSubstring &, const const_IWSubstring &, int &);   // when reading control files

    uint64_t items_allocated () const;
    T * rawdata () const { return _zdata;}
    const uint64_t * row_offset () const { return _row_offset;}
    const T * rawdata_for_row (int) const;

    T maxval () const;
    T minval () const;

    int set (int, int, T);
    int set (const IWString & id1, const IWString & id2, T);

    void increment (int, int, T);

    int do_write (const const_IWSubstring &);
    int do_write (std::ostream &);
    int do_write (IWString_and_File_Descriptor &);
    int write_header (std::ostream &) const;
    int write_header (IWString_and_File_Descriptor &) const;
    int write_data (std::ostream &);
    int write_data (IWString_and_File_Descriptor &);

    T zvalue (int, int) const;
    T zvalue (const IWString &, const IWString &) const;

//  Fast but dangerous version - note no checking!

    T zvalue_i_less_than_j (int i, int j) const { return _zdata[_row_offset[i] + j - i - 1];}

    int neighbours_sorted_by_distance (int, resizable_array<int> &) const;

    int create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<T> &,
                        std::ostream &) const;

    int swap_items(int, int);

//  Four different const combinations
    
    template <typename F> void each (const F & f) const
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        f(_zdata[i]);
      }
    }
    template <typename F> void each (F & f) const
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        f(_zdata[i]);
      }
    }
    template <typename F> void each (const F & f) 
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        f(_zdata[i]);
      }
    }
    template <typename F> void each (F & f) 
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        f(_zdata[i]);
      }
    }
    template <typename F> void each_lambda (F f)  const
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        f(_zdata[i]);
      }
    }

    template <typename F> void change_values (F & f)
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        _zdata[i] = f(_zdata[i]);
      }

      return;
    }

    iterator begin () { return _zdata;}
    const_iterator begin () const { return _zdata;}

//  const_iterator end () const { return _zdata + (_number_molecules * (_number_molecules - 1) / 2);}
    const_iterator end () const { return _end;}

    int remove_items (const int * r);
};

#include "Masquerading_as_Byte.h"

template <typename T>
class IWDistanceMatrixMasquerading_as_Byte : public IWDistanceMatrixBase<unsigned char>,
                                             public Masquerading_as_Byte<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using IWDistanceMatrixBase<unsigned char>::_zdata;
  using Masquerading_as_Byte<T>::_byte_to_T;
  using Masquerading_as_Byte<T>::_minval;
  using Masquerading_as_Byte<T>::_maxval;
  using Masquerading_as_Byte<T>::_string_value;
#endif

  protected:

//  private functions

    int _parse_input_record (const const_IWSubstring &);

    int _read_the_data (iwstring_data_source & input);

    int _create_nn_file (int, const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<int> &,
                        ID_Distance<int> *,
                        std::ostream &) const;
    int _create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<int> &,
                        ID_Distance<int> *,
                        std::ostream &) const;
  public:
    IWDistanceMatrixMasquerading_as_Byte ();
    ~IWDistanceMatrixMasquerading_as_Byte ();

    int debug_print (std::ostream & os) const;
    int echo        (std::ostream & os) const;

    int do_read (const const_IWSubstring &);
    int do_read (iwstring_data_source &);

    int do_write (const const_IWSubstring &);
    int do_write (std::ostream &);
    int do_write (IWString_and_File_Descriptor &);

    int set_minval_and_maxval (T, T);

    int set (int, int, T);

    unsigned char byte_value (int, int) const;
    T zvalue (int, int) const;

    T zvalue_i_less_than_j (int, int) const;

    int create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                        const DM_to_NN_Conditions<int> &,
                        std::ostream &) const;


    template <typename F> void each (const F & f) const
    {
      int n = _number_molecules * (_number_molecules - 1) / 2;

      for (int i = 0; i < n; i++)
      {
        T tmp = _byte_to_T[_zdata[i]];
        f (tmp);
      }
    }
};

class IWDistanceMatrixFloat : public IWDistanceMatrixBase<float>
{
  protected:

  public:
};

class IWDistanceMatrixBaseDouble : public IWDistanceMatrixBase<double>
{
  protected:

  public:
};

class IWDistanceMatrixBaseInt : public IWDistanceMatrixBase<int>
{
  protected:

  public:
};

template <typename T>
class Convert_Similarity_to_Distance
{
  private:
  public:
    void operator () (T & v) const { v = static_cast<T>(1.0) - v;}
};

template <typename T>
class Distance_Matrix_Statistics : public Accumulator <T>
{
  private:
  public:
    void operator () (const T & v) { Accumulator<T>::extra(v);}
};

#ifdef IWDISTANCE_MATRIX_IMPLEMENTATION

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include <fstream>
using std::cerr;
using std::endl;

#include "iwqsort.h"
#include "misc.h"
#include "IWDistanceMatrixBase.h"
#include "iwdmsupport.h"

#include "iwstring_data_source.h"

template <typename T>
IWDistanceMatrixBase<T>::IWDistanceMatrixBase ()
{
  _number_molecules = 0;

  _zdata = NULL;

  _end = NULL;

  _row_offset = NULL;

  _initialiser = static_cast<T>(0);

  _id = NULL;

  return;
}

template <typename T>
IWDistanceMatrixBase<T>::~IWDistanceMatrixBase ()
{
  _free_all_dynamically_allocated_arrays();
}

template <typename T>
void
IWDistanceMatrixBase<T>::_free_all_dynamically_allocated_arrays ()
{
  if (NULL != _zdata)
  {
    delete [] _zdata;
    _zdata = NULL;
  }

  if (NULL != _row_offset)
  {
    delete [] _row_offset;
    _row_offset = NULL;
  }

  if (NULL != _id)
  {
    delete [] _id;
    _id = NULL;
  }

  return;
}

template <typename T>
int
IWDistanceMatrixBase<T>::debug_print (std::ostream & os) const
{
  os << "IWDistanceMatrixBase::debug_print: size " << _number_molecules << endl;

  return os.good();
}

template <typename T>
int
IWDistanceMatrixBase<T>::echo (std::ostream & output) const
{
  IWDistanceMatrixBase<T>::debug_print(output);

  if (0 == _number_molecules)
  {
    cerr << "IWDistanceMatrixBase::echo: empty distance matrix\n";
    return 0;
  }

  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    for (unsigned int j = i + 1; j < _number_molecules; j++)
    {
      unsigned int o = _offset(i, j);

      write_index_and_id(i, _id, output);
      output << ' ';
      write_index_and_id(j, _id, output);
      output << " offset " << o << " value " << _zdata[o] << '\n';
    }
  }

  return output.good();
}

template <typename T>
uint64_t
IWDistanceMatrixBase<T>::items_allocated() const
{
  const uint64_t n64 = static_cast<uint64_t>(_number_molecules);

  return n64 * (n64-1) / 2;
}

template <typename T>
int
IWDistanceMatrixBase<T>::resize(int n)
{
  return resize(n, _initialiser);
}

template <typename T>
void
IWDistanceMatrixBase<T>::_build_row_offset (uint64_t n)
{
  for (uint64_t i = 0; i < n - 1; i++)
  {
    _row_offset[i] = (2 * i * n - i * i - i) / 2;
  }

//#define ECHO_ROW_OFFSET
#ifdef ECHO_ROW_OFFSET
  for (int i = 0; i < n - 1; i++)
  {
    cerr << " i = " << i << " offset " << _row_offset[i] << '\n';
  }
#endif

  return;
}

template <typename T>
int
IWDistanceMatrixBase<T>::resize (int n,
                                 T initialiser)
{
  if (0 == n)
  {
    _free_all_dynamically_allocated_arrays();
    _number_molecules = 0;

    return 1;
  }

  assert (NULL == _zdata);
  assert (n > 0);

  const uint64_t un = static_cast<uint64_t>(n);

  const uint64_t items_needed = un * (un - 1) / 2;

//cerr << "n = " << n << " Allocating " << items_needed << " items of size " << sizeof(T) << endl;

  _zdata = new T[items_needed];
  if (NULL == _zdata)
  {
    cerr << "IWDistanceMatrixBase::resize: cannot allocate array for size " << n << ", number items = " << items_needed << '\n';
    return 0;
  }

  _end = _zdata + items_needed;

  for (unsigned int i = 0; i < items_needed; i++)
  {
    _zdata[i] = initialiser;
  }

//cerr << "Problem size " << n << " needs " << items_needed << " items stored\n";

  _row_offset = new uint64_t[n];
  if (NULL == _row_offset)
  {
    cerr << "IWDistanceMatrixBase::resize: cannot allocate row offset array for size " << n << '\n';
    return 0;
  }

  _build_row_offset(un);

  _id = new IWString[n];

  if (NULL == _id)
  {
    cerr << "IWDistanceMatrixBase::resize:cannot allocate " << n << " ids\n";
    return 0;
  }

  _number_molecules = n;

  return 1;
}

template <typename T>
int
IWDistanceMatrixBase<T>::set_id (int i, const IWString & s)
{
  assert (i >= 0 && static_cast<unsigned int>(i) < _number_molecules);

  _id[i] = s;

//IW_STL_Hash_Map_int::const_iterator f = find(s);

  IW_STL_Hash_Map_int::operator [](s) = i;

  assert (IW_STL_Hash_Map_int::size() > 0);
  assert (IW_STL_Hash_Map_int::size() <= static_cast<unsigned int>(_number_molecules));

  return 1;
}

template <typename T>
uint64_t
IWDistanceMatrixBase<T>::_offset (int i, int j) const
{
  if (i < j)
    return _row_offset[i] + j - i - 1;
  else
    return _row_offset[j] + i - j - 1;
}

template <typename T>
T
IWDistanceMatrixBase<T>::zvalue (int i, int j) const
{
  unsigned int b = _offset(i, j);

//cerr << "Data for " << i << "," << j << " at " << b << '\n';

  return _zdata[b];
}

template <typename T>
T
IWDistanceMatrixBase<T>::zvalue (const IWString & id1,
                              const IWString & id2) const
{
  IW_STL_Hash_Map_int::const_iterator f1 = IW_STL_Hash_Map_int::find(id1);
  if (f1 == IW_STL_Hash_Map_int::end())
  {
    cerr << "IWDistanceMatrixBase::set:no hash value for '" << id1 << "'\n";
    for (IW_STL_Hash_Map_int::const_iterator i = IW_STL_Hash_Map_int::begin(); i != IW_STL_Hash_Map_int::end(); i++)
    {
      cerr << (*i).first << " is " << (*i).second << '\n';
    }
    abort();
    return 0;
  }

  IW_STL_Hash_Map_int::const_iterator f2 = IW_STL_Hash_Map_int::find(id2);
  if (f2 == IW_STL_Hash_Map_int::end())
  {
    cerr << "IWDistanceMatrixBase::set:no hash value for '" << id2 << "'\n";
    abort();
    return 0;
  }

  return zvalue((*f1).second, (*f2).second);
}

template <typename T>
int
IWDistanceMatrixBase<T>::set (int i, int j,
                              T v)
{
  unsigned int o = _offset(i, j);

  _zdata[o] = v;

  return 1;
}

template <typename T>
int
IWDistanceMatrixBase<T>::set (const IWString & id1,
                              const IWString & id2,
                              T v)
{
  IW_STL_Hash_Map_int::const_iterator f1 = IW_STL_Hash_Map_int::find(id1);
  if (f1 == IW_STL_Hash_Map_int::end())
  {
    cerr << "IWDistanceMatrixBase::set:no hash value for '" << id1 << "'\n";
    for (IW_STL_Hash_Map_int::const_iterator i = IW_STL_Hash_Map_int::begin(); i != IW_STL_Hash_Map_int::end(); i++)
    {
      cerr << (*i).first << " is " << (*i).second << '\n';
    }
    return 0;
  }

  IW_STL_Hash_Map_int::const_iterator f2 = IW_STL_Hash_Map_int::find(id2);
  if (f2 == IW_STL_Hash_Map_int::end())
  {
    cerr << "IWDistanceMatrixBase::set:no hash value for '" << id2 << "'\n";
    return 0;
  }

  return set((*f1).second, (*f2).second, v);
}

template <typename T>
void
IWDistanceMatrixBase<T>::increment (int i, int j,
                                    T delta)
{
  const uint64_t o = _offset(i, j);

  _zdata[o] += delta;

  return;
}

template <typename T>
T
IWDistanceMatrixBase<T>::minval () const
{
  T rc = _zdata[0];

  const uint64_t n = items_allocated();

  for (uint64_t i = 1; i < n; i++)
  {
    if (_zdata[i] < rc)
      rc = _zdata[i];
  }

  return rc;
}

template <typename T>
T
IWDistanceMatrixBase<T>::maxval () const
{
  T rc = _zdata[0];

  const uint64_t n = items_allocated();

  for (uint64_t i = 1; i < n; i++)
  {
    if (_zdata[i] > rc)
      rc = _zdata[i];
  }

  return rc;
}

template <typename T>
int
IWDistanceMatrixBase<T>::do_read (const const_IWSubstring & fname,
                                  int ignore_unrecognised_header_records)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "IWDistanceMatrixUse::do_read: cannot open '" << fname << "'\n";
    return 0;
  }

  return do_read(input, ignore_unrecognised_header_records);
}

template <typename T>
void
IWDistanceMatrixBase<T>::_process_size_directive (const const_IWSubstring & svalue,
                                                  int & error)
{
  int n;
  if (! svalue.numeric_value(n) || n < 2)
  {
    cerr << "IWDistanceMatrixBase::_process_size_directive: invalid size directive '" << svalue << "'\n";
    error = 1;
    return;
  }

  if (! resize(n))
    error = 1;

  return;
}

template <typename T>
int
IWDistanceMatrixBase<T>::recognised (const const_IWSubstring & directive,
                                     const const_IWSubstring & svalue,
                                     int & error)
{
  error = 0;

  if ("size" == directive)
  {
    _process_size_directive(svalue, error);
    return 1;
  }

  if (directive.starts_with("ID"))
  {
    _process_id_directive(directive, svalue, error);
    return 1;
  }

  return 0;
}

template <typename T>
int
IWDistanceMatrixBase<T>::do_read (iwstring_data_source & input,
                                  int ignore_unrecognised_header_records)
{
  _number_molecules = 0;
  if (NULL != _zdata)
  {
    delete _zdata;
    _zdata = NULL;
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    buffer.strip_trailing_blanks();

    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    if ('|' == buffer)
      break;

    const_IWSubstring directive, zvalue;
    if (! buffer.split(directive, ' ', zvalue))
    {
      cerr << "IWDistanceMatrixBase::do_read: must have at least two tokens, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }

    int error_occurred = 0;
    if (recognised(directive, zvalue, error_occurred) && 0 == error_occurred)   // GREAT!
      ;
    else if (0 == error_occurred && ignore_unrecognised_header_records)    // Just unrecognised
      ;
    else
    {
      cerr << "error_occurred " << error_occurred << " ignore_unrecognised_header_records " << ignore_unrecognised_header_records << '\n';
      cerr << "IWDistanceMatrixBase::do_read:unrecognised or invalid input record, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return _read_the_data(input);
}

template <typename T>
int
IWDistanceMatrixBase<T>::_read_the_data (iwstring_data_source & input)
{
  assert (_number_molecules > 0);

  if (NULL == _zdata)
  {
    if (! IWDistanceMatrixBase::resize(_number_molecules))
    {
      cerr << "IWDistanceMatrixUse::cannot allocate data for " << _number_molecules << " items\n";
      return 0;
    }
  }

  const uint64_t n64 = static_cast<uint64_t>(_number_molecules);

  const uint64_t items_to_be_read = n64 * (n64 - 1) / 2;   // beware overflow

  if (! input.copy_raw_bytes(_zdata, items_to_be_read * sizeof(T)))
  {
    cerr << "IWDistanceMatrixUse::_read_the_data: cannot read " << items_to_be_read << " items of size " << sizeof(T) << " from input stream\n";
    return 0;
  }

//cerr << "_read_the_data: first values " << _zdata[0] << ' ' << (_zdata[0] - 0.0) << endl;
  if (! iw_little_endian())   // do nothing if big endian
    ;
  else if (1 == sizeof(T))
    ;
  else if (2 == sizeof(T))
    htons_unsigned_short(_zdata, items_to_be_read);
  else if (4 == sizeof(T))
    htonl_unsigned_long(_zdata, items_to_be_read);
  else
    cerr << "Don't know how to byte swap little endian things of size " << sizeof(T) << " bytes\n";

//cerr << "_read_the_data: first values after translate " << _zdata[0] << endl;
//#define ECHO_RAW_VALUES
#ifdef ECHO_RAW_VALUES
  cerr << "Read " << items_to_be_read << " items\n";
  for (unsigned int i = 0; i < items_to_be_read; i++)
  {
    cerr << " i = " << i << " value " << _zdata[i] <<  ' ' << (0.0 == _zdata[i]) << '\n';
  }
#endif

  return 1;
}

template <typename T>
int
IWDistanceMatrixBase<T>::do_write (const const_IWSubstring & fname)
{
  IWString tmp(fname);

  std::ofstream output(tmp.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "IWDistanceMatrixBase::do_write: cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(output);
}

template <typename T>
int
IWDistanceMatrixBase<T>::write_header (std::ostream & output) const
{
  output << "size " << _number_molecules << '\n';

  const uint64_t n = items_allocated();

  const uint64_t bytes_to_write = n * sizeof(T);

  output << "#bytes " << bytes_to_write << '\n';     // write as a comment

  output << '\n';

  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    output << "ID" << i << ' ' << _id[i] << '\n';
  }

  return output.good();
}


template <typename T>
int
IWDistanceMatrixBase<T>::write_header (IWString_and_File_Descriptor & output) const
{
  output << "size " << _number_molecules << '\n';

  const auto n = items_allocated();

  const uint64_t bytes_to_write = n * sizeof(T);

  output << "#bytes " << bytes_to_write << '\n';     // write as a comment

  output << '\n';

  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    output << "ID" << i << ' ' << _id[i] << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return output.good();
}

template <typename T>
int
IWDistanceMatrixBase<T>::do_write (std::ostream & output)
{
  write_header(output);

  output << "|\n";

  return write_data(output);
}

template <typename T>
int
IWDistanceMatrixBase<T>::do_write (IWString_and_File_Descriptor & output)
{
  write_header(output);

  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return write_data(output);
}

template <typename T>
int
IWDistanceMatrixBase<T>::write_data (std::ostream & output)
{
  const uint64_t n = items_allocated();

  uint64_t bytes_to_write = n * sizeof(T);

  int need_to_swap_bytes_back = 1;

//cerr << "Before writing " << _zdata[0] << " endian " << iw_little_endian() << endl;
  if (! iw_little_endian())
    need_to_swap_bytes_back = 0;
  else if (1 == sizeof(T))
    need_to_swap_bytes_back = 0;
  else if (2 == sizeof(T))
    htons_unsigned_short(_zdata, bytes_to_write / 2);
  else if (4 == sizeof(T))
    htonl_unsigned_long(_zdata, bytes_to_write / 4);
  else
  {
    cerr << "Don't know how to byte swap little endian things of size " << sizeof(T) << " bytes\n";
    need_to_swap_bytes_back = 0;
  }

//cerr << "After swapping " << _zdata[0] << endl;

  const uint64_t bsave = bytes_to_write;

//output.write((const char *) _zdata, bytes_to_write);

//cerr << "Writing " << bytes_to_write << " bytes\n";

  const char * w = reinterpret_cast<const char *>(_zdata);

  while (bytes_to_write > 0)
  {
    if (bytes_to_write >= 4096)
    {
      output.write(w, 4096);
      bytes_to_write -= 4096;
      w += 4096;
    }
    else
    {
      output.write(w, bytes_to_write);
      bytes_to_write = 0;
    }
//  cerr << "Still to write " << bytes_to_write << endl;
  }

  if (! need_to_swap_bytes_back)
    ;
  else if (2 == sizeof(T))
    ntohs_unsigned_short(_zdata, bsave / 2);
  else if (4 == sizeof(T))
    ntohl_unsigned_long(_zdata, bsave / 4);

//cerr << "After swapping back  " << _zdata[0] << endl;
  return output.good();
}


template <typename T>
int
IWDistanceMatrixBase<T>::write_data (IWString_and_File_Descriptor & output)
{
  const auto n = items_allocated();

  const uint64_t bytes_to_write = n * sizeof(T);

  int need_to_swap_bytes_back = 1;

//cerr << "write_data:before swap " << _zdata[0] << ", will write " << bytes_to_write << " bytes\n";

  if (! iw_little_endian())
    need_to_swap_bytes_back = 0;
  else if (1 == sizeof(T))
    need_to_swap_bytes_back = 0;
  else if (2 == sizeof(T))
    htons_unsigned_short(_zdata, bytes_to_write / 2);
  else if (4 == sizeof(T))
    htonl_unsigned_long(_zdata, bytes_to_write / 4);
  else
  {
    cerr << "Don't know how to byte swap little endian things of size " << sizeof(T) << " bytes\n";
    need_to_swap_bytes_back = 0;
  }

//cerr << "Swapped form " << _zdata[0] << endl;

  output.write((const char *) _zdata, bytes_to_write);

  if (! need_to_swap_bytes_back)
    ;
  else if (2 == sizeof(T))
    ntohs_unsigned_short(_zdata, bytes_to_write / 2);
  else if (4 == sizeof(T))
    ntohl_unsigned_long(_zdata, bytes_to_write / 4);

//cerr << "Swapped back to " << _zdata[0] << endl;

  return output.good();
}

template <typename T>
IWDistanceMatrixMasquerading_as_Byte<T>::IWDistanceMatrixMasquerading_as_Byte ()
{
  return;
}

template <typename T>
IWDistanceMatrixMasquerading_as_Byte<T>::~IWDistanceMatrixMasquerading_as_Byte ()
{
  return;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::debug_print (std::ostream & os) const
{
  IWDistanceMatrixBase<unsigned char>::debug_print(os);
  Masquerading_as_Byte<T>::debug_print(os);

  return os.good();
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::echo (std::ostream & output) const
{
  if (0 == _number_molecules)
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::echo: empty matrix\n";
    return 0;
  }

  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    for (unsigned int j = i + 1; j < _number_molecules; j++)
    {
      unsigned int o = _offset(i, j);

      unsigned char bvalue = _zdata[o];

      write_index_and_id(i, _id, output);
      output << ' ';
      write_index_and_id(j, _id, output);
      output << " offset " << o << " byte value " << static_cast<unsigned int>(bvalue) << " translated " << this->convert_from_byte(bvalue) << '\n';
    }
  }

  return output.good();
}

template <typename T>
T
IWDistanceMatrixMasquerading_as_Byte<T>::zvalue (int i, int j) const
{
  const uint64_t o = _offset(i, j);

  unsigned char b = _zdata[o];

  return (_byte_to_T[b]);
}

template <typename T>
T
IWDistanceMatrixMasquerading_as_Byte<T>::zvalue_i_less_than_j (int i, int j) const
{
  unsigned char b = _zdata[_row_offset[i] + j - i - 1];

  return _byte_to_T[b];
}

template <typename T>
void
IWDistanceMatrixBase<T>::_process_id_directive (const_IWSubstring directive,    // note local copy
                                   const const_IWSubstring & svalue,
                                   int & error_encountered)
{
  if (0 == _number_molecules)
  {
    cerr << "IWDistanceMatrixBase::recognised: sorry, ID encountered before size\n";
    error_encountered = 1;
    return;
  }
  
  if (NULL == _id)
  {
    _id = new IWString[_number_molecules];
    if (NULL == _id)
    {
      cerr << "IWDistanceMatrixBase::return: sorry, cannot allocate " << _number_molecules << " ID strings\n";
      error_encountered = 1;
      return;
    }
  }

  directive.remove_leading_chars(2);

  int i;
  if (! directive.numeric_value(i) || i < 0 || static_cast<unsigned int>(i) >= _number_molecules)
  {
    cerr << "IWDistanceMatrixBase::recognised: invalid identifier '" << directive << "'\n";
    error_encountered = 1;
    return;
  }

  if (0 != _id[i].length())
    cerr << "Warning duplicate specification of identifier " << i << " now '" << svalue << "' previous '" << _id[i] << "'\n";

  if (contains(svalue))
  {
    cerr << "IWDistanceMatrixBase::_process_id_directive:duplicate id '" << svalue << "'\n";
    error_encountered = 1;
    return;
  }

  _id[i] = svalue;

  IW_STL_Hash_Map_int::operator [] (svalue) = i;
  
  return;
}

template <typename T>
unsigned char
IWDistanceMatrixMasquerading_as_Byte<T>::byte_value (int i, int j) const
{
  unsigned int o = _offset(i, j);

  return _zdata[o];
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::do_read (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::do_read: cannot open '" << fname << "'\n";
    return 0;
  }

  return do_read(input);
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::_parse_input_record (const const_IWSubstring & buffer)
{
  const_IWSubstring directive, svalue;
  if (! buffer.split(directive, ' ', svalue))
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::_parse_input_record: must have at least two tokens\n";
    return 0;
  }

  int error_occurred = 0;
  if (IWDistanceMatrixBase<unsigned char>::recognised(directive, svalue, error_occurred) || error_occurred)
  {
    if (error_occurred)
      return 0;

    return 1;
  }

  if (Masquerading_as_Byte<T>::recognised(directive, svalue, error_occurred))
    return 1;

  cerr << "IWDistanceMatrixMasquerading_as_Byte::_parse_input_record: unrecognised directive '" << directive << "'\n";
  return 0;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::do_read (iwstring_data_source & input)
{
//cerr << "IWDistanceMatrixMasquerading_as_Byte:do_read\n";

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    buffer.strip_trailing_blanks();

    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    if ('|' == buffer)
      break;

//  cerr << "IWDistanceMatrixMasquerading_as_Byte::do_read:examining '" << buffer << "'\n";
    if (! _parse_input_record(buffer))
    {
      cerr << "IWDistanceMatrixMasquerading_as_Byte::do_read: invalid record, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return _read_the_data(input);
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::_read_the_data (iwstring_data_source & input)
{
  if (_maxval < _minval)
  {
    cerr << "IWDistanceMatrixUse::_read_the_data: max " << _maxval << " and minval " << _minval << " inconsistent\n";
    return 0;
  }

  if (! IWDistanceMatrixBase<unsigned char>::_read_the_data(input))
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::_read_the_data: cannot resize/read\n";
    return 0;
  }

  if (! Masquerading_as_Byte<T>::fully_specified())
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::_read_the_data: ranges not fully specified or invalid\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::set (int i, int j,
                                              T v)
{
  unsigned char bvalue;
  if (! Masquerading_as_Byte<T>::convert_to_byte(v, bvalue))
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::set: value for " << i << " " << j << " out of range\n";
    return 0;
  }

  unsigned int o = _offset(i, j);

//cerr << "i = " << i << " j = " << j << " goes to " << o << " value " << v << " transformed to " << static_cast<int>(bvalue) << '\n';

  _zdata[o] = bvalue;

  return 1;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::do_write (const const_IWSubstring & fname)
{
  IWString tmp(fname);

  std::ofstream output(tmp.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::do_write: cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(output);
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::do_write (std::ostream & output)
{
  Masquerading_as_Byte<T>::do_write(output);

  return IWDistanceMatrixBase<unsigned char>::do_write(output);
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::do_write (IWString_and_File_Descriptor & output)
{
  Masquerading_as_Byte<T>::do_write(output);
  return IWDistanceMatrixBase<unsigned char>::do_write(output);
}

template <typename T>
ID_Distance<T>::ID_Distance ()
{
  _id = 0;
  _distance = static_cast<T>(0);

  return;
}

template <typename T>
int
ID_Distance<T>::iwqsortcompare (const ID_Distance<T> & rhs)
{
  if (_distance < rhs._distance)
    return -1;

  if (_distance > rhs._distance)
    return 1;

  return 0;
}

template <typename T>
int
IWDistanceMatrixBase<T>::create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<T> & dmc,
                                      std::ostream & output) const
{
  if (0 == _number_molecules)
  {
    cerr << "IWDistanceMatrixBase::create_nn_file: empty distance matrix\n";
    return 1;
  }

  ID_Distance<T> * tmp = new ID_Distance<T>[_number_molecules];

  int rc = _create_nn_file(id_to_smiles, dmc, tmp, output);

  delete [] tmp;

  return rc;
}

template <typename T>
int
IWDistanceMatrixBase<T>::_create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<T> & dmc,
                                      ID_Distance<T> * tmp,
                                      std::ostream & output) const
{
  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    if (! _create_nn_file(i, id_to_smiles, dmc, tmp, output))
      return 0;
  }

  return output.good();
}

template <typename T>
int
IWDistanceMatrixBase<T>::_neighbours_sorted_by_distance (int target,
                                                        resizable_array<int> & nbr,
                                                        ID_Distance<T> * tmp) const
{
  int rc = _fetch_neighbours(target, tmp);

  for (int i = 0; i < rc; i++)
  {
    nbr.add(tmp[i].id());
  }

  return rc;
}

template <typename T>
int
IWDistanceMatrixBase<T>::neighbours_sorted_by_distance (int target,
                                                        resizable_array<int> & nbr) const
{
  nbr.resize_keep_storage(0);
  nbr.resize (_number_molecules - 1);

  ID_Distance<T> * tmp = new ID_Distance<T>[_number_molecules - 1];

  int rc = _neighbours_sorted_by_distance(target, nbr, tmp);

  delete [] tmp;

  return rc;
}


template <typename T>
int
IWDistanceMatrixBase<T>::_fetch_neighbours (int target,
                                            ID_Distance<T> * tmp) const
{
  unsigned int unsigned_target = static_cast<unsigned int>(target);

  int n = 0;
  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    if (i == unsigned_target)     // no self neighbours
      continue;

    T d = zvalue(target, i);
    tmp[n].set_id(i);
    tmp[n].set_distance(d);
    n++;
  }

  iwqsort(tmp, n);

  return n;
}

template <typename T>
int
IWDistanceMatrixBase<T>::_create_nn_file (int target,
                                      const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<T> & dmc,
                                      ID_Distance<T> * tmp,
                                      std::ostream & output) const
{
  int n = _fetch_neighbours(target, tmp);

  int istart = 0;
  if (dmc.min_neighbours() > 0)
    istart = dmc.min_neighbours();

  if (dmc.max_neighbours() > 0)
    n = dmc.max_neighbours();

  if (dmc.max_distance() > static_cast<T>(0))
  {
    for (int i = istart; i < n; i++)
    {
      if (tmp[i].distance() > dmc.max_distance())
      {
        n = i - 1;
        break;
      }
    }
  }

// Now write the results

  write_smiles_and_id(_id[target], id_to_smiles, output);
  for (int i = 0; i < n; i++)
  {
    const ID_Distance<T> & idi = tmp[i];

    unsigned int id = idi.id();

    if (dmc.dash_ho())
      output << "NBR<" << id << ">\n";
    else
      write_smiles_and_id(_id[id], id_to_smiles, output);

    output << "DIST<" << idi.distance() << ">\n";
  }

  output << "|\n";

  return output.good();
}

template <typename T>
int
IWDistanceMatrixBase<T>::which_item_has_id (const const_IWSubstring & s) const
{
  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    if (s == _id[i])
      return i;
  }

  return -1;
}

template <typename T>
int
IWDistanceMatrixBase<T>::swap_items(int i1,
                                 int i2)
{
  assert (i1 >= 0);
  assert (i2 >= 0);

  if (i1 > i2)
    iwswap(i1, i2);
  else if (i1 < i2)
    ;
  else
  {
    cerr << "IWDistanceMatrixBase::swap_items:identical items " << i1 << " and " << i2 << endl;
    return 0;
  }

  assert (static_cast<unsigned int>(i1) < _number_molecules);
  assert (static_cast<unsigned int>(i2) < _number_molecules);

  iwswap(_id[i1], _id[i2]);

//#define DEBUG_SWAP_ITEMS
#ifdef DEBUG_SWAP_ITEMS
  cerr << "Swapping " << i1 << " and " << i2 << endl;
#endif

// For everything below i1, just swap the I1 and I2 columns

  for (int i = 0; i < i1; i++)
  {
    T * r = _zdata + _row_offset[i];

#ifdef DEBUG_SWAP_ITEMS
    cerr << "Within row " << i << " swap " << (i1 - i - 1) << " with " << (i2 - i - 1) << endl;
#endif
    iwswap(r[i1 - i - 1], r[i2 - i - 1]);
  }

// Now update the I1 row with all the distances from here to I2

  T * i1row = _zdata + _row_offset[i1];

  int ndx = 0;

  for (int i = i1 + 1; i < i2; i++)
  {
    T * r = _zdata + _row_offset[i];

#ifdef DEBUG_SWAP_ITEMS
    cerr << "Within i1 " << i << " swap " << ndx << " with " << (i2 - i - 1) << endl;
#endif
    iwswap(i1row[ndx], r[i2 - i - 1]);

    ndx++;
  }

  ndx++;

  T * i2row = _zdata + _row_offset[i2];

  for (unsigned int i = i2 + 1; i < _number_molecules; i++)
  {
#ifdef DEBUG_SWAP_ITEMS
    cerr << "UPdate " << i << " swap " << ndx << " with " << (i - i2 - 1) << endl;
#endif
    iwswap(i1row[ndx], i2row[i - i2 - 1]);
    ndx++;
  }

  return 1;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<int> & dmc,
                                      std::ostream & output) const
{
  if (0 == _number_molecules)
  {
    cerr << "IWDistanceMatrixMasquerading_as_Byte::create_nn_file: empty distance matrix\n";
    return 1;
  }

  ID_Distance<int> * tmp = new ID_Distance<int>[_number_molecules];

  int rc = _create_nn_file(id_to_smiles, dmc, tmp, output);

  delete [] tmp;

  return rc;
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::_create_nn_file (const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<int> & dmc,
                                      ID_Distance<int> * tmp,
                                      std::ostream & output) const
{
  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    if (! _create_nn_file(i, id_to_smiles, dmc, tmp, output))
      return 0;
  }

  return output.good();
}

template <typename T>
int
IWDistanceMatrixMasquerading_as_Byte<T>::_create_nn_file (int target,
                                      const IW_STL_Hash_Map_String & id_to_smiles,
                                      const DM_to_NN_Conditions<int> & dmc,
                                      ID_Distance<int> * tmp,
                                      std::ostream & output) const
{
  int n = 0;
  for (unsigned int i = 0; i < _number_molecules; i++)
  {
    if (i == static_cast<unsigned int>(target))     // no self neighbours
      continue;

    unsigned char d = byte_value(target, i);
    tmp[n].set_id(i);
    tmp[n].set_distance(d);
    n++;
  }

  iwqsort(tmp, n);

  int istart = 0;
  if (dmc.min_neighbours() > 0)
    istart = dmc.min_neighbours();

  if (dmc.max_neighbours() > 0)
    n = dmc.max_neighbours();

  if (dmc.max_distance() > static_cast<T>(0))
  {
    for (int i = istart; i < n; i++)
    {
      if (tmp[i].distance() > dmc.max_distance())
      {
        n = i - 1;
        break;
      }
    }
  }

// Now write the results

  write_smiles_and_id(_id[target], id_to_smiles, output);
  for (int i = 0; i < n; i++)
  {
    const ID_Distance<int> & idi = tmp[i];

    write_smiles_and_id(_id[idi.id()], id_to_smiles, output);
    unsigned char d = idi.distance();
    output << "DIST<" << _string_value[d] << ">\n";
  }

  output << "|\n";

  return output.good();
}

template <typename T>
int
IWDistanceMatrixBase<T>::remove_items (const int * to_remove)
{
  int dfrom = 0;
  int dto   = 0;
  int x = 0;
  unsigned int y = 1;

  const uint64_t nitems = items_allocated();

  int items_removed = 0;

  for (uint64_t i = 0; i < nitems; i++)
  {
    if (to_remove[x] || to_remove[y])
      items_removed++;
    else
    {
      _zdata[dto] = _zdata[dfrom];
      dto++;
    }

    y++;
    if (y >= _number_molecules)
    {
      x++;
      y = x + 1;
    }

    dfrom++;
  }

  if (0 == items_removed)
  {
    cerr << "IWDistanceMatrixBase::remove_items:nothing removed\n";
    return 0;
  }

// items_removed used for a different purpose now

  items_removed = count_non_zero_occurrences_in_array(to_remove, _number_molecules);

  int new_number_molecules = _number_molecules - items_removed;

  _build_row_offset(new_number_molecules);

  dfrom = 0;
  dto   = 0;

  for (unsigned int i = 0; i < _number_molecules; i++, dfrom++)
  {
    IW_STL_Hash_Map_int::iterator f = this->find(_id[i]);

    if (to_remove[i])
      this->erase(f);
    else
    {
      if (dfrom != dto)
      {
        _id[dfrom] = _id[dto];
        f->second = dto;
      }

      dto++;
    }
  }

  _number_molecules = new_number_molecules;

  return 1;
}

#endif
#endif
