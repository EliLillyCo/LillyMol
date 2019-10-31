#ifndef BIT_SUBSET_V2_H
#define BIT_SUBSET_V2_H

#include <unordered_map>
#include <unordered_set>

#define BIT_SUBSET_HEADER_RECORD "BitSubSet"

#include "iwstring_data_source.h"

#include "sparsefp.h"
#include "dyfp.h"

template <typename T>
class Sparse_Subset
{
  private:
    IWString _tag;

    std::unordered_map<unsigned int, T> _v;

  public:
    Sparse_Subset();

    void set_tag (const IWString & s) { _tag = s;}
    const IWString & tag () const { return _tag;}

    int build (iwstring_data_source &);

    template <typename O> int process (Sparse_Fingerprint & sfp, O op, const int remove_bits_not_mentioned) const; 

    template <typename O> int process_record (const const_IWSubstring & buffer, O op, const int remove_bits_not_mentioned, IWString_and_File_Descriptor & output); 

    template <typename O> int do_write(O &) const;
};

template <typename T>
class Dense_Subset
{
  private:
    IWString _tag;

    std::unordered_map<unsigned int, T> _v;

    int _nbits;  

  public:
    Dense_Subset();
    ~Dense_Subset();

    void set_tag (const IWString & s) { _tag = s;}
    const IWString & tag () const { return _tag;}

    int build (iwstring_data_source &);

    void set_nbits (int s) { _nbits = s;}

    template <typename O> int process_record (const const_IWSubstring & buffer, O op, const int remove_bits_not_mentioned, IWString_and_File_Descriptor & output);

    template <typename O> int do_write(O &) const;
};

template <typename T>
class GFP_Bit_Subset
{
  private:
    int _nfixed;
    Dense_Subset<T> * _fixed;

    int _nsparse;
    Sparse_Subset<T> * _sparse;

    int _single_fingerprint_processing;

//  private functions

    int _build_fingerprint (const IWString & buffer, iwstring_data_source & input);

  public:
    GFP_Bit_Subset();
    ~GFP_Bit_Subset();

    int display_file_format_example (std::ostream &) const;

    int build (const char * fname);
    int build (iwstring_data_source &);
    int build_single_fingerprint (iwstring_data_source & input);

    void set_single_fingerprint_processing (int s) { _single_fingerprint_processing = s;}
    int tag_recognised (const const_IWSubstring & s) const;

    template <typename O> int process_record (const_IWSubstring buffer, O op, const int remove_bits_not_mentioned, IWString_and_File_Descriptor & output) const;

    template <typename O> int do_write(O & os) const;
};

#ifdef GFP_BIT_SUBSET_IMPLEMENTATION

template <typename T>
Dense_Subset<T>::Dense_Subset()
{
  _nbits = 0;

  return;
}

template <typename T>
Dense_Subset<T>::~Dense_Subset()
{
   return;
}

template <typename T>
int
common_build (std::unordered_map<unsigned int, T> & v,
              iwstring_data_source & input)
{
  const_IWSubstring buffer;

  int words_per_record = -1;

  while (input.next_record(buffer))
  {
    if (words_per_record < 0)
      words_per_record = buffer.nwords();

    if ('|' == buffer)
      break;

    unsigned int b;
    T x;

    if (1 == words_per_record)
    {
      if (! buffer.numeric_value(b))
      {
        cerr << "GFP_Bit_Subset::common_build:invalid bit '" << buffer << "'\n";
        return 0;
      }
      x = static_cast<T>(1);
    }
    else
    {
      int i = 0;
      const_IWSubstring token;
      if (! buffer.nextword(token, i) || ! token.numeric_value(b))
      {
        cerr << "GFP_Bit_Subset::comm_build:invalid bit '" << buffer << "'\n";
        return 0;
      }

      if (! buffer.nextword(token, i) || ! token.numeric_value(x))
      {
        cerr << "GFP_Bit_Subset::common_build:invalid numeric qualifier '" << buffer << "'\n";
        return 0;
      }
    }

    v[b] = x;
  }

  return 1;
}

template <typename T>
int
Dense_Subset<T>::build (iwstring_data_source & input)
{
  return common_build (_v, input);
}

template <typename T>
int
Sparse_Subset<T>::build (iwstring_data_source & input)
{
  return common_build(_v, input);
}

template <typename T> template <typename O>
int
Dense_Subset<T>::do_write (O & os) const
{
  os << _tag << '\n';

  for (const auto i : _v)
  {
    os << i.first << ' ' << i.second << '\n';
  }

  os << "|\n";

  return 1;
}

template<typename T> template <typename O>
int
Sparse_Subset<T>::do_write (O & os) const
{
  os << _tag << '\n';

  for (const auto i : _v)
  {
    os << i.first << ' ' << i.second << '\n';
  }

  os << "|\n";

  return 1;
}

template <typename T>
GFP_Bit_Subset<T>::GFP_Bit_Subset()
{
  _nfixed  = 0;
  _fixed = nullptr;

  _nsparse = 0;
  _sparse = nullptr;

  _single_fingerprint_processing = 0;

  return;
}

template <typename T>
GFP_Bit_Subset<T>::~GFP_Bit_Subset ()
{
  if (nullptr != _fixed)
    delete [] _fixed;

  if (nullptr != _sparse)
    delete [] _sparse;

  return;
}

template <typename T>
int
GFP_Bit_Subset<T>::build (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "GFP_Bit_Subset::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build (input);
}

template <typename T>
int
GFP_Bit_Subset<T>::_build_fingerprint (const IWString & buffer,
                                    iwstring_data_source & input)
{
  for (int i = 0; i < _nsparse; ++i)
  {
    if (buffer == _sparse[i].tag())
    {
      if (! _sparse[i].build (input))
      {
        cerr << "GFP_Bit_Subset::_build_fingerprint:cannot build data for '" << buffer << "'\n";
        return 0;
      }

      return 1;
    }
  }


  for (auto i = 0; i < _nfixed; ++i)
  {
    if (buffer != _fixed[i].tag())
      continue;

    if (! _fixed[i].build(input))
    {
      cerr << "GFP_Bit_Subset::_build_fingerprint:cannot process fixed fingerprint '" << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "GFP_Bit_Subset::_build_fingerprint:no match to '" << buffer << "'\n";
  return 0;
}

template <typename T>
int
GFP_Bit_Subset<T>::build (iwstring_data_source & input)
{
  input.set_strip_trailing_blanks(1);

  if (_single_fingerprint_processing)
    return build_single_fingerprint (input);

  const_IWSubstring buffer;

  if (! input.next_record (buffer))
  {
    cerr << "GFP_Bit_Subset::build:cannot read header record\n";
    return 0;
  }

  if (BIT_SUBSET_HEADER_RECORD != buffer)
  {
    cerr << "GFP_Bit_Subset::build:invalid header '" << buffer << "'\n";
    return 0;
  }

  resizable_array_p<IWString> fixed_tag;
  resizable_array<int> fixed_size;
  resizable_array_p<IWString> sparse_tag;

  while (input.next_record(buffer))
  {
    if ('|' == buffer)
      break;

    if (0 == buffer.length() || '#' == buffer[0])
      continue;

    int i = 0;
    const_IWSubstring token1, token2;

    if (! buffer.nextword(token1, i) || ! buffer.nextword(token2, i))
    {
      cerr << "GFP_Bit_Subset:build:record must consist of directive and value, '" << buffer << "' invalid\n";
      return 0;
    }

    if ("fixed:" == token1)
    {
      const_IWSubstring s;
      int nb;
      if (! buffer.nextword(s, i) || ! s.numeric_value(nb) || nb < 1)
      {
        cerr << "GFP_Bit_Subset::build:invalid fixed width fingerprint specification '" << buffer << "'\n";
        return 0;
      }

      fixed_tag.add(new IWString(token1));
      fixed_size.add(nb);
    }
    else if ("sparse:" == token1)
    {
      sparse_tag.add(new IWString (token2));
    }
    else
    {
      cerr << "GFP_Bit_Subset::build:unrecognised directive '" << buffer << "'\n";
      return 0;
    }
  }

  _nfixed = fixed_tag.size();
  _nsparse = sparse_tag.size();

  if (0 == _nfixed && 0 == _nsparse)
  {
    cerr << "GFP_Bit_Subset::build:no fingerprints\n";
    return 0;
  }

  if (_nfixed > 0)
  {
    _fixed = new Dense_Subset<T>[_nfixed];
    for (auto i = 0; i < _nfixed; ++i)
    {
      _fixed[i].set_tag(*fixed_tag[i]);
      _fixed[i].set_nbits(fixed_size[i]);
    }
  }

  if (_nsparse > 0)
  {
    _sparse = new Sparse_Subset<T>[_nsparse];
    for (auto i = 0; i < _nsparse; ++i)
    {
      _sparse[i].set_tag(*sparse_tag[i]);
    }
  }

  IWString tag;
  while (input.next_record(tag))
  {
    if (! _build_fingerprint(tag, input))
    {
      cerr << "GFP_Bit_Subset::build:cannot process '" << tag << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
GFP_Bit_Subset<T>::build_single_fingerprint (iwstring_data_source & input)
{
  assert (nullptr == _sparse);

  _nsparse = 1;
  _sparse = new Sparse_Subset<T>[1];

  return _sparse[0].build(input);
}

template <typename T> template<typename O>
int
GFP_Bit_Subset<T>::do_write (O & os) const
{
  os << BIT_SUBSET_HEADER_RECORD << '\n';

  for (auto i = 0; i < _nsparse; ++i)
  {
    os << "sparse:" << _sparse[i].tag() << '\n';
  }

  for (auto i = 0; i < _nfixed; ++i)
  {
    os << "fixed:" << _fixed[i].tag() << ' ' << _fixed[i].nbits() <<'\n';
  }

  os << "|\n";

  for (auto i = 0; i < _nsparse; ++i)
  {
    _sparse[i].do_write(os);
  }

  for (auto i = 0; i < _nfixed; ++i)
  {
    _fixed[i].do_write(os);
  }

  return 1;
}

template <typename T>
int
GFP_Bit_Subset<T>::tag_recognised (const const_IWSubstring & buffer) const
{
  for (auto i = 0; i < _nsparse; ++i)
  {
    if (buffer.starts_with(_sparse[i].tag()))    // no entirely bullet proof, could have very similar fingerprints
      return 1;
  }

  for (auto i = 0; i < _nfixed; ++i)
  {
    if (buffer.starts_with(_fixed[i].tag()))
      return 1;
  }

  return 0;
}

template <typename T>
Sparse_Subset<T>::Sparse_Subset ()
{
  return;
}

template <typename T> template <typename O>
int
GFP_Bit_Subset<T>::process_record (const_IWSubstring buffer,        // note local copy
                                   O op,
                                   const int remove_bits_not_mentioned,
                                   IWString_and_File_Descriptor & output) const
{
  if (! buffer.ends_with('>'))
  {
    cerr << "GFP_Bit_Subset::process_record:must be TDT form, '" << buffer << "' invalid\n";
    return 0;
  }

  auto openangle = buffer.index('<');

  if (openangle < 1)
  {
    cerr << "GFP_Bit_Subset::process_record:no TDT tag detected '" << buffer << "'\n";
    return 0;
  }

  if (_single_fingerprint_processing && 1 == _nsparse)
    return _sparse[0].process_record(buffer, op, remove_bits_not_mentioned, output);

  const const_IWSubstring tag(buffer.rawchars(), openangle);

  for (auto i = 0; i < _nsparse; ++i)
  {
    if (tag != _sparse[i].tag())
      continue;

    return _sparse[i].process_record(buffer, op, remove_bits_not_mentioned, output);
  }
  for (auto i = 0; i < _nfixed; ++i)
  {
    if (tag !=_fixed[i].tag())
      continue;

    return _fixed[i].process_record(buffer, op, remove_bits_not_mentioned, output);
  }

  cerr << "GFP_Bit_Subset::process_record:unrecognised tag '" << buffer << "'\n";
  return 0;
}

template <typename T> template <typename O>
int
Dense_Subset<T>::process_record (const const_IWSubstring & buffer,
                                 O op,
                                 const int remove_bits_not_mentioned,
                                 IWString_and_File_Descriptor & output)
{
  IWDYFP dyfp;

  if (! dyfp.construct_from_tdt_record(buffer))
  {
    cerr << "Dense_Subset::process_record:invalid fingerprint '" << buffer << "'\n";
    return 0;
  }

  if (0 == _tag.length())
  {
    _tag = buffer;
    _tag.truncate_at_first('<');
  }

  for (const auto i : _v)
  {
    const unsigned int b = i.first;

    int c = op(i.first, dyfp.is_set(b), i.second);

    if (c)
      dyfp.set(b, 1);
    else
      dyfp.set(b, 0);
  }

  IWString tmp;
  dyfp.daylight_ascii_tdt(tmp, _tag);

  output << tmp << "\n";

//output << _tag << '<' << buffer << ">\n";

  return 1;
}

template <typename T> template <typename O>
int
Sparse_Subset<T>::process_record (const const_IWSubstring & buffer,
                                  O op,
                                  const int remove_bits_not_mentioned,
                                  IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint sfp;

  if (! sfp.construct_from_tdt_record(buffer))
  {
    cerr << "Sparse_Subset::process_record:cannot interpret '" << buffer << "'\n";
    return 0;
  }

  if (0 == _tag.length())
  {
    _tag = buffer;
    _tag.truncate_at_first('<');
  }

  int rc = process (sfp, op, remove_bits_not_mentioned);

  if (0 == rc)    // no change
    output << buffer << '\n';
  else
  {
    output << _tag << '<';
    sfp.append_daylight_ascii_form_with_counts_encoded(output);
    output << ">\n";
  }

  return rc;
}

template <typename T> template <typename O>
int
Sparse_Subset<T>::process (Sparse_Fingerprint & sfp,
                           O op,
                           const int remove_bits_not_mentioned) const
{
  int rc = 0;

  for (auto i : _v)
  {
    const auto b = i.first;
    auto c = sfp.count_for_bit(b);

    if (0 == c)
      continue;

#ifdef DEBUG_BS_PROCESS
    cerr << "bit " << b << " old value " << c << " delta " << i.second;
#endif
    
    c = op(b, c, i.second);

#ifdef DEBUG_BS_PROCESS
    cerr << " updat4ed to " << c << endl;
#endif

    if (0 == c)
      sfp.remove_bit(b);
    else
      sfp.set_count(b, c);

    rc++;
  }

  if (! remove_bits_not_mentioned)   // we are done
    return rc;

  std::unordered_set<unsigned int> to_remove;

  int i = 0;
  unsigned int b;
  int c;

  while (sfp.next_bit_set(i, b, c))
  {
    const auto f = _v.find(b);

//  cerr << "What about bit " << b << " " << (f == _v.end()) << endl;
    if (f == _v.end())
      to_remove.insert(b);
  }

  if (0 == to_remove.size())
    return rc;

//cerr << "Need to remove " << to_remove.size() << " bits, sfp starts with " << sfp.nbits() << " bits\n";
  for (const auto i : to_remove)
  {
    sfp.remove_bit(i);
  }

//cerr << "After removing " << to_remove.size() << " bits, sfp has " << sfp.nbits() << " bits\n";

  return rc + to_remove.size();   // really just an arbitrary thing
}

template <typename T>
int
GFP_Bit_Subset<T>::display_file_format_example (std::ostream & output) const
{
  return 0;
}

#endif

#endif
