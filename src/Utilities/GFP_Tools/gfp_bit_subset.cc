#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"

#define SPARSEFP_HAS_CREATE_SUBSET

#include "Utilities/GFP_Tools/gfp.h"

#include "gfp_bit_subset.h"

using std::cerr;
using std::endl;

GFP_Bit_Subset::GFP_Bit_Subset()
{
  _active = 0;

  _number_properties = 0;
  _properties = nullptr;

  _number_fixed = 0;

  _fixed = nullptr;

  _number_sparse = 0;

  _sparse = nullptr;

  return;
}

GFP_Bit_Subset::~GFP_Bit_Subset()
{
  if (nullptr != _properties)
    delete [] _properties;

  if (nullptr != _fixed)
    delete [] _fixed;

  if (nullptr != _sparse)
    delete [] _sparse;

  return;
}

int
GFP_Bit_Subset::do_read (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "GFP_Bit_Subset::do_read:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_read (input);
}

#define SPECIAL_FLAG_TO_INDICATE_ALL_BITS_USED -4

static int
read_mapping (iwstring_data_source & input,
              int nbits,
              int * m)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;

    if ('|' == buffer)
      break;

    buffer.truncate_at_first(' ');

    int b;

    if (! buffer.numeric_value(b))
    {
      cerr << "Invalid bit number '" << buffer << "'\n";
      return 0;
    }

    if (b < 0 || b >= nbits)
    {
      cerr << "read_mapping:invalid bit " << buffer << "'\n";
      return 0;
    }

    m[b] = 1;
  }

//cerr << "set " << count_non_zero_occurrences_in_array(m, nbits) << " of " << nbits << " bits\n";
  if (nbits == count_non_zero_occurrences_in_array(m, nbits))
    m[0] = SPECIAL_FLAG_TO_INDICATE_ALL_BITS_USED;   

  return 1;
}

Fixed_Width_Fingerprint_Subset::Fixed_Width_Fingerprint_Subset()
{
  _active = 0;

  return;
}

int
Fixed_Width_Fingerprint_Subset::do_read(int nb,
                                        iwstring_data_source & input)
{
  int * m = new_int(nb); std::unique_ptr<int> free_m(m);

  if (! read_mapping(input, nb, m))
    return 0;

  if (SPECIAL_FLAG_TO_INDICATE_ALL_BITS_USED == m[0])
  {
    _active = 0;
    _mask.allocate_space_for_bits(nb);   // we need to create the mask so we know how many bits we have. The mask is not used otherwise
    return 1;
  }

  _mask.construct_from_array_of_ints(m, nb);

//cerr << "Mask has " << _mask.nset() << " of " << _mask.nbits() << " bits set\n";

  if (_mask.nset())
    _active = 1;

  return 1;
}

int
Fixed_Width_Fingerprint_Subset::create_subset(IWDYFP & fp) const
{
  if (! _active)
    return 0;

  int initial_nset = fp.nset();

  fp.IW_Bits_Base::iwand(_mask);

  return fp.compute_nset() - initial_nset;   // force nset computation
}

int
Fixed_Width_Fingerprint_Subset::bits_used_in_fixed_width_fingerprint() const
{
  if (! _active)
    return _mask.nbits();

  return _mask.nset();
}

static int
read_sparse_mapping (iwstring_data_source & input,
                     IW_Hash_Set<unsigned int> & m)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;

    if ('|' == buffer)
      break;

    buffer.truncate_at_first(' ');

    unsigned int b;

    if (! buffer.numeric_value(b))
    {
      cerr << "Invalid bit number '" << buffer << "'\n";
      return 0;
    }

    m.insert(b);
  }

  return 1;
}

/*
  Buffer will contain the number of properties present
*/

int
GFP_Bit_Subset::_read_property_mapping(const const_IWSubstring & buffer,
                                       iwstring_data_source & input)
{
  assert (nullptr == _properties);

  const_IWSubstring tmp(buffer);

  tmp.remove_leading_words(1);

  if (! tmp.numeric_value(_number_properties) || _number_properties <= 0)
  {
    cerr << "GFP_Bit_Subset::_read_property_mapping:invalid nproperities '" << buffer << "'\n";
    return 0;
  }

  _properties = new_int(_number_properties, -1);

  return read_mapping (input, _number_properties, _properties);
}

int
GFP_Bit_Subset::_read_fixed_fingerprint (const const_IWSubstring & buffer,
                                         Fixed_Width_Fingerprint_Subset & fp,
                                         iwstring_data_source & input)
{
  const_IWSubstring tmp(buffer);

  tmp.remove_leading_words(1);

  int nb;

  if (! tmp.numeric_value(nb) || nb <= 0)
  {
    cerr << "GFP_Bit_Subset::_read_fixed_fingerprint:invalid nbits '" << buffer << "'\n";
    return 0;
  }

  return fp.do_read(nb, input);
}

int
GFP_Bit_Subset::do_read (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "GFP_Bit_Subset::read_existing_mapping:cannot read header\n";
    return 0;
  }

  buffer.strip_trailing_blanks();

  if (HEADER_RECORD != buffer)
  {
    cerr << "GFP_Bit_Subset::read_existing_mapping:header record mismatch\n";
    cerr << "Expected '" << HEADER_RECORD << "' got '" << buffer << "'\n";
    return 0;
  }

  int properties_present = 0;
  _number_fixed = 0;
  _number_sparse = 0;

  while (input.next_record (buffer))
  {
    if (0 == buffer.length() || buffer.starts_with('#'))
      continue;

    if (buffer.starts_with(COUNT_FIXED))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(_number_fixed) || _number_fixed < 1)
      {
        cerr << "The number of fixed fingerprints must be a whole +ve number '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(COUNT_SPARSE))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(_number_sparse) || _number_sparse < 1)
      {
        cerr << "The number of sparse fingerprints must be a whole +ve number '" << buffer << "'\n";
        return 0;
      }
    }
    else if (PROPERTIES_TAG == buffer)
      properties_present = 1;
    else if ('|' == buffer)   // end of preamble
      break;
    else
    {
      cerr << "Unrecognised record in bit->feature cross reference '" << buffer << "'\n";
      return 0;
    }
  }

// Did we get anything?

  if (properties_present)
    ;
  else if (_number_fixed)
    ;
  else if (_number_sparse)
    ;
  else
  {
    cerr << "GFP_Bit_Subset::no bits present in preamble\n";
    return 0;
  }

  _active = 1;

  if (_number_fixed > 0)
    _fixed = new Fixed_Width_Fingerprint_Subset[_number_fixed];

  if (_number_sparse > 0)
    _sparse = new IW_Hash_Set<unsigned int>[_number_sparse];

  int ndx_fixed = 0;
  int ndx_sparse = 0;
  int ndx_mpr = 0;
  while (input.next_record(buffer))
  {
    int rc = 0;

    if (buffer.starts_with(PROPERTIES_TAG))
    {
      if (ndx_mpr > 0) {
        cerr << "GFP_Bit_Subset::read_existing_mapping:only one set of properties allowed, line " << input.lines_read() << endl;
        return 0;
      }

      rc = _read_property_mapping (buffer, input);
      ndx_mpr++;
    }
    else if (buffer.starts_with(FIXED_TAG))
    {
      if (ndx_fixed >= _number_fixed) {
        cerr << "GFP_Bit_Subset::read_existing_mapping:too many fixed fingerprints\n";
        return 0;
      }

      rc = _read_fixed_fingerprint (buffer, _fixed[ndx_fixed], input);
      ndx_fixed++;
    }
    else if (SPARSE_TAG == buffer)
    {
      if (ndx_sparse >= _number_sparse) {
        cerr << "GFP_Bit_Subset::read_existing_mapping:too many sparse fingerprints\n";
        return 0;
      }

      rc = read_sparse_mapping(input, _sparse[ndx_sparse]);
      ndx_sparse++;
    }
    else
    {
      cerr << "Unrecognised group heading in bit xref file '" << buffer << "'\n";
      return 0;
    }

    if (rc == 0) {
      return 0;
    }
  }

  if (properties_present && 0 == ndx_mpr)
  {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:no properties\n";
    return 0;
  }

  if (ndx_fixed != _number_fixed)
  {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:missing fixed fingerprints, got " << (ndx_fixed + 1) << ", expected " << _number_fixed << "\n";
    return 0;
  }

  if (ndx_sparse != _number_sparse)
  {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:missing sparse fingerprints, got " << (ndx_sparse + 1) << ", expected " << _number_sparse << "\n";
    return 0;
  }

  return 1;
}

int
GFP_Bit_Subset::debug_print (std::ostream & os) const
{
  os << "Mapping contains";
  if (_number_fixed)
    os << ' ' << _number_fixed << " fixed";
  if (_number_sparse)
    os << ' ' << _number_sparse << " sparse";
  os << endl;

  return os.good();
}

int
GFP_Bit_Subset::reduce_to_subset (IW_General_Fingerprint & fp) const
{
  int rc = 0;    // number bits removed

  if (nullptr != _properties)
    rc += _do_property_subset(fp.molecular_properties_integer());

  for (int i = 0; i < _number_fixed; i++)
  {
    rc += _fixed[i].create_subset(fp[i]);
  }

  for (int i = 0; i < _number_sparse; i++)
  {
    rc += _create_sparse_subset(_sparse[i], fp.sparse_fingerprint(i));
  }

  return rc;
}

template <typename T>
void
Molecular_Properties<T>::create_subset(const int * s)
{
  int ndx = 0;

  for (int i = 0; i < _nproperties; i++)
  {
    if (s[i] > 0)
    {
      _property[ndx] = _property[i];
      ndx++;
    }
  }

  _nproperties = ndx;

  return;
}

template void Molecular_Properties<int>::create_subset(const int *);

int
GFP_Bit_Subset::_do_property_subset (Molecular_Properties_Integer & p) const
{
  assert (_number_properties == p.nproperties());

  if (SPECIAL_FLAG_TO_INDICATE_ALL_BITS_USED == _properties[0])
    return 0;

  p.create_subset (_properties);

  return p.nproperties();
}

int
GFP_Bit_Subset::_create_sparse_subset(const IW_Hash_Set<unsigned int> & h,
                                      Sparse_Fingerprint & fp) const
{
  return fp.ReduceToSubset(h);
}

#ifdef NOW_IN_EXTERNAL
int
Sparse_Fingerprint::create_subset(const IW_Hash_Set<unsigned int> & h)
{
  int ndx = 0;

  _nset = 0;
  for (int i = 0; i < _nbits; i++)
  {
    IW_Hash_Set<unsigned int>::const_iterator f = h.find(_bit[i]);

    if (f == h.end())
      continue;

    _bit[ndx] = _bit[i];
    _count[ndx] = _count[i];
    _nset += _count[ndx];
    ndx++;
  }

  int rc = _nbits - ndx;

  assert (_nset >= 0);

  _nbits = ndx;

  return rc;
}
#endif
