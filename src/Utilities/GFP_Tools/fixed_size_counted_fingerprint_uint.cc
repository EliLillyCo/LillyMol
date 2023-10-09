#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwstring/iwstring.h"

#include "fixed_size_counted_fingerprint.h"

/*
  The fingerprint may or may not have a value for nset in it
*/

int
Fixed_Size_Counted_Fingerprint_uint::construct_from_tdt_record (const const_IWSubstring & buffer)
{
  assert (buffer.ends_with ('>'));

  int open_angle_bracket = buffer.index ('<');

  if (open_angle_bracket < 0)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uint::construct_from_tdt_record:no < in TDT record, impossible\n";
    return 0;
  }

  const_IWSubstring mybuffer (buffer);
  mybuffer.remove_leading_chars (open_angle_bracket + 1);
  mybuffer.chop ();

  const_IWSubstring fp, other_info;

  if (! mybuffer.split (fp, ';', other_info))   // just the fingerprint present, no info about nset or nbits
  {
    if (! _construct_from_tdt_record (mybuffer))
      return 0;

    _nbits = _n;
    compute_nset ();

    return 1;
  }

  if (! _construct_from_tdt_record (fp))
    return 0;

  return _parse_nset_and_nbits (other_info);
}

int
Fixed_Size_Counted_Fingerprint_uint::_construct_from_tdt_record (const const_IWSubstring & fp)
{
  IW_Bits_Base b;

  if (! b.construct_from_daylight_ascii_bit_rep (fp.rawchars (), fp.length ()))
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uint::_construct_from_tdt_record:invalid fingerprint data\n";
    cerr << "'" << fp << "'\n";
    return 0;
  }

  if (0 != b.nbits () % IW_BITS_PER_BYTE)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uint::_construct_from_tdt_record:fingerprint representation contains " << b.nbits () << " bits, impossible\n";
    return 0;
  }

  if (! resize (b.nbits () / IW_BITS_PER_BYTE))   // each count stored as an 8 bit byte
    return 0;

  const unsigned char * zbits = b.bits ();

  for (int i = 0; i < _n; i++)
  {
    _count[i] = static_cast<unsigned int> (zbits[i]);
  }

  return _n;
}

int
Fixed_Size_Counted_Fingerprint_uint::write_daylight_ascii_representation (std::ostream & output,
                                             const const_IWSubstring & tag, 
                                             int include_nset)
{
  if (nullptr == _count)
    return 0;

  if (include_nset && _nset <= 0)
    compute_nset ();

  IW_Bits_Base b;

  if (! b.construct_from_array_of_ints (reinterpret_cast<const int *> (_count), _n))
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uint::write_daylight_ascii_representation:cannot allocate bits for n = " << _n << endl;
    return 0;
  }

  IWString ascii;

  b.daylight_ascii_representation (ascii);

  output << tag << ascii;

  if (include_nset)
    output << ';' << _nset;

  output << ">\n";

  return output.good ();
}

int
Fixed_Size_Counted_Fingerprint_uint::bits_in_common (const Fixed_Size_Counted_Fingerprint_uint & rhs) const
{
  int rc = 0;

  for (int i = 0; i < _n; i++)
  {
    if (0 == _count[i])
      ;
    else if (0 == rhs._count[i])
      ;
    else if (_count[i] < rhs._count[i])
      rc += _count[i];
    else
      rc += rhs._count[i];
  }

  return rc;
}

similarity_type_t
Fixed_Size_Counted_Fingerprint_uint::tanimoto (const Fixed_Size_Counted_Fingerprint_uint & rhs) const
{
  int bic = bits_in_common (rhs);

//cerr << "Fixed_Size_Counted_Fingerprint_uint:tanimoto: _nset " << _nset << " rhs._nset " << rhs._nset << " bic " << bic << endl;

  return static_cast<similarity_type_t> (bic) / static_cast<similarity_type_t> (_nset + rhs._nset - bic);
}
