#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "iwbits.h"
#include "iwstring.h"

#define FIXED_SIZE_COUNTED_FINGERPRINT_IMPLEMENTATION 
#include "fixed_size_counted_fingerprint.h"

/*
  The fingerprint may or may not have a value for nset in it
*/

int
Fixed_Size_Counted_Fingerprint_uchar::construct_from_tdt_record (const const_IWSubstring & buffer)
{
  assert (buffer.ends_with ('>'));

  int open_angle_bracket = buffer.index ('<');

  if (open_angle_bracket < 0)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::construct_from_tdt_record:no < in TDT record, impossible\n";
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

// There is at least nset and maybe also nbits

  const_IWSubstring token;
  int i = 0;

  other_info.nextword (token, i, ';');

  if (! token.numeric_value (_nset) || _nset < 0)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::construct_from_tdt_record:invalid nset value '" << token << "'\n";
    return 0;
  }
  else
    compute_nset ();

  if (other_info.nextword (token, i, ';'))    // info on nbits present
  {
    if (! token.numeric_value (_nbits) || _nbits < 0)
    {
      cerr << "Fixed_Size_Counted_Fingerprint_uchar::construct_from_tdt_record:invalid nbits value '" << token << "'\n";
      return 0;
    }
  }

  return _n;
}

int
Fixed_Size_Counted_Fingerprint_uchar::_construct_from_tdt_record (const const_IWSubstring & fp)
{
  IW_Bits_Base b;

  if (! b.construct_from_daylight_ascii_bit_rep (fp.rawchars (), fp.length ()))
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::_construct_from_tdt_record:invalid fingerprint data\n";
    cerr << "'" << fp << "'\n";
    return 0;
  }

  if (0 != b.nbits () % IW_BITS_PER_BYTE)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::_construct_from_tdt_record:fingerprint representation contains " << b.nbits () << " bits, impossible\n";
    return 0;
  }

  if (! resize (b.nbits () / IW_BITS_PER_BYTE))   // each count stored as an 8 bit byte
    return 0;

  memcpy (_count, b.bits (), _n);

  return _n;
}

int
Fixed_Size_Counted_Fingerprint_uchar::write_daylight_ascii_representation (std::ostream & output,
                                             const const_IWSubstring & tag, 
                                             int include_nset)
{
  if (NULL == _count)
    return 0;

  if (include_nset && _nset <= 0)
    compute_nset ();

  IW_Bits_Base b;

  if (! b.construct_from_array_of_bits (_count, _n))
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::write_daylight_ascii_representation:cannot allocate bits for n = " << _n << endl;
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
Fixed_Size_Counted_Fingerprint_uchar::bits_in_common (const Fixed_Size_Counted_Fingerprint_uchar & rhs) const
{
  int rc = 0;

  for (int i = 0; i < _n; i++)
  {
#ifdef DEBUG_FCBIC
    int c = _count[i] * 256 + rhs._count[i];
    cerr << " i = " << i << " lhs " << static_cast<int> (_count[i]) << " rhs " << static_cast<int> (rhs._count[i]) << " index " << c << endl;
    cerr << " In common " << static_cast<int> (bic_table[c]) << endl;
#endif

//  rc += bic_table[_count[i] * 256 + rhs._count[i]];
    if (_count[i] < rhs._count[i])
      rc += _count[i];
    else
      rc += rhs._count[i];
  }

  return rc;
}

/*static int
sum_vector_unsigned_char (const unsigned char * s,
                          int n)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    rc += s[i];
  }

  return rc;
}*/

similarity_type_t
Fixed_Size_Counted_Fingerprint_uchar::tanimoto (const Fixed_Size_Counted_Fingerprint_uchar & rhs) const
{
  int bic = bits_in_common (rhs);

//cerr << "Between fingerprints with " << _nset << " and " << rhs._nset << " bits set, bic " << bic << endl;
//assert (_n == rhs._n);
//assert (bic <= _nset);
//assert (bic <= rhs._nset);

  if (0 == bic)
  {
    if (0 == _nset && 0 == rhs._nset)
      return static_cast<similarity_type_t> (1.0);

    return static_cast<similarity_type_t> (0.0);
  }

  return static_cast<float> (bic) / static_cast<float> (_nset + rhs._nset - bic);
}

int
Fixed_Size_Counted_Fingerprint_uchar::construct_from_array(unsigned char const* b, int nb)
{
  if (! resize(nb))
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uchar::construct_from_array:cannot size for " << nb << " bytes\n";
    return 0;
  }

  _n = nb;
  _nbits = nb;

  _nset = 0;
  for (int i = 0; i < _n; i++)
  {
    _count[i] = b[i];
    if (_count[i])
      _nset++;
  }

  return 1;
}
