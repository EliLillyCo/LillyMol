#include <stdlib.h>

#define MULTICONFORMER_IMPLEMENTATION

#include "gfp.h"
#include "multi_conformer.h"

int
Multiconformer_01::construct_from_tdt_record (const const_IWSubstring & s)
{
  assert (s.ends_with('>'));

  const_IWSubstring tmp(s);

  tmp.chop();
  tmp.remove_up_to_first ('<');

  return construct_from_daylight_ascii_representation (tmp);
}

int
Multiconformer_Fixed_Counted::construct_from_tdt_record (const const_IWSubstring & s)
{
  assert (s.ends_with('>'));

  const_IWSubstring tmp(s);

  tmp.chop();
  tmp.remove_up_to_first ('<');

  return construct_from_daylight_ascii_representation (tmp);
}

int
Multiconformer_Fixed_Counted::construct_from_daylight_ascii_representation(const const_IWSubstring & s)
{
  if (0 == s.length())
  {
    _n = 0;
    _fp = NULL;
    return 1;
  }

  _n = s.nwords(',');

  _fp = new Fixed_Size_Counted_Fingerprint_uchar[_n];

  int i = 0;
  const_IWSubstring token;
  int ndx = 0;

  while (s.nextword(token, i, ','))
  {
    IW_Bits_Base tmp;
    if (! tmp.construct_from_daylight_ascii_representation(token))
    {
      cerr << "Multiconformer_Fixed_Counted::construct_from_daylight_ascii_representation:invalid bit representation\n";
      cerr << "'" << token << "'\n";
      return 0;
    }

    _fp[ndx].construct_from_array(reinterpret_cast<const unsigned char *>(tmp.bits()), tmp.nbits() / IW_BITS_PER_BYTE);

    ndx++;
  }

  return _n;
}

int
Multiconformer_01::construct_from_daylight_ascii_representation(const const_IWSubstring & s)
{
  if (0 == s.length())
  {
    _n = 0;
    _fp = NULL;
    return 1;
  }

  _n = s.nwords(',');

  _fp = new IWDYFP[_n];

  int i = 0;
  const_IWSubstring token;
  int ndx = 0;

  while (s.nextword(token, i, ','))
  {
    if (! _fp[ndx].construct_from_daylight_ascii_bit_rep(token))
    {
      cerr << "Multiconformer_01::construct_from_daylight_ascii_representation:invalid bit representation\n";
      cerr << "'" << token << "'\n";
      return 0;
    }

    ndx++;
  }

  return _n;
}

int
Multiconformer_Sparse::construct_from_daylight_ascii_representation(const const_IWSubstring & s)
{
  int n = s.ccount(',');

  if (0 == s.length())
  {
    _n = 0;
    _fp = NULL;
    return 1;
  }

  _n = n + 1;

  _fp = new Sparse_Fingerprint[_n];

  if (NULL == _fp)
  {
    cerr << "Multiconformer_Sparse::construct_from_daylight_ascii_representation:cannot allocate " << _n << " fingerprints\n";
    return 0;
  }

  n = 0;
  int i = 0;
  const_IWSubstring token;
  while (s.nextword(token, i, ','))
  {
    if (! _fp[n].construct_from_daylight_ascii_representation(token))
    {
      cerr << "Invalid sparse fingerprint specification '" << token << "'\n";
      return 0;
    }

    n++;
  }

  assert (n == _n);

  return _n;
}

int 
Multiconformer_Sparse::construct_from_tdt_record (const const_IWSubstring & s)
{
  assert (s.ends_with('>'));

  const_IWSubstring tmp(s);

  tmp.chop();
  tmp.remove_up_to_first ('<');

  return construct_from_daylight_ascii_representation (tmp);
}

int
Multiconformer_01::debug_print(std::ostream & os) const
{
  os << "Multiconformer_01:contains " << _n << " components";
  if (_n > 0)
    os << ", each of " << _fp[0].nbits() << " bits";
  os << endl;

  return 1;
}

int
Multiconformer_Fixed_Counted::debug_print(std::ostream & os) const
{
  os << "Multiconformer_Fixed_Counted:contains " << _n << " components";
  if (_n > 0)
    os << ", each of " << _fp[0].nbits() << " bits";
  os << endl;

  return 1;
}

int
Multiconformer_Sparse::debug_print (std::ostream & os) const
{
  os << "Multiconformer_Sparse::contains " << _n << " components\n";

  return 1;
}

template class Multiconformer_Base<IWDYFP>;
template class Multiconformer_Base<Fixed_Size_Counted_Fingerprint_uchar>;
template class Multiconformer_Base<Sparse_Fingerprint>;
