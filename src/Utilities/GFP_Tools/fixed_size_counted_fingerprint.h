#ifndef FXDSZSPFP_H
#define FXDSZSPFP_H

class Tversky;
class const_IWSubstring;

/*
  A Fixed_Size_Counted_Fingerprint consists of a set of unsigned bytes holding the count for each 

  Note that the nset value is wierd.  Because these may come from
  compressed forms of sparse fingerprints, _nset may not correspond to
  the number of bits set in this object itself.
*/

template <typename T>
class Fixed_Size_Counted_Fingerprint_Base
{
  protected:
    int _n;
    T * _count;

//  These fingerprints can be generated from sparse fingerprints. If a subset is chosen,
//  then we may need to preserve the original values of nbits and nset. These will be
//  present whenever the fingerprint is read in

    int _nset;
    int _nbits;

//  protected functions

    int _parse_nset_and_nbits (const const_IWSubstring & other_info);

  public:
    Fixed_Size_Counted_Fingerprint_Base ();
    ~Fixed_Size_Counted_Fingerprint_Base ();

    int size () const { return _n;}
    int nset () const { return _nset;}
    int nbits () const { return _nbits;}

    int construct_from_array_of_ints (const int *, int);

//  Since _nset may not be just our own data, it can be set to anything

    void set_nset (int s) { _nset = s;}

    void set (int i, T v) { _count[i] = static_cast<T> (v);}    // no checking, just fast. _nset not changed!

    int resize (int s);

    int compute_nset ();

    int bits_in_common (const Fixed_Size_Counted_Fingerprint_Base<T> &) const;

//  similarity_type_t tanimoto         (const Fixed_Size_Counted_Fingerprint_Base<T> &) const;
    similarity_type_t tversky          (const Fixed_Size_Counted_Fingerprint_Base<T> &, const Tversky &) const;
    similarity_type_t fraction_matched (const Fixed_Size_Counted_Fingerprint_Base<T> &) const;
    similarity_type_t optimistic_distance (const Fixed_Size_Counted_Fingerprint_Base<T> &, const Tversky &);
    similarity_type_t fvb_modified_tanimoto (const Fixed_Size_Counted_Fingerprint_Base<T> &) const;
    similarity_type_t simple_matching (const Fixed_Size_Counted_Fingerprint_Base<T> &) const;
};

class Fixed_Size_Counted_Fingerprint_uchar : public Fixed_Size_Counted_Fingerprint_Base<unsigned char>
{
  private:
    int _construct_from_tdt_record (const const_IWSubstring & fp);

  public:

    int construct_from_array(const unsigned char *, int);

    int construct_from_tdt_record (const const_IWSubstring &);

    int bits_in_common (const Fixed_Size_Counted_Fingerprint_uchar &) const;
    similarity_type_t tanimoto (const Fixed_Size_Counted_Fingerprint_uchar & rhs) const;

    int write_daylight_ascii_representation (std::ostream & os,
                                             const const_IWSubstring & data_item_name, 
                                             int include_nset);
};

class Fixed_Size_Counted_Fingerprint_uint : public Fixed_Size_Counted_Fingerprint_Base<unsigned int>
{
  private:
    int _construct_from_tdt_record (const const_IWSubstring & fp);

  public:
    int construct_from_tdt_record (const const_IWSubstring &);

    int write_daylight_ascii_representation (std::ostream & os,
                                             const const_IWSubstring & data_item_name, 
                                             int include_nset);

    int bits_in_common (const Fixed_Size_Counted_Fingerprint_uint & rhs) const;

    similarity_type_t tanimoto (const Fixed_Size_Counted_Fingerprint_uint & rhs) const;
};

extern const unsigned char bic_table[];


#ifdef FIXED_SIZE_COUNTED_FINGERPRINT_IMPLEMENTATION 

#include "iwstring.h"
#include "misc.h"

#include "various_distance_metrics.h"
#include "tversky.h"
#include "fixed_size_counted_fingerprint.h"

template <typename T>
Fixed_Size_Counted_Fingerprint_Base<T>::Fixed_Size_Counted_Fingerprint_Base ()
{
  _n = 0;
  _count = NULL;
  _nset = 0;
  _nbits = 0;
}

template <typename T>
Fixed_Size_Counted_Fingerprint_Base<T>::~Fixed_Size_Counted_Fingerprint_Base ()
{
  _n = -3;

  if (NULL != _count)
    delete _count;

  return;
}

template <typename T>
int
Fixed_Size_Counted_Fingerprint_Base<T>::compute_nset ()
{
  _nset = 0;

  if (NULL == _count)
    return 0;

  for (int i = 0; i < _n; i++)
  {
    _nset += static_cast<int> (_count[i]);
  }

  return _nset;
}

template <typename T>
int
Fixed_Size_Counted_Fingerprint_Base<T>::resize (int s)
{
  _nset = 0;
  _nbits = 0;

  if (NULL != _count && _n == s)   // already the correct size
    return 1;

  if (NULL != _count)
    delete _count;

  if (0 == s)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_Base::resize: zero size requested, wierd\n";
    _count = NULL;
    _n = 0;
    return 1;
  }
  
  _count = new T[s];

  if (NULL == _count)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_Base::resize: cannot allocate " << s << " items of size " << sizeof (T) << endl;
    return 0;
  }

  _n = s;

  set_vector (_count, _n, static_cast<T> (0));

  return _n;
}

template <typename T>
int 
Fixed_Size_Counted_Fingerprint_Base<T>::bits_in_common (const Fixed_Size_Counted_Fingerprint_Base<T> & rhs) const
{
  int rc = 0;
  for (int i = 0; i < _n; i++)
  {
    if (_count[i] < rhs._count[i])
      rc += _count[i];
    else 
      rc += rhs._count[i];
  }

  return rc;
}

/*template <typename T>
similarity_type_t 
Fixed_Size_Counted_Fingerprint_Base<T>::tanimoto (const Fixed_Size_Counted_Fingerprint_Base<T> & rhs) const
{
  if (0 == _nset || 0 == rhs._nset)
  {
    if (0 == _nset && 0 == rhs._nset)      // otherwise we have identical molecules that will have non-zero distances
      return static_cast<similarity_type_t> (1.0);

    return static_cast<similarity_type_t> (0.0);
  }

  int bic = bits_in_common (rhs);

  return static_cast<float> (bic) / static_cast<float> (_nset + rhs._nset - bic);
}*/

/*
  Great caution needed here because _nset could be larger than _n
*/

template <typename T>
similarity_type_t 
Fixed_Size_Counted_Fingerprint_Base<T>::fvb_modified_tanimoto (const Fixed_Size_Counted_Fingerprint_Base<T> & rhs) const
{
  int n11 = bits_in_common (rhs);

  int n00 = _n - _nset - rhs._nset + n11;    // bits set in neither

  return fligner_verducci_blower (_n, _nset, rhs._nset, n00, n11);
}

template <typename T>
similarity_type_t 
Fixed_Size_Counted_Fingerprint_Base<T>::tversky (const Fixed_Size_Counted_Fingerprint_Base<T> & rhs,
                                                 const Tversky & tv) const
{
  if (0 == _nset || 0 == rhs._nset)
    return tv.tanimoto (_nset, rhs._nset, 0);

  int bic = bits_in_common (rhs);

  return tv.tanimoto (_nset, rhs._nset, bic);
}

/*
  When reading a fingerprint we will optionally have nset and nbits following the fingerprint
*/

template <typename T>
int
Fixed_Size_Counted_Fingerprint_Base<T>::_parse_nset_and_nbits (const const_IWSubstring & other_info)
{
  const_IWSubstring token;
  int i = 0;

  other_info.nextword (token, i, ';');

  if (! token.numeric_value (_nset) || _nset < 0)
  {
    cerr << "Fixed_Size_Counted_Fingerprint_uint::construct_from_tdt_record:invalid nset value '" << token << "'\n";
    return 0;
  }
  else
    compute_nset ();

  if (other_info.nextword (token, i, ';'))    // info on nbits present
  {
    if (! token.numeric_value (_nbits) || _nbits < 0)
    {
      cerr << "Fixed_Size_Counted_Fingerprint_uint::construct_from_tdt_record:invalid nbits value '" << token << "'\n";
      return 0;
    }
  }
  else
    _nbits = _n;

  return _n;
}

template <typename T>
int
Fixed_Size_Counted_Fingerprint_Base<T>::construct_from_array_of_ints (const int * c, int nb)
{
  if (NULL == _count)
    _count = new T[nb];
  else if (NULL != _count && _n == nb)
    ;
  else
  {
    delete [] _count;
    _count = new T[nb];
  }

  _nbits = nb;
  _n = nb;

  _nset = 0;

  for (int i = 0; i < _n; i++)
  {
    _count[i] = static_cast<T>(c[i]);
    _nset += c[i];
  }

  return 1;
}

#endif
#endif
