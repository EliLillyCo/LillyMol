#ifndef MULTICONFORMER_H
#define MULTICONFORMER_H

#include "gfp.h"
#include "fixed_size_counted_fingerprint.h"

template <typename T>
class
Multiconformer_Base
{
  protected:
    int _n;

    T * _fp;

  public:
    Multiconformer_Base();
    ~Multiconformer_Base();

    double tanimoto (const Multiconformer_Base<T> &) const;
    double average_tanimoto (const Multiconformer_Base<T> &) const;
};

class Multiconformer_01 : public Multiconformer_Base<IWDYFP>
{
  private:
  public:
    int construct_from_daylight_ascii_representation (const const_IWSubstring &);
    int construct_from_tdt_record (const const_IWSubstring &);
    int debug_print (std::ostream &) const;

};

class Multiconformer_Fixed_Counted : public Multiconformer_Base<Fixed_Size_Counted_Fingerprint_uchar>
{
  private:
  public:
    int construct_from_daylight_ascii_representation (const const_IWSubstring &);
    int construct_from_tdt_record (const const_IWSubstring &);
    int debug_print (std::ostream &) const;
};

class Multiconformer_Sparse : public Multiconformer_Base<Sparse_Fingerprint>
{
  private:
  public:
    int construct_from_tdt_record (const const_IWSubstring &);
    int construct_from_daylight_ascii_representation (const const_IWSubstring &);
    int debug_print (std::ostream &) const;
};

#ifdef MULTICONFORMER_IMPLEMENTATION

template <typename T>
Multiconformer_Base<T>::Multiconformer_Base()
{
  _n = 0;
  _fp = NULL;

  return;
}

template <typename T>
Multiconformer_Base<T>::~Multiconformer_Base()
{
  assert (-87 != _n);    // do not delete twice

  _n = -87;

  if (NULL != _fp)
    delete [] _fp;

  return;
}

template <typename T>
double
Multiconformer_Base<T>::tanimoto(const Multiconformer_Base<T> & rhs) const
{
  if (0 == _n || 0 == rhs._n)
  {
    if (0 == _n && 0 == rhs._n)
      return 1.0;

    return 0.0;
  }

  double rc = 0.0;

  for (int i = 0; i < _n; i++)
  {
    for (int j = 0; j < rhs._n; j++)
    {
      double tmp = _fp[i].tanimoto(rhs._fp[j]);
//    cerr << "Between component " << i << " and " << j << " value " << tmp << endl;
      if (tmp > rc)
        rc = tmp;
    }
  }

  return rc;
}

template <typename T>
double
Multiconformer_Base<T>::average_tanimoto(const Multiconformer_Base<T> & rhs) const
{
  if (0 == _n || 0 == rhs._n)
  {
    if (0 == _n && 0 == rhs._n)
      return 1.0;

    return 0.0;
  }

  double rc = 0.0;

  for (int i = 0; i < _n; i++)
  {
    for (int j = 0; j < rhs._n; j++)
    {
      double tmp = _fp[i].tanimoto(rhs._fp[j]);
      cerr << " i = " << i << " j = " << j << " value " << tmp << endl;
      rc += tmp;
    }
  }

  return rc / static_cast<double>(_n * rhs._n);
}

#endif
#endif
