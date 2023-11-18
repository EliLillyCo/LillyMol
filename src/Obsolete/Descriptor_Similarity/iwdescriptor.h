#ifndef IWDESCRIPTOR_H
#define IWDESCRIPTOR_H

#include <math.h>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#define DISTANCE_CARTESIAN 1
#define DISTANCE_MANHATTAN 2
#define DISTANCE_TAXONOMIC 3
#define DISTANCE_CANBERRA 4
#define DISTANCE_COSINE 5
#define DISTANCE_MINKINF 6
#define DISTANCE_BRAY_CURTIS 7
#define DISTANCE_CENSORED_EUCLIDIAN 8
#define DISTANCE_IW_RATIO 9
#define DISTANCE_TANIMOTO 10
#define DISTANCE_RUSSEL_RAO 11
#define DISTANCE_CONTINUOUS_TANIMOTO 12
#define DISTANCE_FORBES 13
#define DISTANCE_CONTINUOUS_FORBES 14
#define DISTANCE_CIRCLE_PRODUCT 15
#define DISTANCE_CZEKANOWSKI 16
#define DISTANCE_PEARSON 17
#define DISTANCE_DOT_PRODUCT 18

/*
  There are two template items.

  D is the type of data being stored
  R is the result of a distance computation
*/

template <typename D, typename R>
class IWDescriptors
{
  private:
    IWString _id;

    D * _d;

//  Some metrics require knowledge of the average and variance

    double _average;
    double _pearson_denominator;

//  the continuous tanimoto requires the product of the non-zero elements

    double _product_of_non_zero_items;

//  The tanimoto metric requires the sum of the items

    double _nset;

//  private functions

    R _cartesian_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _manhattan_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _taxonomic_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _canberra_distance  (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _cosine_distance    (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _minkowski_infinite_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _bray_curtiss_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _censored_euclidian_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _iw_ratio_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _tanimoto_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _russel_rao_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _continuous_tanimoto_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _forbes_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _continuous_forbes_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _distance_circle_product (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _distance_czekanowski (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _distance_pearson (const IWDescriptors<D, R> & rhs, int n, const int *) const;
    R _dot_product_distance (const IWDescriptors<D, R> & rhs, int n, const int *) const;

    int _build (const const_IWSubstring & buffer, int i, int number_descriptors, int & fatal);

  public:
    IWDescriptors();
    ~IWDescriptors();

    const IWString & id() const { return _id;}

    int build (const const_IWSubstring &, int, int &);

    R distance (const IWDescriptors<D, R> &, int dtype, int n, const int *) const;

    template <typename C> R distance (const IWDescriptors<D, R> &, const C &, int, const int *) const;

    const D * rawdata() const { return _d;}

//  No allowance for doing a subset of descriptors when doing average and std.
//  this is a serious flaw. Hopefullu it won't happen

    void compute_pearson_data(int);

    void compute_product_of_non_zero_items(int);
    void compute_nset (int);

    double nset () const { return _nset;}

    double average() const { return _average;}
};

extern void set_descriptors_may_contain_big_E (int d);
extern int descriptors_may_contain_big_E();

extern int determine_distance_type (const const_IWSubstring & x, int verbose);

extern int display_distance_metrics (std::ostream &, char);

extern int needs_ave_and_variance(int);

#ifdef IWDESCRIPTOR_IMPLEMENTATION

#include <math.h>

// Does not belong in a header file. TODO ianwatson fix.
//static IWString missing_value ('.');

template <typename D, typename R>
IWDescriptors<D, R>::IWDescriptors()
{
  _product_of_non_zero_items = -1.0;

  _d = NULL;

  return;
}

template <typename D, typename R>
IWDescriptors<D, R>::~IWDescriptors()
{
  if (NULL != _d)
    delete [] _d;

  return;
}

static int
contains_E_in_numeric_fields (const const_IWSubstring & buffer,
                              int istart)
{
  const char * s = buffer.data();

  for (int i = istart; i < buffer.length(); i++)
  {
    if ('E' == s[i])
      return 1;
  }

  return 0;
}

template <typename D, typename R>
int
IWDescriptors<D, R>::build (const const_IWSubstring & buffer,
                    int number_descriptors,
                    int & fatal)
{
  assert (NULL == _d);

  _d = new D[number_descriptors];

  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword (_id, i))
  {
    std::cerr << "IWDescriptors::build: cannot extract id\n";
    fatal = 1;
    return 0;
  }

  if (descriptors_may_contain_big_E() && contains_E_in_numeric_fields(buffer, i))
  {
    IWString tmp (buffer);
    tmp.gsub ('E', 'e');
//  std::cerr << "Transformed to '" << tmp << "'\n";
    return _build (tmp, i, number_descriptors, fatal);
  }

  return _build (buffer, i, number_descriptors, fatal);
}

template <typename D, typename R>
int
IWDescriptors<D, R>::_build (const const_IWSubstring & buffer,
                     int i,
                     int number_descriptors,
                     int & fatal)
{
  int col = 0;

  const_IWSubstring token;

  while (buffer.nextword (token, i))
  {
    if (col >= number_descriptors)
    {
      std::cerr << "IWDescriptors::_build: too many tokens " << number_descriptors << '\n';
      return 0;
    }

    if (! token.numeric_value (_d[col]))
    {
      if ('.' == token)  // missing_value
      {
        std::cerr << "IWDescriptors::_build: skipping missing value with '" << _id << "'\n";
        fatal = 0;
        delete [] _d;
        _d = NULL;
        return 0;
      }

      std::cerr << "IWDescriptors::_build:invalid numeric value '" << token << "'\n";
      fatal = 1;
      return 0;
    }

    col++;
  }

  if (col != number_descriptors)
  {
    std::cerr << "IWDescriptors::_build: column count mismatch, got " << col << " expected " << number_descriptors << '\n';
    return 0;
  }

  return 1;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_cartesian_distance (const IWDescriptors & rhs,
                                  int n,
                                  const int * xref) const
{
  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    rc += (_d[i] - rhs._d[j]) * (_d[i] - rhs._d[j]);
  }

  return static_cast<R>(sqrt (rc));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_manhattan_distance (const IWDescriptors & rhs,
                                  int n,
                                  const int * xref) const
{
  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 > d2)
      rc += d1 - d2;
    else 
      rc += d2 - d1;
  }

  return rc;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_canberra_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 != - d2)
      rc += static_cast<R>(fabs ((d1 - d2) / (d1 + d2)));
  }

  return rc;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_taxonomic_distance (const IWDescriptors & rhs,
                                  int n,
                                  const int * xref) const
{
  R sum1 = static_cast<R>(0.0);
  R sum2 = static_cast<R>(0.0);
  R sum3 = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    sum1 += d1 * d1;
    sum2 += d2 * d2;
    sum3 += d1 * d2;
  }

  return static_cast<R>(sqrt (sum1 + sum2 - sum3 - sum3));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_cosine_distance (const IWDescriptors & rhs,
                               int n,
                               const int * xref) const
{
  R sum1 = static_cast<R>(0.0);
  R sum2 = static_cast<R>(0.0);
  R sum3 = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    sum1 += d1 * d1;
    sum2 += d2 * d2;
    sum3 += d1 * d2;
  }

  if (static_cast<R>(0.0) == sum1 || static_cast<R> (0.0) == sum2)
    return static_cast<R>(1.0);

  return static_cast<R>(1.0) - static_cast<R> (sum3 / sqrt (sum1 * sum2));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_minkowski_infinite_distance (const IWDescriptors & rhs,
                               int n,
                               const int * xref) const
{
  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    D zdiff;
    if (d1 > d2)
      zdiff = d1 - d2;
    else
      zdiff = d2 - d1;

    if (zdiff > rc)
      rc = zdiff;
  }

  return rc;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_bray_curtiss_distance (const IWDescriptors & rhs,
                               int n,
                               const int * xref) const
{
  R sum1 = static_cast<R>(0.0);
  R sum2 = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 > d2)
      sum1 += (d1 - d2);
    else
      sum1 += (d2 - d1);

    sum2 += (d1 + d2);
  }

  return static_cast<R>(sum1 / sum2);
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_censored_euclidian_distance (const IWDescriptors & rhs,
                                  int n,
                                  const int * xref) const
{
  D zero_distance = static_cast<D> (0.0);

  R sum = static_cast<R>(0.0);

  int nonzero = 0;
  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    sum += (_d[i] - rhs._d[j]) * (_d[i] - rhs._d[j]);

    if (zero_distance != _d[i] || zero_distance != rhs._d[j])
      nonzero++;
  }

  return static_cast<R>(sqrt (sum / nonzero));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_iw_ratio_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
  D zero_distance = static_cast<D> (0.0);

  R rc = static_cast<R>(0.0);

  int denominator = 0;

#ifdef DEBUG_IW_RATIO
  std::cerr << "Doing computation over " << n << " descriptors\n";
#endif

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

#ifdef DEBUG_IW_RATIO
    std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (d1 > zero_distance && d2 > zero_distance)
    {
      if (d1 > d2)
        rc += d2 / d1;
      else
        rc += d1 / d2;

      denominator++;
    }
    else if (d1 < zero_distance && d2 < zero_distance)
    {
      if (d1 > d2)
        rc += d1 / d2;
      else
        rc += d2 / d1;

      denominator++;
    }
    else if (zero_distance == d1 && zero_distance == d2)
    {
    }
    else
    {
      rc += 0.0;
      denominator++;
    }
  }

  if (0 == denominator)
    return static_cast<R>(0.0);

#ifdef DEBUG_IW_RATIO
  std::cerr << " rc = " << rc << " denominator " << denominator << '\n';
#endif

  return static_cast<R>(1.0) - static_cast<R> (rc / denominator);
}

//#define DEBUG_TANIMOTO

template <typename D, typename R>
R
IWDescriptors<D, R>::_tanimoto_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
#ifdef DEBUG_TANIMOTO
  std::cerr << "Doing Tanimoto computation over " << n << " descriptors\n";
#endif

//R bic = static_cast<R>(0.0);
  D bic = 0.0;

  for (int i = 0; i < n; i++)
  {
    D d1 = _d[i];
    D d2 = rhs._d[i];

#ifdef DEBUG_TANIMOTO
    std::cerr << " i = " << i << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (d1 >= d2)
      bic += d2;
    else
      bic += d1;
  }

  if (static_cast<D>(0.0) == bic)
  {
    if (0.0 == _nset && 0.0 == rhs._nset)
      return 0.5;         // arbitrarily chosen value

    return 1.0;
  }

  double rc = 1.0 - bic / (_nset + rhs._nset - bic);

#ifdef DEBUG_TANIMOTO
  std::cerr << " bic = " << bic << " na " << _nset <<  " nb " << rhs._nset << " rc = " << rc << '\n';
#endif

  if (rc < 1.0e-05)
    return static_cast<R>(0.0);
  else if (rc > 0.9999)
    return 1.0;
  else
    return rc;
}

/*template <typename D, typename R>
R
IWDescriptors<D, R>::_tanimoto_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
#ifdef DEBUG_TANIMOTO
  std::cerr << "Doing Tanimoto computation over " << n << " descriptors\n";
#endif

  D zero_distance = static_cast<D>(0.0);

  R bic = static_cast<R>(0.0);
  R na  = static_cast<R>(0.0);
  R nb  = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

#ifdef DEBUG_TANIMOTO
    std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (d1 < zero_distance && d2 < zero_distance)   // both negative, switch to +ve
    {
      d1 = - d1;
      d2 = - d2;
    }
    else if (d1 >= zero_distance && d2 >= zero_distance)    // both non-negative, great
      ;
    else       // one positive, the other negative, don't process
    {
      std::cerr << "IWDescriptors::_tanimoto_distance:cannot compute mixed +- " << d1 << " and " << d2 << '\n';
      abort();
    }

    if (d1 >= d2)
      bic += d2;
    else
      bic += d1;

    na += d1;
    nb += d2;
  }

  if (static_cast<R>(0.0) == bic)
  {
    if (0 == na && 0 == nb)
      return 0.5;         // arbitrarily chosen value

    return 1.0;
  }

#ifdef DEBUG_TANIMOTO
  std::cerr << " bic = " << bic << " na " << na <<  " nb " << nb << " rc = " << static_cast<R>(1.0) - bic / (na + nb - bic) << '\n';
#endif

  return static_cast<R>(1.0) - bic / (na + nb - bic);
}*/


#ifdef VERSION_THAT_COMPUTES_EVERYTHING
template <typename D, typename R>
R
IWDescriptors<D, R>::_continuous_tanimoto_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
#ifdef DEBUG_TANIMOTO
  std::cerr << "Doing Tanimoto computation over " << n << " descriptors\n";
#endif

  double bic = 0.0;
  double na  = 0.0;
  double nb  = 0.0;
  int nzero = 0;

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

#ifdef DEBUG_TANIMOTO
    std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (static_cast<D>(0.0) == d1)
    {
      if (static_cast<D>(0.0) == d2)
        nzero++;
      else
        nb += d2 * d2;
    }
    else if (static_cast<D>(0.0) == d2)
      na += d1 * d1;
    else
    {
      bic += d1 * d2;
      na += d1 * d1;
      nb += d2 * d2;
    }
  }

  if (nzero == n)
    return static_cast<R>(0.0);

  double tmp = bic / (na + nb - bic);

// The result is now between -0.333 and 1.0. Scale to 0-1

#ifdef DEBUG_TANIMOTO
  std::cerr << " bic = " << bic << " na " << na <<  " nb " << nb << " tmp = " << tmp << '\n';
#endif

  tmp = (3.0 * tmp + 1.0) / 4.0;

  return static_cast<R>(1.0 - tmp);
}
#endif

template <typename D, typename R>
R
IWDescriptors<D, R>::_continuous_tanimoto_distance (const IWDescriptors & rhs,
                                 int n,
                                 const int * xref) const
{
#ifdef DEBUG_TANIMOTO
  std::cerr << "Doing Tanimoto computation over " << n << " descriptors\n";
#endif

  double bic = 0.0;
  double d12 = 0.0;
  double rhs_d12 = 0.0;

  for (int i = 0; i < n; i++)
  {
    double d1 = _d[i];
    double d2 = rhs._d[i];

#ifdef DEBUG_TANIMOTO
    std::cerr << " i = " << i << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

      bic += d1 * d2;
      d12 += d1 * d1;
      rhs_d12 += d2 * d2;
  }

  if (0.0 == bic)
  {
    if (0.0 == _product_of_non_zero_items && 0.0 == rhs._product_of_non_zero_items)
      return static_cast<R>(0.0);
    else
      return static_cast<R>(1.0);
  }

//std::cerr << _product_of_non_zero_items << " and " << rhs._product_of_non_zero_items << ", bic " << bic << '\n';

// March 2011, change definition

  return 1.0 - bic / (d12 + rhs_d12 - bic);

#ifdef THIS_PART_NOT_NEEDED_ANY_MORE
  double tmp = bic / (_product_of_non_zero_items + rhs._product_of_non_zero_items - bic);

// The result is now between -0.333 and 1.0. Scale to 0-1

#ifdef DEBUG_TANIMOTO
  std::cerr << " bic = " << bic << " na " << _product_of_non_zero_items <<  " nb " << rhs._product_of_non_zero_items << " tmp = " << tmp << '\n';
#endif

  tmp = (3.0 * tmp + 1.0) / 4.0;

  if (tmp < 1.0e-05)
  {
    if (tmp > -0.001)
      tmp = 0.0;
    else
      std::cerr << "_continuous_tanimoto_distance:possible out of range - " << _product_of_non_zero_items << " and " << rhs._product_of_non_zero_items << " bic " << bic << ", result " << tmp << " out of range " << (tmp + 1.0) << '\n';
  }
  else if (tmp > 0.99999)
  {
    if (tmp < 1.001)
      tmp = 1.0;
    else
      std::cerr << "_continuous_tanimoto_distance:possible out of range + " << _product_of_non_zero_items << " and " << rhs._product_of_non_zero_items << " bic " << bic << ", result " << tmp << " out of range " << (tmp - 1.0) << '\n';
  }

  return static_cast<R>(1.0 - tmp);
#endif
}

//#define DEBUG_FORBES

template <typename D, typename R>
R
IWDescriptors<D, R>::_continuous_forbes_distance (const IWDescriptors & rhs,
                               int n,
                               const int * xref) const
{
#ifdef DEBUG_FORBES
  std::cerr << "Doing Forbes computation over " << n << " descriptors\n";
#endif

  double a =  0.0;
  double b  = 0.0;
  double c  = 0.0;

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

#ifdef DEBUG_FORBES
//  std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    a += d1 * d2;
    b += d1 * d1;
    c += d2 * d2;
  }

#ifdef DEBUG_FORBES
  std::cerr << " bic = " << a << " b " << b <<  " c " << c << '\n';
#endif

//return static_cast<D, R> (1.0) - bic / (na + nb - bic);
  return n * a / ((a + b) * (a + c));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_forbes_distance (const IWDescriptors & rhs,
                               int n,
                               const int * xref) const
{
#ifdef DEBUG_FORBES
  std::cerr << "Doing Forbes computation over " << n << " descriptors\n";
#endif

  D zero_distance = static_cast<D> (0.0);

  double a =  0.0;
  double b  = 0.0;
  double c  = 0.0;

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 >= zero_distance && d2 >= zero_distance)
      ;
    else if (d1 < zero_distance && d2 < zero_distance)
    {
      d1 = -d1;
      d2 = -d2;
    }
    else
    {
      std::cerr << "IWDescriptors::_forbes_distance: cannot do mixed +- " << d1 << " and " << d2 << '\n';
      abort();
    }

#ifdef DEBUG_FORBES
//  std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (d1 > d2)
      a += d2;
    else
      a += d1;

    b += d1 * d1;
    c += d2 * d2;
  }

#ifdef DEBUG_FORBES
  std::cerr << " bic = " << a << " b " << b <<  " c " << c << '\n';
#endif

//return static_cast<D, R> (1.0) - bic / (na + nb - bic);
  return n * a / ((a + b) * (a + c));
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_russel_rao_distance (const IWDescriptors & rhs,
                                   int n,
                                   const int * xref) const
{
#ifdef DEBUG_RUSSEL_RAO
  std::cerr << "Doing Russel Rao computation over " << n << " descriptors\n";
#endif

  D zero_distance = static_cast<D> (0.0);

  R bic = static_cast<R>(0.0);
  R nbits  = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

#ifdef DEBUG_TANIMOTO
    std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif

    if (d1 < zero_distance && d2 < zero_distance)   // both negative, switch to +ve
    {
      d1 = - d1;
      d2 = - d2;
    }
    else if (d1 >= zero_distance && d2 >= zero_distance)    // both non-negative, great
      ;
    else       // one positive, the other negative, don't process
      continue;

    if (d1 >= d2)
    {
      bic += d2;
      nbits += d1;
    }
    else
    {
      bic += d1;
      nbits += d2;
    }
  }

#ifdef DEBUG_RUSSEL_RAO
  std::cerr << " bic = " << bic << " nbits " << nbits << '\n';
#endif

  return static_cast<R>(1.0) - bic / nbits;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_distance_circle_product (const IWDescriptors & rhs,
                                               int n,
                                               const int * xref) const
{
  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 < d2)
      rc += d1;
    else
      rc += d2;
  }

  return rc;
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_distance_czekanowski (const IWDescriptors & rhs,
                                               int n,
                                               const int * xref) const
{
  D zero_distance = static_cast<D> (0.0);

  R sum1 = static_cast<R>(0.0);
  R sum2 = static_cast<R>(0.0);
  R sum3 = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    if (d1 >= zero_distance && d2 >= zero_distance)
      ;
    else if (d1 <= zero_distance && d2 <= zero_distance)
    {
      d1 = -d1;
      d2 = -d2;
    }
    else
    {
      std::cerr << "IWDescriptors::_distance_czekanowski:cannot mix positive and negative\n";
      abort();
    }

    if (d1 < d2)
      sum1 += d1;
    else
      sum1 += d2;

    sum2 += d1;
    sum3 += d2;
  }

  return static_cast<D> (0.5) - sum1 / (sum2 + sum3);
}

/*
  Note that this really isn't robust if someone is using just a subset
  of the descriptors because the _average and _pearson_denominator will
  not be correct. But we would lose too much efficiency otherwise
*/

template <typename D, typename R>
R
IWDescriptors<D, R>::_distance_pearson (const IWDescriptors & rhs,
                                        int n,
                                        const int * xref) const
{
//#define DEBUG_DISTANCE_PEARSON
#ifdef DEBUG_DISTANCE_PEARSON
  std::cerr << "Computing pearson distance, ave " << _average << " pearson_denominator  " <<  _pearson_denominator << " rhs.ave " << rhs._average << " pd " << rhs._pearson_denominator << '\n';
#endif

  int c = 0;     // number of columns included in computation
  double prod = 0.0;

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    prod += d1 * d2;
    c++;
  }

  double rc = (prod - _average * rhs._average * static_cast<double>(c)) / sqrt(_pearson_denominator * rhs._pearson_denominator);

#ifdef DEBUG_DISTANCE_PEARSON
  std::cerr << "prod " << prod << ", rc = " << rc << '\n';
#endif

  if (rc > 1.0)
  {
    if (rc > 1.001)
      std::cerr << "IWDescriptors::_distance_pearson:problematic distance " << rc << ", set to 1.0\n";
    rc = 1.0;
  }
  else if (rc < -1.0)
  {
    if (rc < -1.001)
      std::cerr << "IWDescriptors::_distance_pearson:problematic distance " << rc << ", set to -1.0\n";
    rc = -1.0;
  }

#ifdef ECHO_ERROR
  if (rc < -1 || rc > 1.0)
  {
    std::cerr << "Out of range " << rc << ", prod " << prod << ", c = " << c << '\n';
    std::cerr << rc - 1.0 << '\n';
    for (int i = 0; i < n; i++)
    {
      std::cerr << " i = " << i << " xref " << xref[i] << " d1 " << _d[i] << " rhs " << rhs._d[xref[i]];
      if (_d[i] != rhs._d[xref[i]])
        std::cerr << " *";
      std::cerr << '\n';
    }
  }

  assert (rc >= -1.0 && rc <= 1.0);
#endif

  return static_cast<R>(1.0) - rc;   // convert to distance
}

template <typename D, typename R>
R
IWDescriptors<D, R>::_dot_product_distance (const IWDescriptors & rhs,
                                   int n,
                                   const int * xref) const
{
#ifdef DEBUG_DOT_PRODUCT
  std::cerr << "Doing Dot product computation over " << n << " descriptors\n";
#endif

  R rc = static_cast<R>(0.0);

  for (int i = 0; i < n; i++)
  {
    int j = xref[i];
    if (j < 0)
      continue;

    D d1 = _d[i];
    D d2 = rhs._d[j];

    rc += d1 * d2;

#ifdef DEBUG_DOT_PRODUCT
    std::cerr << " i = " << i << " j = " << j << " lhs " << d1 << " rhs " << d2 << '\n';
#endif
  }

  if (rc >= -1.0 && rc <= 1.0)   // great
    return 1.0 - rc;

  if (rc < -1.0)
  {
    if (fabs(rc + 1) < 1.0e-05)
      return 2.0;

    std::cerr << "Invalid dot product " << rc << '\n';
    abort();
  }

  if (rc > 1.0)
  {
    if (fabs(rc - 1.0) < 1.0e-05)
      return 0.0;

    std::cerr << "Invalid dot product " << rc << '\n';
    abort ();
  }

  return 3.0;
}


template <typename D, typename R>
R
IWDescriptors<D, R>::distance (const IWDescriptors & rhs,
                               int distance_computation,
                               int n,
                               const int * xref) const
{
//std::cerr << "Distance metric is " << distance_computation << '\n';

  switch (distance_computation)
  {
    case DISTANCE_CARTESIAN:
      return _cartesian_distance (rhs, n, xref);

    case DISTANCE_MANHATTAN:
      return _manhattan_distance (rhs, n, xref);

    case DISTANCE_TAXONOMIC:
      return _taxonomic_distance (rhs, n, xref);

    case DISTANCE_CANBERRA:
      return _canberra_distance (rhs, n, xref);

    case DISTANCE_COSINE:
      return _cosine_distance (rhs, n, xref);

    case DISTANCE_MINKINF:
      return _minkowski_infinite_distance (rhs, n, xref);

    case DISTANCE_BRAY_CURTIS:
      return _bray_curtiss_distance (rhs, n, xref);

    case DISTANCE_CENSORED_EUCLIDIAN:
      return _censored_euclidian_distance (rhs, n, xref);

    case DISTANCE_IW_RATIO:
      return _iw_ratio_distance (rhs, n, xref);

    case DISTANCE_TANIMOTO:
      return _tanimoto_distance (rhs, n, xref);

    case DISTANCE_RUSSEL_RAO:
      return _russel_rao_distance (rhs, n, xref);

    case DISTANCE_CONTINUOUS_TANIMOTO:
      return _continuous_tanimoto_distance (rhs, n, xref);

    case DISTANCE_FORBES:
      return _forbes_distance (rhs, n, xref);

    case DISTANCE_DOT_PRODUCT:
      return _dot_product_distance (rhs, n, xref);

    case DISTANCE_CONTINUOUS_FORBES:
      return _continuous_forbes_distance (rhs, n, xref);

    case DISTANCE_CIRCLE_PRODUCT:
        return _distance_circle_product (rhs, n, xref);

    case DISTANCE_CZEKANOWSKI:
        return _distance_czekanowski (rhs, n, xref);

    case DISTANCE_PEARSON:
        return _distance_pearson (rhs, n, xref);

    default:
      std::cerr << "What kind of distance am I supposed to be computing " << distance_computation << '\n';
      abort();
  };

  assert (NULL == "IWDescriptors::distance: should not come to here\n");

  return static_cast<R>(0.0);
}

template <typename D, typename R>
void 
IWDescriptors<D, R>::compute_pearson_data(int n)
{
  double sum = 0.0;
  double sum2 = 0.0;

  for (int i = 0; i < n; i++)
  {
    sum += _d[i];
    sum2 += (_d[i] * _d[i]);
  }

  _average = sum / static_cast<double>(n);

  _pearson_denominator = (sum2 - sum * sum / static_cast<double>(n)); // / static_cast<double>(n);

//std::cerr << "IWDescriptors::ave " << _average << " pearson denominator " << _pearson_denominator << '\n';

  return;
}

template <typename D, typename R>
void
IWDescriptors<D, R>::compute_product_of_non_zero_items (int n)
{
  _product_of_non_zero_items = 0.0;

  D zero = static_cast<D>(0.0);

  for (int i = 0; i < n; i++)
  {
    if (zero != _d[i])
      _product_of_non_zero_items += (_d[i] * _d[i]);
  }

  return;
}

template <typename D, typename R>
void
IWDescriptors<D, R>::compute_nset (int n)
{
  _nset = 0.0;

  D zero = static_cast<D>(0.0);

  for (int i = 0; i < n; i++)
  {
    if (zero != _d[i])
      _nset += _d[i];
  }

  return;
}


template <typename D, typename R> template <typename C>
R
IWDescriptors<D, R>::distance (const IWDescriptors<D, R> & rhs,
                               const C & c,
                               int n,
                               const int * xref) const
{
  R rc = static_cast<R>(0.0);
  int valid_pairs = 0;

  for (int i = 0; i < n; i++)
  {
    if (xref[i] < 0)
      continue;

    rc += c(_d[i], rhs._d[xref[i]]);
    valid_pairs++;
  }

  return c.final(rc, n, valid_pairs);
}

#endif

#endif

