#ifndef IW_MIST_H
#define IW_MIST_H 1

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

using std::cerr;
using std::endl;

#ifdef IWCYGWIN
#define IWDIRECTORY_SEPARATOR '\\'
#else
#define IWDIRECTORY_SEPARATOR '/'
#endif

#include "assert.h"

#include "iwaray.h"
#include "iwcrex.h"

/*
  various useful functions, some templates, others just functions
*/

template <typename T>
int
locate_item_in_array (T needle, int n, const T * haystack, int istart = 0)
{
  for (int i = istart; i < n; i++)
  {
    if (needle == haystack[i])
      return i;
  }

  return -1;
}

template <typename T>
int
count_occurrences_of_item_in_array (T needle, int n, const T * haystack)
{
  int rc = 0;
  for (int i = 0; i < n; i++)
  {
    if (needle == haystack[i])
      rc++;
  }

  return rc;
}

template <typename T>
int
count_non_zero_occurrences_in_array (const T * haystack, int n)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    if (static_cast<T> (0) != haystack[i])
      rc++;
  }

  return rc;
}

template <typename T>
T *
array_of (int n, T initialiser)
{
  assert (n > 0);

  T * rc = new T[n];

  for (int i = 0; i < n; i++)
  {
    rc[i] = initialiser;
  }

  return rc;
}

template <typename T>
int
index_of_smallest (int n, const T * items)
{
  assert (NULL != items);

  T minval = items[0];
  int rc = 0;

  for (int i = 1; i < n; i++)
  {
    if (items[i] < minval)
    {
      minval = items[i];
      rc = i;
    }
  }

  return rc;
}

template <typename T>
int
index_of_largest (int n, const T * items)
{
  assert (NULL != items);

  T maxval = items[0];
  int rc = 0;

  for (int i = 1; i < n; i++)
  {
    if (items[i] > maxval)
    {
      maxval = items[i];
      rc = i;
    }
  }

  return rc;
}

template <typename T>
void
set_vector (T * values, int n, T new_value)
{
  for (int i = 0; i < n; i++)
  {
    values[i] = new_value;
  }

  return;
}

template <typename T>
T
sum_vector (const T * values, int n)
{
  T rc = T (0);

  for (int i = 0; i < n; i++)
  {
    rc += values[i];
  }

  return rc;
}

template <typename T>
T
iwmax_of_array (const T * values, int n)
{
  assert (n > 0);

  T rc = values[0];
  for (int i = 1; i < n; i++)
  {
    if (values[i] > rc)
      rc = values[i];
  }

  return rc;
}

template <typename T>
T
iwmin_of_array (const T * values, int n)
{
  assert (n > 0);

  T rc = values[0];
  for (int i = 1; i < n; i++)
  {
    if (values[i] < rc)
      rc = values[i];
  }

  return rc;
}

template <typename T>
void
remove_item_from_array (T * values, int items_in_array, int item_to_remove)
{
  assert (item_to_remove >= 0 && item_to_remove < items_in_array - 1);

  int istop = items_in_array - 1;
  for (int i = item_to_remove; i < istop; i++)
  {
    values[i] = values[i + 1];
  }

  return;
}

template <typename T>
void
copy_vector (T * zto, const T * zfrom, int n)
{
  for (int i = 0; i < n; i++)
  {
    zto[i] = zfrom[i];
  }

  return;
}

/*
  Swap two items using a temporary 
*/

template <typename T>
inline void
iwswap (T & v1, T & v2)
{
  T tmp = v2;
  v2 = v1;
  v1 = tmp;

  return;
}

// max of two numbers
template <typename T>
T
max2 (const T &n1, const T & n2 )
{
  if ( n1 > n2 )
    return n1;
  else
    return n2;
}

// min of two numbers
template <typename T>
T
min2 ( T &n1, T & n2 )
{
  if( n1 < n2 ) return n1;
  else          return n2;
}

extern int * new_int (int, int = 0);
extern float * new_float (int, float = static_cast<float> (0.0));

extern unsigned int * new_unsigned_int (int, unsigned int = 0);

extern double * new_double (int, double = 0.0);

/*
  When creating descriptor files, we often need to write a molecule name
  without any intervening spaces
*/

class IWString;

extern std::ostream & write_space_suppressed_string (const IWString & mname, std::ostream & os, char = '_');
extern std::ostream & write_first_token_of_string (const IWString & mname, std::ostream & os);
extern int append_first_token_of_name(const IWString & , IWString &);

/*
  Many times I want to write out the contents of an array. The character string is
  an optional name of the array
*/

template <typename T> int iw_write_array (const T *, int n, const char *, std::ostream &);

#ifdef IW_WRITE_ARRAY_IMPLEMENTATION

#include <string.h>

template <typename T> int
iw_write_array (const T * a,
                int n,
                const char * array_name,
                std::ostream & output)
{
  for (int i = 0; i < n; i++)
  {
    if (NULL == array_name || 0 == strlen (array_name))
      output << " i = " << i;
    else
      output << array_name << '[' << i << ']';

    output << " = " << a[i] << endl;
  }

  return output.good ();
}

#endif

// Some handy file checking things

#include <sys/types.h>
#include <sys/stat.h>

extern int dash_x (const char * fname);
extern int dash_d (const char * fname);

extern off_t dash_s (const char * fname);
extern int   dash_f (const char * fname);

extern void rick_higgs_byte_swap (int, unsigned int *);
extern void rick_higgs_byte_swap (unsigned int &);

extern int iw_little_endian ();
extern void htons_unsigned_short (void *, unsigned int);
extern void ntohs_unsigned_short (void *, unsigned int);
extern void htonl_unsigned_long  (void *, unsigned int);
extern void ntohl_unsigned_long (void *, unsigned int);

/*
  Various things for my implementation of uuencode
*/

extern int IWuuencode_append (const void * p,
                   int nchars,
                   IWString & destination);

/*
  If you have an encoded form, you need to know how many bytes are needed for decoding
*/

extern int IWuudecode_bytes_needed (int);  // length of string holding encoded form

class IWString;
class const_IWSubstring;

extern int IWuudecode (const unsigned char * encoded, int nchars, unsigned char * destination);
extern int IWuudecode (const IWString & encoded, unsigned char * destination);
extern int IWuudecode (const const_IWSubstring & encoded, unsigned char * destination);

class Int_Comparator_Larger
{
  private:
  public:
    int operator() (int i1, int i2) const;
};

class Int_Comparator_Smaller
{
  private:
  public:
    int operator() (int i1, int i2) const;
};


//  Approximately divide the file into offsets for chunked processing

extern int find_dividing_points (const char * fname, int n,
                                 resizable_array<off_t> & o,
                                 IW_Regular_Expression & rx);
#endif
