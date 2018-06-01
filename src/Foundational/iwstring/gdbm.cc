#include <stdlib.h>

#include "gdbm.h"

#include "iwstring.h"


int
operator == (const datum & lhs, const const_IWSubstring & rhs)
{
  if (lhs.dsize != rhs.nchars ())
    return 0;

  return 0 == rhs.strncmp (lhs.dptr, lhs.dsize);
}

int
operator != (const datum & lhs, const const_IWSubstring & rhs)
{
  if (lhs.dsize != rhs.nchars ())
    return 1;

  return rhs.strncmp (lhs.dptr, lhs.dsize);
}

int
operator == (const const_IWSubstring & lhs, const datum & rhs)
{
  if (rhs.dsize != lhs.nchars ())
    return 0;

  return 0 == lhs.strncmp (rhs.dptr, rhs.dsize);
}

int
operator != (const const_IWSubstring & lhs, const datum & rhs)
{
  if (rhs.dsize != lhs.nchars ())
    return 1;

  return lhs.strncmp (rhs.dptr, rhs.dsize);
}

int
operator == (const datum & lhs, const IWString & rhs)
{
  if (lhs.dsize != rhs.nchars ())
    return 0;

  return 0 == rhs.strncmp (lhs.dptr, lhs.dsize);
}

int
operator != (const datum & lhs, const IWString & rhs)
{
  if (lhs.dsize != rhs.nchars ())
    return 1;

  return rhs.strncmp (lhs.dptr, lhs.dsize);
}

int
operator == (const IWString & lhs, const datum & rhs)
{
  if (rhs.dsize != lhs.nchars ())
    return 0;

  return 0 == lhs.strncmp (rhs.dptr, rhs.dsize);
}

int
operator != (const IWString & lhs, const datum & rhs)
{
  if (rhs.dsize != lhs.nchars ())
    return 1;

  return lhs.strncmp (rhs.dptr, rhs.dsize);
}

ostream &
operator << (ostream & os, const datum & rhs)
{
  return os.write (rhs.dptr, rhs.dsize);
}
