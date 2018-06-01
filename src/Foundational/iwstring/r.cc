#include "stdlib.h"
#include "iwstring.h"

static int
qis_int (const char * s, int nchars, int & result)
{
  int isign = 1;
  if ('-' == *s)
  {
    isign = -1;
    s++;
    nchars--;
  }

  int tmp = 0;
  while (nchars)
  {
    if (*s > '9')
      return 0;
    if (*s < '0')
      return 0;

    tmp = 10 * tmp + *s - '0';
    s++;
    nchars--;
  }

  if (nchars)
    return 0;

  if (-1 == isign)
    tmp = - tmp;

  result = tmp;

  return 1;
}

static int
is_int (const char * s, int nchars, int & result)
{
  int i;
  if ('-' == *s)
    i = 1;
  else
    i = 0;

  int tmp = 0;
  for ( ; i < nchars; i++)
  {
    int d = s[i] - '0';
    if (d > 9 || d < 0)
      return 0;

    tmp = 10 * tmp + d;
  }

  if ('-' == *s)
    tmp = - tmp;

  result = tmp;

  return 1;
}

int
main (int argc, char ** argv)
{
  cerr << "ARGC = " << argc << endl;

  IWString foo = "12345";

  int maxj;

  for (int i = 0; i < 100000; i++)
  {
    int j;
    if (1 == argc)
      (void) foo.is_int (j);
    else
      (void) ::is_int (foo.rawchars (), foo.nchars (), j);
    if (j > maxj)
      maxj = j;
  }

  cerr << "Max j is " << maxj << endl;

  return 0;
}
