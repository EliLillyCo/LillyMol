#include <stdlib.h>
#include <assert.h>

int
getdouble(const char * string, double * value)
{
  assert (nullptr != string);
  assert (nullptr != value);

  char *c;
  double tmp = strtod(string, &c);

// If c is anything but '\0', then we have a problem.

  if (*c)
    return 0;

  *value = tmp;
  return 1;
}

