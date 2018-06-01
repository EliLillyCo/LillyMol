#include <stdlib.h"
#include <iostream>
#include <math.h>
#include <limits.h>

static int
is_prime (long p)
{
  int isqrt = int (sqrt (p)) + 1;
  for (int i = 3; i <= isqrt; i++)
  {
    if (0 == p % i)
      return 0;
  }

  return 1;
}

int
main ()
{
  int primes_needed = 1000;
  int primes_found = 0;

  cout << "static const int nprimes = " << primes_needed << ";\n";

  cout << "static const int primes [] = {\n";
  for (long int i = 3; i < LONG_MAX; i += 2)
  {
    if (is_prime (i))
    {
      cout << "    " << i;
      primes_found++;
      if (primes_found >= primes_needed)
        break;

      cout << ",\n";
    }
  }

  if (primes_found < primes_needed)
    cerr << "Warning, only " << primes_found << " primes found\b";

  cout << "\n};\n";

  return 0;
}

