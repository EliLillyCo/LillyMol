#include <iostream>
#include <memory>

#include "Foundational/iwstring/iwstring.h"
#include "iwbits.h"
#include "dy_fingerprint.h"

using std::ostream;
using std::cerr;

static void
convert_to_byte_range(const int * c,
                      int n,
                      unsigned char * b)
{
  for (int i = 0; i < n; i++)
  {
    if (c[i] < 0)
    {
      cerr << "write_fixed_size_counted_fingerprint:invalid value i = " << i << " v = " << c[i] << " set to zero\n";
      b[i] = static_cast<unsigned char> (0);
    }
    else if (c[i] > 255)
    {
      cerr << "write_fixed_size_counted_fingerprint:value out of range i = " << i << " v = " << c[i] << ", truncated\n";
      b[i] = static_cast<unsigned char> (0);
    }
    else
      b[i] = static_cast<unsigned char> (c[i]);
  }

  return;
}

int
write_fixed_size_counted_fingerprint(const int * c,
                                     int n,
                                     int nset,
                                     int bits_in_original_fingerprint,
                                     ostream & output)
{
  if (0 == n || nullptr == c)
  {
    cerr << "write_fixed_size_counted_fingerprint:invalid input\n";
    return 0;
  }

  unsigned char * b = new unsigned char[n]; std::unique_ptr<unsigned char[]> free_b (b);

  convert_to_byte_range (c, n, b);

  int bytes_allocated;

  char * bytes = du_bin2ascii (&bytes_allocated, n, reinterpret_cast<char *> (b));

  if (nullptr == bytes)
  {
    cerr << "write_fixed_size_counted_fingerprint:cannot allocate fingerprint\n";
    return 0;
  }

  output.write (bytes, bytes_allocated);

  delete bytes;

  if (nset >= 0 || bits_in_original_fingerprint >= 0)    // should check that nset is positive
    output << ';' << nset << ';' << n;

  return output.good ();
}
int
append_fixed_size_counted_fingerprint (const int * c,
                                      int n,
                                      int nset,
                                      int bits_in_original_fingerprint,
                                      IWString & output)
{
  if (0 == n || nullptr == c)
  {
    cerr << "write_fixed_size_counted_fingerprint:invalid input\n";
    return 0;
  }

  unsigned char * b = new unsigned char[n]; std::unique_ptr<unsigned char[]> free_b (b);

  convert_to_byte_range (c, n, b);

  int bytes_allocated;

  char * bytes = du_bin2ascii (&bytes_allocated, n, reinterpret_cast<char *> (b));

  if (nullptr == bytes)
  {
    cerr << "write_fixed_size_counted_fingerprint:cannot allocate fingerprint\n";
    return 0;
  }

  output.strncat (bytes, bytes_allocated);

  delete bytes;

  if (nset >= 0 || bits_in_original_fingerprint >= 0)    // should check that nset is positive
    output << ';' << nset << ';' << n;

  return 1;
}
