/*
  We need to generate the header of precomputed bit counts
*/

#include <stdlib.h>
#include <limits.h>

#include <iostream>
#include <iomanip>

using namespace std;

static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};
static const unsigned short one_bit_16[16] = {32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1};

static int
count_bits (unsigned short s)
{
  int rc = 0;

  for (int i = 0; i < 16; i++)
  {
    if (one_bit_16[i] & s)
      rc++;
  }

  return rc;
}

int
main ()
{
  cout << "static const unsigned char sixteen_bit_count[65536] = {\n";

  for (int i = 0; i <= USHRT_MAX; i++)
  {
    unsigned short s = static_cast<unsigned short> (i);

    int b = count_bits (s);

    cout << "   " << b << ", /* " << s << ' ' << hex << s << dec << " */\n";
  }

  cout << "};\n";
  return 0;
}
