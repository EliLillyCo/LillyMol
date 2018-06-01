/*
  We need to create some static arrays for inclusion.
  Given two byte arguments, what are the
    Number of bits in common
    Number of bits on in (1), and not in (2)
    Number of bits on in (2), and not in (1)
*/

#include <ostream>
#include <iomanip>

/*
*/

static int maxval = 256;

static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};

static const unsigned short one_bit_16[16] = {
  32768, 16384, 8192, 4096, 2048, 1024, 512, 256,
  128, 64, 32, 16, 8, 4, 2, 1};

static int
bits_in_common (ostream & os)
{
  int n = maxval * (maxval - 1) / 2;

  os << "static int _bits_in_common[" << n << "] = {\n";

  int count = 0;

  for (int i = 0; i < maxval; i++)
  {
    for (int j = i + 1; j < maxval; j++)
    {
      unsigned char andij = (unsigned char) (i & j);

      int bic = 0;
      for (int k = 0; k < 8; k++)
      {
        unsigned char tmp = one_bit_8[k] & andij;
        if (tmp)
          bic++;
      }

      os << "    " << bic;
      if (i == maxval - 1 && j == maxval - 1)
        ;
      else
        os << ',';

      os << "      /" << "*   " << count++ << " : " << i << ',' << j << " *" << "/\n";
    }
  }

  os << "\n};\n\n";

  return os.good ();
}

static int
bits_in_a (ostream & os)
{
  os << "static int _bits_in_a[] = {";

  for (int i = 0; i < maxval; i++)
  {
    int iset[8];
    for (int j = 0; j < 8; j++)
      iset[j] = one_bit_8[j] & i;

    for (int j = i + 1; j < maxval; j++)
    {
      unsigned char cj = (unsigned char) j;

      int bits_in_a = 0;
      for (int k = 0; k < 8; k++)
      {
        if (0 == iset[k])   // bit K NOT set in I
          continue;

        if (0 == (one_bit_8[k] & cj))
          bits_in_a++;
      }

      if (j > 1)
        os << ',';

      os << "\n    " << bits_in_a;
    }
  }

  os << "\n};\n\n";

  return os.good ();
}

static int
bits_in_b (ostream & os)
{
  os << "static int _bits_in_b[] = {";

  for (int i = 0; i < maxval; i++)
  {
    int iset[8];
    for (int j = 0; j < 8; j++)
      iset[j] = one_bit_8[j] & i;

    for (int j = i + 1; j < maxval; j++)
    {
      unsigned char cj = (unsigned char) j;

      int bits_in_b = 0;
      for (int k = 0; k < 8; k++)
      {
        if (iset[k])   // bit K set in I
          continue;

        if (one_bit_8[k] & cj)
          bits_in_b++;
      }

      if (j > 1)
        os << ',';

      os << "\n    " << bits_in_b;
    }
  }

  os << "\n};\n\n";

  return os.good ();
}

int
bit_count (ostream & os)
{
  os << "static const unsigned char sixteen_bit_count[65536] = {\n";

  for (unsigned short i = 0; i <= 65535; i++)
  {
    int nb = 0;
    for (int j = 0; j < 16; j++)
    {
      if (one_bit_16[j] & i)
        nb++;
    }

    os << "    " << nb;
    if (i < 65535)
      os << ',';
    else
      os << ' ';

    os << "   /" << "* " << i << ' ' << hex << i << dec << " *" << "/\n";

    if (65535 == i)
      break;
  }

  os << "};\n\n";

  return os.good ();
}

int
main (int argc, char ** argv)
{
//bits_in_common (cout);
//bits_in_a (cout);
//bits_in_b (cout);

  bit_count (cout);

  return 0;
}
