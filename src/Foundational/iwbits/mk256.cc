#include <stdlib.h>
#include <iostream.h>

static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};

int
main ()
{
  cout << "static char * byte_form_space [] = {\n";

  for (unsigned char i = 0; i < 256; i++)
  {
    cout << "         \" ";
    for (int j = 0; j < 8; j++)
    {
      if (one_bit_8[j] & i)
        cout << '1';
      else
        cout << '0';

      cout << ' ';
    }

    cout << "\",\n";

    if (255 == i)
      break;
  }

  cout << "   };\n";
}
