#include <stdlib.h>
#include <iomanip.h>

int
main ()
{

  for (int i = 0; i < 260; i++)
  {
    unsigned char * p = (unsigned char *) &i;
    cout << " i = " << dec << i << " :";
    for (int j = 0; j < sizeof (int); j++)
    {
      cout << ' ' << hex << int (p[j]);
    }
    cout << endl;
  }
}
