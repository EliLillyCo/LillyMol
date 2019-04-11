#include <stdlib.h>

/*
  Very bad idea, a globally accessible variable which determines the min of two numbers in 
  the range 0-255
  The reason for making it global is so we have just one copy of it in the programme
*/

unsigned char bic_table[256*256];

static int
initialise_bic ()
{
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < 256; j++)
    {
      if (i < j)
        bic_table[i * 256 + j] = i;
      else
        bic_table[i * 256 + j] = j;
    }
  }

  return 1;
}

class Just_to_Get_Something_Initialised
{
  private:
  public:
    Just_to_Get_Something_Initialised ();
};

Just_to_Get_Something_Initialised::Just_to_Get_Something_Initialised ()
{
  initialise_bic ();

  return;
}

static Just_to_Get_Something_Initialised notused;


