#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "iwstring.h"

int
main()
{
  IWString_and_File_Descriptor output(1);

  int i = 5;

  output << "Hello " << std::setw(4) << i << " world\n";

  return 0;
}
