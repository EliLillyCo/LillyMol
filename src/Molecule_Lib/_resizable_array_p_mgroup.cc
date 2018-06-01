#include <iostream>
using namespace std;
#include <iostream>
using namespace std;


#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "mgroup.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)
template class resizable_array_p<Molecule_Subset>;
template class resizable_array_base<Molecule_Subset *>;
#else

static void
unused ()
{
  resizable_array_p<Molecule_Subset> m;
  resizable_array_base<Molecule_Subset *> mb;

  mb[0];
  mb.add (NULL);

  m.transfer_in (m);
}

#endif
