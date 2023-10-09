#include <iostream>

#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"

using std::cerr;

static IWString smiles_tag("$SMI<");

int
FP_and_Smiles::construct_from_tdt(IW_TDT & tdt, int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal))
    return 0;

  if (smiles_tag.length())
  {
    if (! tdt.dataitem_value(smiles_tag, _smiles))
    {
      cerr << "NN_Object::construct_from_tdt: no smiles '" << smiles_tag << "'\n";
      return 0;
    }
  }

  return 1;
}
