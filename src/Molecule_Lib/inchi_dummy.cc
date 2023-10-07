#include <stdlib.h>
#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"
#include "rwmolecule.h"

using std::cerr;
using std::endl;

int
inchi_to_inchi_key(const char * inchi,
                   IWString & key)
{
  cerr << "inchi_to_inchi_key not implemented, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
  return 0;
}

int
Molecule::build_from_inchi(const const_IWSubstring & inchi)
{
  cerr << "Molecule::build_from_inchi:not implemented, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
  return 0;
}

int
Molecule::read_molecule_inchi_ds(iwstring_data_source & input)
{
  resize(0);   // clear out anything already present

  const_IWSubstring buffer;

  EXTRA_STRING_RECORD(input, buffer, "read mol inchi");

  return build_from_inchi(buffer);
}

//#define DEBUG_MAKE_INCHI_STRING

int
Molecule::InChI(IWString & inchi_string)      // maybe this could be a const method, depends on what inchi needs...
{
  cerr << "Molecule::InChI:not implemented, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";

  return 0;
}

int
Molecule::write_molecule_inchi(std::ostream & output)
{
  IWString tmp;

  InChI(tmp);   // generate string version of InChI

  output << tmp << ' ' << _molecule_name << '\n';

  return 1;
}
