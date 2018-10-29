#include <stdlib.h>

#include "reaction_duplicate.h"

Reaction_Duplicate::Reaction_Duplicate ()
{
   _molecules_submitted = 0;
   _unique_molecules_encountered = 0;
   _detected_via_standard_smiles = 0;

   return;
}

int
Reaction_Duplicate::is_duplicate(Molecule & m)
{
  _molecules_submitted++;

  const IWString & smi = m.smiles ();

  if (_smiles.contains (smi))
  {
    _detected_via_standard_smiles++;
    return 1;
  }

  _smiles.insert (smi);

  const IWString & usmi = m.unique_smiles ();

  if (_unique_smiles.contains (usmi))
  {
    m.invalidate_smiles ();
    return 1;
  }

  _unique_smiles.insert (usmi);
  m.invalidate_smiles ();

  _unique_molecules_encountered++;

  return 0;
}

int
Reaction_Duplicate::report (std::ostream & output) const
{
  output << "Reaction_Duplicate::report:encountered " << _molecules_submitted << " molecules\n";
  output << _unique_molecules_encountered << " unique structures, " << _detected_via_standard_smiles << " detected via regular smiles\n";

  return output.good ();
}
