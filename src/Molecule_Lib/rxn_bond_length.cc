#include <stdlib.h>

#include "iwreaction.h"

Reaction_Bond_Length::Reaction_Bond_Length ()
{
  _desired_length = 0.0;

  return;
}

void
Reaction_Bond_Length::all_atoms_in_scaffold ()
{
  _atom[0].set_in_scaffold ();
  _atom[1].set_in_scaffold ();

  return;
}

int
Reaction_Bond_Length::construct_from_msi_attribute (const msi_attribute * att)
{
  const_IWSubstring m;
  att->value (m);

  if (3 != m.nwords ())
  {
    cerr << "Reaction_Bond_Length::construct_from_msi_attribute: attribute must have 5 tokens\n";
    cerr << m << endl;
    return 0;
  }

  for (int i = 0; i < 2; i++)
  {
    const_IWSubstring token;
    (void) m.word (i, token);

    if (! _atom[i].construct (token))
    {
      cerr << "Reaction_Bond_Length::construct_from_msi_attribute: cannot parse attribute " << i << " '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }
  }

// The last token is the distance

  const_IWSubstring token = m.word (2);

  if (! token.numeric_value (_desired_length))
  {
    cerr << "Reaction_Bond_Length::construct_from_msi_attribute: invalid distance '" << token << "'\n";
    return 0;
  }

//write_msi(cerr, " ", "built");

  return 1;
}

int
Reaction_Bond_Length::write_msi (std::ostream & os,
                                 const const_IWSubstring & ind,
                                 const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"";
  for (int i = 0; i < 2; i++)
  {
    if (i > 0)
      os << ' ';

    os << _atom[i];
  }

  os << ' ' << _desired_length;

  os << "\")\n";

  return os.good ();
}

int
Reaction_Bond_Length::adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < 2; i++)
  {
    _atom[i].adjust_matched_atoms_in_component (xref);
  }

  return 1;
}

/*int
Reaction_Bond_Length::process (Molecule & m, 
                               const Set_of_Atoms * scaffold_embedding,
                               const int * atoms_in_growing_molecule,
                               const Molecule_and_Embedding ** reagent) const
{
  atom_number_t a[2];

  for (int i = 0; i < 2; i++)
  {
    if (! determine_atom_number (*scaffold_embedding, _atom[i], atoms_in_growing_molecule, reagent, "Reaction_Dihedral_Angle:process:", a[i]))
      return 0;
  }

  cerr << "Setting bond length between atoms " << a[0] << " to " << a[1] << " to " << _desired_length << " Angstroms\n";
  int rc = m.set_bond_length (a[0], a[1], _desired_length);

  if (rc)
    return rc;

  cerr << "Reaction_Bond_Length::process: could not set bond length to " << _desired_length << " Angstroms\n";
  cerr << "Atoms " << a[0] << " and " << a[1] << endl;
  return 0;
}*/
