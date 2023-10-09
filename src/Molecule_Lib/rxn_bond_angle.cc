#include <stdlib.h>
#include <iostream>

#include "iwreaction.h"

using std::cerr;
using std::endl;

Reaction_Bond_Angle::Reaction_Bond_Angle()
{
  _desired_angle = static_cast<angle_t>(0.0);

  return;
}

void
Reaction_Bond_Angle::all_atoms_in_scaffold()
{
  for (int i = 0; i < 3; i++)
  {
    _atom[i].set_in_scaffold();
  }

  return;
}

int
Reaction_Bond_Angle::construct_from_msi_attribute(const msi_attribute * att)
{
  const_IWSubstring m;
  att->value(m);

  if (4 != m.nwords())
  {
    cerr << "Reaction_Bond_Angle::construct_from_msi_attribute: attribute must have 4 tokens\n";
    cerr << m << endl;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  for (int j = 0; j < 3; j++)
  {
    (void) m.nextword(token, i);

    if (! _atom[j].construct(token))
    {
      cerr << "Reaction_Bond_Angle::construct_from_msi_attribute: cannot parse attribute " << j << " '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }
  }

// The last token is the angle

  (void) m.nextword(token, i);

  if (! token.numeric_value(_desired_angle))
  {
    cerr << "Reaction_Bond_Angle::construct_from_msi_attribute: invalid angle '" << token << "'\n";
    return 0;
  }

  _desired_angle = static_cast<angle_t>(_desired_angle * DEG2RAD);

  return 1;
}

int
Reaction_Bond_Angle::write_msi(std::ostream & os,
                               const const_IWSubstring & ind,
                               const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"";
  for (int i = 0; i < 3; i++)
  {
    if (i > 0)
      os << ' ';

    os << _atom[i];
  }

  os << ' ' << RAD2DEG * _desired_angle;

  os << "\")\n";

  return os.good();
}

int
Reaction_Bond_Angle::adjust_matched_atoms_in_component(const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < 3; i++)
  {
    _atom[i].adjust_matched_atoms_in_component(xref);
  }

  return 1;
}

//#define DEBUG_PROCESS_BOND_ANGLE

int
Reaction_Bond_Angle::process(Molecule & m,
                             const Set_of_Atoms * scaffold_embedding,
                             const Enumeration_Temporaries & etmp) const
{
#ifdef DEBUG_PROCESS_BOND_ANGLE
  cerr << "Processing bond angle " << _desired_angle << endl;
  if (nullptr != scaffold_embedding)
    cerr << "Scaffold embedding " << (*scaffold_embedding) << endl;
  cerr << m.smiles() << endl;
#endif

  atom_number_t a[3];

  for (int i = 0; i < 3; i++)
  {
    if (! determine_atom_number(*scaffold_embedding, _atom[i], etmp, "Reaction_Bond_Angle:process:", a[i]))
      return 0;

#ifdef DEBUG_PROCESS_BOND_ANGLE
    cerr << _atom[i] << " is atom " << a[i] << endl;
#endif
  }

// should do bump checking

#ifdef DEBUG_PROCESS_BOND_ANGLE
  cerr << "Atoms " << a[0] << " '" << m.smarts_equivalent_for_atom(a[0]) << "', " << a[1] << " '" << m.smarts_equivalent_for_atom(a[1]) << "', " << a[2] << " '" << m.smarts_equivalent_for_atom(a[2]) << "', angle " << (_desired_angle * RAD2DEG) << "\n";
#endif

  if (! m.are_bonded(a[0], a[1]))
  {
    cerr << "Reaction_Bond_Angle::process::atoms " << a[0] << " and " << a[1] << " not bonded (01)\n";
    return 0;
  }

  if (! m.are_bonded(a[1], a[2]))
  {
    cerr << "Reaction_Bond_Angle::process::atoms " << a[1] << " and " << a[2] << " not bonded (12)\n";
    return 0;
  }

  if (m.set_bond_angle(a[0], a[1], a[2], _desired_angle))
  {
#ifdef DEBUG_PROCESS_BOND_ANGLE
    cerr << "After setting, angle is " << RAD2DEG * m.bond_angle(a[0], a[1], a[2]) << endl;
#endif
    return 1;
  }

  cerr << "Reaction_Bond_Angle::process:set_bond_angle failed '" << m.name() << "'\n";
  cerr << "Atoms " << a[0] << ' ' << a[1] << ' ' << a[2] << endl;
  write_isotopically_labelled_smiles(m, cerr);

  return 0;
}
