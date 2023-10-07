#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "iwreaction.h"

using std::cerr;
using std::endl;

/*
  Jun 2004. Turn off bump checking. The problem is that we may have
  multiple, adjacent torsions, and the dump check won't be valid
  until all have been adjusted.
*/

static distance_t bump_tolerance = 0.0;

Reaction_Dihedral_Angle::Reaction_Dihedral_Angle()
{
  _desired_angle = 0.0;

  return;
}

void
Reaction_Dihedral_Angle::all_atoms_in_scaffold()
{
  for (int i = 0; i < 4; i++)
  {
    _atom[i].set_in_scaffold ();
  }

  return;
}

int
Reaction_Dihedral_Angle::construct_from_msi_attribute(const msi_attribute * att)
{
  const_IWSubstring m;
  att->value(m);

  if (5 != m.nwords())
  {
    cerr << "Reaction_Dihedral_Angle::construct_from_msi_attribute: attribute must have 5 tokens\n";
    cerr << m << endl;
    return 0;
  }

  for (int i = 0; i < 4; i++)
  {
    const_IWSubstring token;
    (void) m.word(i, token);

    if (! _atom[i].construct(token))
    {
      cerr << "Reaction_Dihedral_Angle::construct_from_msi_attribute: cannot parse attribute " << i << " '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }
  }

// The last token is the angle

  const_IWSubstring token = m.word(4);

  if (! token.numeric_value(_desired_angle))
  {
    cerr << "Reaction_Dihedral_Angle::construct_from_msi_attribute: invalid angle '" << token << "'\n";
    return 0;
  }

  _desired_angle = static_cast<angle_t>(_desired_angle * DEG2RAD);

  return 1;
}

int
Reaction_Dihedral_Angle::write_msi(std::ostream & os,
                                   const const_IWSubstring & ind,
                                   const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"";
  for (int i = 0; i < 4; i++)
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
Reaction_Dihedral_Angle::adjust_matched_atoms_in_component(const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < 4; i++)
  {
    _atom[i].adjust_matched_atoms_in_component(xref);
  }

  return 1;
}

//#define DEBUG_PROCESS_DIHEDRAL_ANGLE

int
Reaction_Dihedral_Angle::process(Molecule & m, 
                                 const Set_of_Atoms * scaffold_embedding,
                                 const Enumeration_Temporaries & etmp) const
{
#ifdef DEBUG_PROCESS_DIHEDRAL_ANGLE
  cerr << "Processing dihedral angle " << _desired_angle << endl;
  if (nullptr != scaffold_embedding)
    cerr << "Scaffold embedding " << (*scaffold_embedding) << endl;
#endif

  atom_number_t a[4];

  for (int i = 0; i < 4; i++)
  {
    if (! determine_atom_number(*scaffold_embedding, _atom[i], etmp,
                                "Reaction_Dihedral_Angle:process:", a[i]))
      return 0;

#ifdef DEBUG_PROCESS_DIHEDRAL_ANGLE
    cerr << "embedding member " << _atom[i] << " is atom " << a[i] << '\n';
#endif
  }

  if (! m.are_bonded(a[0], a[1]) || ! m.are_bonded(a[1], a[2]) || ! m.are_bonded(a[2], a[3]))
  {
    cerr << "Reaction_Dihedral_Angle::process: atom do not form a valid dihedral set in product\n";
    if (nullptr != scaffold_embedding)
      cerr << "Scaffold embedding " << (*scaffold_embedding) << endl;

    cerr << "Atoms " << a[0] << ' ' << a[1] << ' ' << a[2] << ' ' << a[3] << endl;
    if (! m.are_bonded(a[0], a[1]))
      cerr << "Atoms " << a[0] << " and " << a[1] << " not bonded\n";
    if (! m.are_bonded(a[1], a[2]))
      cerr << "Atoms " << a[1] << " and " << a[2] << " not bonded\n";
    if (! m.are_bonded(a[2], a[3]))
      cerr << "Atoms " << a[2] << " and " << a[3] << " not bonded\n";

    return 0;
  }

#ifdef DEBUG_PROCESS_DIHEDRAL_ANGLE
  cerr << "Before setting angle " << (m.signed_dihedral_angle(a[0], a[1], a[2], a[3]) * RAD2DEG) << '\n';
#endif

  if (! m.set_dihedral(a[0], a[1], a[2], a[3], _desired_angle))
  {
    cerr << "Reaction_Dihedral_Angle::process:cannot set angle, atoms " << a[0] << ',' << a[1] << ',' << a[2] << ',' << a[3] << endl;
    write_isotopically_labelled_smiles(m, cerr);
    return 0;
  }

#ifdef DEBUG_PROCESS_DIHEDRAL_ANGLE
  cerr << "After setting angle " << (m.signed_dihedral_angle(a[0], a[1], a[2], a[3]) * RAD2DEG) << '\n';
#endif

  if (bump_tolerance > static_cast<distance_t>(0.0))
  {
    if (! m.bump_check(a[0], a[1], a[2], a[3], bump_tolerance))
    {
      cerr << "Reaction_Dihedral_Angle::process: bump check failed tolerance " << bump_tolerance << endl;
      cerr << "Twist " << a[0] << ',' << a[1] << ',' << a[2] << ',' << a[3] << endl;
      return 0;
    }
  }

  return 1;
}
