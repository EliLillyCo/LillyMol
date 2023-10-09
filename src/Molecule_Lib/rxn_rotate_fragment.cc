#include <stdlib.h>

#include "misc.h"

#include "iwreaction.h"

Reaction_Rotate_Fragment::Reaction_Rotate_Fragment()
{
  _desired_angle = static_cast<angle_t> (0.0);

  return;
}

void
Reaction_Rotate_Fragment::all_atoms_in_scaffold ()
{
  for (int i = 0; i < 3; i++)
  {
    _atom[i].set_in_scaffold ();
  }

  return;
}

int
Reaction_Rotate_Fragment::construct_from_msi_attribute (const msi_attribute * att)
{
  const_IWSubstring m;
  att->value (m);

  if (4 != m.nwords ())
  {
    cerr << "Reaction_Rotate_Fragment::construct_from_msi_attribute: attribute must have 4 tokens\n";
    cerr << m << endl;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  for (int j = 0; j < 3; j++)
  {
    (void) m.nextword (token, i);

    if (! _atom[j].construct (token))
    {
      cerr << "Reaction_Rotate_Fragment::construct_from_msi_attribute: cannot parse attribute " << j << " '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }

//  cerr << "Built " << _atom[j] << endl;
  }

// The last token is the angle

  (void) m.nextword (token, i);

  if (! token.numeric_value (_desired_angle))
  {
    cerr << "Reaction_Rotate_Fragment::construct_from_msi_attribute: invalid angle '" << token << "'\n";
    return 0;
  }

  _desired_angle = static_cast<angle_t>(_desired_angle * DEG2RAD);

//write_msi (cerr, " ", "built");

  return 1;
}

int
Reaction_Rotate_Fragment::write_msi (ostream & os,
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

  return os.good ();
}

int
Reaction_Rotate_Fragment::adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < 3; i++)
  {
    _atom[i].adjust_matched_atoms_in_component (xref);
  }

  return 1;
}

//#define DEBUG_PROCESS_ROTATE_FRAGMENT

static int
identify_fragment(Molecule & m,
                  atom_number_t zatom,
                  int * f)

{
  f[zatom] = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 1;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);
    if (f[j])
      continue;

    rc += identify_fragment(m, j, f);
  }

  return rc;
}

int
Reaction_Rotate_Fragment::process (Molecule & m,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp) const
{
#ifdef DEBUG_PROCESS_ROTATE_FRAGMENT
  write_msi (cerr, " ", "process");
  cerr << "Processing bond angle " << _desired_angle << endl;
  if (nullptr != scaffold_embedding)
    cerr << "Scaffold embedding " << (*scaffold_embedding) << endl;
  cerr << m.smiles() << endl;
#endif

  atom_number_t a[3];

  for (int i = 0; i < 3; i++)
  {
    if (! determine_atom_number (*scaffold_embedding, _atom[i], etmp, "Reaction_Rotate_Fragment:process:", a[i]))
      return 0;

#ifdef DEBUG_PROCESS_ROTATE_FRAGMENT
    cerr << _atom[i] << " is atom " << a[i] << endl;
#endif
  }

// should do bump checking

#ifdef DEBUG_PROCESS_ROTATE_FRAGMENT
  cerr << "Atoms " << a[0] << " '" << m.smarts_equivalent_for_atom (a[0]) << "', " << a[1] << " '" << m.smarts_equivalent_for_atom (a[1]) << "', " << a[2] << " '" << m.smarts_equivalent_for_atom (a[2]) << "', angle " << (_desired_angle * RAD2DEG) << "\n";
#endif

  cerr << "Atoms are " << a[0] << " " << a[1] << " and " << a[2] << endl;
  if (m.in_same_ring(a[0], a[1]))
  {
    cerr << "Reaction_Rotate_Fragment::process:atoms in a ring, cannot process\n";
    return 0;
  }

  int matoms = m.natoms();

  int * f = new_int(matoms); std::unique_ptr<int[]> free_f(f);

  f[a[0]] = 1;
  int atoms_being_moved = identify_fragment(m, a[1], f);
  if (1 == atoms_being_moved)
    return 1;

  angle_t current_angle = m.bond_angle(a[0], a[1], a[2]);

  cerr << "Reaction_Rotate_Fragment::process:moving " << atoms_being_moved << " atoms, current_angle " << current_angle << " desired " << _desired_angle << endl;

  if (fabs(current_angle - _desired_angle) < 0.01)
    return 1;

  f[a[0]] = 0;

  const Atom * a0 = m.atomi(a[0]);
  const Atom * a1 = m.atomi(a[1]);
  const Atom * a2 = m.atomi(a[2]);

  Space_Vector<double> v0(a0->x() - a1->x(), a0->y() - a1->y(), a0->z() - a1->z());
  Space_Vector<double> v2(a2->x() - a1->x(), a2->y() - a1->y(), a2->z() - a1->z());

  v0.cross_product(v2);
  v0.normalise();

  coord_t xoffset = a1->x();
  coord_t yoffset = a1->y();
  coord_t zoffset = a1->z();

  Set_of_Atoms s;
  for (int i = 0; i < matoms; i++)
  {
    if (f[i])
    {
      const Atom * ai = m.atomi(i);
      m.setxyz(i, ai->x() - xoffset, ai->y() - yoffset, ai->z() - zoffset);
      s.add(i);
    }
  }

  assert (s.number_elements() == atoms_being_moved);

  double theta = _desired_angle - current_angle;

  cerr << "Around " << v0 << ", atoms " << s << endl;

  m.rotate_atoms(v0, theta, s);

  for (int i = 0; i < matoms; i++)
  {
    if (f[i])
    {
      const Atom * ai = m.atomi(i);
      m.setxyz(i, ai->x() + xoffset, ai->y() + yoffset, ai->z() + zoffset);
    }
  }

  cerr << "After rotation angle is " << RAD2DEG * m.bond_angle(a[0], a[1], a[2]) << endl;

  return 1;
}
