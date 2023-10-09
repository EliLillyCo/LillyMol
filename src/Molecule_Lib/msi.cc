#include <iostream>
#include <ctype.h>
#include <fstream>
using std::cerr;
using std::endl;

#include "Foundational/iwmisc/msi_object.h"
#include "Foundational/iwmisc/misc.h"

#include "molecule.h"
#include "rwmolecule.h"

/*
  The atomic number will be something like "6 C" without the quotes
*/

static int
_determine_atomic_number(const IWString & buffer, int & anum)
{
  if (! buffer.is_int(anum))
  {
    cerr << "determine atomic number: bad form '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

#define NAME_OF_MSI_MOLECULE_OBJECT "Model"
#define NAME_OF_MSI_ATOM_OBJECT "Atom"
#define NAME_OF_MSI_BOND_OBJECT "Bond"

Atom *
_convert_to_atom(const msi_object * msi)
{
  assert(nullptr != msi && msi->ok());
  assert(NAME_OF_MSI_MOLECULE_OBJECT == msi->name());

  const msi_attribute * acl = msi->attribute("ACL");
  if (nullptr == acl)
  {
    cerr << "convert to atom: no ACL type for object\n";
    return nullptr;
  }

  int zz;
  if (! _determine_atomic_number(acl->stringval(), zz))
  {
    cerr << "_convert to atom: cannot convert '" << acl->stringval() << "' to atomic number\n";
    return nullptr;
  }

  Atom * rc = new Atom(zz);
  assert(nullptr != rc);

// get the coordinates.

  const msi_attribute * msi_coords = msi->attribute("XYZ");
  coord_t x, y, z;
  x = 0.0;
  y = 0.0;
  z = 0.0;

  if (nullptr == msi_coords)
  {
    cerr << "process msi atom: no XYZ attribute, coordinates set to 0.0\n";
  }
  else if (3 != msi_coords->number_double_values())
  {
    cerr << *msi_coords;
    cerr << "process msi atom: incorrect number of coordinates, " << msi_coords->number_double_values()
         << " ignored\n";
  }
  else
  {
    x = static_cast<coord_t>( msi_coords->double_multi_value(0) );
    y = static_cast<coord_t>( msi_coords->double_multi_value(1) );
    z = static_cast<coord_t>( msi_coords->double_multi_value(2) );
  }

  rc->setxyz(x, y, z);

  return rc;
}

static int
_process_msi_atom_to_molecule(const msi_object * msi, Molecule & m)
{
  assert(nullptr != msi && msi->ok());
  assert(NAME_OF_MSI_ATOM_OBJECT == msi->name());
  assert(m.ok());

  Atom * a = _convert_to_atom(msi);
  if (nullptr == a)
  {
    cerr << "process atom to molecule: cannot convert to atom\n";
    return 0;
  }

  m.add(a);
  const msi_attribute * chg = msi->attribute("Charge");
  if (nullptr != chg)
  {
    int a = m.natoms() - 1;
    charge_t q;
    if (! chg->value(q))
    {
      cerr << "_process_msi_atom_to_molecule: Bad change value " << *msi;
      return 0;
    }
//  cerr << "Setting charge on atom " << a << " to " << q << endl;
    m.set_charge(a, q);
  }

  assert(m.ok());

  return 1;
}

static int
_convert_msi_atoms_to_molecule(const msi_object & msi_molecule_object,
                                Molecule & m, int * xref)
{
  assert(msi_molecule_object.ok());
  assert(NAME_OF_MSI_MOLECULE_OBJECT == msi_molecule_object.name());
  assert(m.ok());

  int n = msi_molecule_object.number_elements();
  int natoms = 0;
  int rc = 1;
  for (int i = 0; i < n && 1 == rc; i++)
  {
    const msi_object * msi = msi_molecule_object[i];
    assert(nullptr != msi && msi->ok());

    if (NAME_OF_MSI_ATOM_OBJECT == msi->name())
    {
      rc = _process_msi_atom_to_molecule(msi, m);
      xref[msi->object_id()] = natoms;
      natoms++;
    }
  }

  return rc;
}

static int
_process_msi_bond_to_molecule(const msi_object * msi, Molecule & m,
                              const int * xref)
{
  assert(msi->ok());
  assert(NAME_OF_MSI_BOND_OBJECT == msi->name());
  assert(m.ok());
  assert(nullptr != xref);

  const msi_attribute * a1o = msi->attribute("Atom1");
  const msi_attribute * a2o = msi->attribute("Atom2");
  assert(nullptr != a1o);
  assert(nullptr != a2o);

  const msi_attribute * boo = msi->attribute("Order");
  if (nullptr == boo)
    boo = msi->attribute("Type");

  bond_type_t bt;
  int tmp;
  if (nullptr == boo)
    bt = SINGLE_BOND;
  else if (! boo->value(tmp))
  {
    cerr << "_process_msi_bond_to_molecule: illegal value for bond type " << *msi;
    return 0;
  }
  else
    bt = static_cast<bond_type_t>(tmp);

  atom_number_t a1, a2;
  if (! a1o->value(a1) || ! a2o->value(a2))
  {
    cerr << "_process_msi_bond_to_molecule: illegal value for a1 or a2 " << *msi;
    return 0;
  }
  
  return m.add_bond(xref[a1], xref[a2], bt);
}

static int
_convert_msi_bonds_to_molecule(const msi_object & msi_molecule_object,
                               Molecule & m, const int * xref)
{
  assert(NAME_OF_MSI_MOLECULE_OBJECT == msi_molecule_object.name());
  assert(m.ok());
  assert(nullptr != xref);

  int n = msi_molecule_object.number_elements();
  int rc = 1;
  for (int i = 0; i < n && 1 == rc; i++)
  {
    const msi_object * msi = msi_molecule_object[i];
    assert(msi->ok());

    if (NAME_OF_MSI_BOND_OBJECT == msi->name())
      rc = _process_msi_bond_to_molecule(msi, m, xref);
  }

  return rc;
}

static int
_convert_msi_object_to_molecule(const msi_object & msi, Molecule & m,
                                int * xref)
{
  assert(msi.ok());
  assert(m.ok());

  int rc = _convert_msi_atoms_to_molecule(msi, m, xref);

  if (1 == rc)
    rc = _convert_msi_bonds_to_molecule(msi, m, xref);

  if (1 == rc)
  {
    const msi_attribute * msi_name = msi.attribute("Label");
    if (nullptr != msi_name)
      m.set_name(msi_name->stringval());
  }

  return rc;
}

static int
_convert_msi_object_to_molecule(const msi_object & msi, Molecule & m)
{
  assert(msi.ok());
  assert(m.ok());

// We need a temporary array to cross reference msi object numbers
// to atom numbers

   int * xref = new_int(msi.number_elements() + 1, 0);

   int rc = _convert_msi_object_to_molecule(msi, m, xref);

   delete xref;

   return rc;
}

int
Molecule::read_molecule_msi_ds(iwstring_data_source & input)
{
  assert(ok());
  assert(input.good());

  input.set_ignore_pattern("#");   // when regexp works, change to "^#"
  input.set_strip_leading_blanks(1);

  msi_object msi;
  msi.resize(40);    // possibly average for molecules

  if (! msi.read(input) || (! msi.ok()))
  {
    if (input.eof())
      return 0;
    else
      return rwmolecule_error("read mol msi: msi.read failed", input);
  }
  
  int rc = _convert_msi_object_to_molecule(msi, *this);

  check_bonding();

  return rc;
}

int
Molecule::write_molecule_msi (std::ostream & os, const IWString & comments) const
{
  assert (ok());
  assert (os.good());

  os << "# MSI CERIUS2 DataModel File Version 1 0\n";
  os << "(1 Model\n";
  if (name().length())
    os << " (A C Label \"" << name() << "\")\n";

  int object_counter = 0;
  os << "(A I Id " << ++object_counter << ")\n";

  if (comments.length() && comments != _molecule_name)
    os << "(A C Label \"" << comments << "\")\n";

  int matoms = natoms();
  for (int i = 0; i < matoms; i++)
  {
    os << " (" << ++object_counter << " Atom\n";
    os << "  (A I ACL \"" << atomic_number(i) << " " << atomic_symbol(i) <<
          "\")\n";
    if (charge_on_atom(i))
      os << "  (A F Charge " << charge_on_atom(i) << ")\n";

    os.setf(std::ios::showpoint);
    os << "  (A D XYZ (" << x(i) << " " << y(i) << " " << z(i) << "))\n";
    os << " )\n";
  }

/*
  A typical bond record looks like

 (15 Bond
  (A F Order 1)      sometimes we have "Type" instead of Order
  (A O Atom1 2)
  (A O Atom2 3)
 )

 Remember, that MSI atom objects start with object 1 (the molecule object
 is assigned object 1).
*/

  int medges = nedges();
  for (int i = 0; i < medges; i++)
  {
    const Bond * b = bondi(i);

    os << " (" << ++object_counter << " Bond\n";
    os << "  (A F Order " << b->btype() << ")\n";
    os << "  (A O Atom1 " << b->a1() + 2 << ")\n";
    os << "  (A O Atom2 " << b->a2() + 2 << ")\n";
    os << " )\n";
  }

  os << ")\n";

  return 1;
}
