#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "iwreaction.h"

using std::cerr;
using std::endl;

Reaction_3D_Replace::Reaction_3D_Replace()
{
  _n = 0;
  _weight = nullptr;

  _a1 = nullptr;
  _a2 = nullptr;

  return;
}

Reaction_3D_Replace::~Reaction_3D_Replace()
{
  if (nullptr != _weight)
    delete [] _weight;

  if (nullptr != _a1)
    delete [] _a1;

  if (nullptr != _a2)
    delete [] _a2;

  return;
}

int
Reaction_3D_Replace::construct_from_msi_attribute(const msi_attribute * att)
{
  const_IWSubstring m;
  att->value(m);

  _n = m.nwords();

  if (0 == _n || (_n != (_n / 2) * 2))
  {
    cerr << "Reaction_3D_Replace::construct_from_msi_attribute: attribute must have an even number of tokens\n";
    cerr << m << endl;
    return 0;
  }

  _n = _n / 2;

  _a1 = new Matched_Atom_in_Component[_n];
  _a2 = new Matched_Atom_in_Component[_n];

  int i = 0;
  const_IWSubstring token;

  _n = 0;

  while (m.nextword(token, i))
  {
//  cerr << "Reaction_3D_Replace::construct_from_msi_attribute:token '" << token << "'\n";
    if (! _a1[_n].construct(token))
    {
      cerr << "Reaction_3D_Replace::construct_from_msi_attribute: cannot parse '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }

    m.nextword(token, i);

    if (! _a2[_n].construct(token))
    {
      cerr << "Reaction_3D_Replace::construct_from_msi_attribute: cannot parse '" << token << "'\n";
      cerr << m << endl;
      return 0;
    }

    _n++;
  }

  _weight = new double[_n];

  _weight[0] = 1.0;

  for (int i = 1; i < _n; i++)
  {
    _weight[i] = 0.1;
  }

//write_msi (cerr, " ", "built");

  return 1;
}

int
Reaction_3D_Replace::write_msi (std::ostream & os,
                                const const_IWSubstring & ind,
                                const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"";

  for (int i = 0; i < _n; i++)
  {
    if (i > 0)
      os << ' ';

    os << _a1[i] << ' ' << _a2[i];
  }

  os << "\")\n";

  return os.good();
}


int
Reaction_3D_Replace::adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < _n; i++)
  {
    _a1[i].adjust_matched_atoms_in_component(xref);
    _a2[i].adjust_matched_atoms_in_component(xref);
  }

  return 1;
}

static int
identify_moving_atoms_by_atom_number (const Molecule & m,
                                      atom_number_t only_greater_than,
                                      atom_number_t astart,
                                      int * moving)
{
  moving[astart] = 1;

  int rc = 1;

  const Atom * a = m.atomi(astart);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(astart, i);

    if (moving[j])
      continue;

    if (j < only_greater_than)
      continue;

    rc += identify_moving_atoms_by_atom_number(m, only_greater_than, j, moving);
  }

  return rc;
}

static int
identify_moving_atoms (const Molecule & m,
                       atom_number_t astart,
                       int * moving)
{
  moving[astart] = 1;

  int rc = 1;

  const Atom * a = m.atomi(astart);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(astart, i);

    if (moving[j])
      continue;

    rc += identify_moving_atoms(m, j, moving);
  }

  return rc;
}

extern "C" void u3b_(const double * w, double * c1, double * c2, const int * n, const int * mode, double *rms, double * u, double * t, int  * ier);

int
Reaction_3D_Replace::process(Molecule & m,
                             const Set_of_Atoms * scaffold_embedding,
                             const Enumeration_Temporaries & etmp) const
{
//#define DEBUG_PROCESS_3D_REPLACE
#ifdef DEBUG_PROCESS_3D_REPLACE
  write_msi(cerr, " ", "process");
  if (nullptr != scaffold_embedding)
    cerr << "Scaffold embedding " << (*scaffold_embedding) << endl;
  cerr << m.smiles() << endl;
  atom_number_t first_fixed_atom;
#endif

  double * c1 = new double[_n * 3 * 2]; std::unique_ptr<double[]> free_c1(c1);
  double * c2 = c1 + (_n * 3);

  Set_of_Atoms moving_atoms;

  for (int i = 0; i < _n; i++)
  {
    atom_number_t j;
    if (! determine_atom_number(*scaffold_embedding, _a1[i], etmp, "Reaction_Rotate_Fragment:process:", j))
      return 0;

    const Atom * a = m.atomi(j);
#ifdef DEBUG_PROCESS_3D_REPLACE
    if (0 == i)
     first_fixed_atom = j;
#endif

    c1[3 * i] = a->x();
    c1[3 * i + 1] = a->y();
    c1[3 * i + 2] = a->z();

    if (! determine_atom_number(*scaffold_embedding, _a2[i], etmp, "Reaction_Rotate_Fragment:process:", j))
      return 0;

    a = m.atomi(j);

    c2[3 * i] = a->x();
    c2[3 * i + 1] = a->y();
    c2[3 * i + 2] = a->z();

    moving_atoms.add(j);
  }

  const int matoms = m.natoms();

  // Which atoms will be moved.
  int * moving = new_int(matoms); std::unique_ptr<int[]> free_moving(moving);

  // u3b_ will fail if given just one atom. In that case we can do a simple
  // translation.
  if (moving_atoms.size() == 1) {
    identify_moving_atoms_by_atom_number(m, moving_atoms[0], moving_atoms[0], moving);
    return DoTranslation(m, c1, c2, moving);
  }

  int mode = 1;
  double u[9];
  double rms;
  double t[3];
  int ier = 0;

  u3b_(_weight, c1, c2, &_n, &mode, &rms, u, t, &ier);

  if (0 == ier)
    ;
  else if (-1 == ier)
    cerr << "Reaction_Rotate_Fragment::process:superposition not unique, but optimal\n";
  else
  {
    cerr << "Reaction_Rotate_Fragment::process:u3b failed\n";
    return 0;
  }

#ifdef DEBUG_PROCESS_3D_REPLACE
  cerr << "RMS " << rms << endl;

  for (int i = 0; i < 3; i++)
  {
    cerr << "t[" << i << "] = " << t[i] << endl;
  }
  cerr << "First fixed atom " << m.x(first_fixed_atom) << ',' << m.y(first_fixed_atom) << ',' << m.z(first_fixed_atom) << endl;
  cerr << "First moving atom " << m.x(moving_atoms[0]) << ',' << m.y(moving_atoms[0]) << ',' << m.z(moving_atoms[0]) << endl;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      cerr << ' ' << u[i * 3 + j];
    }
    cerr << endl;
  }
#endif

  // Identify the bond that separates the moving from non-moving part of
  // the molecule.  We assume that it is attached to the first of the
  // moving atoms If there is just one atom in the alignment, we assume
  // that everything connected and with a higher atom number is moving

#ifdef DEBUG_PROCESS_3D_REPLACE
  cerr << "There are " << moving_atoms.number_elements() << " moving atoms\n";
#endif

  if (1 == moving_atoms.number_elements())
    identify_moving_atoms_by_atom_number(m, moving_atoms[0], moving_atoms[0], moving);
  else
  {
    moving[moving_atoms[0]] = 1;
    for (int i = 1; i < moving_atoms.number_elements(); i++)
    {
      identify_moving_atoms(m, moving_atoms[i], moving);
    }
  }

  double rotmat11 = u[0];
  double rotmat12 = u[1];
  double rotmat13 = u[2];
  double rotmat21 = u[3];
  double rotmat22 = u[4];
  double rotmat23 = u[5];
  double rotmat31 = u[6];
  double rotmat32 = u[7];
  double rotmat33 = u[8];

  for (int i = 0; i < matoms; i++)
  {
    if (0 == moving[i])
      continue;

#ifdef DEBUG_PROCESS_3D_REPLACE
    cerr << "Atom " << i << " '" << m.smarts_equivalent_for_atom(i) << "' is moving\n";
#endif

    const Atom * a = m.atomi(i);

    double x0 = a->x() - t[0];
    double y0 = a->y() - t[1];
    double z0 = a->z() - t[2];

    double xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    double yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    double zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    m.setxyz( i, static_cast<coord_t> (xx), static_cast<coord_t> (yy), static_cast<coord_t> (zz) );
  }

  return 1;
}

// There was only one matched atom specified in the 3d replacement
// so all we need do is translate the molecule.
// The coodinates of two positions are in `initial` and `destination`.
// For those atoms that are set in `moving`, translate those
// atoms according to the diff btw `initial` and `destination`.
int
Reaction_3D_Replace::DoTranslation(Molecule& m,
                const double * initial,
                const double * destination,
                const int* moving) const {
  Coordinates initial_c(initial[0], initial[1], initial[2]);
  Coordinates destination_c(destination[0], destination[1], destination[2]);

  Coordinates delta = destination_c - initial_c;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (! moving[i]) {
      continue;
    }
    m.setx(i, m.x(i) - delta.x());
    m.sety(i, m.y(i) - delta.y());
    m.setz(i, m.z(i) - delta.z());
  }

  return 1;
}
