/*
  Stuff for dealing with mopac files.

Here is an example

 MNDO ESP 1SCF NOINTER NOMM
 O=C(C1=CC=CC=C1)NC1=CC=C2NC=C(C3CCN(C)CC3)C2=C1 PBCHM22611554
 R(1) 329511
   XX     0.000   0     0.00000   0     0.00000   0       0   0   0
   XX     1.000   0     0.00000   0     0.00000   0       1   0   0
   XX     1.000   0    90.00000   0     0.00000   0       2   1   0
   C      4.252   0    90.00000   0   -90.85146   0       3   2   1
   C      1.510   1    88.07703   0   -15.56192   0       4   3   2
   C      1.540   1   109.94521   1   -81.79817   0       4   5   3
   C      1.539   1   109.94521   1   122.98030   1       4   5   6
   C      1.340   1   126.32602   1  -118.50986   1       5   4   7
   C      1.480   1   126.32602   1   180.00000   0       5   4   8


  Looks like the first three lines define the x and y axes. Also looks
  as if this will preserve the coordinates of the molecule.
*/

#include "molecule.h"
#include "mopac.h"

class Atom_and_Number
{
  private:
    const Atom * _atom;
    int          _n;

  public:
    int n () const { return _n;}
    const Atom * a () const { return _atom;}

    int set_n (int i) { return _n = i;}
    void set_atom (const Atom * q) { _atom = q;}
};

static int
output_record (ofstream & output,  
               const Mopac_Output_Control_Object & mopac_output_control,
               const Atom_and_Number & a1,
               const Atom_and_Number & a2,
               const Atom_and_Number & a3,
               const Atom_and_Number & a4)
{
  const Atom * a = a1.a ();
  const char * s = a->atomic_symbol ();
  output << " " << s;
  if (1 == strlen (s))
    output << " ";

  distance tmp = a->distance (*(a2.a ()));
  output << endl;

  return 1;
};

int
Molecule::_write_molecule_mop (ofstream & output, 
                               const Mopac_Output_Control_Object & mopac_output_control,
                               int * ordering) const
{
  assert (ok ());
  assert (output.good ());
  assert (NULL != ordering);

  atom_number_t start_atom = 0;
  if (INVALID_ATOM_NUMBER != mopac_output_control.start_atom ())
    start_atom = mopac_output_control.start_atom ();

  return 1;
}

/*
  This variant preserves native ordering
*/

int
Molecule::_write_molecule_mop_native (ofstream & output, 
                      const Mopac_Output_Control_Object & mopac_output_control) const
{
  assert (ok ());
  assert (output.good ());

// The first atom is treated specially.

  Atom atom0 (0);
  Atom atom1 (0);
  Atom atom2 (0);

  atom0.set_xyz (0.0, 0.0, 0.0);
  atom1.set_xyz (1.0, 0.0, 0.0);
  atom2.set_xyz (1.0, 1.0, 0.0);

  int matoms = natoms ();

  if (0 == matoms)
    return 1;

  if (! output_record (output, mopac_output_control, _things[0], &atom2, &atom1, &atom0);
    return 0;

  if (1 == matoms)
    return 1;

  if (! output_record (output, mopac_output_control, _things[1], _things[0], &atom2, &atom1))
    return 0;

  if (2 == matoms)
    return 1;

  if (! output_record (output, mopac_output_control, _things[2], _things[1], _things[0], &atom2))
    return 0;

  if (3 == matoms)
    return 1;

  for (int i = 3; i < matoms; i++)
  {
    if (! output_record (output,  mopac_output_control, _things[i], _things[i - 1],
                                  _things[i - 2], things[i - 3]))
      return 0;
  }

  return 1;
}

int
Molecule::write_molecule_mop (ofstream & output, 
                              const Mopac_Output_Control_Object & mopac_output_control) const
{
  assert (ok ());
  assert (output.good ());

  output << mopac_output_control.keywords () << endl;
  output << "Blank\n";
  output << name () << endl;
  output << " XX     0.000   0     0.00000   0     0.00000   0       0   0   0\n";
  output << " XX     1.000   0     0.00000   0     0.00000   0       1   0   0\n";
  output << " XX     1.000   0    90.00000   0     0.00000   0       2   1   0\n";

  int rc = _write_molecule_mop (output, mopac_output_control);

  output << "   0      0.000   0     0.00000   0     0.00000   0       0   0   0\n";

  return rc;
}

int
Molecule::_write_molecule_mop (ofstream & output, 
                               const Mopac_Output_Control_Object & mopac_output_control) const
{
  assert (output.good ());

  if (mopac_output_control.native_ordering ())
    return _write_molecule_mop_native (output, mopac_output_control);

  int * ordering = new_int (_number_elements);

  int rc = _write_molecule_mop_order (output, mopac_output_control, ordering);

  delete ordering;

  return rc;
}

int
Molecule::write_molecule_mop (const char * fname, 
                              const Mopac_Output_Control_Object & mopac_output_control) const
{

  ofstream output (fname, ios::out);
  if (! output.good ())
  {
    cerr << "Molecule::write_molecule_mop: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_molecule_mop (output, mopac_output_control);
}
