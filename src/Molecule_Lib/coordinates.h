#ifndef ATOM_COORDINATES_H
#define ATOM_COORDINATES_H

#include "space_vector.h"
#include "iwmtypes.h"

class Atom;

class Coordinates : public Space_Vector<coord_t>
{
  private:

  public:
    Coordinates () {};
    Coordinates (const Space_Vector<coord_t> &);

    Coordinates (coord_t cx, coord_t cy, coord_t cz) : Space_Vector<coord_t> (cx, cy, cz) {};
    Coordinates (const Coordinates & r) : Space_Vector<coord_t> (r.x (), r.y (), r.z ()) {};
    Coordinates (const Atom & a);
};

class Coordinates_double : public Space_Vector<double>
{
  private:

  public:
    Coordinates_double () {};
    Coordinates_double (const Space_Vector<double> &);

    Coordinates_double (double cx, double cy, double cz) : Space_Vector<double> (cx, cy, cz) {};
    Coordinates_double (const Coordinates & r) : Space_Vector<double> (r.x (), r.y (), r.z ()) {};
    Coordinates_double (const Atom & a);
};

#endif
