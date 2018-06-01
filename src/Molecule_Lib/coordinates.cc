#include <stdlib.h>

#include "coordinates.h"
#include "atom.h"

Coordinates::Coordinates (const Atom & a) : Space_Vector<coord_t> (a._x, a._y, a._z)
{
}

Coordinates::Coordinates (const Space_Vector<coord_t> & a) : Space_Vector<coord_t> (a)
{
}
