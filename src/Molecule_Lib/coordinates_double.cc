#include <stdlib.h>

#include "coordinates.h"
#include "atom.h"

Coordinates_double::Coordinates_double (const Atom & a) : Space_Vector<double> (a.x (), a.y (), a.z ())
{
}

Coordinates_double::Coordinates_double (const Space_Vector<double> & a) : Space_Vector<double> (a)
{
}
