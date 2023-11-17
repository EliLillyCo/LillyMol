#include <stdlib.h>
#include <iostream>

#include "nn_specification.h"
#include "nn_stuff.h"

using std::cerr;
using std::endl;

/*
  Common code for parsing neighbour conditions

  Typical specifications look like:

  0.22
  0.18-3N
  0.22,W=linear
  0.18-3N,W=step0.5
  3N,W=def,S1=0.4
*/

static int
parse_special_class_specification (const const_IWSubstring & token,
                                   int & special_class,
                                   double & special_class_probability)
{
  assert (token.starts_with ("S:"));

  const_IWSubstring mytoken (token);
  mytoken.remove_leading_chars (2);    // get rid of the 'S:'

  const_IWSubstring sclass, sthreshold;
  if (! mytoken.split (sclass, '=', sthreshold))
  {
    cerr << "parse_special_class_specification: invalid special class '" << token << "'\n";
    return 0;
  }

  int c;     // store in temporary variable first

  if (! sclass.numeric_value (c) || c < 0)
  {
    cerr << "parse_special_class_specification: invalid special class value '" << token << "', '" << sclass << "'\n";
    return 0;
  }

  if (! sthreshold.numeric_value (special_class_probability) || special_class_probability <= 0.0 || special_class_probability > 1.0)
  {
    cerr << "parse_special_class_specification: invalid probability '" << token << "'\n";
    return 0;
  }

  special_class = c;

  return 1;
}

static int
is_distance (const const_IWSubstring & d, double & distance)
{
  double tmp;

  if (d.numeric_value (tmp) && tmp >= 0.0 && tmp <= maximum_possible_distance ())
  {
    distance = tmp;
    return 1;
  }

  return 0;
}

static int
is_neighbours (const const_IWSubstring & n,
               int & number_neighbours)
{
  if (! n.ends_with ('N'))
    return 0;

  const_IWSubstring myn = n;
  myn.chop (1);

  if (myn.numeric_value (number_neighbours) && number_neighbours > 0)
    return 1;

  return 0;
}

/*
  Parse the Neighbours and Distance part
  0.20
  4N
  0.1-3N
*/

static int
parse_neighbour_specification (const const_IWSubstring & token,
                               int & neighbours,
                               double & distance)
{
  const_IWSubstring c1, c2;

  if (token.split (c1, '-', c2))
  {
    if (0 == c1.length () || 0 == c2.length ())
    {
      cerr << "parse_neighbour_specification: invalid specification '" << token << "'\n";
      return 0;
    }

    if (is_distance (c1, distance) && is_neighbours (c2, neighbours))
      return 1;

    cerr << "parse_model_specifications: cannot process both '" << c1 << "' and '" << c2 << "'\n";
    return 0;
  }

//cerr << "Before checks " << neighbours << " dist " << distance << " token '" << token << "'\n";

  if (is_distance (token, distance))
    ;
  else if (is_neighbours (token, neighbours))
    ;
  else
  {
    cerr << "parse_model_specifications: cannot process '" << token << "'\n";
    return 0;
  }

//cerr << "After checks " << neighbours << " dist " << distance << endl;

  return 1;
}

int
parse_model_specifications (const const_IWSubstring & c,
                            const_IWSubstring & weight_function,
                            int & special_class,
                            double & special_class_probability,
                            double & distance,
                            int & number_neighbours)
{
  int i = 0;
  const_IWSubstring token;
  while (c.nextword (token, i, MODEL_COMPONENT_SEPARATOR))
  {
    if (token.starts_with ("W:"))
    {
      weight_function = token;
      weight_function.remove_leading_chars (2);
    }
    else if (token.starts_with ("S:"))
    {
      if (! parse_special_class_specification (token, special_class, special_class_probability))
      {
        cerr << "parse_model_specifications: invalid special class '" << token << "'\n";
        return 0;
      }
    }
    else if (! parse_neighbour_specification (token, number_neighbours, distance))
    {
      cerr << "Invalid neighbour specification '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}
