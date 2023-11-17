/*
  There is a lot of code common to nn_training and gfp_predictions
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "nn_stuff.h"
#include "prediction_weight_function.h"

using std::cerr;
using std::endl;

int
display_standard_nn_stuff (std::ostream & os)
{
  os << " -n n1,n2,n3..    number of neighbours to use\n";
  os << " -n n1-n2         a range of neighbours to use\n";
  os << " -t delta=dx      step size for testing distance thresholds\n";
  os << " -t min=x.x       minimum distance threshold to test\n";
  os << " -t max=x.x       maximum distance threshold to test\n";
  os << " -t require=nn    neighbours to require with each -t combination\n";
  os << " -W ...           weight function options, enter '-W help' for details\n";

  return os.good ();
}

static int
last_step_n_qualifier (resizable_array<int> & neighbours_to_use,
                       const const_IWSubstring & n)
{
  int j;

  if (! n.numeric_value (j) || j < 1)
  {
    cerr << "Neighbour specifications must be whole positive numbers\n";
    return 0;
  }

  neighbours_to_use.add (j);

  return 1;
}

static int
parse_n_qualifier_range (resizable_array<int> & neighbours_to_use,
                         const const_IWSubstring & n)
{
  const_IWSubstring n1, n2;

  if (! n.split (n1, '-', n2) || 0 == n1.length () || 0 == n2.length ())
  {
    cerr << "Invalid range specifier '" << n << "'\n";
    return 0;
  }

  int r1;
  if (! n1.numeric_value (r1) || r1 < 1)
  {
    cerr << "Invalid first range component '" << n1 << "' from '" << n << "'\n";
    return 0;
  }

  int r2;
  if (! n2.numeric_value (r2) || r2 < r1)
  {
    cerr << "Invalid first range component '" << n2 << "' from '" << n << "'\n";
    return 0;
  }

  for (int i = r1; i <= r2; i++)
  {
    neighbours_to_use.add (i);
  }

  return 1;
}

static int
parse_n_qualifier (resizable_array<int> & neighbours_to_use,
                   const const_IWSubstring & n)
{
  if ("none" == n)
    return 1;

  int i = 0;
  const_IWSubstring token;

  while (n.nextword (token, i, ','))
  {
    if (token.contains ('-'))
    {
      if (! parse_n_qualifier_range (neighbours_to_use, token))
      {
        cerr << "INvalid range specifier '" << token << "'\n";
        return 0;
      }
    }
    else if ( ! last_step_n_qualifier (neighbours_to_use, token))    // just a single number
    {
      cerr << "Invalid number '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
parse_dash_n (Command_Line & cl,
              char flag,
              resizable_array<int> & neighbours_to_use,
              int verbose)
{
  if (cl.option_present (flag))
  {
    int i = 0;
    const_IWSubstring n;
    while (cl.value ('n', n, i++))
    {
      if (! parse_n_qualifier (neighbours_to_use, n))
      {
        cerr << "INvalid -n qualifier '" << n << "'\n";
        return 0;
      }
    }

    if (verbose)
    {
      cerr << "Will base estimates on ";
      for (int i = 0; i < neighbours_to_use.number_elements (); i++)
      {
        cerr << ' ' << neighbours_to_use[i];
      }

      cerr << " neighbours\n";
    }
  }

  return 1;
}

int 
parse_dash_t (Command_Line & cl, 
              char flag,
              resizable_array<float> & distances_to_use,
              resizable_array<int> & min_neighbours_required,
              double max_distance,
              int verbose)
{
  if (cl.option_present (flag))
  {
    double mind = 0.0;
    double maxd = 0.0;
    double dx = 0.0;

    int i = 0;
    const_IWSubstring t;
    while (cl.value (flag, t, i++))
    {
      if (t.starts_with ("max="))
      {
        t.remove_leading_chars (4);
        if (! t.numeric_value (maxd) || maxd <= 0.0)
        {
          cerr << "INvalid maximum distance specifier 'maxd=" << t << "'\n";
          return 0;
        }
      }
      else if (t.starts_with ("min="))
      {
        t.remove_leading_chars (4);
        if (! t.numeric_value (mind) || mind <= 0.0)
        {
          cerr << "Invalid minimum distance specifier 'mind=" << t << "'\n";
          return 0;
        }
      }
      else if (t.starts_with ("delta="))
      {
        t.remove_leading_chars (6);
        if (! t.numeric_value (dx) || dx <= 0.0)
        {
          cerr << "Invalid distance delta specification 'delta=" << t << "'\n";
          return 0;
        }
      }
      else if (t.starts_with ("require="))
      {
        t.remove_leading_chars (8);
        if (! parse_n_qualifier (min_neighbours_required, t))
        {
          cerr << "Invalid min neighbour specification 'require=" << t << "'\n";
          return 0;
        }
      }
      else if ("none" == t)
      {
        if (verbose)
          cerr << "No distances studies performed\n";
        return 1;
      }
      else
      {
        cerr << "Unrecognised -t qualifier '" << t << "'\n";
        return 0;
      }
    }

    if (mind > 0.0 && 0.0 == dx && 0.0 == maxd)    // OK to enter just a minimum number
      ;
    else if (maxd <= mind || dx <= 0.0 || mind > maxd)
    {
      cerr << "Invalid or incomplete specification for distances: max " << maxd << " dx " << dx << endl;
      return 0;
    }

//  When doing tests setting either the distance or the number of neighbours, we may have just one distance

    if (mind + dx > maxd)
    {
      distances_to_use.add (static_cast<float> (mind));
      if (verbose)
        cerr << "Will test threshold " << mind << endl;

      return 1;
    }

    if (0.0 == mind)
      mind = dx;

    for (double x = mind; x <= (maxd + 1.0e-04); x+= dx)    // allow for some numeric fuzz
    {
      if (max_distance > 0.0 && x - dx > max_distance)
      {
        cerr << "Largest distance in input is " << max_distance << " cannot do " << x << endl;
        break;
      }

      distances_to_use.add (static_cast<float> (x));
      if (verbose)
        cerr << " will test threshold " << x << endl;
    }
  }

  return 1;
}

/*
  Fingerprint distances are between 0 and 1.0 by default
*/

static similarity_type_t _maximum_possible_distance = 1.0;

similarity_type_t
maximum_possible_distance ()
{
  return _maximum_possible_distance;
}

void
set_maximum_possible_distance (similarity_type_t d)
{
  assert (d > 0.0);

  _maximum_possible_distance = d;

  return;
}
