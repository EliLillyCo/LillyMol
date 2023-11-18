#include <stdlib.h>
#include <math.h>

#include "re2/re2.h"

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "prediction_weight_function.h"
#include "nn_stuff.h"

using std::cerr;
using std::endl;

int
display_weight_function_options (std::ostream & os, char c)
{
  os << " -" << c << " default       distance sensitive weighting function from Dave Cummins\n";
  os << " -" << c << " linear        linear distance function\n";
  os << " -" << c << " flat          flat distance function\n";
  os << " -" << c << " stair         each distance constant ratio of previous\n";
  os << " -" << c << " stair=ratio   specify ratio for stair function\n";
  os << " -" << c << " test          test distance function\n";
  os << " -" << c << " file=fname    read weights from 'fname' - first token on each line\n";
  os << " -" << c << " echo          echo values after assignment\n";

  return os.good();
}

static double maximum_distance = 1.0;
static double distance_delta = 0.01;
static int number_points_in_weight_array = 101;

double
maximum_distance_possible()
{
  return maximum_distance;
}

int
set_max_distance_and_resolution (double d,
                                 int n)
{
  assert (d > 0.0);

  maximum_distance = d;

  assert (n > 1);

  number_points_in_weight_array = n;

  distance_delta = d / static_cast<double> (n - 1);

  return 1;
}

int 
compute_weight_function_index (double d, int & ndx)
{
  if (d < 0.0 || d > maximum_distance)
  {
    cerr << "Invalid distance " << d << " must be between 0 and " << maximum_distance << endl;
    return 0;
  }

  ndx = static_cast<int> (d / maximum_distance * static_cast<double> (number_points_in_weight_array - 1) + 0.499);

  return 1;
}

int
initialise_distances (Command_Line & cl,
                      char flag,
                      int verbose)
{
  if (cl.option_present (flag))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value (flag, d, i++))
    {
      if (d.starts_with ("max="))
      {
        d.remove_leading_chars (4);
        if (! d.numeric_value (maximum_distance) || maximum_distance <= 0.0)
        {
          cerr << "Invalid maximum distance '" << d << "'\n";
          return 0;
        }

        if (verbose)
          cerr << "maximum distance " << maximum_distance << endl;
      }
      else if (d.starts_with ("nwt="))
      {
        d.remove_leading_chars (4);

        if (! d.numeric_value (number_points_in_weight_array) || number_points_in_weight_array < 2)
        {
          cerr << "Invalid number of points in weight arrays '" << d << "'\n";
          return 0;
        }

        if (verbose)
          cerr << number_points_in_weight_array << " points in each weight array\n";
      }
      else
      {
        cerr << "Unrecognised distance specification '" << d << "'\n";
        return 0;
      }
    }
  }

  return 1;
}

/*
  Near neighbour estimations use a weighting function
*/

NN_Weight_Function::NN_Weight_Function()
{
  _weight = nullptr;

  return;
}

NN_Weight_Function::~NN_Weight_Function()
{
  if (NULL != _weight)
    delete [] _weight;

  return;
}

void
NN_Weight_Function::_allocate_weight_array_if_needed()
{
  if (NULL != _weight)
    delete [] _weight;

  _weight = new weight_t[number_points_in_weight_array];

  return;
}

int
NN_Weight_Function::construct_from_command_line (Command_Line & cl,
                                                 char flag,
                                                 int verbose)
{
  if (! cl.option_present (flag))
  {
    return initialise_default_weight_function();
  }

  int echo = 0;

  int i = 0;
  IWString w;
  while (cl.value (flag, w, i++))
  {
    if (! construct_from_command_line_token (w, flag))
    {
      cerr << "NN_Weight_Function::construct_from_command_line: unrecognised '" << w << "'\n";
      return 0;
    }
  }

  if (NULL == _weight)
  {
    initialise_default_weight_function();
  }

  if (echo)
    echo_weight_function (cerr);

  return 1;
}

int
NN_Weight_Function::construct_from_command_line_token (const const_IWSubstring & w,
                              char flag)
{
  if (w.starts_with("maxd="))
  {
    const_IWSubstring tmp(w);
    tmp.remove_leading_chars(5);
    float d;
    if (! tmp.numeric_value(d) || d <= 0.0)
    {
      cerr << "Invalid maximum distance directive '" << w << "'\n";
      return 0;
    }
  }
  else if ("linear" == w)
  {
    initialise_linear_weight_function();
  }
  else if ("flat" == w)
  {
    initialise_flat_weight_function();
  }
  else if ("test" == w)
  {
    initialise_test_weight_function();
  }
  else if (w.starts_with ("stair=") || w.starts_with ("step="))
  {
    const_IWSubstring myw (w);

    myw.remove_up_to_first ('=');
    weight_t step;
    if (! myw.numeric_value (step) || step <= 0.0)
    {
      cerr << "The stair ratio must be positive, '" << w << "' is invalid\n";
      return 8;
    }

    initialise_stair_weight_function (step);
  }
  else if (w.starts_with ("stair"))
  {
    initialise_stair_weight_function(0.8);
  }
  else if (w.starts_with ("slope="))
  {
    const_IWSubstring myw (w);

    myw.remove_up_to_first ('=');
    weight_t slope;
    if (! myw.numeric_value (slope) || slope <= 0.0)
    {
      cerr << "The stair ratio must be positive, '" << w << "' is invalid\n";
      return 8;
    }

    initialise_sloping_weight_function (slope);
  }
  else if (w.starts_with ("def:"))
  {
    if (! initialise_default_weight_function (w))
    {
      cerr << "Invalid weight function specifier '" << w << "'\n";
      return 0;
    }
  }
  else if ("def" == w)
  {
    initialise_default_weight_function();
  }
  else if (w.starts_with ("exp-"))
  {
    initialise_exponential_weight_function (w);
  }
  else if (w.starts_with ("file="))
  {
    const_IWSubstring fname (w);

    fname.remove_leading_chars (5);
    if (! _read_weights_from_file (fname))
    {
      cerr << "Cannot read weights from '" << w << "'\n";
      return 0;
    }
  }
  else if ("help" == w)
  {
    display_weight_function_options (cerr, flag);
  }
  else
  {
    cerr << "Unrecognised weight function directive '" << w << "'\n";
    return 0;
  }

  return 1;
}

void
NN_Weight_Function::echo_weight_function (std::ostream & os) const
{
  os << "Weight function " << _type << endl;

  for (int i = 0; i < number_points_in_weight_array; i++)
  {
    os << (i * distance_delta) << ' ' << _weight[i] << endl;
  }

  return;
}

static double switchover = 0.30;

int
NN_Weight_Function::initialise_default_weight_function()
{
  _allocate_weight_array_if_needed();

  double md = maximum_possible_distance();

  double delta = md / static_cast<double> (number_points_in_weight_array - 1);

  for (int i = 0; i < number_points_in_weight_array; i++)
  {
    double distvec = i * delta;

    double p = (switchover - distvec) / switchover;
    if (p < 0.0)
      p = 0.0;

    if (distvec < md)
      _weight[i] = p * pow (1.0 - distvec, 1.0 / 6.0) + (1.0 - p) * pow (1.0 - distvec, 3.0);
    else
      _weight[i] = 1.0e-08;
  }

  _type = "def";

  return 1;
}

/*
  The weight function specifier must look like

  def:a=xx,b=yy
*/

int
NN_Weight_Function::initialise_default_weight_function (const const_IWSubstring & w)
{
#ifdef USE_REGEXP
  std::unique_ptr<re2::RE2> rx = std::make_unique<re2::RE2>("^def:a=(.+),b=(.+)");

#else
  const_IWSubstring myw (w);

  myw.remove_leading_chars (6);    // get rid of 'def:a='

  const_IWSubstring sa, sb;
  myw.split (sa, ',', sb);

  assert (sb.starts_with ("b="));
  sb.remove_leading_chars (2);

  double a;
  if (! sa.numeric_value (a))
  {
    cerr << "NN_Weight_Function::initialise_default_weight_function: invalid 'a' specifier '" << w << "' '" << sa << "'\n";
    return 0;
  }

  if (a >= 0.0)
    cerr << "Very strange, 'a' values are usually negative\n";

  sb.strip_trailing_blanks();

  double b;
  if (! sb.numeric_value (b))
  {
    cerr << "NN_Weight_Function::initialise_default_weight_function: invalid 'b' specifier '" << w << "' '" << sb << "'\n";
    return 0;
  }
#endif

  _allocate_weight_array_if_needed();

  _weight[0] = 1.0;

  double log10b = log10 (b);

  double md = maximum_possible_distance();

  double delta = md / static_cast<double> (number_points_in_weight_array - 1);

  for (int i = 1; i < number_points_in_weight_array; i++)
  {
    double d = i * delta;

    double tmp = a * (log10b - log10 (d));

    _weight[i] = 1.0 / (1.0 + pow (10.0, tmp));
  }

  _type = w;

  return 1;
}

/*
  Must look like exp-nn
*/

int
NN_Weight_Function::initialise_exponential_weight_function (const const_IWSubstring & w)
{
  assert (w.starts_with ("exp-"));

  const_IWSubstring myw (w);
  myw.remove_leading_chars (4);

  double efactor;

  if (! myw.numeric_value (efactor))
  {
    cerr << "NN_Weight_Function::initialise_exponential_weight_function: invalid exponent factor '" << w << "'\n";
    return 0;
  }

  _allocate_weight_array_if_needed();

  double md = maximum_possible_distance();

  double delta = md / static_cast<double> (number_points_in_weight_array - 1);

  for (int i = 0; i < number_points_in_weight_array; i++)
  {
    double d = i * delta;

    _weight[i] = exp (- d * efactor);
  }

  _type = w;

  return 1;
}

int
NN_Weight_Function::initialise_linear_weight_function()
{
  _allocate_weight_array_if_needed();

  double md = maximum_possible_distance();

  double delta = md / static_cast<double> (number_points_in_weight_array - 1);

  for (int i = 0; i < number_points_in_weight_array; i++)
  {
    _weight[i] = 1.0 - i * delta;
  }

  _type = "linear";

  return 1;
}

int
NN_Weight_Function::initialise_flat_weight_function()
{
  _allocate_weight_array_if_needed();

  set_vector (_weight, number_points_in_weight_array, static_cast<weight_t> (0.5));

  _type = "flat";

  return 1;
}

int
NN_Weight_Function::initialise_test_weight_function()
{
  _allocate_weight_array_if_needed();

  _weight[0] = 1.0;
  for (int i = 1; i < number_points_in_weight_array; i++)
  {
    _weight[i] = 0.9 * _weight[i - 1];
  }

  return 1;
}

int
NN_Weight_Function::initialise_stair_weight_function (weight_t ratio)
{
  _allocate_weight_array_if_needed();

  _weight[0] = 1.0;
  for (int i = 1; i < number_points_in_weight_array; i++)
  {
    _weight[i] = ratio * _weight[i - 1];
  }

  _type = "stair=";
  _type << ratio;

  return 1;
}

int
NN_Weight_Function::_read_weights_from_file (const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.ok())
  {
    cerr << "NN_Weight_Function::_read_weights_from_file: cannot open '" << fname << "'\n";
    return 0;
  }

  if (! _read_weights_from_file (input))
    return 0;

  _type << "File=" << fname;

  return 1;
}

int
NN_Weight_Function::_read_weights_from_file (iwstring_data_source & input)
{
  if (NULL == _weight)
    _weight = new weight_t[number_points_in_weight_array];

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    weight_t w;
    if (! buffer.numeric_value (w))
    {
      cerr << "Invalid weight value '" << buffer << "' on line " << input.lines_read() << endl;
      return 0;
    }

    _weight[input.lines_read()] = w;

    if (input.lines_read() >= number_points_in_weight_array)
      break;
  }

// All the other weights are 0.0

  for (int i = input.lines_read(); i < number_points_in_weight_array; i++)
  {
    _weight[i] = 0.0;
  }

  return 1;
}

int
NN_Weight_Function::initialise_sloping_weight_function (weight_t slope)
{
  _allocate_weight_array_if_needed();

  if (slope < 0.0)
    slope = -slope;

  weight_t delta = slope / static_cast<weight_t> (number_points_in_weight_array - 1);

  _weight[0] = 1.0;

  _type =  "slope=";
  _type << slope;

  for (int i = 1; i < number_points_in_weight_array; i++)
  {
    _weight[i] = 1.0 - static_cast<weight_t> (i) * delta;

    if (_weight[i] > 0.0)
      continue;

//  We don't allow zero or negative weights

    for (int j = i; j < number_points_in_weight_array; j++)
    {
      _weight[j] = _weight[i - 1];
    }
    return 1;
  }

  return 1;
}
