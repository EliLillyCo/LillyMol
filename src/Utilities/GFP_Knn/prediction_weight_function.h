#ifndef PREDWF_H
#define PREDWF_H 2

#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"

#include "nn_stuff.h"

/*
  Since we may be dealing with distances from arbitrary sources, we need
  a means of specifying the largest distance possible. By default, it is 1.0
*/

extern double maximum_distance_possible ();

/*
  Once we have initialised all the distance conditions, we can compute
  the index into a weight function array of any distance
*/

extern int compute_weight_function_index (double, int &);

extern int initialise_distances (Command_Line & cl, char flag, int verbose);

/*
  We impute activity to each target based on:
    the activity of each neighbout
    the distance from that neighbour
    a weighting function

  We compute the weight function at 0.01 increments only. That way, it can be of
  arbitrary complexity
*/

class NN_Weight_Function
{
  private:
    weight_t * _weight;

    IWString _type;

// private functions
    
    int _read_weights_from_file (iwstring_data_source &);
    int _read_weights_from_file (const const_IWSubstring &);
    void _allocate_weight_array_if_needed ();

  public:
    NN_Weight_Function ();
    ~NN_Weight_Function ();

    int construct_from_command_line (Command_Line & cl, char, int);
    int construct_from_command_line_token (const const_IWSubstring &, char);

    const IWString & weight_function_type () const { return _type;}

    void echo_weight_function (std::ostream & os) const;

    int initialise_default_weight_function ();
    int initialise_default_weight_function (const const_IWSubstring &);
    int initialise_exponential_weight_function (const const_IWSubstring &);
    int initialise_linear_weight_function ();
    int initialise_flat_weight_function ();
    int initialise_test_weight_function ();
    int initialise_stair_weight_function (weight_t);
    int initialise_sloping_weight_function (weight_t);

    weight_t weight (similarity_type_t t) const
      {
        int j = static_cast<int> (t * 100.0);    // the index into the weight array - no rounding here

        return _weight[j];
      }

    weight_t weight (int i) const { return _weight[i];}
};

extern int display_weight_function_options (std::ostream & os, char = 'W');

#endif
