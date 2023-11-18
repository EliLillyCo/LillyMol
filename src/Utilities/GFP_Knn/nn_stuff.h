#ifndef NN_STUFF_H
#define NN_STUFF_H

#include <iostream>

class Command_Line;

#include "Foundational/iwbits/iwbits.h"

extern int
parse_dash_n (Command_Line & cl,
              char flag,
              resizable_array<int> & neighbours_to_use,
              int verbose);

extern int 
parse_dash_t (Command_Line & cl, 
              char flag,
              resizable_array<float> & distances_to_use,
              resizable_array<int> & min_neighbours_required,
              double max_distance,
              int verbose);

extern int display_standard_nn_stuff (std::ostream &);

extern similarity_type_t maximum_possible_distance ();
extern void set_maximum_possible_distance (similarity_type_t);

typedef float weight_t;

#endif
