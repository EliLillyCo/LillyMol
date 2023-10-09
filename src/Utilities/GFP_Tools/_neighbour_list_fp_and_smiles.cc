#include <stdlib.h>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "gfp.h"

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

/*
  We can save time by pre-computing the string representations of distances
*/

#ifdef NO_LONGER_USED____
static IWString * string_distance = nullptr;

int
initialise_string_distances ()
{
  string_distance = new IWString[1001];
  if (nullptr == string_distance)
  {
    cerr << "Yipes, cannot allocate 1000 strings for distances\n";
    return 0;
  }

  for (int i = 0; i < 1001; i++)
  {
    double f = 0.001 * i;

    string_distance[i].append_number (f, 3);
  }

  string_distance[1000] = "1";

  return 1;
}

static void
append_string_distance (similarity_type_t d,
                        IWString & output_buffer)
{
  int i = static_cast<int> (d * static_cast<float> (1000.0) + static_cast<float> (0.499999));

  output_buffer << string_distance[i];

  return;
}
#endif

#define NEIGHBOUR_LIST_IMPLEMENTATION

#include "neighbour_list.h"

template class Neighbour_List<float, FP_and_Smiles, FP_and_Smiles>;
