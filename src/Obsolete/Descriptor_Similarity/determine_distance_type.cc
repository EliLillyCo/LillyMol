#include <stdlib.h>
#include <iostream>

#include "iwdescriptor.h"

using std::cerr;

int
display_distance_metrics (std::ostream & output,
                          char c)
{
  output << " -" << c << " <type>        type of distance computation\n";
  output << " -" << c << " cartesian     cartesian distance\n";
  output << " -" << c << " manhattan     manhattan distance\n";
  output << " -" << c << " taxonomic     taxonomic distance\n";
  output << " -" << c << " canberra      canberra  distance\n";
  output << " -" << c << " cosine        cosine    distance\n";
  output << " -" << c << " iwratio       IW ratio  distance\n";
  output << " -" << c << " braycurtiss   braycurtiss distance\n";
  output << " -" << c << " minkinf       minkowski Infinite distance\n";
  output << " -" << c << " censeuc       censored euclidian istance\n";
  output << " -" << c << " tanimoto      Tanimoto - adapted to continuous descriptors\n";
  output << " -" << c << " ctan          continuous tanimoto\n";
  output << " -" << c << " rsrao         Russel-Rao metric\n";
  output << " -" << c << " circle        Circle Product distance (just sum of mins)\n";
  output << " -" << c << " czekan        Czekanowski metric\n";
  output << " -" << c << " pearson       Pearson correlation (use -O 1)\n";
  output << " -" << c << " dot           Dot product\n";

  return output.good ();
}
int
determine_distance_type (const const_IWSubstring & x,
                         int verbose)
{
  int distance_computation = 0;

  if (x.starts_with ("cart"))
  {
    distance_computation = DISTANCE_CARTESIAN;
    if (verbose)
      cerr << "The cartesian distance will be used\n";
  }
  else if (x.starts_with ("manh"))
  {
    distance_computation = DISTANCE_MANHATTAN;
    if (verbose)
      cerr << "The manhattan distance will be used\n";
  }
  else if (x.starts_with ("taxon"))
  {
    distance_computation = DISTANCE_TAXONOMIC;
    if (verbose)
      cerr << "The taxonomic distance will be used\n";
  }
  else if (x.starts_with ("canbe"))
  {
    distance_computation = DISTANCE_CANBERRA;
    if (verbose)
      cerr << "The canberra distance will be used\n";
  }
  else if (x.starts_with ("cosine"))
  {
    distance_computation = DISTANCE_COSINE;
    if (verbose)
      cerr << "The cosine distance will be used\n";
  }
  else if (x.starts_with ("minkinf"))
  {
    distance_computation = DISTANCE_MINKINF;
    if (verbose)
      cerr << "The minkowski infinite distance will be used\n";
  }
  else if (x.starts_with ("braycurtis"))
  {
    distance_computation = DISTANCE_BRAY_CURTIS;
    if (verbose)
      cerr << "The braycurtiss infinite distance will be used\n";
  }
  else if (x.starts_with ("censeuc"))
  {
    distance_computation = DISTANCE_CENSORED_EUCLIDIAN;
    if (verbose)
      cerr << "The censored euclidian distance will be used\n";
  }
  else if (x.starts_with ("iwratio"))
  {
    distance_computation = DISTANCE_IW_RATIO;
    if (verbose)
      cerr << "The iw ratio distance will be used\n";
  }
  else if (x.starts_with ("tanimoto"))
  {
    distance_computation = DISTANCE_TANIMOTO;
    if (verbose)
      cerr << "The Tanimoto distance will be used\n";
  }
  else if (x.starts_with ("rrao") || x.starts_with ("rsrao"))
  {
    distance_computation = DISTANCE_RUSSEL_RAO;
    if (verbose)
      cerr << "The Russel Rao distance will be used\n";
  }
  else if (x.starts_with ("ctan"))
  {
    distance_computation = DISTANCE_CONTINUOUS_TANIMOTO;
    if (verbose)
      cerr << "The continuous Tanimoto distance will be used\n";

//    negative_distances_allowed = 1;
  }
  else if (x.starts_with ("forbes"))
  {
    distance_computation = DISTANCE_FORBES;
    if (verbose)
      cerr << "The Forbes distance will be used\n";
  }
  else if (x.starts_with ("cforbes"))
  {
    distance_computation = DISTANCE_CONTINUOUS_FORBES;
    if (verbose)
      cerr << "The Forbes continuous distance will be used\n";
  }
  else if (x.starts_with ("circle"))
  {
    distance_computation = DISTANCE_CIRCLE_PRODUCT;
    if (verbose)
      cerr << "The circle product distance will be used\n";
  }
  else if (x.starts_with ("czekan"))
  {
    distance_computation = DISTANCE_CZEKANOWSKI;
    if (verbose)
      cerr << "The Czekanowski distance will be used\n";
  }
  else if (x.starts_with ("pearson"))
  {
    distance_computation = DISTANCE_PEARSON;
    if (verbose)
      cerr << "The Pearson correlaton will be used\n";
  }
  else if (x.starts_with ("dot"))
  {
    distance_computation = DISTANCE_DOT_PRODUCT;
    if (verbose)
      cerr << "The dot product will be used\n";
  }
  else
  {
    cerr << "Unrecognised distance computation '" << x << "'\n";
    return 0;
  }

  return distance_computation;
}


int
needs_ave_and_variance(int d)
{
  if (DISTANCE_PEARSON == d)
    return 1;

  return 0;
}
