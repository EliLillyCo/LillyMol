#ifndef DISTANCE_SCALING_H
#define DISTANCE_SCALING_H 1

#include "Foundational/data_source/iwstring_data_source.h"

/*
  We want to be able to convert from some other fingerprint
  distances back to a familiar scale.
*/

class Distance_Scaling
{
  private:
    float * _conv;

// private functions

    int _parse_scaling_record (const const_IWSubstring & buffer);

  public:
    Distance_Scaling();
    ~Distance_Scaling();

    int build (const char *);
    int build (iwstring_data_source &);

    int active() const { return nullptr != _conv;}

    float convert (float) const;
};

#endif
