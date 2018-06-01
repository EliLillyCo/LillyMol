#ifndef IWRANDOM_H
#define IWRANDOM_H

#include <stdlib.h>

#include "mtrand.h"

typedef unsigned long int random_number_seed_t;
typedef double random_number_t;

class Random_Number_Working_Storage : public MTRand_int32
{
  private:

  public:
    Random_Number_Working_Storage ();
    Random_Number_Working_Storage (random_number_seed_t seed) : MTRand_int32(seed) {}

//  int initialised () const { return _initialised;}

    random_number_seed_t choose_random_seed ();

//  int debug_print (ostream & = cerr) const;

    void set_seed (random_number_seed_t s) { MTRand_int32::seed(s);}

//  unsigned short * buffer () const { return (unsigned short *) & _buffer;}

    int random_one_or_zero () { return static_cast<int>(MTRand_int32::closed_open() + 0.50);}

    random_number_t random_number () { return MTRand_int32::closed_closed();}

    int intbtwij (int low, int high) { return static_cast<int>(low + MTRand_int32::closed_open() * static_cast<double>(high - low + 1));}
};


extern random_number_t iwrandom(void);
extern int intbtwij (int, int);
extern void iw_set_rnum_seed(random_number_seed_t);
extern random_number_seed_t iw_random_seed ();    // also sets the global default
extern random_number_seed_t random_seed_based_on_time_and_pid();
extern random_number_seed_t random_seed_from_dev_random ();

template <typename T> T random_number_between (const T & low, const T & high);

class IW_Box_Muller_Normally_Distributed : public MTRand_int32
{
  private:
    double _n1;
    int    _n2_valid;
    double _n2;
  public:
    IW_Box_Muller_Normally_Distributed();
    IW_Box_Muller_Normally_Distributed(random_number_seed_t);

    double operator () ();
};

class IW_Inverse_Normal_Normally_Distributed : public MTRand_int32
{
  private:
  public:
    IW_Inverse_Normal_Normally_Distributed();

    double operator () ();
};

#ifdef RANDOM_NUMBER_BETWEEN_IMPLEMENTATION

template <typename T>
T
random_number_between (const T & low, const T & high)
{
  double range = static_cast<double> (high - low + 1);

  T delta = static_cast<T>(range * iwrandom());

  return low + delta;
}

#endif

#endif
