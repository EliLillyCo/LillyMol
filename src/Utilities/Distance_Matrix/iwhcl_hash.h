#ifndef IWHCL_HASH_H
#define IWHCL_HASH_H

#include <unordered_map>
#include "iw_stl_hash_map.h"

/*
  due to massive problems with 8 byte int's we kludge things
*/

typedef unsigned long long iwuint64;
//typedef unsigned int iwuint64;

class Pair_of_Uints
{
  private:
    unsigned int _multiplier;
    unsigned int _n1;
    unsigned int _n2;

  public:
    Pair_of_Uints(unsigned int m) : _multiplier(m) {};

    void set (unsigned int u1, unsigned int u2) { _n1 = u1; _n2 = u2;}

    unsigned int n1 () const { return _n1;}
    unsigned int n2 () const { return _n2;}
    unsigned int multiplier() const { return _multiplier;}

    int operator==(const Pair_of_Uints &) const;
};

class Pair_of_Uints_Hash_Fn
{
  private:

  public:
    size_t operator() (const Pair_of_Uints &) const;
};

class IW64bithash
{
  private:
  public:

#if defined (__INTEL_COMPILER)
    const size_t bucket_size = 4;
    const size_t min_buckets = 8;
    bool  operator () (const iwuint64 &, const iwuint64 &) const;
#endif

    size_t operator () (const iwuint64 &) const;
};

//typedef hash_map<iwuint64, unsigned char> Dist_Hash;

int
iwhcl_hash (const unsigned int n,
       const unsigned int multiplier,
       const int iopt,
       int * ia,
       int * ib,
       float * crit,
       float * membr,
       unsigned int * nn,
       unsigned char * disnn,
       int * flag,
       std::unordered_map<Pair_of_Uints, unsigned char, Pair_of_Uints_Hash_Fn> & diss);

#endif
