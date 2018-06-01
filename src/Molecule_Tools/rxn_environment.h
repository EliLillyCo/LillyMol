#ifndef RXN_ENVIRONMENT_H
#define RXN_ENVIRONMENT_H

#include <unordered_map>

#include "db_cxx.h"
#include "sys/types.h"
#include "sys/stat.h"

#include "accumulator.h"

#include "rxn_file.h"

class RXN_Environment_Assessment
{
  private:
    Accumulator<double> _acc_fraction;

    float _fraction_found;    // fraction of bonds examined that are in the set

    int _found;
    int _notfound;

//  private functions

    void _default_values();

  public:
    RXN_Environment_Assessment();

    void extra_fraction(const float s) { _acc_fraction.extra(s);}

    void set_found(const int f, const int n) { _found = f; _notfound = n;}

    void reset();

    template <typename T> int do_print(T & output) const;
};

class RXN_Environment
{
  private:
    int _max_radius;
    int _nreactions;

    std::unordered_map<uint64_t, int> * _found;

//  private functions

    int _do_store(Db & db, const int radius, const std::pair<uint64_t, int> & f) const;

    int _fetch_count(const uint64_t x, int & r) const;

    uint64_t _form_hash(Molecule & m, const Bond * b, const int * atype, const int * changed, int & d)  const;

    template <typename T> int _do_store_value(Db & db, const char * key, const T & v) const;

  public:
    RXN_Environment();
    ~RXN_Environment();

    int summary(std::ostream & output) const;

    int initialise(const int s);

    int nreactions() const { return _nreactions;}

    int gather(ISIS_RXN_FILE_Molecule & m, const int * atype, const int * changed);

    int do_store(Db & db) const;

    int do_read(Db & db);

    int assess(Molecule & m, const int * atype, const int * changed,    RXN_Environment_Assessment & rxnenva) const;
    int assess(Molecule & m, const int * atype, const Set_of_Atoms * e, RXN_Environment_Assessment & rxnenva) const;
};

#endif
