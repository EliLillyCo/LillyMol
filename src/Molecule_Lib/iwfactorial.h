#ifndef MOLECULE_LIB_IWFACTORIAL_H_
#define MOLECULE_LIB_IWFACTORIAL_H_

#include "iwmtypes.h"

template <typename T>
class IW_Factorial
{
  private:
    int _size;
    T * _f;

//  private functions

    void _default_values ();

  public:
    IW_Factorial ();
    IW_Factorial (int);

    ~IW_Factorial ();

    int resize (int);

    int max_index () const { return _size + 1;}

    inline T operator [] (int i) const { return _f[i];}
};

#endif  // MOLECULE_LIB_IWFACTORIAL_H_
