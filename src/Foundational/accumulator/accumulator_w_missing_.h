#ifndef ACCUM_WMISSING_IMP
#define ACCUM_WMISSING_IMP

#include "accumulator.h"

template <typename T>
Accumulator_with_Missing_Values<T>::Accumulator_with_Missing_Values ()
{
  _nmissing = 0;

  return;
}

template <typename T>
void
Accumulator_with_Missing_Values<T>::reset ()
{
  Accumulator<T>::reset ();

  _nmissing = 0;

  return;
}

#endif
