#include <stdlib.h>
#include <iostream>

#include "Foundational/iwmisc/misc.h"

#include "molecule.h"

using std::cerr;
using std::endl;

Symmetry_Class_and_Canonical_Rank::Symmetry_Class_and_Canonical_Rank()
{
  _canonical_rank = nullptr;
  _symmetry_class = nullptr;

  return;
}

Symmetry_Class_and_Canonical_Rank::~Symmetry_Class_and_Canonical_Rank()
{
  invalidate();

  return;
}

int
Symmetry_Class_and_Canonical_Rank::invalidate()
{
  DELETE_IF_NOT_NULL_ARRAY(_canonical_rank);
  DELETE_IF_NOT_NULL_ARRAY(_symmetry_class);

  return 1;
}

static int
common_allocate_array(int*& p, int s)
{
  if (nullptr == p)
    delete[] p;

  p = new int[s];

  return (nullptr != p);
}

int
Symmetry_Class_and_Canonical_Rank::allocate_arrays(int s)
{
  assert(s > 0);

  int rc = common_allocate_array(_canonical_rank, s);
  if (0 != rc)
    rc = common_allocate_array(_symmetry_class, s);

  return rc;
}

static int
common_copy(int*& lhs, const int* rhs, int n)
{
  if (nullptr == lhs && nullptr == rhs)
    return 1;

  if (nullptr != lhs && nullptr == rhs)
  {
    delete lhs;
    lhs = nullptr;
    return 1;
  }

  // At this stage rhs is NOT null

  assert(nullptr != rhs);

  if (nullptr == lhs)
  {
    lhs = new int[n];
    if (nullptr == lhs)
    {
      cerr << "common_copy:memory failure, cannot allocate " << n << " items\n";
      return 0;
    }
  }

  copy_vector(lhs, rhs, n);

  return 1;
}

int
Symmetry_Class_and_Canonical_Rank::store_values_from(const Symmetry_Class_and_Canonical_Rank& rhs,
                                                     int natoms)
{
  common_copy(_symmetry_class, rhs._symmetry_class, natoms);
  common_copy(_canonical_rank, rhs._canonical_rank, natoms);

  return 1;
}
