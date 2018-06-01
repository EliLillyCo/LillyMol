#include <stdlib.h>

#include "misc.h"

#include "molecule.h"


Symmetry_Class_and_Canonical_Rank::Symmetry_Class_and_Canonical_Rank ()
{
  _canonical_rank = NULL;
  _symmetry_class = NULL;

  return;
}

Symmetry_Class_and_Canonical_Rank::~Symmetry_Class_and_Canonical_Rank ()
{
  invalidate();

  return;
}

int
Symmetry_Class_and_Canonical_Rank::invalidate ()
{
  DELETE_IF_NOT_NULL_ARRAY(_canonical_rank);
  DELETE_IF_NOT_NULL_ARRAY(_symmetry_class);

  return 1;
}

static int
common_allocate_array (int * & p,
                       int s)
{
  if (NULL == p)
    delete [] p;

  p = new int[s];

  return (NULL != p);
}

int
Symmetry_Class_and_Canonical_Rank::allocate_arrays (int s)
{
  assert (s > 0);

  int rc = common_allocate_array (_canonical_rank, s);
  if (0 != rc)
    rc = common_allocate_array (_symmetry_class, s);

  return rc;
}

static int
common_copy (int * & lhs,
             const int * rhs,
             int n)
{
  if (NULL == lhs && NULL == rhs)
    return 1;

  if (NULL != lhs && NULL == rhs)
  {
    delete lhs;
    lhs = NULL;
    return 1;
  }

// At this stage rhs is NOT null

  assert (NULL != rhs);

  if (NULL == lhs)
  {
    lhs = new int[n];
    if (NULL == lhs)
    {
      cerr << "common_copy:memory failure, cannot allocate " << n << " items\n";
      return 0;
    }
  }

  copy_vector (lhs, rhs, n);

  return 1;
}

int
Symmetry_Class_and_Canonical_Rank::store_values_from (const Symmetry_Class_and_Canonical_Rank & rhs,
                                                      int natoms)
{
  common_copy (_symmetry_class, rhs._symmetry_class, natoms);
  common_copy (_canonical_rank, rhs._canonical_rank, natoms);

  return 1;
}
