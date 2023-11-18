#ifndef IWREAD_POOL_H
#define IWREAD_POOL_H

#include "iwstring_data_source.h"

template <typename P>
int
read_pool (P * pool,
           int pool_size,
           iwstring_data_source & input,
           int number_descriptors)
{
  const_IWSubstring buffer;
  int items_in_pool = 0;
  while (input.next_record (buffer))
  {
    if (items_in_pool == pool_size)
    {
      cerr << "Pool is full " << pool_size << endl;
      return 1;
    }

    P & di = pool[items_in_pool];
    int fatal;
    if (! di.build (buffer, number_descriptors, fatal))
    {
      cerr << "Invalid descriptor record '" << buffer << "'\n";
      return 0;
    }

    items_in_pool++;
  }

  pool_size = items_in_pool;

  return items_in_pool;
}

static int
read_pool (P * & pool,
           int & pool_size,
           iwstring_data_source & input,
           IWString & header)
{
  if (! input.next_record (header))
  {
    cerr << "Cannot read header record from input file\n";
    return 0;
  }

  header.remove_leading_words (1);

  if (0 == pool_size)
  {
    pool_size = input.records_remaining ();
    if (pool_size <= 0)
    {
      cerr << "Possibly empty target descriptor file\n";
      return 0;
    }

    pool = new Descriptors_and_Neighbours[pool_size];

    assert (NULL != pool);

    if (verbose)
      cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return read_pool (input, header.nwords ());
}


static int
read_pool (IWString & fname,
           P * & pool,
           int & pool_size,
           IWString & header)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_pool (input, header);
}

#endif
