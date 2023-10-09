#include <iostream>

/*
  All these programmes need to build an in-memory pool

  Note that this only works properly reading one pool file. Sometime I
  may fix it to allow the pool to be built from multiple files
*/

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"

using std::cerr;

int
build_pool (iwstring_data_source & input,
            IW_General_Fingerprint * pool,
            int max_pool_size,
            int & pool_size)
{
  assert (max_pool_size > 0);
  assert (pool_size >= 0);

  if (pool_size >= max_pool_size)
  {
    cerr << "Build pool, pool already full, max size " << max_pool_size << " current size " << pool_size << '\n';
    return 0;
  }

  off_t offset = input.tellg ();

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    if (! pool[pool_size].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg ();
      continue;
    }

    pool[pool_size].set_offset (offset);

    pool_size++;

    if (pool_size >= max_pool_size)
    {
      cerr << "Pool is full, max " << max_pool_size << '\n';
      return 1;
    }

    offset = input.tellg ();
  }

  return 1;
}

int
build_pool (const const_IWSubstring & fname,
            IW_General_Fingerprint * pool,
            int max_pool_size,
            int & pool_size,
            const IWString & identifier_tag)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open pool file '" << fname << "'\n";
    return 0;
  }

  if (0 == max_pool_size)  // Need to determine size
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    re2::StringPiece tmpsp(tmp.data(), tmp.length());
    RE2 pcn(tmp);
    max_pool_size = input.grep(pcn);

    if (0 == max_pool_size)
    {
      cerr << "No occurrences of " << pcn.pattern() << "' in input\n";
      return 0;
    }

    pool = new IW_General_Fingerprint[max_pool_size];
    if (nullptr == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << max_pool_size << '\n';
      return 62;
    }

    cerr << "Pool automatically sized to " << max_pool_size << '\n';
  }

  return build_pool(input, pool, max_pool_size, pool_size);
}
