#include <stdlib.h>

#include "re2/re2.h"

#include "Foundational/data_source/iwstring_data_source.h"

#include "misc.h"

int
find_dividing_points(const char * fname, 
                     int n,
                     resizable_array<off_t> & offset,
                     re2::RE2& rx)
{
  assert (n > 1);

  offset.resize(0);

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "find_dividing_points:cannot open '" << fname << "'\n";
    return 0;
  }

  off_t file_size = input.file_size();

  if (file_size < n)
  {
    cerr << "find_dividing_points:cannot divide file '" << fname << "' of sie " << file_size << " into " << n << " sections\n";
    return 0;
  }

  off_t bytes_per_chunk = file_size / n;

  offset.add(static_cast<off_t>(0));

  for (int i = 1; i < n; i++)
  {
    off_t o = i * bytes_per_chunk;
    if (! input.seekg(o))
    {
      cerr << "find_dividing_points:cannot seek to " << o << " in '" << fname << "'\n";
      return 0;
    }

    const_IWSubstring buffer;
    input.next_record(buffer);
    o = input.tellg();

    while (input.next_record(buffer))
    {
      const re2::StringPiece tmp(buffer.data(), buffer.length());
      if (RE2::PartialMatch(tmp, rx))
      {
        offset.add(o);
        break;
      }

      o = input.tellg();
    }

    if (offset.number_elements() != i + 1)
    {
      cerr << "find_dividing_points:eof\n";
      return 0;
    }
  }

  return n;
}
