#ifndef IWSET_OF_DESCRIPTORS_IMPLEMENTATION_H
#define IWSET_OF_DESCRIPTORS_IMPLEMENTATION_H

#include "set_of_descriptors.h"

template <typename P>
Set_of_Descriptors<P>::Set_of_Descriptors ()
{
  return;
}

template <typename P>
Set_of_Descriptors<P>::~Set_of_Descriptors ()
{
  return;
}

template <typename P>
int
Set_of_Descriptors<P>::build (const IWString & fname,
                              IWString & header)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Set_of_Descriptors::build: cannot open '" << fname << "'\n";
    return 0;
  }

  return build (input, header);
}

template <typename P>
int
Set_of_Descriptors<P>::build (iwstring_data_source & input,
                              IWString & header)
{
  if (! input.next_record (header))
  {
    cerr << "Cannot read header record from input file\n";
    return 0;
  }

  header.remove_leading_words (1);

  if (0 == _number_elements)
  {
    int tmp = input.records_remaining ();
    if (tmp <= 0)
    {
      cerr << "Possibly empty target descriptor file\n";
      return 0;
    }

    if (! iwaray<P>::resize (tmp))
    {
      cerr << "Set_of_Descriptors::build: cannot allocate pool for " << tmp << " descriptors\n";
      return 0;
    }
  }

  return _build (input, header.nwords ());
}

template <typename P>
int
Set_of_Descriptors<P>::_build (iwstring_data_source & input,
                               int number_descriptors)
{
  const_IWSubstring buffer;
  int items_in_pool = 0;
  while (input.next_record (buffer))
  {
    if (items_in_pool == _number_elements)
    {
      cerr << "Pool is full " << _number_elements << endl;
      return 1;
    }

    P & di = _things[items_in_pool];
    int fatal;
    if (! di.build (buffer, number_descriptors, fatal))
    {
      cerr << "Invalid descriptor record '" << buffer << "'\n";
      return 0;
    }

    items_in_pool++;
  }

  _number_elements = items_in_pool;

  return items_in_pool;
}


#endif
