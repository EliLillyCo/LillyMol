#ifndef IW_SET_OF_DESCRIPTORS_H
#define IW_SET_OF_DESCRIPTORS_H

// A set of descriptors is as read from a file

#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"

template <typename D>
class Set_of_Descriptors : public iwaray<D>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using iwaray<D>::_number_elements;
    using iwaray<D>::_things;
#endif

  protected:

//  private functions

    int _build (iwstring_data_source &, int);

  public:
    Set_of_Descriptors();
    ~Set_of_Descriptors();

    int build (iwstring_data_source & fname, IWString & header); 
    int build (const IWString & fname, IWString & header);

    template <typename T> int each (const T &);
};

#ifdef SET_OF_DESCRIPTORS_IMPLEMENTATION

template <typename P>
Set_of_Descriptors<P>::Set_of_Descriptors()
{
  return;
}

template <typename P>
Set_of_Descriptors<P>::~Set_of_Descriptors()
{
  return;
}

template <typename P>
int
Set_of_Descriptors<P>::build (const IWString & fname,
                              IWString & header)
{
  iwstring_data_source input (fname);
  if (! input.ok())
  {
    std::cerr << "Set_of_Descriptors::build: cannot open '" << fname << "'\n";
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
    std::cerr << "Cannot read header record from input file\n";
    return 0;
  }

  header.remove_leading_words(1);

  if (0 == _number_elements)
  {
    int tmp = input.records_remaining();
    if (tmp <= 0)
    {
      std::cerr << "Possibly empty target descriptor file\n";
      return 0;
    }

    if (! iwaray<P>::resize(tmp))
    {
      std::cerr << "Set_of_Descriptors::build: cannot allocate pool for " << tmp << " descriptors\n";
      return 0;
    }
  }

  return _build (input, header.nwords());
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
      std::cerr << "Pool is full " << _number_elements << '\n';
      return 1;
    }

    P & di = _things[items_in_pool];
    int fatal;
    if (! di.build (buffer, number_descriptors, fatal))
    {
      std::cerr << "Invalid descriptor record '" << buffer << "'\n";
      if (fatal)
        return 0;
      continue;
    }

    items_in_pool++;
  }

  _number_elements = items_in_pool;

  return items_in_pool;
}

template <typename P> template <typename T>
int
Set_of_Descriptors<P>::each (const T & t)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    rc += t(_things[i]);
  }

  return rc;

}
#endif

#endif
