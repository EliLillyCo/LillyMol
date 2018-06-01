#ifndef COLLECTION_TEMPLATE_H
#define COLLECTION_TEMPLATE_H

/*
  There are a couple of instances where we need a collection of
  objects as well as a text description of what the items are.
  The first two examples are partial charges, and atom types
*/

#include "iwstring.h"

template <typename T>
class Collection_Template : public resizable_array<T>
{
  private:
    IWString _type;

  public:

    Collection_Template<T> & operator = (const Collection_Template<T> & rhs)
      {
        resizable_array<T>::operator= (rhs);
        _type = rhs._type;

        return *this;
      }

    IWString & ztype () { return _type;}
    const IWString & ztype () const { return _type;}

    void set_type (const IWString & t) { _type = t;}
};

#endif
