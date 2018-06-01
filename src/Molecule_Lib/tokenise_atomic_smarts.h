#ifndef ASMARTS_COMPONENT_H
#define ASMARTS_COMPONENT_H

#include "iwstring.h"

class Atomic_Smarts_Component : public const_IWSubstring
{
  private:
    int _unary_operator;
    int _op;              // the operator following the token
    Atomic_Smarts_Component * _next;

//  private functions

    int _parse (const_IWSubstring &);

  public:
    Atomic_Smarts_Component ();
    ~Atomic_Smarts_Component ();

    int ok () const;
    int debug_print (std::ostream &) const;

    Atomic_Smarts_Component * next () const { return _next;}

    int op () const { return _op;}
    int unary_operator () const { return _unary_operator;}

    int parse (const_IWSubstring);
};

std::ostream & operator << (std::ostream &, const Atomic_Smarts_Component &);


#endif
