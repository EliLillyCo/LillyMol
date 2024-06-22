#ifndef MOLECULE_LIB_TOKENISE_ATOMIC_SMARTS_H_
#define MOLECULE_LIB_TOKENISE_ATOMIC_SMARTS_H_

#include "Foundational/iwstring/iwstring.h"

class Atomic_Smarts_Component : public IWString
{
  private:
    int _unary_operator;
    int _op;              // the operator following the token
    Atomic_Smarts_Component * _next;

    // If this has been converted from a primitive to an environment,
    // this flag will be set.
    int _from_primitive;

//  private functions

    int _parse(const_IWSubstring &);

  public:
    Atomic_Smarts_Component();
    ~Atomic_Smarts_Component();

    int ok() const;
    int debug_print(std::ostream &) const;

    Atomic_Smarts_Component * next() const { return _next;}

    int op() const { return _op;}
    int unary_operator() const { return _unary_operator;}

    int from_primitive() const {
      return _from_primitive;
    }

    int parse(const_IWSubstring);

    // Parsing a mixed smarts [N,$(n)], conver all tokens to $( form.
    int ConvertToEnvironment();
};

extern std::ostream & operator << (std::ostream &, const Atomic_Smarts_Component &);

#endif  // MOLECULE_LIB_TOKENISE_ATOMIC_SMARTS_H_
