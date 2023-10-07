#ifndef QRY_AND_DEMERIT_H
#define QRY_AND_DEMERIT_H


#include "Molecule_Lib/substructure.h"

class IWString;
class Demerit;
class Molecule_to_Match;

class Query_and_Demerit_Value: public Substructure_Query
{
  private:
    IWString _description;
    int _reject;
    int _demerit;
    int _demerit_each_occurrence;

    int _molecules_examined;
    int _molecules_demerited;

//  private functions

    void _default_values ();

  public:
    Query_and_Demerit_Value ();
    Query_and_Demerit_Value (const const_IWSubstring &);

    int debug_print (std::ostream &) const;
    int ok () const;

    int set_rejection (int);
    int determine_action (const IWString &);

    int evaluate (Molecule_to_Match &, Demerit &);
};

#endif
