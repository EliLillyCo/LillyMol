#ifndef SMI_ID_D_H
#define SMI_ID_D_H

#include "gfp.h"

/*
  Each molecule will have a list of Smiles_ID_Dist items that are its neighbours
*/

class Smiles_ID_Dist
{
  protected:
    IWString _smiles;
    IWString _id;
    similarity_type_t _distance;

  public:
    Smiles_ID_Dist ();
    Smiles_ID_Dist (const IWString &, const IWString &, similarity_type_t);
    Smiles_ID_Dist (const const_IWSubstring &, const const_IWSubstring &, similarity_type_t);

    Smiles_ID_Dist & operator = (const Smiles_ID_Dist &);

    int build (iwstring_data_source & input, int & fatal);

    const IWString & smiles () const { return _smiles;}
    void  set_smiles        (const IWString & s)  { _smiles   = s;}

    const IWString & id () const { return _id;}
    void  set_id        (const IWString & s)  { _id       = s;}

    similarity_type_t distance () const { return _distance;}
    void set_distance (similarity_type_t d) { _distance = d;}
};

#endif

