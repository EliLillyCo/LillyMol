#ifndef MOLECULE_LIB_QRY_WCHARGE_H_
#define MOLECULE_LIB_QRY_WCHARGE_H_

#include "accumulator.h"
#include "iwmtypes.h"
#include "substructure.h"

/*
  A Query_and_Charge_Stats object is a substructure_query and details on
  charges found on matched atoms. It is used primarily in iwdescr.
*/

class Charge_Distribution : public Accumulator<charge_t>
{
  private:
  
  public:
};

class Query_and_Charge_Stats: public Substructure_Query
{
  private:
    resizable_array_p<Charge_Distribution> _charge_distributions;
    int _embeddings;

//  The maximum number of atoms in the query

    int _max_atoms;

// private functions

    void _default_values ();
    int  _allocate_charge_distributions ();
    int  _determine_max_atoms ();

  public:
//  Query_and_Charge_Stats (const char *);
    Query_and_Charge_Stats (const const_IWSubstring &);
    Query_and_Charge_Stats ();
    ~Query_and_Charge_Stats ();

    int tally_embedding (const Molecule &, Substructure_Results &);

    int report (int, std::ostream &) const;
};

#endif  // MOLECULE_LIB_QRY_WCHARGE_H_
