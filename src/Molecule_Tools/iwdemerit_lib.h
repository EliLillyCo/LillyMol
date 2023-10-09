#ifndef MOLECULE_TOOLS_IWDEMERIT_LIB_H
#define MOLECULE_TOOLS_IWDEMERIT_LIB_H
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/qry_wstats.h"

#include "demerit.h"

namespace lilly_medchem_rules {

class MCDemerit {
  private:
    int _verbose = 0;

    int _keep_going_after_rejection = 0;

    int _hard_lower_atom_count_cutoff = 0;
    int _soft_lower_atom_count_cutoff = 0;
    int _soft_upper_atom_count_cutoff = 0;
    int _hard_upper_atom_count_cutoff = 0;

    int _lower_atom_count_demerit = 100;
    int _upper_atom_count_demerit = 100;

    int _reduce_to_largest_fragment = 0;

    // If reduction to largest fragment is requested, queries are
    // subdivided into those that are applied to the whole molecule
    // and those that look only at the largest fragment.
    resizable_array_p<Substructure_Hit_Statistics> _whole_molecule;
    resizable_array_p<Substructure_Hit_Statistics> _largest_fragment;

    // private functions.
    void SeparateDependingOnFragmentMatch(resizable_array_p<Substructure_Hit_Statistics> & queries);

    // Run each of the queries in `queries`, and assign demerits based
    // on the demerit value in `demerit_numeric_value_index`.
    void RunASetOfQueries(Molecule_to_Match & target,
                 resizable_array_p<Substructure_Hit_Statistics> & queries,
                 int demerit_numeric_value_index,
                 Demerit & demerit) const;

    void DoAtomCountDemerits(Molecule & m, Demerit & demerit) const;

  public:
    MCDemerit();

    int Build(Command_Line& cl);

    // The function that actually does work.
    // `demerit_numeric_value_index` is the numeric value index within the
    // arrays of numeric values in the queries.
    Demerit Process(Molecule& m, int demerit_numeric_value_index);
};

}  // namespace lilly_medchem_rules

#endif // MOLECULE_TOOLS_IWDEMERIT_LIB_H
