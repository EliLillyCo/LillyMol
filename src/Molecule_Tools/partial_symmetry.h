#ifndef MOLECULE_TOOLS_PARTIAL_SYMMETRY_H
#define MOLECULE_TOOLS_PARTIAL_SYMMETRY_H

#include "Molecule_Lib/molecule.h"

namespace partial_symmetry {

struct InLayer {
  atom_number_t atom;
  int score;
};

class PartialSymmetry {
  private:
    Molecule& _m;
    int _matoms;
    InLayer* _score;
    int * _symmetric_at_radius;
    Set_of_Atoms _current_shell;
    Set_of_Atoms _next_shell;

    bool _computation_done;

  // private functions.

    int Expand(atom_number_t starting_atom, int radius);
    int Expand(atom_number_t starting_atom1, atom_number_t starting_atom2, int radius);
    int SortAndAssign(int radius);
    void AdjustForRings(const int * dm);
    void AllScoresZero();

  public:
    PartialSymmetry(Molecule& m);
    ~PartialSymmetry();

    const int * SymmetricAtRadius();
};

}  // namespace partial_symmetry
#endif
