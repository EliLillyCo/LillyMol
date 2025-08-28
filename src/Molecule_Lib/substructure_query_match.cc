#include "Foundational/iwaray/iwaray.h"

#include "molecule.h"
#include "substructure.h"
#include "target.h"

namespace lillymol {

// Return true if all the queries in `queries` match `m`.
int
AllQueriesMatch(Molecule& m, resizable_array_p<Substructure_Query> & queries) {
  Molecule_to_Match target(&m);

  return AllQueriesMatch(target, queries);
}

int
AllQueriesMatch(Molecule_to_Match& target, resizable_array_p<Substructure_Query>& queries) {
  for (Substructure_Query* q : queries) {
    if (! q->substructure_search(target)) {
      return 0;
    }
  }

  return 1;
}

// Return true if any of the queries in `queries` matches `m`.
int
AnyQueryMatches(Molecule& m, resizable_array_p<Substructure_Query> & queries) {
  Molecule_to_Match target(&m);

  return AnyQueryMatches(target, queries);
}

int
AnyQueryMatches(Molecule_to_Match& target, resizable_array_p<Substructure_Query>& queries ) {

  for (Substructure_Query* q : queries) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

} // namespace lillymol
