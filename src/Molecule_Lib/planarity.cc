#include <iostream>
#include <memory>
#include <type_traits>

#include "graphLib.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Lib/planarity.h"

namespace iwplanarity {

using std::cerr;

namespace {

int
BuildEapsGraph(const Molecule& m, graphP g) {
  for (const Bond* b : m.bond_list()) {
    const int a1 = b->a1();
    const int a2 = b->a2();

    if (int rc = gp_AddEdge(g, a1 + 1, 0, a2 + 1, 0); rc != OK) {
      cerr << "BuildEapsGraph:failed to add bond " << a1 << ' ' << a2 << " rc " << rc << '\n';
      return 0;
    }
  }

  return 1;
}

static void
FillObstructionBonds(graphP g, PlanarityResult& result) {
#ifdef NO_LONGER_GEING_COMPUTED
  Newer versions of the graph library do not have these functions.
  disable for now...
  TODO:ianwatson investigate newer versions
  for (int e = gp_LowerBoundEdges(g); e < gp_UpperBoundEdges(g); e += 2) {
    if (gp_EdgeNotInUse(g, e)) {
      continue;
    }

    const int a1 = gp_GetNeighbor(g, gp_GetTwin(g, e));
    const int a2 = gp_GetNeighbor(g, e);

    if (a1 >= 1 && a2 >= 1 && a1 != a2) [[likely]] {
      result.obstruction_bonds.emplace_back(a1 - 1, a2 - 1);
    }
  }
#endif
}

// Custom deleter for a graphP.
auto FreeGraph = [](graphP g) {
  if (g != nullptr) {
    gp_Free(&g);
  }
};

using EapsGraph = std::unique_ptr<std::remove_pointer_t<graphP>, decltype(FreeGraph)>;

}  // namespace

PlanarityResult
Planarity(const Molecule& m) {
  PlanarityResult result;

  EapsGraph g(gp_New(), FreeGraph);
  if (g == nullptr) {
    cerr << "Planarity:cannot allocate graph\n";
    return result;
  }

  if (gp_InitGraph(g.get(), m.natoms()) != OK) {
    cerr << "Planarity:cannot initialise graph\n";
    result.status = iwplanarity::PlanarityStatus::kError;
    return result;
  }

  if (! BuildEapsGraph(m, g.get())) {
    cerr << "Planarity:cannot build EAPS Graph\n";
    result.status = iwplanarity::PlanarityStatus::kError;
    return result;
  }

  const int rc = gp_Embed(g.get(), EMBEDFLAGS_OUTERPLANAR);
  // cerr << "gp_Embed rc " << rc << '\n';

  if (rc == OK) {
    result.status = PlanarityStatus::kPlanar;
    return result;
  }

  if (rc != NONEMBEDDABLE) {
    result.status = PlanarityStatus::kError;
    return result;
  }

  result.status = PlanarityStatus::kNonPlanar;
  FillObstructionBonds(g.get(), result);

  return result;
}

} // namespace iwplanarity
