#include "Molecule_Tools/ring_system_shape.h"

#include <algorithm>
#include <vector>

namespace ring_system_shape {
namespace {

Set_of_Atoms
RingSystemAtoms(const Molecule& m, const int* ring_system_membership,
                int ring_system_id) {
  Set_of_Atoms result;
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (ring_system_membership[i] == ring_system_id) {
      result.add(i);
    }
  }

  return result;
}

Set_of_Atoms
UniqueAttachmentRingAtoms(const std::vector<RingSystemAttachment>& attachments) {
  Set_of_Atoms result;

  for (const RingSystemAttachment& attachment : attachments) {
    result.add_if_not_already_present(attachment.ring_atom);
  }

  return result;
}

struct BfsWorkspace {
  explicit BfsWorkspace(int matoms) : separation(matoms, -1) {
    queue.reserve(matoms);
  }

  std::vector<int> separation;
  std::vector<atom_number_t> queue;
};

bool
IsTerminalExocyclicMultipleBond(const Molecule& m, const Bond& bond,
                                atom_number_t outside_atom) {
  return !bond.is_single_bond() && m.ncon(outside_atom) == 1;
}

std::vector<RingSystemAttachment>
IdentifyAttachments(const Molecule& m, const int* ring_system_membership,
                    int ring_system_id, isotope_t only_process_isotope) {
  std::vector<RingSystemAttachment> result;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_system_membership[i] != ring_system_id) {
      continue;
    }
    if (only_process_isotope > 0 && m.isotope(i) != only_process_isotope) {
      continue;
    }

    const Atom& atom = m[i];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(i);
      if (ring_system_membership[other] == ring_system_id) {
        continue;
      }

      // Terminal multiply bonded exocyclic atoms, like =O, are not exit points.
      if (IsTerminalExocyclicMultipleBond(m, *bond, other)) {
        continue;
      }

      result.push_back(RingSystemAttachment{static_cast<atom_number_t>(i), other});
    }
  }

  return result;
}

bool
BfsFrom(const Molecule& m, const int* ring_system_membership, int ring_system_id,
        const Set_of_Atoms& ring_system_atoms, atom_number_t start,
        BfsWorkspace& workspace, RingSystemSpan& result) {
  for (int i = 0; i < ring_system_atoms.number_elements(); ++i) {
    workspace.separation[ring_system_atoms[i]] = -1;
  }
  workspace.queue.clear();

  workspace.separation[start] = 0;
  workspace.queue.push_back(start);

  for (size_t queue_index = 0; queue_index < workspace.queue.size(); ++queue_index) {
    const atom_number_t current = workspace.queue[queue_index];

    const Atom& atom = m[current];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(current);
      if (ring_system_membership[other] != ring_system_id) {
        continue;
      }
      if (workspace.separation[other] >= 0) {
        continue;
      }

      workspace.separation[other] = workspace.separation[current] + 1;
      workspace.queue.push_back(other);
    }
  }

  if (workspace.queue.size() != ring_system_atoms.size()) {
    return false;
  }

  result.from = start;
  result.max_separation = 0;
  result.farthest_atoms.clear();

  for (atom_number_t atom : ring_system_atoms) {
    if (workspace.separation[atom] > result.max_separation) {
      result.max_separation = workspace.separation[atom];
      result.farthest_atoms.clear();
      result.farthest_atoms.push_back(atom);
    } else if (workspace.separation[atom] == result.max_separation && atom != start) {
      result.farthest_atoms.push_back(atom);
    }
  }

  return true;
}

int
ShortestPathBetweenExitPoints(const Molecule& m, const int* ring_system_membership,
                              int ring_system_id, atom_number_t start,
                              atom_number_t destination, BfsWorkspace& workspace) {
  workspace.queue.clear();

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_system_membership[i] == ring_system_id) {
      workspace.separation[i] = -1;
    }
  }

  workspace.separation[start] = 0;
  workspace.queue.push_back(start);

  for (size_t queue_index = 0; queue_index < workspace.queue.size(); ++queue_index) {
    const atom_number_t current = workspace.queue[queue_index];
    if (current == destination) {
      return workspace.separation[current];
    }

    const Atom& atom = m[current];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(current);
      if (ring_system_membership[other] != ring_system_id) {
        continue;
      }
      if (workspace.separation[other] >= 0) {
        continue;
      }

      workspace.separation[other] = workspace.separation[current] + 1;
      workspace.queue.push_back(other);
    }
  }

  return -1;
}

const RingSystemSpan*
FindRingSystemSpan(const std::vector<RingSystemSpan>& ring_system_spans,
                   atom_number_t from) {
  for (const RingSystemSpan& ring_system_span : ring_system_spans) {
    if (ring_system_span.from == from) {
      return &ring_system_span;
    }
  }

  return nullptr;
}

bool
IsTerminalSingleAtomSubstituent(const Molecule& m, const Bond& bond,
                                atom_number_t outside_atom) {
  return bond.is_single_bond() && m.ncon(outside_atom) == 1;
}

bool
ContainsAtom(const std::vector<atom_number_t>& atoms, atom_number_t needle) {
  return std::find(atoms.begin(), atoms.end(), needle) != atoms.end();
}

void
SetRodDeficit(const Molecule& m, const int* ring_system_membership, int ring_system_id,
              const Set_of_Atoms& exit_points,
              const std::vector<RingSystemSpan>& ring_system_spans,
              BfsWorkspace& workspace, RingSystemShape& result) {
  if (exit_points.size() != 2) {
    return;
  }

  const atom_number_t a1 = exit_points[0];
  const atom_number_t a2 = exit_points[1];

  result.observed_separation = ShortestPathBetweenExitPoints(
      m, ring_system_membership, ring_system_id, a1, a2, workspace);
  if (result.observed_separation < 0) {
    return;
  }

  const RingSystemSpan* span1 = FindRingSystemSpan(ring_system_spans, a1);
  const RingSystemSpan* span2 = FindRingSystemSpan(ring_system_spans, a2);
  if (span1 == nullptr || span2 == nullptr) {
    return;
  }

  const int deficit1 = span1->max_separation - result.observed_separation;
  const int deficit2 = span2->max_separation - result.observed_separation;
  result.rod_deficit = std::max(deficit1, deficit2);
}

RingSystemShapeClass
Classify(const Set_of_Atoms& exit_points,
         const std::vector<RingSystemSpan>& ring_system_spans) {
  if (exit_points.size() <= 1) {
    return RingSystemShapeClass::kNotApplicable;
  }

  if (exit_points.size() > 2) {
    return RingSystemShapeClass::kMultiSubstituted;
  }

  const atom_number_t a1 = exit_points[0];
  const atom_number_t a2 = exit_points[1];

  const RingSystemSpan* span1 = FindRingSystemSpan(ring_system_spans, a1);
  const RingSystemSpan* span2 = FindRingSystemSpan(ring_system_spans, a2);
  if (span1 == nullptr || span2 == nullptr) {
    return RingSystemShapeClass::kInvalid;
  }

  if (ContainsAtom(span1->farthest_atoms, a2) &&
      ContainsAtom(span2->farthest_atoms, a1)) {
    return RingSystemShapeClass::kRodLike;
  }

  return RingSystemShapeClass::kNotRodLike;
}

}  // namespace

int
NonRingBranchPointCount(const Molecule& m, const int* ring_system_membership) {
  if (ring_system_membership == nullptr) {
    return 0;
  }

  int result = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    // This is the cheapest and most important exclusion. Ring-system shape is
    // handled separately.
    if (ring_system_membership[i] > 0) {
      continue;
    }

    if (m.ncon(i) <= 2) {
      continue;
    }

    int meaningful_neighbours = 0;
    bool has_terminal_multiple_bond = false;
    const Atom& atom = m[i];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(i);

      if (IsTerminalExocyclicMultipleBond(m, *bond, other)) {
        has_terminal_multiple_bond = true;
        break;
      }

      if (IsTerminalSingleAtomSubstituent(m, *bond, other)) {
        continue;
      }

      ++meaningful_neighbours;
    }

    if (has_terminal_multiple_bond) {
      continue;
    }

    if (meaningful_neighbours > 2) {
      ++result;
    }
  }

  return result;
}

int
RingAtomBranchPointCount(const Molecule& m, const int* ring_system_membership) {
  return RingAtomBranchPointCount(m, ring_system_membership, 0);
}

int
RingAtomBranchPointCount(const Molecule& m, const int* ring_system_membership,
                         isotope_t only_process_isotope) {
  if (ring_system_membership == nullptr) {
    return 0;
  }

  int result = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_system_membership[i] <= 0) {
      continue;
    }
    if (only_process_isotope > 0 && m.isotope(i) != only_process_isotope) {
      continue;
    }

    if (m.ncon(i) < 4) {
      continue;
    }

    int substantial_outside_connections = 0;
    const Atom& atom = m[i];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(i);
      if (ring_system_membership[other] == ring_system_membership[i]) {
        continue;
      }

      if (!bond->is_single_bond()) {
        continue;
      }

      if (ring_system_membership[other] > 0) {
        ++substantial_outside_connections;
      } else if (m.ncon(other) > 1) {
        ++substantial_outside_connections;
      }
    }

    if (substantial_outside_connections >= 2) {
      ++result;
    }
  }

  return result;
}

bool
AnalyseRingSystemShape(const Molecule& m, const int* ring_system_membership,
                       int ring_system_id, RingSystemShape& result) {
  return AnalyseRingSystemShape(m, ring_system_membership, ring_system_id, 0, result);
}

bool
AnalyseRingSystemShape(const Molecule& m, const int* ring_system_membership,
                       int ring_system_id, isotope_t only_process_isotope,
                       RingSystemShape& result) {
  result = RingSystemShape();

  if (ring_system_membership == nullptr || ring_system_id <= 0) {
    return false;
  }

  const Set_of_Atoms ring_system_atoms =
      RingSystemAtoms(m, ring_system_membership, ring_system_id);
  if (ring_system_atoms.empty()) {
    return false;
  }

  result.attachments = IdentifyAttachments(m, ring_system_membership, ring_system_id,
                                          only_process_isotope);
  const Set_of_Atoms exit_points = UniqueAttachmentRingAtoms(result.attachments);

  BfsWorkspace workspace(m.natoms());
  result.ring_system_spans.reserve(exit_points.size());
  for (atom_number_t exit_point : exit_points) {
    RingSystemSpan ring_system_span;
    if (!BfsFrom(m, ring_system_membership, ring_system_id, ring_system_atoms, exit_point,
                 workspace, ring_system_span)) {
      result.shape_class = RingSystemShapeClass::kInvalid;
      return false;
    }
    result.ring_system_spans.push_back(std::move(ring_system_span));
  }

  SetRodDeficit(m, ring_system_membership, ring_system_id, exit_points,
                result.ring_system_spans, workspace, result);
  result.shape_class = Classify(exit_points, result.ring_system_spans);
  return true;
}

}  // namespace ring_system_shape
