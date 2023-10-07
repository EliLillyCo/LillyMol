// Identify an extended pi system starting with a set of atoms.

#include <utility>

#include "Foundational/iwmisc/misc.h"

#include "extended_pi.h"

namespace radius_pi_extension {

// Algorithm will fail if anyone uses this number for a marker in `in_system`.
constexpr int kNextShell = 321511992;

bool
BothHavePiElectrons(Molecule& m,
                    atom_number_t a1,
                    atom_number_t a2)
{
  if (m.is_aromatic(a1) && m.is_aromatic(a2))
    return true;

  int pi1, pi2;
  m.pi_electrons(a1, pi1);
  m.pi_electrons(a2, pi2);
  return pi1 && pi2;
}

bool
IsConjugatedExtension(Molecule& m, const Bond * b) {
  if (b->is_aromatic())
    return true;
  if (! b->is_single_bond())
    return true;

  return BothHavePiElectrons(m, b->a1(), b->a2());
}

void
ExtendConjugated(Molecule& m,
                 int * in_system,
                 const Set_of_Atoms& outer_shell,
                 const RadiusPiExtensionParams& params)
{
  m.compute_aromaticity_if_needed();

  Set_of_Atoms atoms_in_next_shell;
  atoms_in_next_shell.resize_keep_storage(0);

  for (atom_number_t i : outer_shell) {
    const Atom & a = m.atom(i);
    for (const Bond * b : a) {
      const atom_number_t j = b->other(i);
      if (in_system[j] > 0) {
        continue;
      }

      if (! IsConjugatedExtension(m, b))
        continue;

      in_system[j] = kNextShell;
      atoms_in_next_shell.add(j);
    }
  }

  if (atoms_in_next_shell.empty())
    return;

  for (int r = 0; r < params.conjugated_bond_extend; ++r) {
    Set_of_Atoms next_next_shell;
    next_next_shell.reserve(m.natoms());
    for (auto i : atoms_in_next_shell) {
      in_system[i] = params.aromatic_extension_marker;
      const Atom& a = m.atom(i);
      for (const Bond * b : a) {
        if (! IsConjugatedExtension(m, b))
          continue;
        atom_number_t j = b->other(i);
        if (in_system[j])
          continue;

        in_system[j] = kNextShell;
        next_next_shell.add(j);
      }
    }
    if (next_next_shell.empty())
      return;
    atoms_in_next_shell = std::move(next_next_shell);
  }

  atoms_in_next_shell.set_vector(in_system, params.aromatic_extension_marker);
}

void
ExpandToRadius(const Molecule & m,
               int * in_system,
               const RadiusPiExtensionParams& params,
               Set_of_Atoms& atoms_in_next_shell) {

  const int matoms = m.natoms();
  atoms_in_next_shell.resize_keep_storage(0);

  for (int i = 0; i < matoms; ++i) {
    if (in_system[i] == 0 || in_system[i] == kNextShell)
      continue;

    const Atom & a = m.atom(i);
    for (const Bond * b : a) {
      const atom_number_t j = b->other(i);
      if (in_system[j] > 0) {
        continue;
      }

      in_system[j] = kNextShell;
      atoms_in_next_shell.add(j);
    }
  }

  for (int i = 0; i < params.radius - 1; ++i) {
    Set_of_Atoms next_next_shell;
    next_next_shell.reserve((i + 1) * 4);   // Heuristic.
    for (const auto j : atoms_in_next_shell) {
      in_system[j] = params.in_radius_marker + i * params.in_radius_step;
      const Atom& a = m.atom(j);
      for (const Bond * b : a) {
        const atom_number_t k = b->other(j);
        if (in_system[k] > 0) {
          continue;
        }
        in_system[k] = kNextShell;
        next_next_shell.add(k);
      }
    }
    atoms_in_next_shell = std::move(next_next_shell);
    if (atoms_in_next_shell.empty()) {  // Must be done after the move.
      break;
    }
  }

  atoms_in_next_shell.set_vector(in_system, params.in_radius_marker + params.radius * params.in_radius_step);
}

void
ExtendAromatic(Molecule& m,
               int * in_system,
               Set_of_Atoms& outer_shell,
               const RadiusPiExtensionParams& params)
{
  m.compute_aromaticity_if_needed();

  Set_of_Atoms next_shell;
  next_shell.reserve(m.natoms());

  for (const atom_number_t i : outer_shell) {
    const Atom& a = m.atom(i);
    for (const Bond * b : a) {
      if (! b->is_aromatic())
        continue;

      const atom_number_t j = b->other(i);
      if (in_system[j])
        continue;

      in_system[j] = kNextShell;
      next_shell.add(j);
    }
  }

  if (next_shell.empty())
    return;

  for (int r = 0; r < params.aromatic_bond_extend - 1; ++r) {
    Set_of_Atoms next_next_shell;
    next_next_shell.reserve(m.natoms());
    for (auto i : next_shell) {
      in_system[i] = params.aromatic_extension_marker;
      const Atom& a = m.atom(i);
      for (const Bond * b : a) {
        if (! b->is_aromatic())
          continue;
        atom_number_t j = b->other(i);
        if (in_system[j])
          continue;

        in_system[j] = kNextShell;
        next_next_shell.add(j);
      }
    }
    if (next_next_shell.empty())
      return;
    next_shell = std::move(next_next_shell);
  }

  for (const atom_number_t i : next_shell) {
    in_system[i] = params.aromatic_extension_marker;
  }
}

int
RadiusPiExtension(Molecule& m, int * in_system, const RadiusPiExtensionParams& params)
{
  const int matoms = m.natoms();

  Set_of_Atoms atoms_in_next_shell;
  if (params.radius > 0) {
    ExpandToRadius(m, in_system, params, atoms_in_next_shell);
  } else {
    for (int i = 0; i < matoms; ++i) {
      if (in_system[i])
        atoms_in_next_shell.add(i);
    }
  }

  if (params.aromatic_bond_extend > 0)
    ExtendAromatic(m, in_system, atoms_in_next_shell, params);

  if (params.conjugated_bond_extend > 0)
    ExtendConjugated(m, in_system, atoms_in_next_shell, params);

  return count_non_zero_occurrences_in_array (in_system, matoms);
}

};  // namespace radius_pi_extension
