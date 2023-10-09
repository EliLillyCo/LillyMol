// Functions relating to conversion from Molecule to Graph form.

#include <iostream>
#include <memory>

#define COMPILING_MOL2GRAPH_CC 1

#include "Foundational/iwmisc/misc.h"

#include "chiral_centre.h"
#include "molecule.h"

using std::cerr;

static int exclude_triple_bonds_from_graph_reduction = 0;

void
set_exclude_triple_bonds_from_graph_reduction(int s)
{
  exclude_triple_bonds_from_graph_reduction = s;
}

Mol2Graph::Mol2Graph()
{
  _exclude_triple_bonds_from_graph_reduction = 0;
  _revert_all_directional_bonds_to_non_directional = 1;
  _preserve_cc_double_bonds_saturated = 0;
  _preserve_cc_double_bonds_no_heteroatoms = 0;

  _remove_chiral_centres = 1;

  _append_molecular_formula = 0;

  _aromatic_distinguishing_formula = 0;

  _active = false;

  return;
}

int
Mol2Graph::debug_print(std::ostream & output) const
{
  output << "Mol2Graph:debug_print:\n";
  output << " _remove chiral centres " << _remove_chiral_centres << '\n';
  output << " _exclude_triple_bonds_from_graph_reduction "
         << _exclude_triple_bonds_from_graph_reduction << '\n';
  output << " _revert_all_directional_bonds_to_non_directional "
         << _revert_all_directional_bonds_to_non_directional << '\n';
  output << " _preserve_cc_double_bonds_saturated " << _preserve_cc_double_bonds_saturated << '\n';
  output << " _preserve_cc_double_bonds_no_heteroatoms " << _preserve_cc_double_bonds_no_heteroatoms
         << '\n';
  output << " _append_molecular_formula " << _append_molecular_formula << '\n';
  output << " _aromatic_distinguishing_formula " << _aromatic_distinguishing_formula << '\n';

  return output.good();
}

void
Mol2Graph::TurnOnMostUsefulOptions() {
  _remove_chiral_centres = 1;
  _exclude_triple_bonds_from_graph_reduction = 1;
  _revert_all_directional_bonds_to_non_directional = 1;
  _preserve_cc_double_bonds_no_heteroatoms = 1;
  _preserve_cc_double_bonds_saturated = 1;

  return;
}

int
Mol2Graph::construct(Command_Line & cl, const char flag, const int verbose)
{
  const_IWSubstring h;
  for (int i = 0; cl.value(flag, h, i); ++i)
  {
    if ('s' == h)
    {
      _preserve_cc_double_bonds_saturated = 1;
      if (verbose)
        cerr
            << "During graph reduction, will preserve double bonds adjacent to fully saturated atoms\n";
    }
    else if ('c' == h)
    {
      _preserve_cc_double_bonds_no_heteroatoms = 1;
      if (verbose)
        cerr
            << "During graph reduction, will preserve double bonds adjacent to only carbon atoms\n";
    }
    else if ("kptriple" == h) {
      _exclude_triple_bonds_from_graph_reduction = true;
      if (verbose)
        cerr << "Triple bonds not changed during graph formation\n";
    }
    else if ("chiral" == h)
    {
      _remove_chiral_centres = 0;
      if (verbose)
        cerr << "Will NOT remove chiral centres during graph reduction\n";
    }
    else if ("rmchiral" == h)
    {
      _remove_chiral_centres = 1;
      if (verbose)
        cerr << "Will remove chiral centres during graph reduction\n";
    }
    else if (h == "def") {
      _active = true;
    }
    else if ("help" == h)
    {
      // clang-format off
      cerr << "Mol2Graph::construct:the following options are recognised\n";
      cerr << " -" << flag << " def        activate with all defaults\n";
      cerr << " -" << flag << " s          preserve C=C bonds that are adjacent to fully saturated atoms\n";
      cerr << " -" << flag << " c          preserve C=C bonds that are adjacent to only carbon atoms\n";
      cerr << " -" << flag << " kptriple   preserve triple bonds during graph formation\n";
      cerr << " -" << flag << " chiral     preserve chiral centres during graph reduction\n";
      cerr << " -" << flag << " rmchiral   remove chiral centres during graph reduction\n";
      // clang-format on
      exit(0);
    }
    else
    {
      cerr << "Mol2Graph::construct:unrecognised -" << flag << " qualifier '" << h << "'\n";
      return 0;
    }
  }

  _active = true;

  return 1;
}


int
Molecule::change_to_graph_form()
{
  int changes = 0;    // changes which require recomputation of nbonds & implicit hydrogens

  // We remove all explicit hydrogens

  Set_of_Atoms atoms_to_be_removed;

  for (Bond * b : _bond_list)
  {
    if (b->is_single_bond())
      continue;

    if (b->is_triple_bond() && exclude_triple_bonds_from_graph_reduction)
      continue;

    b->set_bond_type(SINGLE_BOND);
    b->set_permanent_aromatic(0);
    changes++;
  }

  if (_chiral_centres.number_elements())
  {
    _chiral_centres.resize(0);
    changes++;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (1 == a->atomic_number())
    {
      atoms_to_be_removed.add(i);
      changes++;
      continue;
    }

    //  Remove charges unless it would create valence problems

    if (a->formal_charge())
    {
      if (7 == a->atomic_number() && 4 == a->ncon())
        ;
      else
      {
        a->set_formal_charge(0);
        changes++;
      }
    }

    if (a->implicit_hydrogens_known())
    {
      a->set_implicit_hydrogens_known(0);
      changes++;
    }

    if (a->is_isotope())
    {
      a->set_isotope(0);
    }
  }

  if (atoms_to_be_removed.number_elements())
    remove_atoms(atoms_to_be_removed);

  if (revert_all_directional_bonds_to_non_directional())
    changes++;

  if (changes)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->recompute_nbonds();
      int unused;
      _things[i]->recompute_implicit_hydrogens(unused);
    }
  }

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  return 1;
}

/*
  Are all the connections from atom ZATOM fully saturated
*/

int
Molecule::_all_connections_saturated(const atom_number_t zatom, const atom_number_t ignore) const
{
  const Atom * a = _things[zatom];

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (j == ignore)
      continue;

    const Atom * aj = _things[j];

    if (aj->ncon() < aj->nbonds())    // unsaturated
      return 0;
  }

  return 1;
}

int
Molecule::_double_bond_needs_changing_for_graph_form(const Bond & b,
                                                     const Mol2Graph & mol2graph) const
{
  assert(b.is_double_bond());

  if (! mol2graph.some_kind_of_double_bond_preservation_active())
    return 1;

  // Now we have non-aromatic double bonds

  const Atom * a1 = _things[b.a1()];
  if (6 != a1->atomic_number())    // must be changed
    return 1;

  const Atom * a2 = _things[b.a2()];
  if (6 != a2->atomic_number())    // must be changed
    return 1;

  if (mol2graph.preserve_cc_double_bonds_saturated())
  {
    if (! _all_connections_saturated(b.a1(), b.a2()))    // must be changed
      return 1;

    if (! _all_connections_saturated(b.a2(), b.a1()))    // must be changed
      return 1;

    return 0;    // all neighbours fully saturated, does not need changing
  }

  if (mol2graph.preserve_cc_double_bonds_no_heteroatoms())
  {
    if (! _all_connections_saturated(b.a1(), b.a2()))    // must be changed
      return 1;

    if (! _all_connections_saturated(b.a2(), b.a1()))    // must be changed
      return 1;

    if (attached_heteroatom_count(b.a1() > 0) || attached_heteroatom_count(b.a2() > 0))
      return 1;

    return 0;
  }

  return 1;
}

int
Molecule::change_to_graph_form(const Mol2Graph & mol2graph)
{
  int changes = 0;    // changes which require recomputation of nbonds & implicit hydrogens

  const int nb = _bond_list.number_elements();

  int * change_to_single = new_int(nb, 1);
  std::unique_ptr<int[]> free_change_to_single(change_to_single);

  if (mol2graph.some_kind_of_double_bond_preservation_active())
    compute_aromaticity_if_needed();

  for (int i = 0; i < nb; ++i)
  {
    Bond * b = _bond_list[i];

    if (b->is_single_bond())    // no need to change
      change_to_single[i] = 0;
    else if (b->is_aromatic())    // non single, aromatic bond, definitely change
      b->set_permanent_aromatic(0);
    else if (b->is_triple_bond())
      change_to_single[i] = ! mol2graph.exclude_triple_bonds_from_graph_reduction();
    else if (b->is_double_bond() && ! mol2graph.some_kind_of_double_bond_preservation_active())
      ;
    else if (! _double_bond_needs_changing_for_graph_form(*b, mol2graph))
      change_to_single[i] = 0;
  }

  for (int i = 0; i < nb; i++)
  {
    if (! change_to_single[i])
      continue;

    Bond * b = _bond_list[i];

    b->set_bond_type(SINGLE_BOND);
    b->set_permanent_aromatic(0);
    changes++;
  }

  if (mol2graph.remove_chiral_centres() && _chiral_centres.number_elements())
  {
    _chiral_centres.resize(0);
    changes++;
  }

  Set_of_Atoms atoms_to_be_removed;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (1 == a->atomic_number())
    {
      atoms_to_be_removed.add(i);
      changes++;
      continue;
    }

    //  Remove charges unless it would create valence problems

    if (a->formal_charge())
    {
      if (7 == a->atomic_number() && 4 == a->ncon())
        ;
      else
      {
        a->set_formal_charge(0);
        changes++;
      }
    }

    if (a->implicit_hydrogens_known())
    {
      a->set_implicit_hydrogens_known(0);
      changes++;
    }

    if (a->is_isotope())
      a->set_isotope(0);
  }

  if (atoms_to_be_removed.number_elements())
    remove_atoms(atoms_to_be_removed);

  if (mol2graph.revert_all_directional_bonds_to_non_directional())
  {
    if (revert_all_directional_bonds_to_non_directional())
      changes++;
  }

  if (changes)
  {
    int unused;
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->recompute_nbonds();
      _things[i]->recompute_implicit_hydrogens(unused);
    }
  }

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  return 1;
}
