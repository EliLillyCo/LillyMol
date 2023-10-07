#include "topological_sterimol.h"

#include <stdlib.h>

#include <iomanip>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"

using std::cerr;

Topological_Sterimol::Topological_Sterimol()
{
  _do_partial_charge_descriptors = 0;

  return;
}

int
Topological_Sterimol::write_descriptors(std::ostream& os) const
{
  for (int i = 0; i < NTOPOLOGICAL_STERIMOL; i++) {
    os << ' ' << _descriptor[i];
  }

  if (_do_partial_charge_descriptors) {
    for (int i = 0; i < NTOPO_STERIMOL_CHARGE; i++) {
      os << ' ' << _partial_charge_descriptor[i];
    }
  }

  return os.good();
}

static int
identify_sterimol_atoms(Molecule& m, int* in_fragment, atom_number_t zatom)
{
  assert(in_fragment[zatom]);

#ifdef DEBUG_IDENTIFY_STERIMOL_ATOMS
  cerr << "Adding atom " << zatom << " to fragment\n";
#endif

  int rc = 1;  // the number of atoms we add

  const Atom* a = m.atomi(zatom);

  int acon = a->ncon();
  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(zatom, i);

    if (in_fragment[j]) {
      continue;
    }

    in_fragment[j] = 1;
    rc += identify_sterimol_atoms(m, in_fragment, j);
  }

  return rc;
}

int
Topological_Sterimol::compute_descriptors(Molecule& m, const int* don_acc,
                                          int* in_fragment, const atom_number_t a1,
                                          const atom_number_t a2)
{
  assert(m.ok_2_atoms(a1, a2));

  if (m.in_same_ring(a1, a2)) {
    cerr << "Topological_Sterimol::compute_descriptors: ring atoms " << a1 << " and "
         << a2 << '\n';
    return 0;
  }

  Accumulator_Int<int> acc_rad, acc_arom, acc_heteroatom;

  int atoms_at_distance[4];
  set_vector(atoms_at_distance, 4, 0);

  const int matoms = m.natoms();

  set_vector(in_fragment, matoms, 0);

  in_fragment[a1] = 9;
  in_fragment[a2] = 1;

#ifdef DEBUG_IDENTIFY_STERIMOL_ATOMS
  cerr << "Starting fragment, from atom " << a1 << " to " << a2 << '\n';
#endif

  int n = identify_sterimol_atoms(m, in_fragment, a2);

  // cerr << n << " atoms included in fragment\n";

  if (0 == n || matoms == n) {
    cerr << "Topological_Sterimol::compute_descriptors: failed to identify fragment\n";
    cerr << "matoms = " << matoms << " found " << n << " atoms in fragment\n";
    return 0;
  }

  set_vector(_descriptor, NTOPOLOGICAL_STERIMOL, 0);

  set_vector(_partial_charge_descriptor, NTOPO_STERIMOL_CHARGE,
             static_cast<charge_t>(0.0));

  _descriptor[TOPOLOGICAL_STERIMOL_NATOMS] = n;

  if (_do_partial_charge_descriptors) {
    _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_Q0] = m.charge_on_atom(a1);
    _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_Q1] = m.charge_on_atom(a2);
  }

  charge_t max_partial_charge, min_partial_charge;
  if (_do_partial_charge_descriptors) {
    charge_t q = m.charge_on_atom(a2);
    min_partial_charge = q;
    max_partial_charge = q;
  } else {
    // initialised just to keep the compiler quiet
    min_partial_charge = static_cast<charge_t>(0.0);
    max_partial_charge = static_cast<charge_t>(0.0);
  }

  charge_t total_absolute_charge = static_cast<charge_t>(0.0);

  for (int i = 0; i < matoms; i++) {
    if (1 != in_fragment[i]) {  // skip atom A1
      continue;
    }

    if (m.is_ring_atom(i)) {
      _descriptor[TOPOLOGICAL_STERIMOL_NRINGS]++;
    }

    if (m.is_aromatic(i)) {
      _descriptor[TOPOLOGICAL_STERIMOL_AROM]++;
    }

    Atom* a = const_cast<Atom*>(m.atomi(i));  // nbonds is non-const

    if (0 == a->formal_charge()) {
      ;
    } else if (a->formal_charge() < 0) {
      _descriptor[TOPOLOGICAL_STERIMOL_NEG]++;
    } else {
      _descriptor[TOPOLOGICAL_STERIMOL_POS]++;
    }

    if (1 == don_acc[i] || 2 == don_acc[i]) {  // acceptor
      _descriptor[TOPOLOGICAL_STERIMOL_ACCEPT]++;
    }
    if (2 == don_acc[i] || 3 == don_acc[i]) {  // acceptor
      _descriptor[TOPOLOGICAL_STERIMOL_DONOR]++;
    }

    if (a->ncon() < a->nbonds()) {
      _descriptor[TOPOLOGICAL_STERIMOL_UNSAT]++;
    }

    const Element* e = a->element();

    if (6 != e->atomic_number() && 1 != e->atomic_number()) {
      _descriptor[TOPOLOGICAL_STERIMOL_HETERO]++;
    }

    if (_do_partial_charge_descriptors) {
      charge_t q = m.charge_on_atom(i);
      if (q > max_partial_charge) {
        max_partial_charge = q;
      } else if (q < min_partial_charge) {
        min_partial_charge = q;
      }

      if (q >= static_cast<charge_t>(0.0)) {
        total_absolute_charge += q;
      } else {
        total_absolute_charge += -q;
      }
    }

    if (i != a2) {
      const auto d = m.bonds_between(a2, i);

      if (d < 4) {
        atoms_at_distance[d]++;
      }

      acc_rad.extra(d);
      if (m.is_aromatic(i)) {
        acc_arom.extra(d);
      }
      if (6 != e->atomic_number()) {
        acc_heteroatom.extra(d);
      }

      if (d > _descriptor[TOPOLOGICAL_STERIMOL_MXDIST]) {
        _descriptor[TOPOLOGICAL_STERIMOL_MXDIST] = d;
        if (_do_partial_charge_descriptors) {
          charge_t q = m.charge_on_atom(i);

          if (q < static_cast<charge_t>(0.0)) {
            q = -q;
          }

          if (q > _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_QEXT]) {
            _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_QEXT] = q;
          }
        }
      }

      if ((1 == don_acc[i] || 2 == don_acc[i]) &&
          d > _descriptor[TOPOLOGICAL_STERIMOL_MXACCDIST]) {
        _descriptor[TOPOLOGICAL_STERIMOL_MXACCDIST] = d;
      }
      if ((2 == don_acc[i] || 3 == don_acc[i]) &&
          d > _descriptor[TOPOLOGICAL_STERIMOL_MXDONDIST]) {
        _descriptor[TOPOLOGICAL_STERIMOL_MXDONDIST] = d;
      }
    }
  }

  _descriptor[TOPOLOGICAL_STERIMOL_CIGAR] = _descriptor[TOPOLOGICAL_STERIMOL_MXDIST] *
                                            100 /
                                            _descriptor[TOPOLOGICAL_STERIMOL_NATOMS];

  if (acc_rad.n() > 0) {
    _descriptor[TOPOLOGICAL_STERIMOL_AVE_DIST_10] =
        static_cast<int>(acc_rad.average() * 10.0 + 0.499999);
  }

  _descriptor[TOPOLOGICAL_STERIMOL_ATOMS_AT_1] = atoms_at_distance[1];
  _descriptor[TOPOLOGICAL_STERIMOL_ATOMS_AT_2] = atoms_at_distance[2];
  _descriptor[TOPOLOGICAL_STERIMOL_ATOMS_AT_3] = atoms_at_distance[3];

  if (acc_arom.n() > 0) {
    _descriptor[TOPOLOGICAL_STERIMOL_MIN_AROM] = acc_arom.minval();
    _descriptor[TOPOLOGICAL_STERIMOL_MAX_AROM] = acc_arom.maxval();
  }

  if (acc_heteroatom.n() > 0) {
    _descriptor[TOPOLOGICAL_STERIMOL_MIN_HTRO] = acc_heteroatom.minval();
    _descriptor[TOPOLOGICAL_STERIMOL_MAX_HTRO] = acc_heteroatom.maxval();
  }

  if (_do_partial_charge_descriptors) {
    _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_MAXQ] = max_partial_charge;
    _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_MINQ] = min_partial_charge;
    _partial_charge_descriptor[TOPOLOGICAL_STERIMOL_TABS] = total_absolute_charge;
  }

  return 1;
}
