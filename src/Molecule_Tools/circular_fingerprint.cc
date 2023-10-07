#include <stdlib.h>

#include "Foundational/iwmisc/misc.h"
#include "circular_fingerprint.h"

#define PROCESSING_FINISHED 1
#define READY_TO_PROCESS 2
#define NEXT_TIME 3

static int default_circular_fingerprint_additive = 1;

static int all_bonds_same_type = 0;

void
set_default_circular_fingerprint_additive (int s)
{
  default_circular_fingerprint_additive = s;

  return;
}

Circular_Fingerprint_Generator::Circular_Fingerprint_Generator ()
{
  _additive = default_circular_fingerprint_additive;
  _matoms = 0;
  _allocated = 0;
  _min_radius = 0;
  _max_radius = 3;
  _atype = nullptr;
  _atom = nullptr;
  _processing_status = nullptr;

  return;
}

Circular_Fingerprint_Generator::~Circular_Fingerprint_Generator ()
{
  if (_allocated > 0)
    _free_dynamic_arrays();

  return;
}

void
Circular_Fingerprint_Generator::_free_dynamic_arrays ()
{
  if (nullptr != _atype)
    delete [] _atype;

  if (nullptr != _atom)
    delete [] _atom;

  if (nullptr != _processing_status)
    delete [] _processing_status;

  _allocated = 0;

  return;
}

int
Circular_Fingerprint_Generator::_initialise_molecule (const Molecule & m,
                                                const int * atype,
                                                const int * include_these_atoms)
{
  if (m.natoms() > _allocated)
    _free_dynamic_arrays();

  if (0 == _allocated)
  {
    _matoms = m.natoms();
    _allocated = _matoms + _matoms;
    _atype = new int[_allocated];
    _atom = new const Atom * [_allocated];
    _processing_status = new int[_allocated];
  }

  _matoms = m.natoms();
  m.atoms((const Atom **) _atom);

  for (auto i = 0; i < _matoms; ++i)
  {
    if (include_these_atoms[i])
      _atype[i] = atype[i];
    else
      _atype[i] = 0;     // we assume that an atom type of 0 is invalid
  }

  return 1;
}

int
Circular_Fingerprint_Generator::generate_fingerprint(const Molecule & m,
                                                     const int * atype,
                                                     const int * include_these_atoms,
                                                     Sparse_Fingerprint_Creator & sfc)
{
  if (! _initialise_molecule (m, atype, include_these_atoms))
    return 0;

  for (int i = 0; i < _matoms; i++)
  {
    unsigned int e = _atype[i];

    if (0 == e)
      continue;

#ifdef DEBUG_ECFP_BITS
    if (0 == min_shell_radius)
      cerr << "Starting with atom " << i << " bit " << e << endl;
#endif

    sfc.hit_bit(e);

    set_vector(_processing_status, _matoms, 0);

    _processing_status[i] = PROCESSING_FINISHED;

    const Atom * ai = _atom[i];

    const int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = ai->other(i, j);

      _processing_status[k] = READY_TO_PROCESS;
    }

    _generate_shells(0, e, sfc);
  }

  return sfc.nset();
}

void
Circular_Fingerprint_Generator::_increment(unsigned int & sum_so_far,
          int bc,
          int atom_constant) const
{
  if (_additive)
    sum_so_far += bc * atom_constant;
  else
    sum_so_far *= bc * atom_constant;

  return;
}

static int 
bond_constant (const Bond * bondi)
{
  if (all_bonds_same_type)
    return 1;

  if (bondi->is_aromatic())
    return 11;
  if (bondi->is_triple_bond())
    return 7;
  if (bondi->is_double_bond())
    return 5;
  
  return 3;
}

int
Circular_Fingerprint_Generator::_generate_shells (int radius,
                                                  unsigned int sum_so_far,
                                                  Sparse_Fingerprint_Creator & sfc)
{
  radius++;

  if (_additive)
    sum_so_far *= 7879;   // an arbitrary prime number

//#define DEBUG_ECFP_BITS 1
#ifdef DEBUG_ECFP_BITS
  cerr << "On entry sum_so_far " << sum_so_far << " radius " << radius << endl;
#endif

// Check for tail addition outside the loop

  for (int i = 0; i < _matoms; i++)
  {
    if (READY_TO_PROCESS != _processing_status[i])
      continue;

    const Atom * ai = _atom[i];

    const int acon = ai->ncon();

    const Bond * const * bonds = ai->rawdata();    // for efficiency

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = bonds[j];

      const atom_number_t k = b->other(i);

      if (0 == _atype[k])      // not included in fingerprint
        continue;

      if (PROCESSING_FINISHED == _processing_status[k])   // we are extending the shell
      {
        const int bc = bond_constant(b);
        _increment(sum_so_far, bc, _atype[i]);
//      cerr << "at level " << radius << " added atom " << i << endl;
      }
      else if (READY_TO_PROCESS == _processing_status[k])
        ;
      else
        _processing_status[k] = NEXT_TIME;
    }
  }

  if (radius >= _min_radius)
  {
#ifdef DEBUG_ECFP_BITS
    cerr << "Hit bit " << sum_so_far << " at radius " << radius << endl;
#endif
    sfc.hit_bit(sum_so_far);
  }

  if (_max_radius > 0 && radius >= _max_radius)
    return 1;

  int continue_processing = 0;

  for (int i = 0; i < _matoms; i++)
  {
    if (READY_TO_PROCESS == _processing_status[i])
      _processing_status[i] = PROCESSING_FINISHED;
    else if (NEXT_TIME == _processing_status[i])
    {
      _processing_status[i] = READY_TO_PROCESS;
      continue_processing = 1;
    }
  }

  if (! continue_processing)
    return 1;

  return _generate_shells(radius, sum_so_far, sfc);
}
