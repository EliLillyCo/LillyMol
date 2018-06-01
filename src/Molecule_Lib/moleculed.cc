#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
using std::cerr;
using std::endl;

#include "iwminmax.h"
#include "misc.h"

// get the private functions

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#define COMPILING_MOLECULED
#define COMPILING_SMILES_CC

class Smiles_First_Atom;
class CRDM_args;

#include "molecule.h"
#include "pearlman.h"
#include "path.h"
#include "misc2.h"

static int full_distance_matrix = 1;

void
set_full_distance_matrix(const int s)
{
  full_distance_matrix = s;
}

int
Molecule::_initialise_distance_matrix()
{
  assert(NULL == _distance_matrix);
  assert(_number_elements > 0);

  _distance_matrix = new_int(_number_elements * _number_elements);

  return 1;
}

/*
  We store the distance matrix in "mostly" upper triangular form,
  so this central routine assumes a2 > a1
*/

int
Molecule::_bonds_between (atom_number_t a1, atom_number_t a2)
{
  assert(a1 < a2);

  if (NULL == _distance_matrix)
    _initialise_distance_matrix();

  int * row = &_distance_matrix[_number_elements * a1];

//#define DEBUG_BONDS_BETWEEN
#ifdef DEBUG_BONDS_BETWEEN
  cerr << "Qbonds_between: between " << a1 << " and " << a2 << endl;

  int precision = 2;
  if (_number_elements > 9)
    precision = 3;
  else
    precision = 4;

  cerr << "MX is ";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << setw(precision) << row[i];
  }
  cerr << endl;
#endif

// If it is already known, or is on the leading edge, grab it.

//cerr << "Atom " << a1 << " to " << a2 << " row = " << row[a2] << endl;

  if (row[a2] > 0)
    return row[a2];
  else if (row[a2] < 0)
    return - row[a2];

//cerr << "Yipes, DM incomplete, atoms " << a1 << " and " << a2 << ", d = " << row[a2] << " continuing...\n";

// The distance has not been computed. Work it out. Identify the
// most positive, negative distance along the row (if present).

  const int invalid_dist_value = - (nedges() + 1);   // longer than longest path in molecule

  int dist = invalid_dist_value;
  for (int i = 0; i < _number_elements; i++)
  {
    if (row[i] < 0 && row[i] > dist)
      dist = row[i];
  }

// If there were no negative numbers along the row, initialise some

  if (invalid_dist_value == dist)
  {
    const Atom * a = _things[a1];

    int a1con = a->ncon();
    for (int i = 0; i < a1con; i++)
    {
      atom_number_t j = a->other(a1, i);
      row[j] = -1;
    }
    dist = -1;

    if (-1 == row[a2])
      return 1;
  }

  int nb = nedges();

  while (1)
  {
    int positive_dist = - dist;
    for (int i = 0; i < _number_elements; i++)
    {
      if (row[i] != dist)
        continue;

#ifdef DEBUG_BONDS_BETWEEN
      cerr << "At distance " << dist << " processing " << i << endl;
#endif

      const Atom * a = _things[i];

      int icon = a->ncon();
      for (int j = 0; j < icon; j++)
      {
        atom_number_t k = a->other(i, j);

#ifdef DEBUG_BONDS_BETWEEN
        cerr << "  Attached to atom " << j << " current = " << row[k] << endl;
#endif       
        if (0 == row[k])
          row[k] = dist - 1;
        else if (row[k] > positive_dist + 1)
          row[k] = dist - 1;
        else if (row[k] < dist - 1)
          row[k] = dist - 1;

        if (k == a2)
          return positive_dist + 1;
      }

      row[i] = positive_dist;    // only change it when all connections processed.
    }

    dist--;
    if (- dist > nb)
    {
      cerr << "Fatal, error, cannot find dist " << a1 << " to " << a2 << ' ' << smiles() << endl;
      cerr << "Fragments " << fragment_membership(a1) << " and " << fragment_membership(a2) << endl;
      debug_print(cerr);
      iwabort();
    }
  }
}

int
Molecule::bonds_between (atom_number_t a1, atom_number_t a2)
{
  if (a1 == a2)
    return 0;

//cerr << "Molecule::bonds_between: atoms " << a1 << " and " << a2 << " dm = " << _distance_matrix << endl;

  if (NULL == _distance_matrix)
    _initialise_distance_matrix();
  else
  {
    if (a1 > a2)
      return _bonds_between(a2, a1);
    else
      return _bonds_between(a1, a2);
  }

// The atoms must be in the same fragment

  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

//assert(_fragment_information.fragment_membership (a1) == _fragment_information.fragment_membership (a2));
  if( _fragment_information.fragment_membership(a1) != _fragment_information.fragment_membership(a2) )
    return ATOMS_NOT_BONDED;

  if (a1 > a2)
    return _bonds_between(a2, a1);
  else
    return _bonds_between(a1, a2);
}

//#define DEBUG_ATOMS_BETWEEN

int
Molecule::atoms_between (atom_number_t a1,
                         atom_number_t a2,
                         Set_of_Atoms & s)
{
  int d = bonds_between(a1, a2);

  if (1 == d)
  {
    s.resize(0);
    return 0;
  }

  if (s.number_elements())
    s.resize_keep_storage(0);
  else
    s.resize(d - 1);

#ifdef DEBUG_ATOMS_BETWEEN
  cerr << "Molecule::atoms_between:atoms " << a1 << " '" << smarts_equivalent_for_atom(a1) << "' and " << a2 << " '" << smarts_equivalent_for_atom(a2) << " are " << d << " bonds apart\n";
#endif

  return _atoms_between(a1, a2, d - 1, s);
}

int
Molecule::_atoms_between (atom_number_t a1,
                          atom_number_t a2,
                          int distance_needed,
                          Set_of_Atoms & s)
{
#ifdef DEBUG_ATOMS_BETWEEN
  cerr << "Molecule::_atoms_between:contine to atom " << a1 << " '" << smarts_equivalent_for_atom(a1) << "' and " << a2 << " '" << smarts_equivalent_for_atom(a2) << "' d = " << distance_needed << endl;
#endif

  assert(NULL != _distance_matrix);

  const Atom * a = _things[a1];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a1, i);

    if (j == a2)     // done, we got to A2
      return 1;

#ifdef DEBUG_ATOMS_BETWEEN
    cerr << "Distance between " << j << " and " << a2 << " is " << _distance_matrix[j * _number_elements + a2] << endl;
#endif

    int d;
    if (j < a2)
      d = _bonds_between(j, a2);
    else
      d = _bonds_between(a2, j);

    if (d != distance_needed)
      continue;

    s.add(j);

    return 1 + _atoms_between(j, a2, distance_needed - 1, s);
  }

  cerr << "Molecule::_atoms_between:yipes, from " << a1 << " '" << smarts_equivalent_for_atom(a1) << "' nothing " << distance_needed << " bonds to " << a2 << " '" << smarts_equivalent_for_atom(a2) << "'\n";
  iwabort();

  return 0;
}

int
Molecule::longest_path()
{
  iwmax<int> rc(0);

  for (int i = 0; i < _number_elements; i++)
  {
    for (int j = i + 1; j < _number_elements; j++)
    {
      rc.extra(_bonds_between(i, j));
    }
  }

  return rc.maxval();
}

//#define DEBUG_COMPUTE_ROW_DM

/*
  We have two versions of _compute_row_of_distance_matrix.
  This one here is designed for speed, and does no ring perception
  or anything like that.
*/

/*void
Molecule::_compute_row_of_distance_matrix (int * row_of_distance_matrix,
                                 atom_number_t current_atom,
                                 int distance)
{
#ifdef DEBUG_COMPUTE_ROW_DM
  cerr << "_compute_row_of_distance_matrix: atom " << current_atom << " (" << atomic_symbol (current_atom) << ") at distance " << distance << endl;
#endif

  row_of_distance_matrix[current_atom] = distance;

  const Atom * a = _things[current_atom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (current_atom, i);

    if (row_of_distance_matrix[j] < 0 || row_of_distance_matrix[j] > distance)
      _compute_row_of_distance_matrix (row_of_distance_matrix, j, distance + 1);
  }

  return;
}*/

/*
  this version is relatively slow compared to current version
*/

#ifdef SLOW_BUT_AVOIDS_CT_ARRAY
void
Molecule::_compute_distance_matrix()
{
  int nb = _bond_list.number_elements();

  set_vector(_distance_matrix, _number_elements * _number_elements, _number_elements + _number_elements);

  const Bond * const * allbonds = _bond_list.rawdata();

  for (int i = 0; i < _number_elements; i++)   // diagonals are zero
  {
    _distance_matrix[i * _number_elements + i] = 0;
  }

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = allbonds[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    int * row1 = _distance_matrix + (a1 * _number_elements);
    int * row2 = _distance_matrix + (a2 * _number_elements);

    row1[a2] = 1;
    row2[a1] = 1;
  }

  while (1)
  {
    int keep_going = 0;

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = allbonds[i];

      atom_number_t a1 = b->a1();
      atom_number_t a2 = b->a2();

      int * row1 = _distance_matrix + (a1 * _number_elements);
      int * row2 = _distance_matrix + (a2 * _number_elements);

      int tmp;
      for (int j = 0; j < _number_elements; j++)
      {
        tmp = row2[j] + 1;
        if (tmp < row1[j])
        {
          row1[j] = tmp;
          keep_going = 1;
        }
        tmp = row1[j] + 1;
        if (tmp < row2[j])
        {
          row2[j] = tmp;
          keep_going = 1;
        }
      }
    }

#ifdef DEBUG_DISTANCE_MATRIX
    cerr << "End of cycle, kg " << keep_going << endl;
    for (int j = 0; j < _number_elements; j++)
    {
      cerr << " Atom " << j << ':';
      const int * r = _distance_matrix + (j * _number_elements);

      for (int k = 0; k < _number_elements; k++)
      {
        cerr << ' ' << r[k];
      }
      cerr << endl;
    }
#endif

    if (0 == keep_going)
      break;
  }

  return;

//when I was writing this I wanted a check with new and old versions
#ifdef CHECK_DM_VS_NEW

  int * tmp = new int[_number_elements*_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);
  copy_vector(tmp, _distance_matrix, _number_elements*_number_elements);
  _compute_distance_matrix_large_mem();
  for (int i = 0; i < _number_elements; ++i)
  {
    for (int j = 0; j < _number_elements; ++j)
    {
      if (tmp[i * _number_elements + j] != _distance_matrix[i * _number_elements + j])
      {
        cerr << "Distance matrix mismatch, atoms " << i << " and " << j << " old " << tmp[i * _number_elements + j] << " new " << _distance_matrix[i * _number_elements + j] << endl;
      }
    }
  }
#endif

}
#endif

/*
  Faster larger memory version
*/

//#define PREVIOUS_BEST_VERSION
#ifdef PREVIOUS_BEST_VERSION
void
Molecule::_compute_distance_matrix()
{
  const int distance_not_set = _number_elements + _number_elements;

  set_vector(_distance_matrix, _number_elements * _number_elements, distance_not_set);

  int maxcon = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    _distance_matrix[i * _number_elements + i] = 0;    // diagonals are zero

    if (_things[i]->ncon() > maxcon)
      maxcon = _things[i]->ncon();
  }

  maxcon++;

  int * ct = new int[_number_elements * maxcon]; std::unique_ptr<int[]> free_ct(ct);

  for (int i = 0; i < _number_elements; ++i)
  {
    ct[i * maxcon] = _things[i]->connections(i, ct + i * maxcon + 1);
  }

  for (auto i = 0; i < _number_elements; ++i)
  {
//  _distance_matrix[i * _number_elements + i] = 0;    already done above

    const auto icon = ct[i * maxcon];
    const int * cti = ct + i * maxcon + 1;

    for (auto j = 0; j < icon; ++j)
    {
      const auto k = cti[j];

      _distance_matrix[i * _number_elements + k] = 1;
      _distance_matrix[k * _number_elements + i] = 1;
    }
  }

  for (auto d = 2; d < _number_elements; ++d)
  {
    bool keep_going = false;

    for (auto i = 0; i < _number_elements; ++i)
    {
      int * distance_to_i  = _distance_matrix + i * _number_elements;

      for (auto j = 0; j < _number_elements; ++j)
      {
        if (j == i)
          continue;

        if (distance_to_i[j] > _number_elements)    // obviously not set
          continue;

        const auto jcon = ct[j * maxcon];
        const int * ctj = ct + j * maxcon + 1;

        for (auto k = 0; k < jcon; ++k)
        {
          const auto l = ctj[k];

          if (distance_to_i[j] + 1 < distance_to_i[l])
          {
            distance_to_i[l] = distance_to_i[j] + 1;
            _distance_matrix[l * _number_elements + i] = distance_to_i[l];
            keep_going = true;
          }
        }
      }
    }

    if (! keep_going)
      break;
  }

  return;
}
#endif

#define VERSION_WITH_FRONTIER
#ifdef VERSION_WITH_FRONTIER
static void
print_frontier(const int zatom,
               const int d,
               const int * frontier,
               const int atoms_on_frontier)
{
  cerr << "Frontier information relative to atom " << zatom << ", dist " << d << " " << atoms_on_frontier << " on frontier\n";

  for (int i = 0; i < atoms_on_frontier; ++i)
  {
    cerr << ' ' << frontier[i];
  }

  cerr << endl;
  
  return;
}

void
Molecule::_compute_distance_matrix()
{
  const int distance_not_set = _number_elements + _number_elements;

  set_vector(_distance_matrix, _number_elements * _number_elements, distance_not_set);

  int maxcon = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    _distance_matrix[i * _number_elements + i] = 0;    // diagonals are zero

    if (_things[i]->ncon() > maxcon)
      maxcon = _things[i]->ncon();
  }

  maxcon++;

  int * ct = new int[_number_elements * maxcon]; std::unique_ptr<int[]> free_ct(ct);

  for (int i = 0; i < _number_elements; ++i)
  {
    ct[i * maxcon] = _things[i]->connections(i, ct + i * maxcon + 1);
  }

  int * frontier = new int[2 * _number_elements]; std::unique_ptr<int[]> free_frontier(frontier);
  int * next_frontier = frontier + _number_elements;

// rather than copy data from arrays, we just swap pointers...

  int * p_frontier = frontier;
  int * p_next_frontier = next_frontier;

  for (int i = 0; i < _number_elements; ++i)
  {
    int * distance_to_i = _distance_matrix + i * _number_elements;

    p_frontier[0] = i;
    int atoms_on_current_frontier = 1;

    for (int d = 1; d < _number_elements; ++d)
    {
      int atoms_on_next_frontier = 0;
      for (int j = atoms_on_current_frontier - 1; j >= 0; --j)
      {
        const atom_number_t k = p_frontier[j];

        const int * ctk = ct + k * maxcon + 1;

        for (auto l = ct[k*maxcon] - 1; l >= 0; --l)
        {
          const auto x = ctk[l];

          if (distance_to_i[k] + 1 < distance_to_i[x])
          {
            distance_to_i[x] = d;

            p_next_frontier[atoms_on_next_frontier] = x;
            atoms_on_next_frontier++;
          }
        }
      }

      if (0 == atoms_on_next_frontier)
        break;

      std::swap(p_frontier, p_next_frontier);
      atoms_on_current_frontier = atoms_on_next_frontier;
    }
  }

  return;
}
#endif

//#define THIS_IS_WAY_SLOWER_NOT_SURE_WHY
#ifdef THIS_IS_WAY_SLOWER_NOT_SURE_WHY

void
Molecule::_compute_distance_matrix()
{
  const int distance_not_set = _number_elements + _number_elements;

  std::fill_n(_distance_matrix, _number_elements * _number_elements, distance_not_set);

  for (int i = 0; i < _number_elements; ++i)
  {
    _distance_matrix[i * _number_elements + i] = 0;    // diagonals are zero
  }

  int * to_check = new_int(_number_elements + _number_elements); std::unique_ptr<int[]> free_to_check(to_check);

  const int ne = nedges();

  int * b1 = new int[ne + ne]; std::unique_ptr<int[]> free_b1(b1);
  int * b2 = b1 + ne;

  for (int i = 0; i < ne; ++i)
  {
    const Bond * bi = _bond_list[i];

    const atom_number_t a1 = bi->a1();
    const atom_number_t a2 = bi->a2();

    if (a1 < a2)
    {
      b1[i] = a1;
      b2[i] = a2;
    }
    else
    {
      b1[i] = a2;
      b2[i] = a1;
    }

    to_check[a1] = 1;
    to_check[a2] = 1;

    _distance_matrix[a1 * _number_elements + a2] = 1;
    _distance_matrix[a2 * _number_elements + a1] = 1;
  }

  int * next_to_check = to_check + _number_elements;

  for (auto d = 2; d < _number_elements; ++d)
  {
    cerr << " d = " << d << endl;

    bool keep_going = false;

    for (int i = 0; i < ne; ++i)
    {
      const atom_number_t a1 = b1[i];
      const atom_number_t a2 = b2[i];

      for (int j = 0; j < _number_elements; ++j)
      {
        if (d - 1 == _distance_matrix[a1 * _number_elements + j])
        {
          if (_distance_matrix[j * _number_elements + a2] > d)
          {
            _distance_matrix[j * _number_elements + a2] = d;
            _distance_matrix[a2 * _number_elements + j] = d;
            keep_going = true;
          }
        }
        if (d - 1 == _distance_matrix[j * _number_elements + a2])
        {
          if (_distance_matrix[a1 * _number_elements + j] > d)
          {
            _distance_matrix[j * _number_elements + a1] = d;
            _distance_matrix[a1 * _number_elements + j] = d;
            keep_going = true;
          }
        }
      }
    }

    if (! keep_going)
      break;
  }

  return;
}
#endif


/*
  _compute_row_of_distance_matrix has identified a ring closure.
  It is between the top of stack and atom ASTOP.
*/

/*static void
found_ring (atom_number_t astop,
            const int * atom_stack,
            int stack_ptr,
            int * ring_atom,
            int * row_of_distance_matrix)
{
  int distance = row_of_distance_matrix[astop] + 1;

  while (stack_ptr)
  {
    atom_number_t a = atom_stack[stack_ptr--];

    ring_atom[a]++;

    if (a == astop)
      return;

    if (row_of_distance_matrix[a] > distance)
      row_of_distance_matrix[a] = distance;

    distance++;

  }

  return;
}*/

#define DEBUG_CRDM

/*
  Fundamental algorithm. We do a depth first traversal through the
  molecule in order to identify ring atoms.

  To avoid passing a lot of arguments, we put them in a class
*/

class CRDM_args
{
//friend
//  void Molecule::_compute_row_of_distance_matrix (CRDM_args &,
//             int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &));
  private:
    int * _row_of_distance_matrix;
    int   _distance;
    int * _atom_stack;
    int   _stack_ptr;
    int   _smiles_order;

    resizable_array_p<Ring> _raw_rings_found;

//  private functions

    int _is_sssr_ring (const Ring * r, const int * ring_membership);

  public:
    CRDM_args (int);
    ~CRDM_args();

    int debug_print (std::ostream &) const;

    void push (atom_number_t a) { _atom_stack[_stack_ptr++] = a;}
    void pop() { _stack_ptr--;}

    atom_number_t current_atom() const { return _atom_stack[_stack_ptr - 1];}

    int smiles_order();

    void set_row (int * r) { _row_of_distance_matrix = r;}

    void found_ring (atom_number_t, int *);

    int  look_for_sssr_rings (int * ring_membership,
             resizable_array_p<Ring> & sssr_rings,
             resizable_array_p<Ring> & raw_rings);
};

CRDM_args::CRDM_args (int matoms)
{
  _distance = 0;
  _atom_stack = new int[matoms];
  _stack_ptr = 0;
  _smiles_order = 0;
}

CRDM_args::~CRDM_args()
{
  DELETE_IF_NOT_NULL(_atom_stack);
}

int
CRDM_args::debug_print (std::ostream & os) const
{
  os << "details on CRDM_args object, smiles order " << _smiles_order << " distance " << _distance << endl;

  if (0 == _stack_ptr)
    os << "Stack empty\n";
  for (int i = 0; i < _stack_ptr; i++)
  {
    os << "Stack level " << i << " atom " << _atom_stack[i] << endl;
  }

  return os.good();
}

int
CRDM_args::smiles_order()
{
  int rc = _smiles_order;

  _smiles_order++;

  return rc;
}

#define DEBUG_FOUND_RING

void
CRDM_args::found_ring (atom_number_t astop,
                       int * ring_membership)
{
  int distance = _row_of_distance_matrix[astop] + 1;

  Ring * s = new Ring;
  s->resize(8);

  _raw_rings_found.add(s);

#ifdef DEBUG_FOUND_RING
  cerr << "Found new ring\n";
#endif

  for (int i = _stack_ptr - 1; i >= 0; i--)
  {
    atom_number_t a = _atom_stack[i];

    ring_membership[a]++;
    if (s->array_is_full())
      s->resize(s->number_elements() + 5);

    s->add(a);

    if (a == astop)
      return;

    if (_row_of_distance_matrix[a] > distance)
      _row_of_distance_matrix[a] = distance;

#ifdef DEBUG_FOUND_RING
    cerr << "i = " << i << " ring atom " << a << " distance set to " << _row_of_distance_matrix[a] << endl;
#endif

    distance++;
  }

  return;
}

int
CRDM_args::_is_sssr_ring (const Ring * r, const int * ring_membership)
{
  int ring_size = r->number_elements();
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r->item(i);

    assert(ring_membership[j] > 0);

    if (ring_membership[j] > 1)
      return 0;
  }

  return 1;
}

int
CRDM_args::look_for_sssr_rings (int * ring_membership,
             resizable_array_p<Ring> & sssr_rings,
             resizable_array_p<Ring> & raw_rings)
{
  int rc = 0;

  int nr = _raw_rings_found.number_elements();
  for (int i = nr - 1; i >= 0; i--)
  {
    Ring * r = _raw_rings_found[i];
    if (_is_sssr_ring(r, ring_membership))
    {
      _raw_rings_found.remove_no_delete(i);
      sssr_rings.add(r);
      rc++;
    }
  }

  for (int i = 0; i < _raw_rings_found.number_elements(); i++)
  {
    const Ring * r = _raw_rings_found[i];
    r->set_vector(ring_membership, IW_RING_MEMBERSHIP_IS_A_RING_ATOM);
  }

  raw_rings.transfer_in(_raw_rings_found);

  return rc;
}

/*void
Molecule::_compute_row_of_distance_matrix (CRDM_args & crdm,
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &))
{
  atom_number_t current_atom = crdm.current_atom();

  int & distance = crdm._distance;

  int * row_of_distance_matrix = crdm._row_of_distance_matrix;

  row_of_distance_matrix[current_atom] = distance;

  Atom * a = _things[current_atom];

//_smiles_order[current_atom] = crdm.smiles_order();

#ifdef DEBUG_CRDM
  cerr << "_compute_row_of_distance_matrix: atom " << current_atom << " (" << atomic_symbol (current_atom) << ") at distance " << distance << endl;
#endif

  int acon = a->ncon();
  if (1 == acon)
  {
    crdm.pop();
    return;
  }

  int initial_distance = distance;    // save our distance

  distance++;

// First check for any ring closures

  int unprocessed_connections = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (current_atom, i);

#ifdef DEBUG_CRDM
    cerr << "From atom " << current_atom << " try atom " << j << " dist = " << row_of_distance_matrix[j] << endl;
#endif

    if (row_of_distance_matrix[j] < initial_distance - 1)
    {
      crdm.found_ring (j, _ring_membership);
      distance = row_of_distance_matrix[j] + 1;
//    _ring_closure_bonds.add (current_atom, j);
    }
    else if (row_of_distance_matrix[j] > _number_elements)
      unprocessed_connections++;
  }

  if (0 == unprocessed_connections)    // all done with this atom
    return;

// Now continue building the smiles order

  atom_number_t next_atom;
  while ((this->*identify_next_atom) (_smiles_order, current_atom, next_atom))
  {
//  _bonds_in_fragment[_number_fragments]++;

#ifdef DEBUG_CRDM
    cerr << "From atom " << current_atom << " to atom " << next_atom << " dist = " << row_of_distance_matrix[next_atom] << endl;
#endif

    crdm.push (next_atom);
    crdm.debug_print (cerr);
    _compute_row_of_distance_matrix (crdm, identify_next_atom);
  }

  distance--;

  crdm.pop();

  return;
}*/

/*int
Molecule::_recompute_distance_matrix (int (Molecule::*identify_first_atom) (const int *, atom_number_t &),
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &))
{
  cerr << "Entering _recompute_distance_matrix\n";
  if (NULL == _smiles_order)
    _smiles_order = new_int (_number_elements, -1);
  else
    set_vector (_smiles_order, _number_elements, -1);

  if (NULL == _ring_membership)
    _ring_membership = new_int (_number_elements);
  else
    set_vector (_ring_membership, _number_elements, 0);

  if (NULL == _distance_matrix)
    _distance_matrix = new_int (_number_elements * _number_elements, _number_elements + 9);
  else
    set_vector (_distance_matrix, _number_elements * _number_elements, _number_elements + 9);

  CRDM_args crdm_args (_number_elements);

  atom_number_t a;
  while ((this->*identify_first_atom) (_smiles_order, a))
  {
    _smiles_order[a] = crdm_args.smiles_order();

    crdm_args.push (a);
    crdm_args.debug_print (cerr);
    crdm_args.set_row (_distance_matrix + a * _number_elements);

    _compute_row_of_distance_matrix (crdm_args, identify_next_atom);
  }

  crdm_args.look_for_sssr_rings (_ring_membership, _sssr_rings, _raw_rings);

// Now process all remaining rows of the distance matrix the fast way

  int singly_connected_atoms = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    int * row_of_distance_matrix = _distance_matrix + i * _number_elements;
    if (row_of_distance_matrix[0] || row_of_distance_matrix[1])
      continue;

    if (1 == _things[i]->ncon())
    {
      singly_connected_atoms++;
      continue;
    }

    int d = 0;   // passed by reference

    _compute_row_of_distance_matrix (row_of_distance_matrix, i, d);
  }

  if (singly_connected_atoms)
  {
    for (int i = 0; i < _number_elements && singly_connected_atoms; i++)
    {
      const Atom * a = _things[i];
      if (1 != a->ncon())
        continue;

      int * ri = _distance_matrix + _number_elements * i;
      if (ri[0] || ri[1])
        continue;

      atom_number_t j = a->other (i, 0);

      int * rj = _distance_matrix + _number_elements * j;

      for (int k = 0; k < _number_elements; k++)
      {
        ri[k] = rj[k] + 1;
      }

      ri[i] = 0;
      singly_connected_atoms--;
    }
  }

  return 1;
}*/

//#define CHECK_DISTANCE_MATRIX

int
Molecule::recompute_distance_matrix()
{
  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

  if (NULL == _distance_matrix)
    _distance_matrix = new_int(_number_elements * _number_elements, _number_elements + 9);
  else
    set_vector(_distance_matrix, _number_elements * _number_elements, _number_elements + 9);

  if (1 == _number_elements)
    return 1;

  if (2 == _number_elements)
  {
    _distance_matrix[1] = 1;
    _distance_matrix[3] = 1;

    return 1;
  }

//cerr << "_compute_distance_matrix....\n";

  _compute_distance_matrix();

  int rc = 1;

#ifdef CHECK_DISTANCE_MATRIX
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 != _distance_matrix[_number_elements * i + i])
    {
      cerr << "Non zero diagonal on distance matrix, i = " << i << " value = " << _distance_matrix[_number_elements * i + i] << endl;
      rc = 0;
    }

    for (int j = i + 1; j < _number_elements; j++)
    {
      if (_distance_matrix[_number_elements * i + j] != _distance_matrix[_number_elements * j + i])
      {
        cerr << "Distance matrix error, from " << i << " to " << j << " is " << _distance_matrix[_number_elements * i + j] << endl;
        cerr << "Distance matrix error, from " << j << " to " << i << " is " << _distance_matrix[_number_elements * j + i] << endl;
        rc = 0;
      }
      else if (0 == _distance_matrix[_number_elements * i + j])
      {
        cerr << "Zero distance matrix entry " << i << ", " << j << endl;
        rc = 0;
      }
      else if (_distance_matrix[_number_elements * i + j] < 0)
      {
        cerr << "Incomplete distance matrix determinations, atoms " << i << " and " << j << endl;
      }
      else
      {
//      cerr << "Bonds between " << i << " and " << j << " is " << _distance_matrix[_number_elements * i + j] << endl;
      }
    }
  }

  if (0 == rc)
    iwabort();
#endif

  return rc;
}

// arch-tag: c8079359-e517-4d71-8303-1697d2ab7377
