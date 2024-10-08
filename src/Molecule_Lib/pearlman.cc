/*
  Ring finder just like Bob Pearlman's
*/

#include <iostream>
#include <memory>
#include <iomanip>

using std::cerr;
using std::endl;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "Foundational/iwmisc/misc.h"
#include "assert.h"

#define COMPILING_MOLECULER_CC

#include "molecule.h"
#include "path.h"
#include "misc2.h"

#include "pearlman.h"

static int global_setting_perceive_sssr_rings = 1;

int
perceive_sssr_rings()
{
  return global_setting_perceive_sssr_rings;
}

void
set_perceive_sssr_rings(int s)
{
  global_setting_perceive_sssr_rings = s;
}

static int file_scope_accumulate_non_sssr_rings = 1;

void
set_accumulate_non_sssr_rings(int s)
{
  file_scope_accumulate_non_sssr_rings = 1;
}

int
accumulate_non_sssr_rings()
{
  return file_scope_accumulate_non_sssr_rings;
}

void
Beep::_default_values()
{
  _pi_electrons = 0;

  _atomic_number_score = 0;

  _single_bond_count = 0;
}

Beep::Beep()
{
  _default_values();

  return;
}

Beep::Beep(int bonds_in_molecule) : IW_Bits_Base(bonds_in_molecule)
{
  _default_values();

  return;
}

int
Beep::ok() const
{
  return IW_Bits_Base::ok();
}

int
Beep::debug_print(std::ostream & os) const
{
  os << "Info on Beep with " << nbits() << " bits, " << _pi_electrons << " pi electrons\n";

  if (! ok())
    os << "Warning, OK function fails\n";

  for (int i = 0; i < nbits(); i++)
  {
    if (is_set(i))
      os << " Bond " << i << " set\n";
  }

  return 1;
}

/*
  The _default_values function needs some care because of the two
  constructors. See which variables are set in each case, and make
  sure we don't overwrite any values set in the constructor.
*/

void
Path_Message::_default_values()
{
  _to_be_removed = 0;

  _last_atom = _start_atom;    // to ensure this message is not sent back

  return;
}

Path_Message::Path_Message(atom_number_t a,
                           const int peosa,
                           const int saans,
                           const int ee,
                           const int bonds_in_molecule) :
              Beep(bonds_in_molecule),
              _start_atom(a),
              _pi_electrons_on_start_atom(peosa),
              _start_atom_atomic_number_score(saans),
              _start_edge(ee),
              _bonds_in_molecule(bonds_in_molecule)
{
  _default_values();

  _pi_electrons = _pi_electrons_on_start_atom;

  _atomic_number_score = saans;

//cerr << "Path_Message::Path_Message:_start_atom_atomic_number_score " << _start_atom_atomic_number_score << endl;

  set(_start_edge, 1);

  return;
}

Path_Message::Path_Message(const Path_Message * p) :
                         Beep(p->_bonds_in_molecule),
                         _start_atom(p->_start_atom),
                         _pi_electrons_on_start_atom(p->_pi_electrons_on_start_atom),
                         _start_atom_atomic_number_score(p->_start_atom_atomic_number_score),
                         _start_edge(p->_start_edge),
                         _bonds_in_molecule(p->_bonds_in_molecule)
{
  _default_values();

  _pi_electrons = p->_pi_electrons;
  _atomic_number_score = p->_atomic_number_score;

  int bytes_to_copy = _whole_bytes;
  if (_extra_bits)
    bytes_to_copy++;

  memcpy(_bits, p->_bits, bytes_to_copy);

//cerr << "Path_Message::copy constructor:_start_atom_atomic_number_score " << _start_atom_atomic_number_score << endl;

  return;
}

int
Path_Message::debug_print(std::ostream & os)
{
  assert(os.good());

  os << "Info on Path message starting with atom " << _start_atom << " (score " << _start_atom_atomic_number_score << ") edge " <<
        _start_edge << ", " << nset() << " bits set, score " << _atomic_number_score << endl;

  (void) Beep::debug_print(os);

  os << "Last atom is " << _last_atom << endl;

  return os.good();
}

void
Path_Message::include_in_path(const atom_number_t a, const int zbond,
                               const int pe, const int h)
{
  assert(! is_set(zbond));
  assert(pe >= 0);

  _last_atom = a;
  _pi_electrons += pe;
  _atomic_number_score += h;

//cerr << "_atomic_number_score incremented by " << h << " to " << _atomic_number_score << endl;

  set(zbond, 1);

  return;
}

int
Path_Message::node_collision(const Path_Message * p) const
{
  if (_start_atom == p->_start_atom && _start_edge != p->_start_edge)
    ;
  else
    return 0;

  return 0 == bits_in_common(*p);
}

int
Path_Message::inverse_edge_collision(const Path_Message * p) const
{
  if (_start_atom != p->_start_atom && _start_edge == p->_start_edge)
    ;
  else
    return 0;

  return 1 == bits_in_common(*p);
}

//#define DEBUG_PATH_MESSAGES

Rings_Found::Rings_Found(int nr, int nb) :
                                    _expected_nrings(nr),
                                    _bonds_in_molecule(nb)
{
  _matrix_of_beeps = new Beep *[_bonds_in_molecule];
  for (int i = 0; i < _bonds_in_molecule; i++)
  {
    _matrix_of_beeps[i] = nullptr;
  }

// Make these large to avoid resizing

  _inverse_edge_collision_rings.resize(nb);
  _node_collision_rings.resize(nb);
  _sssr_rings_perceived.resize(nb);

  return;
}

Rings_Found::~Rings_Found()
{
  for (int i = 0; i < _bonds_in_molecule; i++)
  {
    if (nullptr != _matrix_of_beeps[i])
      delete _matrix_of_beeps[i];
  }

  delete [] _matrix_of_beeps;

  return;
}

void
Rings_Found::initialise_single_bond_count(const Molecule & m)
{
  _is_single_bond.allocate_space_for_bits(_bonds_in_molecule);

  for (int i = 0; i < _bonds_in_molecule; i++)
  {
    if (m.bondi(i)->is_single_bond())
      _is_single_bond.set(i);
  }

  return;
}

static int
is_duplicate(const resizable_array_p<Beep> & beeps,
             const Beep * b)
{
  int nb = beeps.number_elements();
  for (int i = 0; i < nb; i++)
  {
    if (*(beeps[i]) == *b)
      return 1;
  }

  return 0;
}

/*
  Tnodes send newly found rings to this entry.
  
  If we can reject the new ring, we do so.
  If we can identify it as definitely a new ring, we do so.
  Otherwise, we just place it in the receive buffer for
  subsequent analysis.
*/

//#define DEBUG_RINGS_FOUND

void
Rings_Found::node_collision_ring(Beep * r)
{
#ifdef DEBUG_RINGS_FOUND
  cerr << "Rings_Found::node_collision_ring checking\n";
  r->printon(cerr);
  cerr << endl;
#endif

// If there are other Beep's in the receive buffer, check for this
// being a duplicate

  if (_node_collision_rings.number_elements() && is_duplicate(_node_collision_rings, r))
  {
#ifdef DEBUG_RINGS_FOUND
    cerr << "Duplicate (incoming) ring discarded\n";
#endif

    delete r;
    return;
  }

  _node_collision_rings.add(r);

  return;
}

void
Rings_Found::inverse_edge_collision_ring(Beep * r)
{
#ifdef DEBUG_RINGS_FOUND
  cerr << "Rings_Found::inverse_edge_collision_ring checking\n";
  r->printon(cerr);
  cerr << endl;
#endif

  int nr = _inverse_edge_collision_rings.number_elements();

// If there are other Beep's in the receive buffer, check for this
// being a duplicate

  if (nr && is_duplicate(_inverse_edge_collision_rings, r))
  {
#ifdef DEBUG_RINGS_FOUND
    cerr << "Duplicate (incoming) ring discarded\n";
#endif

    delete r;
    return;
  }

  _inverse_edge_collision_rings.add(r);

  return;
}

/*
  When sorting beeps, we want to favour those which can most likely form
  an aromatic system

  Beeps with zero pi electrons are automatically disfavoured
*/

static int
beep_pi_electron_comparitor(Beep * const * bp1, Beep * const * bp2)
{
  const Beep * b1 = *bp1;
  const Beep * b2 = *bp2;

//#define DEBUG_BEEP_PI_ELECTRON_COMPARITOR
#ifdef DEBUG_BEEP_PI_ELECTRON_COMPARITOR
  b1->printon(cerr);
  cerr << endl;
  b2->printon(cerr);
  cerr << endl;
  cerr << "Comparing pi electrons " << b1->pi_electrons() << " (" << b1->nset() << " bonds) and " << b2->pi_electrons() << " (" << b2->nset() << " bonds)\n";
  b1->debug_print(cerr);
  b2->debug_print(cerr);
#endif

  if (b1->pi_electrons() == b2->pi_electrons())
    return 0;

  if (0 == b1->pi_electrons())
    return 1;

  if (0 == b2->pi_electrons())
    return -1;

  int p1 = b1->pi_electrons() % 4;
  int p2 = b2->pi_electrons() % 4;

  if (p1 == p2)     // favour the one with the most pi electrons
  {
    if (b1->pi_electrons() < b2->pi_electrons())
      return 1;
    else
      return -1;
  }

// favour those with an even number of pi electrons

  int even1 = (p1 < (p1 ^ 1));
  int even2 = (p2 < (p2 ^ 1));

//cerr << "Comparing b1 with " << b1->pi_electrons() << " p1 = " << p1 << " even1 = " << even1 << endl;
//cerr << "With      b2 with " << b2->pi_electrons() << " p2 = " << p2 << " even2 = " << even2 << endl;
//cerr << "p1 ^ 1 is " << (p1 ^ 1) << " p2 ^ 1 is " << (p2 ^ 1) << endl;

  if (even1 && even2)
  {
    if (b1->pi_electrons() < b2->pi_electrons())
      return 1;
    else
      return -1;
  }

  if (even1)
    return -1;

  if (even2)
    return 1;

  if (b1->pi_electrons() < b2->pi_electrons())
    return 1;
  else
    return -1;
}

//#define DEBUG_BEEP_COMPARATOR

static int
beep_comparator(Beep * const * bp1, Beep * const * bp2)
{
  int pe = beep_pi_electron_comparitor(bp1, bp2);

#ifdef DEBUG_BEEP_COMPARATOR
  cerr << "Comparing beeps with " << (*bp1)->nset() << " and " << (*bp2)->nset() << " bonds, pe = " << pe << endl;
#endif

  if (0 != pe)
    return pe;

  const Beep * b1 = *bp1;
  const Beep * b2 = *bp2;

#ifdef DEBUG_BEEP_COMPARATOR
  cerr << "Comparing atomic_number_score " << b1->atomic_number_score() << " and " << b2->atomic_number_score() << endl;
#endif

  if (b1->atomic_number_score() < b2->atomic_number_score())
    return 1;
  else if (b1->atomic_number_score() > b2->atomic_number_score())
    return -1;

  if (b1->single_bond_count() < b2->single_bond_count())
    return 1;
  else if (b1->single_bond_count() > b2->single_bond_count())
    return -1;

  return 0;
}

/*
  When adding non sssr rings, we need to avoid duplicates
*/

int
Rings_Found::_beep_is_unique_over_non_sssr_beeps(const Beep * new_beep) const
{
  int nq = _non_sssr_rings.number_elements();
  if (0 == nq)
    return 1;

#ifdef DEBUG_BEEP_IS_UNIQUE_OVER_NON_SSSR_BEEPS
  cerr << "Checking " << nq << " non sssr beeps\n";
  cerr << "New beep is ";
  new_beep->printon(cerr);
  cerr << endl;
#endif

  for (int i = 0; i < nq; i++)
  {
    const Beep * b = _non_sssr_rings[i];
    if (*new_beep == *b)
      return 0;
  }

  return 1;
}

/*
  This function just checks whether or not a beep is unique. It does
  not do any updating of the internal arrays.
*/

int
Rings_Found::_beep_is_unique(const Beep * b) const
{
  if (! global_setting_perceive_sssr_rings)
    return 1;

  Beep tmp = *b;

  while (1)
  {
    int f = tmp.first_bit();

    assert(f >= 0);

    const Beep * bf = _matrix_of_beeps[f];
//  if (nullptr == bf)
//    cerr << "Matrix entry " << f << " is empty, new ring accepted\n";

    if (nullptr == bf)
      return 1;

    tmp.iwxor(*bf);
//  cerr << "After XOR ";
//  tmp.printon (cerr);
//  cerr << endl;
//  cerr << "nset = " << tmp.nset() << endl;

    if (0 == tmp.nset())
      return 0;
  }
}

//#define DEBUG_FIND_UNIQUE_RINGS

/*
  Determine whether a new ring is a linear combination of existing rings.
  If it is, we add B to the list of perceived rings, and we expand
  the matrix of beeps
*/

int
Rings_Found::_is_sssr_ring(const Beep * b)
{
  Beep * tmp = new Beep;
  *tmp = *b;
//tmp->make_copy(b);

#ifdef DEBUG_FIND_UNIQUE_RINGS
  cerr << "SSSR: Determining uniqueness of beep ";
  b->printon(cerr);
  cerr << ' ' << b->nset() << endl;
#endif

  while (1)
  {
    int f = tmp->first_bit();

    assert(f >= 0);

    const Beep * bf = _matrix_of_beeps[f];

    if (nullptr == bf)
    {
      _matrix_of_beeps[f] = tmp;
      return 1;
    }

    tmp->iwxor(*bf);
//  cerr << "After XOR ";
//  tmp->printon(cerr);
//  cerr << endl;
//  cerr << "nset = " << tmp->nset() << endl;

    if (0 == tmp->nset())
    {
      delete tmp;
      return 0;
    }
  }
}

/*
  Feb 2007. Significant re-think of this for esssr rings.
  Here I don't keep a beep for each bit. The test is to search through all
  the previously found beeps (of smaller size) and see if every bit in this new
  beep is already set in one of the already found rings
*/

int
Rings_Found::_is_esssr_ring(const Beep * b)
{
  Beep tmp = (*b);

#ifdef DEBUG_FIND_UNIQUE_RINGS
  cerr << "ESSSR: Determining uniqueness of beep ";
  b->printon(cerr);
  cerr << ' ' << b->nset() << endl;
#endif

  int nset = b->nset();

  int n = _sssr_rings_perceived.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Beep * bi = _sssr_rings_perceived[i];

    if (bi->nset() == nset)    // only check smaller rings
      continue;

    tmp.unset_bits_in_rhs(*bi);

    if (! tmp.any_bits_set())   // all bits found in existing rings
      return 0;
  }

  return 1;   // There are bonds in B that aren't in any smaller rings.
}

#ifdef NOT_USED
static void
remove_duplicate_beeps(resizable_array_p<Beep> & beeps)
{
  int n = beeps.number_elements();

  cerr << "At start, there are " << n << " beeps\n";

  for (int i = 0; i < n; i++)
  {
    const Beep & bi = *(beeps[i]);

    for (int j = n - 1; j > i; j--)
    {
      if (*(beeps[j]) == bi)
      {
        beeps.remove_item(j);
        j--;
        n--;
      }
    }
  }

  cerr << "At end, there are " << n << " beeps\n";
  return;
}
#endif

void
Rings_Found::_assign_single_bond_counts(Beep * beep)
{
  beep->set_single_bond_count(beep->bits_in_common(_is_single_bond));

  return;
}

void
Rings_Found::_assign_single_bond_counts(resizable_array_p<Beep> & beeps)
{
  int n = beeps.number_elements();

  for (int i = 0; i < n; i++)
  {
    _assign_single_bond_counts(beeps[i]);
  }

  return;
}

//#define DEBUG_DETERMINE_UNIQUENESS

/*
  The array BEEPS will be either the incoming array of inverse edge rings,
  or the incoming array of node collision rings

  We do this as a two pass process in order to be able to save those rings
  which are unique, but then rejected due to SSSR considerations - the extra
  ring in cubane for example.
*/

void
Rings_Found::_determine_uniqueness(resizable_array_p<Beep> & beeps)
{
  int nb = beeps.number_elements();

#ifdef DEBUG_DETERMINE_UNIQUENESS
  cerr << "Rings_Found::_determine_uniqueness: examining " << nb << " beeps\n";
#endif

  for (int i = 0; i < nb; i++)
  {
    const Beep * b = beeps[i];
    if (! _beep_is_unique(b))
    {
//    cerr << "Beep " << i << " not unique, removed\n";
      beeps.remove_item(i);
      i--;
      nb--;
    }
  }

#ifdef DEBUG_DETERMINE_UNIQUENESS
  cerr << "Processing " << nb << " new unique beeps\n";
#endif

  if (0 == nb)
    return;

// Sort these beeps by the pi electron count (possibly aromatic first)

  if (nb > 1)
  {
    if (_expected_nrings > 2)
      _assign_single_bond_counts(beeps);

    beeps.sort(beep_comparator);

#ifdef DEBUG_DETERMINE_UNIQUENESS
    cerr << "After sorting by pi electrons we have\n";
    for (int i = 0; i < nb; i++)
    {
      const Beep * b = beeps[i];
      cerr << "i = " << i << " beep has " << b->nset() << " bits set, " << b->pi_electrons() << " pi electrons, " << b->atomic_number_score() << " heteroatoms\n";
    }
#endif
  }
#ifdef DEBUG_DETERMINE_UNIQUENESS
  else
  {
    const Beep * b = beeps[0];

    cerr << "Single beep has " << b->nset() << " bits set, " << b->pi_electrons() << " pi electrons, " << b->atomic_number_score() << " heteroatoms\n";
  }
  cerr << "Uniqueness checking " << nb << " beeps being added to molecule's rings\n";
#endif

  for (int i = 0; i < nb; i++)
  {
    Beep * b = beeps[i];

    if (perceive_sssr_rings())
    {
      if (_sssr_rings_perceived.number_elements() < _expected_nrings && _is_sssr_ring(b))
        _sssr_rings_perceived.add(b);
      else if (file_scope_accumulate_non_sssr_rings && _beep_is_unique_over_non_sssr_beeps(b))
        _non_sssr_rings.add(b);
      else
        delete b;
    }
    else
    {
      if (_is_esssr_ring(b))
        _sssr_rings_perceived.add(b);
      else if (file_scope_accumulate_non_sssr_rings)
        _non_sssr_rings.add(b);
    }
  }

#ifdef DEBUG_DETERMINE_UNIQUENESS
  cerr << "_determine_uniqueness:_sssr_rings_perceived now contains " << _sssr_rings_perceived.number_elements() << " rings, non sssr " << _non_sssr_rings.number_elements() << endl;
#endif

  beeps.resize_no_delete(0);    // we have grabbed all we want

  return;
}

/*
  At the end of each transmit cycle through the Tnodes, this
  function is called to analyse all the rings encountered so far.
  The only thing we know is that _receive_buffer contains no
  duplicates
*/

void
Rings_Found::process_new_rings(const Molecule & m)
{
  const int nri = _inverse_edge_collision_rings.number_elements();
  const int nrn = _node_collision_rings.number_elements();

#ifdef DEBUG_RINGS_FOUND
  cerr << "Rings found has " << _sssr_rings_perceived.number_elements() << 
          " rings perceived, expecting " << _expected_nrings << endl;
  cerr << nri << " inverse edge collision rings and " <<
          nrn << " node collision rings incoming\n";
#endif

  if (0 == nri && 0 == nrn)     // no new rings found this iteration
    return;  

// Process the inverse_edge_collision (shorter) rings first. Check them
// for uniqueness.

  if (nri)
    _determine_uniqueness(_inverse_edge_collision_rings);

#ifdef DEBUG_FIND_UNIQUE_RINGS
  cerr << "After uniqueness " << _inverse_edge_collision_rings.number_elements() << " inverse edge collision rings\n";
#endif

  if (0 == nrn)   // none to worry about
    ;
  else if (global_setting_perceive_sssr_rings)
  {
    if (_sssr_rings_perceived.number_elements() < _expected_nrings)
      _determine_uniqueness(_node_collision_rings);
  }
  else
    _determine_uniqueness(_node_collision_rings);

#ifdef DEBUG_FIND_UNIQUE_RINGS
  cerr << "After uniqueness determination, have " << _sssr_rings_perceived.number_elements() <<
           " rings, and " << _non_sssr_rings.number_elements() << " non sssr\n";

  for (int i = 0; i < _bonds_in_molecule; i++)
  {
    cerr << " Row " << setw(2) << i << " ";
    const Beep * b = _matrix_of_beeps[i];
    if (nullptr == b)
      cerr << "NULL";
    else
      b->printon(cerr);
    cerr << endl;
  }
#endif       

  return;
}

// Define this symbol to get all the private functions here

//#define DEBUG_TNODES

Tnode::Tnode(atom_number_t a, int acon, int pi, int bonds_in_molecule, int h) :
                  _a(a), _acon(acon), 
                  _bonds_in_molecule( bonds_in_molecule),
                  _pi_electrons(pi),
                  _atomic_number_score(h)
{
  assert (_acon > 0 && _acon <= _bonds_in_molecule);

  _con = new int[acon];
  _bond = new int[acon];

  _nsend = 0;

  _receive_buffer.resize(bonds_in_molecule);

  assert(ok());
}

Tnode::~Tnode()
{
  if (-9 == _acon)
    cerr << "Deleting an already deleted tnode\n";

  delete [] _con;
  delete [] _bond;

  _acon = -9;

  return;
}

int
Molecule::_initialise_tnode(atom_number_t zatom,
                            const int * process_these, int id,
                            Tnode ** tnodes,
                            const int * pi,
                            const int * sac)
{
  Tnode * tn = tnodes[zatom];

  assert(tn->acon() == _things[zatom]->ncon());

  int start_atom_atomic_number_score = sac[zatom];

#ifdef DEBUG_TNODES
  cerr << "Tnode " << zatom << " initialised, score " << start_atom_atomic_number_score << ", pi " << pi[zatom] << endl;
#endif

  int rc = 0;    // the number of bonds we process

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->involves(zatom))
      continue;

    atom_number_t j = b->other(zatom);

    if (id != process_these[j])
      continue;

    tn->is_connected_to(rc, j, i);

    Path_Message * tmp = new Path_Message(zatom, pi[zatom], start_atom_atomic_number_score, i, nb);
    tnodes[j]->receive(tmp);
    rc++;

#ifdef DEBUG_TNODES
    cerr << "  Bond to atom " << j << " is bond " << *_bond_list[rc] << endl;
    cerr << "Path message has " << tmp->pi_electrons() << " pi electrons\n";
#endif
  }

  return rc;
}

int
Tnode::ok() const
{
  if (_acon <= 0)
    return 0;

  if (_acon > _bonds_in_molecule)
    return 0;

  if (nullptr == _con)
    return 0;

  return 1;
}

int
Tnode::debug_print(std::ostream & os) const
{
  os << "Info on Tnode for atom " << _a << " acon = " << _acon << 
        ' ' << _pi_electrons << " pi electrons\n";
  os << "Completed " << _nsend << " send operations\n";

  if (! ok())
    os << "Warning, OK fails\n";

  if (_receive_buffer.number_elements())
    os << "Receive buffer contains " << _receive_buffer.number_elements() << 
          " path messages\n";

  return 1;
}

void 
Tnode::is_connected_to(int ndx,
                       atom_number_t a,
                       int bond_number)
{
  _con[ndx] = a;
  _bond[ndx] = bond_number;

  _acon = ndx + 1;

  return;
}

int
Tnode::print_paths(std::ostream & os) const
{
  assert(os.good());

  os << _receive_buffer.number_elements() << " paths currently at tnode for atom " << 
        _a << endl;

  int np = _receive_buffer.number_elements();
  if (0 == np)
    os << "  No paths in receive buffer\n";
  else
  {
    for (int i = 0; i < np; i++)
    {
      _receive_buffer[i]->debug_print(os);
    }
  }

  np = _send_buffer.number_elements();
  if (0 == np)
    os << "  No paths in outgoing buffer\n";
  else
  {
    for (int i = 0; i < np; i++)
    {
      _send_buffer[i]->debug_print(os);
    }
  }

  return os.good();
}

/*
  A send operation involves sending the contents of the send_buffer
  to neighbours.
*/

void
Tnode::send(Tnode ** tnodes)
{
#ifdef DEBUG_TNODES
  cerr << "Tnode for atom " << _a << " sending " << 
          _send_buffer.number_elements() << " path messages\n";
#endif

  int nr = _send_buffer.number_elements();

  if (0 == nr)
    return;

  Path_Message ** sb = _send_buffer.rawdata();    // for efficiency

#ifdef DEBUG_TNODES
  for (int i = 0; i < nr; i++)
  {
    const Path_Message * pi = sb[i];
  }
#endif

  for (int i = 0; i < nr; i++)
  {
    Path_Message * p = sb[i];
    for (int j = 0; j < _acon; j++)
    {
      atom_number_t k = _con[j];

      assert(nullptr != tnodes[k]);

      if (k != p->last_atom() && ! p->is_set(_bond[j]))
      {
        Path_Message * tmp = new Path_Message(p);
        tmp->include_in_path(_a, _bond[j], _pi_electrons, _atomic_number_score);
#ifdef DEBUG_TNODES
        cerr << "  Sending message " << i << " to " << k << " bond " << _bond[j] << endl;
#endif              
        tnodes[k]->receive(tmp);
      }
    }
  }

  _send_buffer.resize_keep_storage(0);    // Delete all the original Path_Message objects
#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced(stderr);
#endif

  _nsend++;

  return;
}

static void
remove_beeps_which_started_here(resizable_array_p<Path_Message> & path_messages,
                                atom_number_t a)
{
  int np = path_messages.number_elements();

  Path_Message ** pm = path_messages.rawdata();

  for (int i = 0; i < np; i++)
  {
    Path_Message * p = pm[i];

    if (a == p->start_atom())
    {
      path_messages.remove_item(i);
      i--;
      np--;
    }
  }

  return;
}

/*
  Check a group of Path_Messages for inverse edge collisions 

  Jan 2001. 
    Need to be a little careful about counting heteroatoms and pi
    electrons in an inverse edge collision. An inverse edge collision
    is when two Path_Message's that started at opposite ends of a common
    bond collide.

    So we don't count the properties of the two atoms at either end of
    this common bond, we subtract out the properties of the atoms that
    started the Path_Message's.
*/

void
Tnode::_check_inverse_edge_collisions(Rings_Found & rings_found)
{
  int nm = _receive_buffer.number_elements();

  Path_Message ** pm = _receive_buffer.rawdata();     // for efficiency

  for (int i = 0; i < nm; i++)
  {
    Path_Message * pi = pm[i];
    for (int j = i + 1; j < nm; j++)
    {
      Path_Message * pj = pm[j];
      if (pi->inverse_edge_collision(pj))
      {
        pi->set_to_be_removed(1);
        pj->set_to_be_removed(1);
  
#ifdef  DEBUG_TNODES
        cerr << "Found inverse edge collision, score " << _atomic_number_score << " between\n";
        pi->debug_print(cerr);
        cerr << " AND\n";
        pj->debug_print(cerr);
        cerr << endl;
#endif              

        Beep * b = new Beep;
        *b = *pi;
        b->iwor(*pj);

//      These Path_Messages both end at the same atom. It's pi count is not yet added. Don't double count the two atoms that are shared between the two Path_Messages

        b->set_pi_electrons(pi->pi_electrons() + pj->pi_electrons() + _pi_electrons - pi->pi_electrons_on_start_atom() - pj->pi_electrons_on_start_atom());
        b->set_atomic_number_score(pi->atomic_number_score() + pj->atomic_number_score() + _atomic_number_score - pi->start_atom_atomic_number_score() - pj->start_atom_atomic_number_score());
        assert(b->atomic_number_score() > 0);

#ifdef DEBUG_TNODES
        cerr << "Inverse edge collision: combining beeps with " << pi->pi_electrons() << " (" << pi->nset() << " bonds) and " << pj->pi_electrons() << " pi electrons\n";
        cerr << pj->nset() << endl;
        cerr << "I " << pi->atomic_number_score() << " first " << pi->start_atom_atomic_number_score() << endl;
        cerr << "J " << pj->atomic_number_score() << " first " << pj->start_atom_atomic_number_score() << endl;
        cerr << "Resulting beep contains " << b->pi_electrons() << " pi electrons, score " << b->atomic_number_score() << endl;
#endif

        rings_found.inverse_edge_collision_ring(b);
      }
    }
  }

  return;
}

/*
  Check a group of Path_Messages for node collisions. The Path_Messages have
  started at the same atom and are ending at the same atom
*/

void
Tnode::_check_node_collisions(Rings_Found & rings_found)
{
  int nm = _receive_buffer.number_elements();

  Path_Message ** pm = _receive_buffer.rawdata();

  for (int i = 0; i < nm; i++)
  {
    Path_Message * pi = pm[i];
    for (int j = i + 1; j < nm; j++)
    {
      Path_Message * pj = pm[j];
      if (pi->node_collision(pj))
      {
        pi->set_to_be_removed(1);
        pj->set_to_be_removed(1);

#ifdef  DEBUG_TNODES
        cerr << "Found node collision, score " << _atomic_number_score << " between\n";
        pi->debug_print(cerr);
        cerr << " AND\n";
        pj->debug_print(cerr);
        cerr << endl;
#endif              
  
        Beep * b = new Beep;
        *b = *pi;
        b->iwor(*pj);

//      These Path_Messages started on the same atom, so subtract out it's pi electrons to avoid double counting

        b->set_pi_electrons(pi->pi_electrons() + pj->pi_electrons() -
                            pi->pi_electrons_on_start_atom() + _pi_electrons);
        b->set_atomic_number_score(pi->atomic_number_score() + pj->atomic_number_score() + _atomic_number_score - pi->start_atom_atomic_number_score());
        assert(b->atomic_number_score() > 0);

#ifdef  DEBUG_TNODES
        cerr << "Node collision: combining beeps with " << pi->pi_electrons() << " (" << pi->nset() << " bonds) and " << pj->pi_electrons() << " pi electrons (" << pj->nset() << " bonds)\n";
        cerr << "I " << pi->atomic_number_score() << " first " << pi->start_atom_atomic_number_score() << endl;
        cerr << "J " << pj->atomic_number_score() << " first " << pj->start_atom_atomic_number_score() << endl;
        cerr << "Resulting beep contains " << b->pi_electrons() << " pi electrons, score " << b->atomic_number_score() << endl;
#endif

        rings_found.node_collision_ring(b);
      }
    }
  }

  return;
}

/*
  Check a group of Path_Messages for direct edge collisions
*/

void
Tnode::_check_direct_edge_collisions()
{
  int nm = _receive_buffer.number_elements();

  Path_Message ** pm = _receive_buffer.rawdata();

  for (int i = 0; i < nm; i++)
  {
    Path_Message * pi = pm[i];
    for (int j = i + 1; j < nm; j++)
    {
      Path_Message * pj = pm[j];
      if (pi->direct_edge_collision(pj))
      {
        pi->set_to_be_removed(1);
      }
    }
  }

  return;
}

/*
  A number of Path_Messages have accumulated in the _receive_buffer
  of a Tnode. See if there are any new rings. Newly perceived rings
  are handed to RINGS_FOUND for resolution.

  Initially, _receive_buffer has some number of Path_Messages. These
  need to be analysed, and the non-coliding messages placed in the
  _send_buffer.

*/

int
Tnode::check_for_rings(Rings_Found & rings_found)
{
  assert(_send_buffer.empty());

  remove_beeps_which_started_here(_receive_buffer, _a);

// Make sure inverse edge collisions are done first, as they
// produce rings of smaller size

  _check_inverse_edge_collisions(rings_found);
  _check_node_collisions        (rings_found);
  _check_direct_edge_collisions();

  int nr = _receive_buffer.number_elements();
  Path_Message ** rb = _receive_buffer.rawdata();

  for (int i = 0; i < nr; i++)
  {
    Path_Message * mi = rb[i];

    if (mi->to_be_removed())
      delete mi;
    else
      _send_buffer.add(mi);
  }

  _receive_buffer.resize_no_delete(0);

  return 1;
}

static int
discern_fused_neighbours(const resizable_array_p<Beep> & beeps,
                         const resizable_array_p<Ring> & rings,
                         int fid)
{
  int nr = beeps.number_elements();
  assert(nr == rings.number_elements());

  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    const Beep * bi = beeps[i];

    for (int j = i + 1; j < nr; j++)
    {
      const Beep * bj = beeps[j];
      int bic = bi->bits_in_common(*bj);

      if (bic)
      {
        Ring * rj = rings[j];
        ri->set_fused_to(rj, bic);
        rj->set_fused_to(ri, bic);
      }
    }
  }

// Now that fused neighbours are known, give them fused system identifiers

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    if (ri->is_fused())    // fused system identifier already assigned - done
      continue;

    if (0 == ri->fused_ring_neighbours())     // must be a spiro fusion
      continue;

    rc += ri->propagate_fused_system_identifier(fid);
    if (rc == nr)
      break;

    fid++;
  }

  return 1;
}

#ifdef BONDS_KNOW_RING_MEMBERSHIP

/*int
BondList::assign_ring_membership_to_bonds (const resizable_array_p<Beep> & beeps)
{
  int nb = beeps.number_elements();
  for (int i = 0; i < nb; i++)
  {
    const Beep * b = beeps[i];
    for (int j = 0; j < _number_elements; j++)
    {
      if (b->is_set(j))
        _things[j]->in_another_ring();
    }
  }

  return 0;
}*/


#endif

int
Molecule::_make_ring(const Beep * bp,
                     Ring * r,
                     int * atom_in_ring)
{
  assert (r->empty());

  const int * fragment_membership = _fragment_information.fragment_membership();

  int bonds_in_ring = 0;

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    if (bp->is_set(i))
    {
      Bond * b = _bond_list[i];

      bonds_in_ring++;

      b->in_another_ring();

      if (r->empty())
      {
        r->add(b->a1());
        r->add(b->a2());
        r->set_fragment_membership(fragment_membership[b->a1()]);
      }
      else
      {
        atom_in_ring[b->a1()] = 1;     // to be added to the ring later
        atom_in_ring[b->a2()] = 1;
      }
    }
  }

  if (r->elements_allocated() < bonds_in_ring)
    r->resize(bonds_in_ring);

// We need to build up the ring with the atoms in order

  bonds_in_ring--;
  atom_number_t first_atom = r->item(1);

  r->set_vector(atom_in_ring, 0);   // the first two atoms will have been set by other bonds

  while (bonds_in_ring--)
  {
    const Atom * a = _things[first_atom];
    int acon = a->ncon();
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(first_atom, j);
      if (atom_in_ring[k])
      {
        r->add(k);
        atom_in_ring[k] = 0;
        first_atom = k;
        break;
      }
    }
  }

  return 1;
}

//#define DEBUG_MAKE_RINGS

int
Molecule::_make_rings(const resizable_array_p<Beep> & beeps,
                      resizable_array_p<Ring> & rings,
                      int * tmp)
{
  int nr = beeps.number_elements();

#ifdef DEBUG_MAKE_RINGS
  cerr << "Molecule::_make_rings:processing " << nr << " beeps\n";
#endif

  if (rings.elements_allocated() < nr)
    rings.resize(nr);

  for (int i = 0; i < nr; i++)
  {
    Ring * r = new Ring;
    if (! _make_ring(beeps[i], r, tmp))
    {
      delete r;
      return 0;
    }

    rings.add(r);
  }

  return nr;
}

int
Molecule::_make_rings(Rings_Found & rings_found,
                      resizable_array_p<Ring> & sssr_rings,
                      resizable_array_p<Ring> & non_sssr_rings,
                      int fid,
                      int * tmp)
{
  if (! _make_rings(rings_found.sssr_beeps(), sssr_rings, tmp))
    return 0;

// ring membership is based on sssr rings only

  int nr = sssr_rings.number_elements();

#ifdef DEBUG_MAKE_RINGS
  cerr << "Molecule::_make_rings: found " << nr << " sssr rings\n";
#endif

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = sssr_rings[i];

#ifdef DEBUG_MAKE_RINGS
    cerr << "Ring " << i << " has " << ri->number_elements() << " members\n";
#endif

    const int rs = ri->number_elements();
    for (int j = 0; j < rs; j++)
    {
      atom_number_t k = ri->item(j);

      _ring_membership[k]++;
//    _things[k]->in_another_ring();
    }
  }

  if (rings_found.non_sssr_beeps().number_elements())
  {
    if (! _make_rings(rings_found.non_sssr_beeps(), non_sssr_rings, tmp))
      return 0;
  }

  discern_fused_neighbours(rings_found.sssr_beeps(), sssr_rings, fid);

  return 1;
}

int
Molecule::_unused_fused_system_identifier() const
{
  int ns = _sssr_rings.number_elements();

  int rc = -99;
  for (int i = 0; i < ns; i++)
  {
    const Ring * ri = _sssr_rings[i];
    //if (! ri->is_fused())
    //  continue;

    if (ri->fused_system_identifier() > rc)
      rc = ri->fused_system_identifier();
  }

  if (rc < 0)
    return 0;

  return rc + 1;     // one beyond the largest one found
}

/*
  Some rings have been found, create the molecule's rings from these
  raw forms.

  First step is to find a fused system identifier
*/

int
Molecule::_make_rings(Rings_Found & rings_found,
                      resizable_array_p<Ring> & sssr_rings,
                      resizable_array_p<Ring> & non_sssr_rings)
{
  int fid = _unused_fused_system_identifier();

  int * tmp = new_int(_number_elements); std::unique_ptr<int[]> free_tmp(tmp);

  return _make_rings(rings_found, sssr_rings, non_sssr_rings, fid, tmp);
}

/*
  Unfortunately the path lenght comparitor function is based on Path
  objects rather than Ring objects. Need to cast it
*/

#define RING_SORT_FN ( int (*) (Ring * const *, Ring * const *))

int
Molecule::_pearlman_sssr(const int * process_these, int id,
                         Tnode ** tnodes, int expected_nrings)
{
// Although not necessary, we keep an iteration counter to guard against infinite loops

//cerr << "Upon entry to _pearlman_sssr\n";
//_bond_list.debug_print(cerr);

  int iterations = 0;

  Rings_Found rings_found(expected_nrings, _bond_list.number_elements());
  if (expected_nrings > 2)
    rings_found.initialise_single_bond_count(*this);

  while (rings_found.rings_found() < expected_nrings)
  {
#ifdef DEBUG_TNODES
    cerr << "Iteration " << iterations << " found " <<  rings_found.rings_found() << " expected " << expected_nrings << " SSSR " << _sssr_rings.number_elements() << " NON SSSR " << _non_sssr_rings.number_elements() << endl;
    cerr << "Tnode send iteration " << iterations << endl;
#endif

    if (iterations)      // first send is done by initialise
    {
      for (int i = 0; i < _number_elements; i++)
      {
        if (id == process_these[i])
          tnodes[i]->send(tnodes);
      }
    }
    
#ifdef DEBUG_TNODES
    cerr << "After send process\n";
    for (int i = 0; i < _number_elements; i++)
    {
      if (id == process_these[i])
        tnodes[i]->print_paths(cerr);
    }
    cerr << "Now check for rings\n";
#endif

    for (int i = 0; i < _number_elements; i++)
    {
      if (id == process_these[i])
        tnodes[i]->check_for_rings(rings_found);
    }

    rings_found.process_new_rings(*this);

    iterations++;
    if (iterations < _number_elements)   // keep looking
      continue;

    if (! perceive_sssr_rings())   // transfer last found rings to sssr set
      break;

    cerr << "Molecule::_pearlman_sssr: too many iterations, beware incomplete SSSR ring determination\n";
    cerr << "Expected " << expected_nrings << " rings, found " << rings_found.rings_found() << ", likely a very complex fused system\n";
    break;
  }

  const int rings_found_here = rings_found.rings_found();

// Mark all these atoms as being in no rings - will increment later

  for (int i = 0; i < _number_elements; i++)
  {
    if (id == process_these[i] && _ring_membership[i] < 0)
      _ring_membership[i] = 0;
  }
 
// We need to give make_rings a fused ring system identifier it can use

  const int ns = _sssr_rings.number_elements();
  int fid = -99;
  for (int i = 0; i < ns; i++)
  {
    const Ring * ri = _sssr_rings[i];
    if (ri->is_fused() && ri->fused_system_identifier() > fid)
      fid = ri->fused_system_identifier();
  }

  if (fid < 0)
    fid = 0;
  else
    fid++;

  resizable_array_p<Ring> sssr_rings;
  resizable_array_p<Ring> non_sssr_rings;

  if (! _make_rings(rings_found, sssr_rings, non_sssr_rings))
    return 0;

// If there was an incomplete determination, try to transfer in some random non sssr rings
// Right now, we take the first item in the array, likely a small ring

  while (rings_found_here < expected_nrings && non_sssr_rings.number_elements())
  {
    _sssr_rings.transfer_in(non_sssr_rings, 0);
//  _sssr_rings.add(non_sssr_rings.pop());
  }

  sssr_rings.sort(RING_SORT_FN path_length_comparitor_longer);

#ifdef DEBUG_TNODES
  int nr = sssr_rings.number_elements();
  cerr << "Found " << nr << " rings\n";
  for (int i = 0; i < nr;i++)
  {
    const Ring * r = sssr_rings[i];
    cerr << "Ring " << i << " is " << *r << endl;
  }
#endif

  _sssr_rings.transfer_in(sssr_rings);
  _sssr_rings.sort(RING_SORT_FN path_length_comparitor_longer);

  int nsr = _sssr_rings.number_elements();
  for (int i = 0; i < nsr; i++)
  {
    _sssr_rings[i]->set_ring_number(i);
  }

  if (file_scope_accumulate_non_sssr_rings && non_sssr_rings.number_elements())
  {
    _non_sssr_rings.transfer_in(non_sssr_rings);
    _non_sssr_rings.sort(RING_SORT_FN path_length_comparitor_longer);
  }

  return 1;
}

int
Molecule::_pearlman_sssr(const int * process_these, int id,
                         Tnode ** tnodes,
                         const int * pi,
                         const int * sac)
{
  int bonds_in_molecule = _bond_list.number_elements();

  int atoms_being_processed = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (id != process_these[i])
      continue;

    const Atom * a = _things[i];

    tnodes[i] = new Tnode(i, a->ncon(), pi[i], bonds_in_molecule, sac[i]);
    atoms_being_processed++;

//  cerr << "Atom " << i << " has " << pi[i] << " pi electrons\n";
  }

  if (atoms_being_processed < 3)
  {
    cerr << "Molecule::_pearlman_sssr: only " << atoms_being_processed << " atoms\n";
    return 1;
  }

  int bonds_being_processed = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (id == process_these[i])
      bonds_being_processed += _initialise_tnode(i, process_these, id, tnodes, pi, sac);
  }

  assert(0 == bonds_being_processed % 2);
  bonds_being_processed = bonds_being_processed >> 1;    // divide by 2

  int expected_nrings = bonds_being_processed - atoms_being_processed + 1;

  assert(expected_nrings >= 0);
  if (0 == expected_nrings)
    return 1;

//cerr << "bonds = " << bonds_being_processed << " atoms = " << atoms_being_processed << endl;
//cerr << "Expecting to find " << expected_nrings << " rings\n";

// Huge problems with the large ring system. This just fails. Therefore
// pre-emtively turn on accumulation of non_sssr_rings if the system looks too complex.
// This is a massive kludge, but I'm not even sure what the right set of rings for that molecule
// would even be - somehow one of the redundant 3 membered rings would need to be kept presumably...

  int restore_file_scope_accumulate_non_sssr_rings = file_scope_accumulate_non_sssr_rings;

  if (expected_nrings > 4 && ! file_scope_accumulate_non_sssr_rings)
  {
    file_scope_accumulate_non_sssr_rings = 1;
    cerr << "Temporarily perceiving non sssr rings\n";
  }

  int rc = _pearlman_sssr(process_these, id, tnodes, expected_nrings);

  file_scope_accumulate_non_sssr_rings = restore_file_scope_accumulate_non_sssr_rings;

  return rc;
}

/*
  Complexity. Consider
FC1=CN2C(O)N(C2=O)C1=O PBCHM68787938
  there are two possible 6 membered rings. Which to choose? If we exclude pi electrons 
  associated with exocyclic double bonds, these two rings are identical. 
  
  If bonded to an exocyclic double bond then it has zero pi electrons to share
*/

#ifdef NEWER_VERSION_IS_CLEANER
int
Molecule::_pi_electrons_in_ring(const atom_number_t zatom, int & result) const
{
  Atom * a = const_cast<Atom *>(_things[zatom]);
//return a->pi_electrons(result);     can turn this on or off. Works slightly better turned off

  const int acon = a->ncon();
  if (2 == acon)
    return a->pi_electrons(result);

  if (acon == a->nbonds())
    return a->pi_electrons(result);

  if (6 != a->atomic_number())
    return a->pi_electrons(result);

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);
    if (b->is_single_bond()) 
      continue;

    const atom_number_t j = b->other(zatom);

    const Atom * aj = _things[j];
    if (6 != aj->atomic_number() && 1 == aj->ncon())
    {
      result = 0;
      return 1;
    }
  }

  return a->pi_electrons(result);
}
#endif

int
Molecule::_pi_electrons_in_ring(const atom_number_t zatom, int & result) const
{
  Atom * a = const_cast<Atom *>(_things[zatom]);
//return a->pi_electrons(result);     can turn this on or off. Works slightly better turned off

  if (! a->pi_electrons(result)) {
    return 0;
  }
  if (result > 2) {
    result = 2;
  }

  // An obviously exocyclic double bond to a heteroatom contributes zero.
  for (const Bond * b : *a) {
    if (b->is_single_bond()) {
      continue;
    }
    const Atom * aj = _things[b->other(zatom)];
    if (aj->ncon() == 1 && aj->atomic_number() != 6) {
      result = 0;
      return 1;
    }
  }

  return 1;
}


int
Molecule::_pearlman_sssr(const int * process_these, int id,
                         Tnode ** tnodes)
{
  int * tmp = new int[_number_elements + _number_elements]; std::unique_ptr<int[]> free_tmp(tmp);
  int * sac = tmp + _number_elements;
  std::fill_n(sac, _number_elements, 0);

  for (int i = 0; i < _number_elements; i++)
  {
    if (nullptr != _aromaticity && AROMATIC == _aromaticity[i])   // oct 05. No this is wrong... 2017
      tmp[i] = 2;    // not really pi count, just to indicate favoured
    else if (! _pi_electrons_in_ring(i, tmp[i]))    // could not be computed
      ;
//  else if (tmp[i] > 2)
//    tmp[i] = 2;

//  cerr << "atom " << i << " " << _things[i]->atomic_symbol() << " ncon " << _things[i]->ncon() << " nbonds " << _things[i]->nbonds() << " pi " << tmp[i] << endl;

    Atom * a = const_cast<Atom *>(_things[i]);
    const int acon = a->ncon();

    sac[i] += 900 * a->element()->unique_id() + 5 * acon + a->implicit_hydrogens();
    assert (sac[i] > 0);

    if (acon < 4)
      ;
    else if (sac[i] > 40)   // disfavour highly connected things in rings
      sac[i] -= 40;
    else
      sac[i] = sac[i] / 2;

    continue;
//  const atomic_number_t z = a->atomic_number();

//  if (6 == z || 1 == z)
//    continue;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      sac[k]++;
    }

    assert(sac[i] > 0);
  }

  return _pearlman_sssr(process_these, id, tnodes, tmp, sac);
}

int
Molecule::_pearlman_sssr(const int * process_these, int id)
{
//cerr << "Molecule::_pearlman_sssr\n";
//debug_print(cerr);

  if (_number_elements <= 2)
    return 1;

  assert(nullptr != _ring_membership);

//debug_print(cerr);

  Tnode ** tnodes = new Tnode * [_number_elements]; std::unique_ptr<Tnode *[]> free_tnodes(tnodes);

  int rc = _pearlman_sssr(process_these, id, tnodes);

  for (int i = 0; i < _number_elements; i++)
  {
    if (id == process_these[i])
      delete tnodes[i];
  }

  return rc;
}
