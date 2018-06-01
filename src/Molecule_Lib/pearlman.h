#ifndef PEARLMAN_SSSR_H
#define PEARLMAN_SSSR_H

#include "iwmtypes.h"

#include "iwbits.h"

/*
  We simultaneously discover members of the SSSR, and we need some means
  of consistently preferentially choosing one ring over another. The
  variable _atomic_number_score is used for that purpose.
*/

class Beep: public IW_Bits_Base
{
  protected:
    int _pi_electrons;       // for deciding between rings in SSSR decisions. We favour possibly aromatic rings

//  int _heteroatoms;

    int _atomic_number_score;

//  We can break ties with this measure N1=C2C3=C(C2=CC=C3)C=C1 PBCHM21760248

    int _single_bond_count;

//  private functions

    void _default_values ();

  public:
    Beep ();
    Beep (int);

    int ok () const;
    int debug_print (std::ostream &) const;

    void set_pi_electrons (int p) { _pi_electrons = p;}
    int pi_electrons () const { return _pi_electrons;}

//  int heteroatoms_encountered () const { return _heteroatoms;}
//  void set_heteroatoms_encountered (int h) { _heteroatoms = h;}

    int atomic_number_score () const { return _atomic_number_score;}
    void set_atomic_number_score (int h) { _atomic_number_score = h;}

    int  single_bond_count () const { return _single_bond_count;}
    void set_single_bond_count (int s) { _single_bond_count = s;}
};

typedef int edge_type_t;

class Rings_Found;

class Path_Message : public Beep
{
  private:
    const atom_number_t _start_atom;
    const int           _pi_electrons_on_start_atom;
    const int           _start_atom_atomic_number_score;
    const int           _start_edge;
    const int           _bonds_in_molecule;

    atom_number_t _last_atom;

//  When a Tnode examines its Path_Messages, it can mark each one as
//  needing to be removed.

    int _to_be_removed;

//  private functions

    void _default_values ();

  public:
    Path_Message (atom_number_t, int, int, int, int);
    Path_Message (const Path_Message *);

    int ok () const;
    int debug_print (std::ostream &);

    atom_number_t start_atom () const { return _start_atom;}
    atom_number_t last_atom () const { return _last_atom;}
    void set_last_atom (atom_number_t a) { _last_atom = a;}

    edge_type_t   start_edge () const { return _start_edge;}

    void set_to_be_removed (int i) { _to_be_removed = i;}
    int to_be_removed () const { return _to_be_removed;}

    void include_in_path (atom_number_t, int, int, int);

    int no_collision (const Path_Message * p) const { return
         _start_atom != p->_start_atom && _start_edge != p->_start_edge;}

    int node_collision (const Path_Message * p) const;

    int inverse_edge_collision (const Path_Message * p) const;

    int direct_edge_collision (const Path_Message * p) const { return
      _start_atom == p->_start_atom && _start_edge == p->_start_edge;}

    int pi_electrons_on_start_atom () const { return _pi_electrons_on_start_atom;}
    int start_atom_atomic_number_score () const { return _start_atom_atomic_number_score;}
};

/*
  As tnodes detect collisions in their in-buffer's, we need a global
  entity which accepts reports of new rings.
*/

class Molecule;
class Ring;

class Rings_Found
{
  private:
    const int _expected_nrings;

    resizable_array_p<Beep> _inverse_edge_collision_rings;
    resizable_array_p<Beep> _node_collision_rings;
    resizable_array_p<Beep> _sssr_rings_perceived;

//  We also keep track of those rings which are small rings, but not part
//  of the SSSR set - the extra ring in cubane for example

    resizable_array_p<Beep> _non_sssr_rings;

//  The matrix is used for determining linear independence

    Beep **                 _matrix_of_beeps;

    const int _bonds_in_molecule;

//  we want to know which bonds are single bonds. Can save some time
//  by precomputing that info here

    IW_Bits_Base _is_single_bond;

//  Beep _bonds_covered;

//  private functions

    int  _beep_is_unique_over_non_sssr_beeps (const Beep *) const;
    int  _beep_is_unique (const Beep *) const;
    int  _is_sssr_ring   (const Beep *);
    int  _is_esssr_ring  (const Beep * b);

    void _assign_single_bond_counts (Beep * beep);
    void _assign_single_bond_counts (resizable_array_p<Beep> & beeps);

    void _determine_uniqueness (resizable_array_p<Beep> &);

  public:
    Rings_Found (int, int);
    ~Rings_Found ();

    void initialise_single_bond_count (const Molecule &);

    int rings_found () const { return _sssr_rings_perceived.number_elements ();}
    int all_rings_found () const { return _expected_nrings == _sssr_rings_perceived.number_elements ();}

    void node_collision_ring (Beep * b);
    void inverse_edge_collision_ring (Beep * b);

    void process_new_rings (const Molecule &);

    const resizable_array_p<Beep> & sssr_beeps () const { return _sssr_rings_perceived;}
    const resizable_array_p<Beep> & non_sssr_beeps () const { return _non_sssr_rings;}

};

class Tnode
{
  private:
    const atom_number_t _a;
    int                 _acon;
    const int           _bonds_in_molecule;
    resizable_array_p<Path_Message> _receive_buffer;
    resizable_array_p<Path_Message> _send_buffer;
    atom_number_t       * _con;
    int                 * _bond;
    int                   _nsend;
    int                   _pi_electrons;
    int                   _atomic_number_score;

//  private functions

    void _check_inverse_edge_collisions (Rings_Found & rings_found);
    void _check_node_collisions (Rings_Found & rings_found);
    void _check_direct_edge_collisions ();

  public:
    Tnode (atom_number_t, int, int, int, int);
    ~Tnode ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int acon () const { return _acon;}

    void is_connected_to (int, atom_number_t, int);

    int print_paths (std::ostream &) const;

    void receive (Path_Message * p) { _receive_buffer.add (p);}

    void send (Tnode **);

    int check_for_rings (Rings_Found &);
};


/*
  We can perceive either SSSR rings, or ESSSR rings
*/

extern int  perceive_sssr_rings();
extern void set_perceive_sssr_rings(int s);

#endif
