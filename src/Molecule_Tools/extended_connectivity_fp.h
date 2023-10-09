#ifndef IW_EC_GENERATOR_H
#define IW_EC_GENERATOR_H

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"

/*
  To avoid passing around a lot of arguments, we just put them all in a struct
*/

class Gather_Single_Bit;

template <typename T>
class ECFP_Args
{
  private:
    const int * _include_atom;
    int   _flag;

    const T * _atype;

    int * _processing_status;

  public:
    ECFP_Args(const int natoms);
    ~ECFP_Args();

    void set_include_atom(const int * s) { _include_atom = s;}
    void set_flag(const int s) { _flag = s;}
    void set_atype(const T * s) { _atype = s;}

    int include_atom(const atom_number_t a) const { return _flag == _include_atom[a];}

    const T * atype() const { return _atype;}
    int * processing_status() const { return _processing_status;}

    T atype(const int s) const { return _atype[s];}
    int processing_status(const int s) const { return _processing_status[s];}
};

class EC_Fingerprint_Generator
{
  private:
    int _min_radius;
    int _max_radius;
    int _precise_fingerprints;
    int _bit_increment;

    int _fingerprint_ring_closures; 
    int _add_tails;
    int _differentiate_ring_and_chain_bonds;

//  private functions

    template <typename T> int _prepare_first_shell (const Molecule & m, const atom_number_t zatom, Set_of_Atoms & first_shell, const ECFP_Args<T> & ecfp_args) const;
    int _bond_constant (const Bond * b) const;
    template <typename S, typename T>
    void _expand_shell (const Molecule & m, const Set_of_Atoms & expand_from, const int r, const ECFP_Args<T> & ecfp_args,
                                         unsigned int sum_so_far, const int each_shell_gets_different_fingerprint, S * sfc) const;
    template <typename T> void _generate_fingerprint (Molecule & m, const atom_number_t zatom, const T * atype,
                              const int * include_atom, const int flag, const int each_shell_gets_different_fingerprint, Gather_Single_Bit * bits) const;
    template <typename T> int _generate_fingerprint (Molecule & m,
                              const int each_shell_gets_different_fingerprint,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              Sparse_Fingerprint_Creator * sfc) const;

    template <typename S, typename T>
    void _do_add_tails (const Molecule & m, const Set_of_Atoms & expand_from, const int radius,
                                         const ECFP_Args<T> & ecfp_args, unsigned int sum_so_far, const int each_shell_gets_different_fingerprint,
                                         S * sfc) const;
  public:
    EC_Fingerprint_Generator ();


    void set_min_radius (int s) { _min_radius = s;}
    void set_max_radius (int s) { _max_radius = s;}
    void set_precise_fingerprints(int s) { _precise_fingerprints = s;}
    void set_bit_increment(const int s) { _bit_increment = s;}
    void set_fingerprint_ring_closures(const int s) { _fingerprint_ring_closures = s;}
    void set_differentiate_ring_and_chain_bonds(const int s) { _differentiate_ring_and_chain_bonds = s;}

    int max_radius () const { return _max_radius;}


#define ECFP_IMPL_INT ;
#define ECFP_IMPL_VOID ;


    template <typename T> int generate_fingerprint (Molecule & m, const T * atype,
                              const int * include_atom, const int flag, Sparse_Fingerprint_Creator & sfc) const ECFP_IMPL_INT

    template <typename T> int generate_fingerprint (Molecule & m, const T * atype,
                              const int * include_atom, const int flag, Sparse_Fingerprint_Creator * sfc) const ECFP_IMPL_INT   // each shell generated separately

//  get bit for a given radius, starting at a given atom

    template <typename T> void generate_fingerprint (Molecule & m, const atom_number_t zatom, const T * atype,
                              const int * include_atom, const int flag, resizable_array<unsigned int> & bits) const ECFP_IMPL_VOID
};

#undef ECFP_IMPL_VOID
#undef ECFP_IMPL_INT

#endif
