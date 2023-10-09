#include <stdint.h>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "extended_connectivity_fp.h"

using std::cerr;
using std::endl;

class Gather_Single_Bit
{
  private:
    resizable_array<unsigned int> _b;

  public:
    Gather_Single_Bit ();

    void hit_bit (const unsigned int s) { _b.add(s);}
    void hit_bit (const unsigned int s, const int notused) { _b.add(s);}

    int nbits() const {return _b.number_elements();}

    const resizable_array<unsigned int> & get () const { return _b;}
};

Gather_Single_Bit::Gather_Single_Bit ()
{
  _b = 0;

  return;
}

EC_Fingerprint_Generator::EC_Fingerprint_Generator()
{
  _min_radius = 0;
  _max_radius = 3;
  _precise_fingerprints = 0;
  _bit_increment = 1;

  _fingerprint_ring_closures = 0;   // hurts svmfp performance slightly

  _add_tails = 0;     // very slightly helps svmfp performance

  _differentiate_ring_and_chain_bonds = 0;

  return;
}

#define EC_FINISHED -9
#define EC_NEXT_TIME -1
#define EC_CURRENT_SHELL -2

template <typename T>
int
EC_Fingerprint_Generator::generate_fingerprint(Molecule & m,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              Sparse_Fingerprint_Creator & sfc) const
{
  return _generate_fingerprint(m, 0, atype, include_atom, flag, &sfc);
}

template <typename T>
int
EC_Fingerprint_Generator::generate_fingerprint(Molecule & m,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              Sparse_Fingerprint_Creator * sfc) const
{
  return _generate_fingerprint(m, 1, atype, include_atom, flag, sfc);
}

//#define DEBUG_EXTENDED_CONNECTIVITY
#ifdef DEBUG_EXTENDED_CONNECTIVITY
static void
write_labelled_by_inclusion(const Molecule & m,
                            const int * include_atom,
                            const int flag,
                            const char * s,
                            std::ostream & output)
{
  Molecule mcopy(m);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (flag == include_atom[i])
    {
      mcopy.set_isotope(i, i);
      if (0 == i)
        mcopy.set_isotope(0, 1);
    }
  }

  output << mcopy.smiles() << ' ' << s << endl;

  return;
}
#endif

template <typename T>
int
EC_Fingerprint_Generator::_generate_fingerprint(Molecule & m,
                              const int each_shell_gets_different_fingerprint,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              Sparse_Fingerprint_Creator * sfc) const
{
  assert (0 == each_shell_gets_different_fingerprint || 1 == each_shell_gets_different_fingerprint);

  m.compute_aromaticity_if_needed();    // because bond typing includes aromaticity

  const int matoms = m.natoms();

  if (0 == matoms)
    return 1;

#ifdef DEBUG_EXTENDED_CONNECTIVITY
  write_labelled_by_inclusion(m, include_atom, flag, "for fingerprints SFC", cerr);
#endif

  ECFP_Args<T> ecfp_args(matoms);
  ecfp_args.set_include_atom(include_atom);
  ecfp_args.set_flag(flag);
  ecfp_args.set_atype(atype);
  int * processing_status = ecfp_args.processing_status();

  Set_of_Atoms first_shell;
  first_shell.resize(4);

  for (int i = 0; i < matoms; ++i)
  {
    if (flag != include_atom[i])
      continue;

    if (0 == _min_radius)
    {
      sfc[0].hit_bit(atype[i], _bit_increment);

      if (0 == _max_radius)
        continue;
    }

    for (int j = 0; j < matoms; ++j)
    {
      if (flag == include_atom[j])
        processing_status[j] = 0;
      else
        processing_status[j] = EC_FINISHED;
    }

    first_shell.resize_keep_storage(0);

    const auto b = 7879 * atype[i] + _prepare_first_shell(m, i, first_shell, ecfp_args);

//  cerr << "From atom " << m.smarts_equivalent_for_atom(i) << " type " << atype[i] << " b " << b << endl;
    if (_min_radius <= 1)
    {
      sfc[each_shell_gets_different_fingerprint].hit_bit(static_cast<unsigned int>(b), _bit_increment);
      if (1 == _max_radius)
        continue;
    }

#ifdef DEBUG_EXTENDED_CONNECTIVITY
    cerr << "shell expansion beginning with atom " << i << " type " << m.smarts_equivalent_for_atom(i) << endl;
#endif

    ecfp_args.processing_status()[i] = EC_FINISHED;

    _expand_shell(m, first_shell, 2, ecfp_args, b, each_shell_gets_different_fingerprint, sfc + 2 * each_shell_gets_different_fingerprint);   // depends on each_shell_gets_different_fingerprint being 1
  }

  return 1;
}

int
EC_Fingerprint_Generator::_bond_constant(const Bond * b) const
{
  if (b->is_aromatic())
    return 11;

  int rc;
  if (b->is_single_bond())
    rc = 3;
  else if (b->is_double_bond())
    rc = 5;
  else
    rc = 7;

  if (_differentiate_ring_and_chain_bonds && b->nrings())
    rc *= 17;

  return rc;
}

template <typename T>
int
EC_Fingerprint_Generator::_prepare_first_shell(const Molecule & m,
                                   const atom_number_t zatom,
                                   Set_of_Atoms & first_shell,
                                   const ECFP_Args<T> & ecfp_args) const
{
//unsigned int rc = ecfp_args.atype(zatom);

  unsigned int rc = 0;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  first_shell.resize_keep_storage(0);

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const auto j = b->other(zatom);

    if (! ecfp_args.include_atom(j))
      continue;

    const auto bcst = _bond_constant(b);
//  cerr << "   BC " << bcst << " to type " << ecfp_args.atype(j) << endl;

    if (_precise_fingerprints >= 1)
      rc += 15 * ecfp_args.atype(zatom) + bcst * ecfp_args.atype(j);
    else
      rc += bcst * ecfp_args.atype(j);

    ecfp_args.processing_status()[j] = EC_NEXT_TIME;

    first_shell.add(j);
  }

#ifdef DEBUG_EXTENDED_CONNECTIVITY
  cerr << "_prepare_first_shell: " << m.smarts_equivalent_for_atom(zatom) << " returning " << rc << endl;
#endif

  return rc;
}

//#define DEBUG_EXPAND_SHELL

template <typename S, typename T>
void
EC_Fingerprint_Generator::_expand_shell(const Molecule & m,
                                        const Set_of_Atoms & expand_from,
                                        const int radius,
                                        const ECFP_Args<T> & ecfp_args,
                                        unsigned int sum_so_far,
                                        const int each_shell_gets_different_fingerprint,
                                        S * sfc) const
{
  if (_add_tails)
    _do_add_tails(m, expand_from, radius, ecfp_args, sum_so_far, each_shell_gets_different_fingerprint, sfc);

#ifdef DEBUG_EXPAND_SHELL
  cerr << "Begin radius " << radius << " sum_so_far " << sum_so_far << " atoms " << expand_from << endl;
  Molecule mcopy(m);
  write_isotopically_labelled_smiles(mcopy, false, cerr);
#endif

//sum_so_far = radius * sum_so_far + s;
  sum_so_far = sum_so_far * 7879;

  unsigned int s = 0;

  Set_of_Atoms next_shell;
  next_shell.resize(10);
  int rings_encountered = 0;   // will only get even sized rings. Implement something sometime maybe...

  expand_from.set_vector(ecfp_args.processing_status(), EC_CURRENT_SHELL);

  const int n = expand_from.number_elements();

  for (auto i = 0; i < n; ++i)
  {
    const atom_number_t j = expand_from[i];

#ifdef DEBUG_EXPAND_SHELL
    for (int k = 0; k < radius; ++k) { cerr << ' '; } cerr << radius << " atom " << j << " type " << ecfp_args.atype(j) << endl;
#endif

    const Atom * a = m.atomi(j);

    const auto acon = a->ncon();

    for (auto k = 0; k < acon; ++k)
    {
      const Bond * b = a->item(k);

      const auto l = b->other(j);

//    cerr << " from " << j << " to " << l << ' ' << m.smarts_equivalent_for_atom(l) << " status " << ecfp_args.processing_status(l) << endl;
      if (EC_FINISHED == ecfp_args.processing_status(l))
        continue;

      const int l_in_expand = (EC_CURRENT_SHELL == ecfp_args.processing_status(l));

      if (l_in_expand)        // we have found an odd sized ring
      {
        if (! _fingerprint_ring_closures)
          continue;

        if (l < j)      // if we are going to fingerprint this connection, only do it once
          continue;
      }

      if (l_in_expand)    // keep it out of the next shell
        ;
      else if (0 == next_shell.add_if_not_already_present(l))
        rings_encountered++;

#ifdef DEBUG_EXPAND_SHELL
      cerr << " _expand_shell went from atom " << j << " to atom " << l << ' ' << m.smarts_equivalent_for_atom(l) << endl;
#endif

      if (_precise_fingerprints)
        s += ecfp_args.atype(k) * (_bond_constant(b) + 17 * ecfp_args.atype(j));
      else
        s += _bond_constant(b) * ecfp_args.atype(l); 
    }
  }

#ifdef DEBUG_EXPAND_SHELL
  for (int j = 0; j < radius; ++j) { cerr << ' '; }
  cerr << radius << " sum " << s << " max radius " << _max_radius << endl;
#endif

  if (rings_encountered && _fingerprint_ring_closures)
    sfc[0].hit_bit(sum_so_far + 19 * radius + rings_encountered);

  sum_so_far += s;

  if (radius == _max_radius || 0 == s)
  {
    if (s > 0)
      sfc[0].hit_bit(sum_so_far, _bit_increment);
#ifdef DEBUG_EXPAND_SHELL
    cerr << "Hit bit " << sum_so_far << " at radius " << radius << " (s = " << s << ") aggregator has " << sfc[0].nbits() << endl;
#endif
    return;
  }
  else if (radius < _min_radius)  // do not set here
    ;
  else
    sfc[0].hit_bit(sum_so_far, _bit_increment);

  expand_from.set_vector(ecfp_args.processing_status(), EC_FINISHED);

  _expand_shell(m, next_shell, radius+1, ecfp_args, sum_so_far, each_shell_gets_different_fingerprint, sfc + each_shell_gets_different_fingerprint);

  return;
}

template <typename T>
void
EC_Fingerprint_Generator::generate_fingerprint(Molecule & m,
                              const atom_number_t zatom,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              resizable_array<unsigned int> & bits) const
{
#ifdef DEBUG_GENERATE_FINGERPRINT_ABITS
  cerr << "EC_Fingerprint_Generator::generate_fingerprint:molecule has " << m.natoms() << " atoms\n";
  for (int i = 0; i < m.natoms(); ++i)
  {
    cerr << " atom " << i << " type " << m.atomic_symbol(i) << " include " << (flag == include_atom[i]) << endl;
  }
#endif

  if (0 == _max_radius)
  {
    bits.add(atype[zatom]);
    return;
  }

  Gather_Single_Bit gsb;

  _generate_fingerprint(m, zatom, atype, include_atom, flag, 0, &gsb);

  bits = gsb.get();

#ifdef DEBUG_GENERATE_FINGERPRINT_ABITS
  cerr  << "Bits array contains " << bits.size() << " items\n";
  for (int i = 0; i < bits.number_elements(); ++i)
  {
    cerr << " " << bits[i] << endl;
  }
#endif

  return;
}

template <typename T>
void
EC_Fingerprint_Generator::_generate_fingerprint(Molecule & m,
                              const atom_number_t zatom,
                              const T * atype,
                              const int * include_atom,
                              const int flag,
                              const int each_shell_gets_different_fingerprint,
                              Gather_Single_Bit * gsb) const
{
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

#ifdef DEBUG_EXTENDED_CONNECTIVITY
  write_labelled_by_inclusion(m, include_atom, flag, "for fingerprints GSB", cerr);
  cerr << "Default radius 0 bit " << atype[zatom] << endl;
#endif

  gsb[0].hit_bit(atype[zatom]);

  if (0 == _max_radius)
    return;

  ECFP_Args<T> ecfp_args(matoms);
  ecfp_args.set_include_atom(include_atom);
  ecfp_args.set_flag(flag);
  ecfp_args.set_atype(atype);

  int * processing_status = ecfp_args.processing_status();
  for (int i = 0; i < matoms; ++i)
  {
    if (flag != include_atom[i])
      processing_status[i] = EC_FINISHED;
    else
      processing_status[i] = 0;
  }

  Set_of_Atoms first_shell;

  auto b = _prepare_first_shell(m, zatom, first_shell, ecfp_args);

//cerr << "Completed first shell, _max_radius " << _max_radius << endl;
  if (1 == _max_radius)
  {
    if (b > 0)
      gsb[each_shell_gets_different_fingerprint].hit_bit(b);
    return;
  }
  else if (_min_radius > 1)
    ;
  else
  {
    b = 7879 * atype[zatom] + b;
    gsb[each_shell_gets_different_fingerprint].hit_bit(b);
  }

#ifdef DEBUG_EXTENDED_CONNECTIVITY
    cerr << "shell expansion beginning with atom " << zatom << " type " << m.smarts_equivalent_for_atom(zatom) << " sum " << b << endl;
#endif

  ecfp_args.processing_status()[zatom] = EC_FINISHED;

  _expand_shell(m, first_shell, 2, ecfp_args, b, each_shell_gets_different_fingerprint, gsb + 2 * each_shell_gets_different_fingerprint);

  return;
}

template <typename T>
ECFP_Args<T>::ECFP_Args(const int matoms)
{
  _processing_status = new int[matoms];

  return;
}

template <typename T>
ECFP_Args<T>::~ECFP_Args()
{
  delete [] _processing_status;
}

template <typename S, typename T>
void
EC_Fingerprint_Generator::_do_add_tails(const Molecule & m,
                                        const Set_of_Atoms & expand_from,
                                        const int radius,
                                        const ECFP_Args<T> & ecfp_args,
                                        unsigned int sum_so_far,
                                        const int each_shell_gets_different_fingerprint,
                                        S * sfc) const
{
  const int n = expand_from.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = expand_from[i];

    const Atom * a = m.atomi(j);

    const int acon = a->ncon();

    for (int k = 0; k < acon; ++k)
    {
      const Bond * b = a->item(k);

      const atom_number_t l = b->other(j);

      if (EC_FINISHED == ecfp_args.processing_status(l))
        continue;

      const auto bc = _bond_constant(b);

      sfc[0].hit_bit(7879 * sum_so_far + bc * ecfp_args.atype(l));
    }
  }

  return;
}
template int EC_Fingerprint_Generator::generate_fingerprint<unsigned int>(Molecule&, unsigned int const*, int const*, int, Sparse_Fingerprint_Creator&) const;
template void EC_Fingerprint_Generator::generate_fingerprint<unsigned int>(Molecule&, int, unsigned int const*, int const*, int, resizable_array<unsigned int>&) const;
