#define COMPILING_RXN_FILE

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "rxn_file.h"
#include "atom_typing.h"

static void
form_individual_atom_difference_bit(int reagent_atom_type, int product_atom_type,
                                    Sparse_Fingerprint_Creator & sfc)
{
  if (reagent_atom_type == product_atom_type)
  {
    sfc.hit_bit(reagent_atom_type);    // is this correct?
    return;
  }

  sfc.hit_bit(2 * reagent_atom_type);
  sfc.hit_bit(2 * product_atom_type);

  if (reagent_atom_type > product_atom_type)
    std::swap(reagent_atom_type, product_atom_type);

  sfc.hit_bit(9 * product_atom_type - 2 * reagent_atom_type);    // arbitrary numbers
}

static int
bond_constant(const Bond & b)
{
  if (b.is_aromatic())
    return 102;
  if (b.is_single_bond())
    return 82;
  if (b.is_double_bond())
    return 39;

  return 14;
}

/*
  Do these atoms exhibit any changed bonding
*/

template <typename T>
void
form_possible_lost_or_changed_bond_fingerprint(T reagent_atom_type_1, const Bond * b1, T reagent_atom_type_2,
                                      int product_atom_type_1, const Bond * b2, int product_atom_type_2,
                                      Sparse_Fingerprint_Creator & sfc)
{
  if (NULL == b1 && NULL == b2)    // not bonded in either reagent or in product
    return;

  if (reagent_atom_type_1 < reagent_atom_type_2)
    std::swap(reagent_atom_type_1, reagent_atom_type_2);

  if (product_atom_type_1 < product_atom_type_2)
    std::swap(product_atom_type_1, product_atom_type_2);

  if (NULL == b2)   // bond destroyed, atoms not connected in product
  {
    sfc.hit_bit(64 * reagent_atom_type_1 - 6 * reagent_atom_type_2 + bond_constant(*b1));
    sfc.hit_bit(64 * reagent_atom_type_1 - 9 * reagent_atom_type_2 + bond_constant(*b1) + 7 * product_atom_type_1 - 3 * product_atom_type_2);

    return;
  }

  if (NULL == b1)    // new bond formed
  {
    sfc.hit_bit(83 * product_atom_type_1 - 19 * product_atom_type_2 + bond_constant(*b2));
    sfc.hit_bit(83 * product_atom_type_1 - 19 * product_atom_type_2 + bond_constant(*b2) + 11 * reagent_atom_type_1 - 2 * reagent_atom_type_2);

    return;
  }

  if (b1->same_bond_type(*b2))
    return;

  const auto b1c = bond_constant(*b1);
  const auto b2c = bond_constant(*b2);

  sfc.hit_bit(101 * reagent_atom_type_1 - 15 * reagent_atom_type_2 + b1c);
  sfc.hit_bit(901 * reagent_atom_type_1 - 15 * reagent_atom_type_2 + 73 * b1c - 21 * b2c);
  sfc.hit_bit(287 * product_atom_type_1 - 53 * product_atom_type_2 + b2c);
  sfc.hit_bit(2877 * product_atom_type_1 - 9 * product_atom_type_2 + 62 * b1c - 9 * b2c);
}

template <typename T>
int
RXN_File::reaction_fingerprint(const int * changing_atoms,
                               const T * reagent_atom_type,
                               const Product_Atom_Types & product_atom_type,
                               const int expand,
                               Sparse_Fingerprint_Creator & sfc)
{
  const ISIS_RXN_FILE_Molecule & r = _reagent[0];

  const int matoms = r.natoms();
//cerr << "Reagent has " << matoms << " atoms, np " << _np << endl;

// In this intiial loop, we look for changes from reagent to product. But we also need to then look back the other way.
// Some temporary structures to hold the atoms we need to examine

  Set_of_Atoms * changing_in_reagent = new Set_of_Atoms[_np]; std::unique_ptr<Set_of_Atoms[]> free_changing_in_reagent(changing_in_reagent);
  Set_of_Atoms * changing_in_product = new Set_of_Atoms[_np]; std::unique_ptr<Set_of_Atoms[]> free_changing_in_product(changing_in_product);

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == changing_atoms[i])
      continue;

    const int mi = r.atom_map(i);

    const int p = _product_locator[mi];
    if (p < 0)
      continue;

    const atom_number_t j = _product[p].which_is_mapped_atom(mi);
    if (j < 0)           // should not happen
      continue;

    changing_in_reagent[p].add(i);
    changing_in_product[p].add(j);

    form_individual_atom_difference_bit(reagent_atom_type[i], product_atom_type.atom_type(p, j), sfc);

    const Atom * a = r.atomi(i);

    const int icon = a->ncon();

    for (int k = 0; k < icon; ++k)
    {
      const Bond * b1 = a->item(k);

      const atom_number_t l = b1->other(i);    // bonded to I in the reagent

      if (0 == changing_atoms[l])
        continue;

      const int ml = r.atom_map(l);
      if (ml <= 0)
        continue;

      const int pl = _product_locator[ml];
      if (pl < 0)    // should not happen
        continue;

      const atom_number_t o = _product[pl].which_is_mapped_atom(ml);
      if (0 < 0)    // should not happen
        continue;

      const Bond * b2 = _product[pl].bond_between_atoms_if_present(j, o);
      form_possible_lost_or_changed_bond_fingerprint(reagent_atom_type[i], b1, reagent_atom_type[l],
                                      product_atom_type.atom_type(p, j), b2, product_atom_type.atom_type(pl, o),
                                      sfc);
    }
  }

  for (int i = 0; i < _np; ++i)
  {
//  cerr << " in product " << changing_in_product[i] << " in eagent " << changing_in_reagent[i] << endl;
    _reaction_fingerprint_rev(i, reagent_atom_type, product_atom_type, changing_in_reagent[i], changing_in_product[i], sfc);
  }

  if (0 == expand)
    return 1;

  return 1;
}

template <typename T>
int
RXN_File::_reaction_fingerprint_rev(const int p,       // processing product molecule P
                                    const T * reagent_atom_type,
                                    const Product_Atom_Types & product_atom_type,
                                    const Set_of_Atoms & changing_in_reagent,
                                    const Set_of_Atoms & changing_in_product,
                                    Sparse_Fingerprint_Creator & sfc)
{
  assert (changing_in_reagent.size() == changing_in_product.size());

  const int n = changing_in_product.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t a11 = changing_in_reagent[i];
    const atom_number_t a21 = changing_in_product[i];

    form_individual_atom_difference_bit(reagent_atom_type[i], product_atom_type.atom_type(p, a21), sfc);

    for (int j = i + 1; j < n; ++j)
    {
      const atom_number_t a12 = changing_in_reagent[j];
      const atom_number_t a22 = changing_in_product[j];

      const Bond * b1 = _reagent[0].bond_between_atoms_if_present(a11, a12);
      const Bond * b2 = _product[p].bond_between_atoms_if_present(a21, a22);

      form_possible_lost_or_changed_bond_fingerprint(reagent_atom_type[a11], b1, reagent_atom_type[a12],
                                      product_atom_type.atom_type(p, a21), b2, product_atom_type.atom_type(p, a22),
                                      sfc);
    }
  }

  return 1;
}

template int RXN_File::reaction_fingerprint<unsigned int>(int const*, unsigned int const*, Product_Atom_Types const&, int, Sparse_Fingerprint_Creator&);
