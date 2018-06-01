#ifndef IW_SET_OF_ATOMS_H
#define IW_SET_OF_ATOMS_H

#include "iwmtypes.h"
#include "iwaray.h"

class Molecule;

class Set_of_Atoms : public resizable_array<atom_number_t>
{
  private:
  public:
    Set_of_Atoms ();
    Set_of_Atoms (int);
    Set_of_Atoms (const Set_of_Atoms &);

    Set_of_Atoms & operator = (const Set_of_Atoms &);

    int write (std::ostream &) const;

    int increment_vector (int *, int = 1) const;

    template <typename T> int set_vector(T * v, const T) const;

    int any_members_set_in_array (const int * haystack, int needle) const;
    int any_members_set_in_array (const int *) const;    // looks for non-zero items

    int count_members_set_in_array (const int * haystack, int needle) const;

    int offset_atom_numbers (int);
    std::ostream & write_atom_numbers (std::ostream &, int = 0, int = 0) const;

    int adjust_for_loss_of_atom (atom_number_t, int = 0);
    int all_members_set_in_array (const int * v, int target) const;

    int all_members_non_zero_in_array (const int * v) const;
    int number_members_non_zero_in_array (const int * v) const;

    int any_members_in_common (const Set_of_Atoms &) const;
    atom_number_t first_member_in_common (const Set_of_Atoms &) const;
    int members_in_common (const Set_of_Atoms &) const;

    int write_as_mdl_v30_collection_block (const const_IWSubstring & zname, const const_IWSubstring & subname, std::ostream &) const;

    int contains_atoms (const atom_number_t, const atom_number_t) const;

    template <typename T> void each (Molecule &, T &) const;
    template <typename V, typename O> int each (const V *, O) const;
    template <typename V, typename O> void each (V *, O) const;
};

template <typename T>
void
Set_of_Atoms::each (Molecule & m, T & o) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    o (m, _things[i]);
  }

  return;
}

template <typename V, typename O>
int
Set_of_Atoms::each (const V * v, O o) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (! o(v[_things[i]]))
      return 0;
  }

  return 1;
}

template <typename V, typename O>
void
Set_of_Atoms::each (V * v, O o) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    o(v[_things[i]]);
  }

  return;
}

std::ostream & operator << (std::ostream &, const Set_of_Atoms &);
//std::ostream & operator << (std::ostream &, const Set_of_Atoms *);    bad idea


template <typename T>
int
Set_of_Atoms::set_vector (T * v, const T x) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i];
    if (j >= 0)
      v[j] = x;
  }

  return _number_elements;
}

#endif
