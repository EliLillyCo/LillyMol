#include <stdlib.h>
#include <assert.h>

#include "molecule.h"
#include "mdl.h"

Set_of_Atoms::Set_of_Atoms ()
{
//cerr << "Set_of_Atoms::Set_of_Atoms called\n";
}

Set_of_Atoms::Set_of_Atoms(int initial_size) : resizable_array<atom_number_t>(initial_size)
{
}

Set_of_Atoms::Set_of_Atoms(const Set_of_Atoms & rhs) : resizable_array<atom_number_t>(rhs)
{
//cerr << "Set_of_Atoms::Set_of_Atoms(const Set_of_Atoms &) called\n";
//resizable_array<atom_number_t>::operator =(rhs);
}

int
Set_of_Atoms::write(std::ostream & os) const
{
  assert (os.good());
  os << "Atoms";
  for (int i = 0; i < _number_elements; i++)
  {
    os << " " << _things[i];
  }

  os << endl;
  return os.good();
}

std::ostream &
operator << (std::ostream & os, const Set_of_Atoms & s)
{
  int ns = s.number_elements();

  os << "Atoms";
  for (int i = 0; i < ns; i++)
  {
    os << ' ' << s.item(i);
  }

  return os;
}

/*ostream &
Set_of_Atoms::write_atom_numbers(ostream & os, int offset, int width) const
{

  for (int i = 0; i < _number_elements; i++)
  {
    if (width > 0)
      os.width(width);

    os << (_things[i] + offset);
  }

  return os;
}*/

std::ostream &
operator << (std::ostream & os, const Set_of_Atoms * s)
{
  return os << (*s);
}

int
Set_of_Atoms::increment_vector(int * v, int increment) const
{
  for (int i = 0; i < _number_elements; i++)    // should we worry about INVALID_ATOM_NUMBER
  {
    const atom_number_t a = _things[i];
    v[a] += increment;
  }

  return _number_elements;
}

template int Set_of_Atoms::set_vector(int *, const int) const;
template int Set_of_Atoms::set_vector(float *, const float) const;

/*
  Sometimes we need a fast way of knowing whether or not
  all the members of the set are also set in some array -
  first application in kekule determinations
*/

int
Set_of_Atoms::all_members_set_in_array(const int * v, int target) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];
    if (j < 0)
      continue;

    if (target != v[j])
      return 0;
  }

  return 1;
}

int
Set_of_Atoms::all_members_non_zero_in_array(const int * v) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];
    if (j < 0)
      continue;

    if (0 == v[j])
      return 0;
  }

  return 1;
}

int
Set_of_Atoms::number_members_non_zero_in_array(const int * v) const
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];
    if (j < 0)
      continue;

    if (0 != v[j])
      rc++;
  }

  return rc;
}

int
Set_of_Atoms::any_members_set_in_array(const int * haystack, int needle) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (j < 0)
      continue;

    if (needle == haystack[j])
      return 1;
  }

  return 0;
}

int
Set_of_Atoms::count_members_set_in_array(const int * haystack, int needle) const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (j < 0)
      continue;

    if (needle == haystack[j])
      rc++;
  }

  return rc;
}

int
Set_of_Atoms::any_members_set_in_array(const int * haystack) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (j < 0)
      continue;

    if (haystack[j])
      return 1;
  }

  return 0;
}

/*
  An atom is being lost from the owning molecule. Adjust any atom
  numbers.
*/

int
Set_of_Atoms::adjust_for_loss_of_atom(atom_number_t lost_atom,
                                       int keep_lost_atom_in_list)
{
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    atom_number_t a = _things[i];
    if (a > lost_atom)
      _things[i]--;
    else if (a < lost_atom)    // no change
      ;
    else if (keep_lost_atom_in_list)
      _things[i] = INVALID_ATOM_NUMBER;
    else
      remove_item(i);
  }

  return 1;
}


int
Set_of_Atoms::offset_atom_numbers (int offset)
{
  for (int i = 0; i < _number_elements; i++) {    // should we worry about invalid_atom_number
  {
    _things[i] += offset;
  }

    assert (_things[i] >= 0);

  }

  return 1;
}

Set_of_Atoms &
Set_of_Atoms::operator = (const Set_of_Atoms & rhs)
{
  resizable_array<atom_number_t>::operator = (rhs);

  return *this;
}

int
Set_of_Atoms::write_as_mdl_v30_collection_block(const const_IWSubstring & zname,
                                                 const const_IWSubstring & subname,
                                                 std::ostream & output) const
{
  output << "M  V30 BEGIN COLLECTION\n";
  output << "M  V30 " << zname;
  if (subname.length() > 0)
    output << '/' << subname;
  output << endl;

  IWString output_buffer;
  output_buffer.resize(200);

  output_buffer << "M  V30 ATOMS=" << _number_elements;

  for (int i = 0; i < _number_elements; i++)
  {
    output_buffer << ' ' << (_things[i] + 1);
  }

  write_v30_record(output_buffer, output);

  output << "M  V30 END COLLECTION\n";

  return output.good();
}

int
Set_of_Atoms::any_members_in_common(const Set_of_Atoms & rhs) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];
    if (j < 0)
      continue;

    if (rhs.contains(j))
      return 1;
  }

  return 0;
}

atom_number_t
Set_of_Atoms::first_member_in_common(const Set_of_Atoms & rhs) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (rhs.contains(j))
      return j;
  }
  return INVALID_ATOM_NUMBER;
}


int
Set_of_Atoms::members_in_common(const Set_of_Atoms & rhs) const
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (rhs.contains(j))
      rc++;
  }

  return rc;
}

int
Set_of_Atoms::contains_atoms(const atom_number_t a1, const atom_number_t a2) const
{
  int nfound = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if (a1 == _things[i])
      nfound++;
    else if (a2 == _things[i])
      nfound++;
    else
      continue;

    if (2 == nfound)
      return 1;
  }

  return 0;
}
