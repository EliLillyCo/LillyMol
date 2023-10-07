#include <stdlib.h>
#include <iostream>

#include <iostream>
#include <sstream>

#include "Foundational/iwmisc/misc.h"
#include "jw_path_with_values.h"

using std::cerr;
using std::endl;

JW_Path::JW_Path (unsigned a1, unsigned int a2, int total_size)
{
  _length = 2;
  _a1 = a1;
  _a2 = a2;

  _max_path_length = total_size;
  _numbers_in_path.allocate_space_for_bits (_max_path_length);
  _numbers_in_path.set (a1);
  _numbers_in_path.set (a2);
  return;
}

JW_Path::JW_Path (unsigned a1, unsigned int a2, unsigned int a3, int total_size)
{
  _length = 3;
  _a1 = a1;
  _a2 = a2;

  _max_path_length = total_size;
  _numbers_in_path.allocate_space_for_bits (_max_path_length);
  _numbers_in_path.set (a1);
  _numbers_in_path.set (a2);
  _numbers_in_path.set (a3);
  return;
}

JW_Path::JW_Path (JW_Path & path, unsigned a1, unsigned int a2)
{
  _length = path.length () +2;
  _a1 = a1;
  _a2 = a2;

  _max_path_length = path.max_path_length();

  _numbers_in_path.allocate_space_for_bits (_max_path_length);

  for (int i=0; i<_max_path_length; i++)
    if (path._numbers_in_path.is_set (i))
      _numbers_in_path.set (i);
  
  _numbers_in_path.set (a1);
  _numbers_in_path.set (a2);

  return;
}

void JW_Path::debug_print_path ()
{
  int n = _max_path_length;
  for (int i=0; i<n; i++)
    if (_numbers_in_path.is_set(i))
      cerr<<"\t"<<i;
  cerr<<endl;
};


int Path_with_values::_number_values = 1;

Path_with_values::Path_with_values (atom_number_t a1, atom_number_t a2, double * values, int n_atoms)
{
  _length = 2;
  _a1 = a1;
  _a2 = a2;
  
  _max_path_length = n_atoms;
  _atom_numbers_in_path = new_int (_max_path_length);

  _atom_numbers_in_path [a1] = 1;
  _atom_numbers_in_path [a2] = 1;
  
  _values = new double [_number_values];
  
  copy_vector (_values, values, _number_values);
  
  return;
};
  
Path_with_values::Path_with_values (atom_number_t a1, atom_number_t a2, atom_number_t middle_atom, double * values, int n_atoms)
{
  _length = 3;
  _a1 = a1;
  _a2 = a2;
  
  _max_path_length = n_atoms;
  _atom_numbers_in_path = new_int (_max_path_length);
  
  _values = new double [_number_values];

  copy_vector (_values, values, _number_values);
  
  _atom_numbers_in_path [a1] = 1;
  _atom_numbers_in_path [a2] = 1;
  _atom_numbers_in_path [middle_atom] = 1;
};


Path_with_values::Path_with_values (const Path_with_values &p, atom_number_t a1, atom_number_t a2, double * values)
{
  _length = p._length + 2;
  
  _a1 = a1;
  _a2 = a2;
  
  _max_path_length = p._max_path_length;
  _atom_numbers_in_path = new int [_max_path_length];
  
  copy_vector (_atom_numbers_in_path, p._atom_numbers_in_path, _max_path_length);
  
  _atom_numbers_in_path[a1] = 1;
  _atom_numbers_in_path[a2] = 1;
  
  _values = new double[_number_values];
  copy_vector (_values, values, _number_values);
};

Path_with_values::~Path_with_values () 
{
  delete [] _atom_numbers_in_path;
  delete [] _values;
};

void Path_with_values::set_number_values ( int number_values) { _number_values = number_values;};
int Path_with_values::number_values () { return _number_values;};

int Path_with_values::value (int i, double & value) 
{
  if ((i<_number_values) && (i>=0))
    {
      value = _values [i];
      return 1;
    }
  return 0;
};

void Path_with_values::debug_print_path ()
{
  int n = _max_path_length;
  for (int i=0; i<n; i++)
    if (_atom_numbers_in_path [i])
      cerr<<"\t"<<i;
  cerr<<endl;
};

Path_with_unsigned_int_values::Path_with_unsigned_int_values (atom_number_t a, int n_atoms, unsigned int value1, 
						 unsigned int value2, unsigned int value3)
{
  _length = 0;
  _max_path_length = n_atoms;
  _value1 = value1;
  _value2 = value2;
  _value3 = value3;
  
  _outmost_layer = new_int (_max_path_length);

  _outmost_layer[a] = 1;
  
  _atom_numbers_in_path = new_int (_max_path_length);

  _atom_numbers_in_path [a] = 1;
};

Path_with_unsigned_int_values::Path_with_unsigned_int_values (const Set_of_Atoms & a_set,
                                                 int n_atoms, unsigned int value1, 
						 unsigned int value2, unsigned int value3)
{
  _length = 0;
  _max_path_length = n_atoms;
  _value1 = value1;
  _value2 = value2;
  _value3 = value3;
  
  _outmost_layer = new_int (_max_path_length);

  a_set.set_vector (_outmost_layer, 1);
  
  _atom_numbers_in_path = new int [_max_path_length];

  copy_vector (_atom_numbers_in_path, _outmost_layer, _max_path_length);

  return;
};

Path_with_unsigned_int_values::Path_with_unsigned_int_values (const Path_with_unsigned_int_values &p,
                                                 const Set_of_Atoms & a_set, unsigned int value1, 
						 unsigned int value2, unsigned int value3)
{
  _value1 = value1;
  _value2 = value2;
  _value3 = value3;
  
  _length = p.length() +1;
  
  _max_path_length = p.max_path_length();
  
  _outmost_layer = new_int (_max_path_length);

  a_set.set_vector (_outmost_layer, 1);
  
  _atom_numbers_in_path = new int [_max_path_length];

  copy_vector (_atom_numbers_in_path, p._atom_numbers_in_path, _max_path_length);

  a_set.set_vector (_atom_numbers_in_path, 1);

  return;
};

Path_with_unsigned_int_values::~Path_with_unsigned_int_values ()
{
  delete [] _outmost_layer;
  delete [] _atom_numbers_in_path;
}

unsigned int Path_with_unsigned_int_values::value1 () { return _value1; };
unsigned int Path_with_unsigned_int_values::value2 () { return _value2; };
unsigned int Path_with_unsigned_int_values::value3 () { return _value3; };

atom_number_t Path_with_unsigned_int_values::outmost_layer_contains (int i) { return _outmost_layer [i];};

void Path_with_unsigned_int_values::set_atom_in_path (atom_number_t i)
{
  if ((i>-1) && (i<_max_path_length))
    _atom_numbers_in_path [i] = 1;
};

void Path_with_unsigned_int_values::set_outpost_layer (atom_number_t i)
{
  if ((i>-1) && (i<_max_path_length))
    {
      _outmost_layer [i] = 1;
      _atom_numbers_in_path [i] = 1;
    }
};

int Path_with_unsigned_int_values::contains (atom_number_t atom) { return _atom_numbers_in_path [atom];};
