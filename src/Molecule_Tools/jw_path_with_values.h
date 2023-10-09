#include "Foundational/iwbits/iwbits.h"
#include "Molecule_Lib/molecule.h"

class JW_Path 
{
 private:
  int _length;
  unsigned int _a1;
  unsigned int _a2;

  int _max_path_length;
  IW_Bits_Base _numbers_in_path;
 public:
  JW_Path (unsigned int, unsigned int, int);
  JW_Path (unsigned int, unsigned int, unsigned int, int);
  JW_Path (JW_Path &, unsigned int, unsigned int);

  ~JW_Path () {};

  int max_path_length () const { return _max_path_length;};
  int length () const { return _length;}
  unsigned int a1() const { return _a1;}
  unsigned int a2() const {return _a2;}

  static void set_number_values (int);
  static int number_values ();
  
  int contains (unsigned int  a) const { return _numbers_in_path.is_set (a);}
  void debug_print_path ();
};

class Path_with_values
{
 private:
  int _length;
  atom_number_t _a1;
  atom_number_t _a2;
  
  int _max_path_length;
  int * _atom_numbers_in_path;
  double * _values;
  
  static int _number_values;

 public:
  Path_with_values (atom_number_t, atom_number_t, double *, int);
  
  Path_with_values (atom_number_t, atom_number_t, atom_number_t, double *, int);

  Path_with_values (const Path_with_values &, atom_number_t, atom_number_t, double *);

  ~Path_with_values ();

  int max_path_length () const { return _max_path_length;}
  int length () const { return _length;}
  atom_number_t a1 () const { return _a1;}
  atom_number_t a2 () const { return _a2;}

  static void set_number_values (int);
  static int number_values ();

  int value (int, double &);

  const double * raw_values () const { return _values;}

  int contains (atom_number_t a) const { return _atom_numbers_in_path[a];}
  
  void debug_print_path ();
};

class Path_with_unsigned_int_values
{
private:
  int _length;
  int _max_path_length;
  unsigned int _value1;
  unsigned int _value2;
  unsigned int _value3;

  int * _outmost_layer;
  int * _atom_numbers_in_path;
public:
  Path_with_unsigned_int_values (atom_number_t, int, unsigned int, unsigned int, unsigned int);
  Path_with_unsigned_int_values (const Set_of_Atoms &, int, unsigned int, unsigned int, unsigned int);
  Path_with_unsigned_int_values (const Path_with_unsigned_int_values &, const Set_of_Atoms &, unsigned int, unsigned int, unsigned int);
  ~Path_with_unsigned_int_values ();

  int max_path_length () const { return _max_path_length;}
  int length () const { return _length;}
 
  unsigned int value1 ();
  unsigned int value2 ();
  unsigned int value3 ();

  atom_number_t outmost_layer_contains (int);

  void set_outpost_layer (atom_number_t);
  void set_atom_in_path (atom_number_t);

  int contains (atom_number_t);
};

