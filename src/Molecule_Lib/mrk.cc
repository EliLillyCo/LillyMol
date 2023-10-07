/*
  Stuff for reading and writing the Merck Force field tupe
*/

#include <iostream>
#include <iomanip>
using std::cerr;
using std::endl;
#ifdef UNIX
#include <unistd.h>
#endif
#include <time.h>

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"
#include "rwmolecule.h"

/*
  Here's an example of the top of the file

5112_1                                                                9410315351
MOL eb90162  O  E =      92.7560   G =  3.39E-01  MMFF94 
   47    51
  -14.1063    13.6992    42.8776     6 63 0     1 C1     1LIGA -0.3316      MMFF
  -14.7019    12.9288    41.9075     6 64 0     2 C2     1LIGA -0.0540      MMFF

*/

int
Molecule::read_molecule_mrk_ds(iwstring_data_source & input)
{
  if (! input.next_record(_molecule_name))
  {
    cerr << "read mol eof mrk\n";
    return 0;
  }

  if ('1' != _molecule_name[79])
  {
    cerr << "Molecule::read_molecule_mrk_ds: column 80 not '1', proceeding...\n";
  }

  if (_molecule_name.length() > 70)
  {
    _molecule_name.iwtruncate(70);
    _molecule_name.strip_trailing_blanks();
  }

  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Molecule::read_molecule_mrk_ds: cannot read second record\n";
    return 0;
  }

  if (! buffer.starts_with("MOL "))
  {
    cerr << "molecule::read_molecule_mrk_ds: second record does not start with 'MOL ', line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "Molecule::read_molecule_mrk_ds: cannot read third record\n";
    return 0;
  }

  if (buffer.nwords() < 2)
  {
    cerr << "Molecule::read_molecule_mrk_ds: invalid atoms and bonds record, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

// Read any format for atoms and bonds

  const_IWSubstring token;
  int i = 0;

  (void) buffer.nextword(token, i);

  int na;

  if (! token.numeric_value(na) || na < 0)
  {
    cerr << "Molecule::read_molecule_mrk_ds: invalid number of atoms, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  (void) buffer.nextword(token, i);

  int nb;
  if (! token.numeric_value(nb) || nb < 0)
  {
    cerr << "Molecule::read_molecule_mrk_ds: invalid number of bonds, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  resize(na);

  for (int i = 0; i < na; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Molecule::read_molecule_mrk_ds: premature eof\n";
      return 0;
    }

    if (! _read_mrk_atom_record(buffer))
    {
      cerr << "Molecule::read_molecule_mrk_ds: invalid record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  while (_bond_list.number_elements() < nb)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Molecule::read_molecule_mrk_ds: premature eof reading bond records\n";
      return 0;
    }

    if (! _read_mrk_bond_record(buffer, na))
    {
      cerr << "Molecule::read_molecule_mrk_ds: invalid bond record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

int
Molecule::_read_mrk_bond_record(const const_IWSubstring & buffer, int na)
{
  const_IWSubstring token;
  int i = 0;

  while (buffer.nextword(token, i))
  {
    atom_number_t a1, a2;
    if (! token.numeric_value(a1) || a1 < 1 || a1 > na)
    {
      cerr << "Molecule::_read_mrk_bond_record: invalid a1 specification, " << na << " atoms in molecule\n";
      cerr << buffer << endl;
      return 0;
    }
    if (! buffer.nextword(token, i) || ! token.numeric_value(a2) || a2 < 1 || a2 > na || a2 == a1)
    {
      cerr << "Molecule::_read_mrk_bond_record: invalid a2 specification, " << na << " atoms in molecule\n";
      cerr << buffer << endl;
      return 0;
    }

    if (! buffer.nextword(token, i))
    {
      cerr << "Molecule::_read_mrk_bond_record: invalid bond specification\n";
      cerr << buffer << endl;
      return 0;
    }

    bond_type_t bt;
    if ('1' == token)
      bt = SINGLE_BOND;
    else if ('2' == token)
      bt = DOUBLE_BOND;
    else if ('3' == token)
      bt = TRIPLE_BOND;
    else
    {
      cerr << "Molecule::_read_mrk_bond_record: unrecognised bond type '" << token << "'\n";
      return 0;
    }

    a1--;
    a2--;

    add_bond(a1, a2, bt, 1);     // 1 means partially built molecule
  }

  return 1;
}

static const formal_charge_t charge_code_to_formal_charge [] = {0, 1, -1, 0, 2, -2, 3, -3, 4, -4};

int
Molecule::_read_mrk_atom_record(const const_IWSubstring & buffer)
{
  if (buffer.nwords() < 6)
  {
    cerr << "Molecule::_read_mrk_atom_record: invalid atom record\n";
    return 0;
  }

  const_IWSubstring token;
  int i = 0;

  coord_t c[3];
  for (int j = 0; j < 3; j++)
  {
    (void) buffer.nextword(token, i);
    if (! token.numeric_value(c[j]))
    {
      cerr << "Molecule::_read_mrk_atom_record: invalid coordinate '" << token << "'\n";
      return 0;
    }
  }

  (void) buffer.nextword(token, i);

  atomic_number_t z;
  if (! token.numeric_value(z) || z < 0 || z > HIGHEST_ATOMIC_NUMBER)
  {
    cerr << "Molecule::_read_mrk_atom_record: invalid atomic number '" << token << "'\n";
    return 0;
  }
  
  Atom * a = new Atom(z);
  a->setxyz(c[0], c[1], c[2]);

// We have cases where the atom type column may be absent, so we need to work with strict column numbers

  buffer.from_to(42, 42, token);    // charge code

  int charge_code;
  if (! token.numeric_value(charge_code) || charge_code < 0 || charge_code > 9)
  {
    cerr << "Molecule::_read_mrk_atom_record: invalid charge code '" << token << "'\n";
    return 0;
  }

  if (3 == charge_code)    // radical
  {
    a->set_implicit_hydrogens(0, 1);
    a->set_implicit_hydrogens_known(1);
  }
  else if (0 == charge_code)
    ;
  else
  {
    formal_charge_t fc = charge_code_to_formal_charge[charge_code];
    a->set_formal_charge(fc);
  }

  add(a);

// The partial charge is in columns 63-70
  
  if (buffer.length() < 70)
  {
    if (nullptr != _charges)    // yipes, charges present but nothing for this record
    {
      cerr << "Molecule::_read_mrk_atom_record: missing partial charge\n";
      return 0;
    }

    return 1;    // ok, charges not present
  }

  buffer.from_to(62, 69, token);
  token.strip_leading_blanks();

  charge_t q;
  if (! token.numeric_value(q) || ! reasonable_atomic_partial_charge_value(q))
  {
    cerr << "Molecule::_read_mrk_atom_record: invalid or unreasonable charge value '" << token << "'\n";
    return 0;
  }

  if (static_cast<charge_t>(0.0) == q && nullptr == _charges)    // no non-zero charges encountered yet
    return 1;

  if (nullptr == _charges)
    allocate_charges();

  _charges->seti(_number_elements - 1, q);

  return 1;
}

static IWString * digit5 = nullptr;
static int max_digit5 = 0;

static IWString digit1[4];

static void
initialise_digit5(int m)
{
  assert(m > 0);

  max_digit5 = m;

  if (nullptr != digit5)
    delete [] digit5;

  digit5 = new IWString[max_digit5];

  for (int i = 0; i < max_digit5; i++)
  {
    digit5[i] << i;

    digit5[i].shift(5 - digit5[i].length(), ' ');
  }

  digit1[0] = "0";
  digit1[1] = "1";
  digit1[2] = "2";
  digit1[3] = "3";

  return;
}

static void
append_with_leading_zero(IWString & buffer, 
                         int i,
                         int nspaces)
{
  if (i < 10)
    buffer << '0';

  if (3 == nspaces && i < 100)
    buffer << '0';

  buffer << i;
}

static void
write_partial_charge(std::ostream & os,
                     charge_t q)
{
#if defined(__GNUG__) || defined (IW_INTEL_COMPILER)
  const std::ios::fmtflags old_flags = os.flags(std::ios::fixed);
#else
  const long old_flags = os.flags(std::ios::fixed);
#endif

  int old_precision  = os.precision(4);

  os << std::setw(8) << q;

  os.precision(old_precision);
  os.flags(old_flags);

  return;
}

int
Molecule::write_molecule_mrk(std::ostream & output)
{
  if (nullptr == digit5)
    initialise_digit5(200);
  else if (_number_elements >= max_digit5)
    initialise_digit5(_number_elements + 20);

  IWString title(_molecule_name);

  if (title.length() > 70)
    title.iwtruncate(70);
  else if (title.length() < 70)
    title.extend(70, ' ');

  time_t t = time(NULL);
  
#ifdef _WIN32
  struct tm * tm = NULL;
#else
  struct tm * tm = localtime(&t);
#endif

  append_with_leading_zero(title, tm->tm_year - 100, 2);    // y3k bug!
  append_with_leading_zero(title, tm->tm_yday + 1, 3);
  append_with_leading_zero(title, tm->tm_hour, 2);
  append_with_leading_zero(title, tm->tm_min, 2);

  title << '1';

  output << title << '\n';

  output << "MOL ";
#ifdef _WIN32
  output << "        ";
#else
  IWString usr = getlogin();
  if (usr.length() > 8)
    usr.iwtruncate(8);
  else if (usr.length() < 8)
    usr.extend(8 - usr.length(), ' ');

  output << usr;
#endif

  output << " F \n";

  int nb = _bond_list.number_elements();

  output << digit5[_number_elements] << ' ' << digit5[nb] << '\n';

  IWString buffer;
  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    a->write_coordinates(output, 1);     // last arg 1 means include space between 10 column groups

    if (kInvalidAtomicNumber == a->atomic_number())
      buffer << digit5[0];
    else
      buffer << digit5[a->atomic_number()];

    buffer << ' ';

    buffer << " 0";      // atom type not read on input

    buffer << ' ';
    formal_charge_t fc = a->formal_charge();

    if (0 == fc)
      buffer << '0';
    else if (1 == fc)
      buffer << '1';
    else if (-1 == fc)
      buffer << '2';
    else if (2 == fc)
      buffer << '4';
    else if (-2 == fc)
      buffer << '5';
    else if (3 == fc)
      buffer << '6';
    else if (-3 == fc)
      buffer << '7';
    else if (4 == fc)
      buffer << '8';
    else if (-4 == fc)
      buffer << '9';
    else
      buffer << '*';

    buffer << ' ' << digit5[i + 1];    // atom sequence number not read on input

    IWString atom_name;
    atom_name << a->atomic_symbol();
    atom_name << (i + 1);
    if (atom_name.length() > 4)
      atom_name.iwtruncate(4);
    else
      atom_name.extend(4, ' ');

    buffer << ' ' << atom_name;

    buffer << "   1";    // residue name
    buffer << "LIGA";    // residue type

//  Write the buffer if we will use the C++ library for writing charges

    if (nullptr != _charges)
    {
      output << buffer;
      write_partial_charge(output, _charges->item(i));
      buffer.resize_keep_storage(0);
    }
    else
      buffer << "  0.0000";

    buffer << "      LIGA\n";

    output << buffer;

    buffer.resize_keep_storage(0);
  }

  buffer.resize(80);

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    atom_number_t a1 = b->a1() + 1;
    atom_number_t a2 = b->a2() + 1;

    if (buffer.length())
      buffer << "  ";

    buffer << digit5[a1] << ' ' << digit5[a2] << ' ';
    if (b->is_single_bond())
      buffer << " 1";
    else if (b->is_double_bond())
      buffer << " 2";
    else if (b->is_triple_bond())
      buffer << " 3";
    else
    {
      cerr << "Molecule::write_molecule_mrk: what kind of bond between atoms " << a1 << " and " << a2 << endl;
      buffer << '0';
    }

    if (4 == i % 5)
    {
      output << buffer << '\n';
      buffer.resize_keep_storage(0);
    }
  }

  if (buffer.length())
  {
    buffer += '\n';
    output << buffer;
  }

  return output.good();
}
