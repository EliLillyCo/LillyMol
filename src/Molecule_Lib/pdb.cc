#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <assert.h>
#include <time.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "iwstring_data_source.h"
#include "misc.h"

#include "molecule.h"
#include "rwmolecule.h"

/*
  Jun 04. Need to preserve atom names when reading and writing
  pdb files. The right way to do that would be to include an
  atom name attribute with each Atom, but I don't want to do
  that for efficiency reasons. Therefore we have a kludge...

  Note that Molecule object operations that remove atoms,
  swap atoms and things like that will not change this array!
*/

static int store_pdb_atom_information = 0;

void
set_store_pdb_atom_information (int s)
{
  store_pdb_atom_information = s;
}

PDB_Stored_Atom_Information::PDB_Stored_Atom_Information()
{
  _atom_number = -1;

  return;
}

void
PDB_Stored_Atom_Information::_do_copy (const PDB_Stored_Atom_Information & rhs)
{
  _atom_number = rhs._atom_number;
  _atom_name = rhs._atom_name;
  _residue_name = rhs._residue_name;
  _occupancy_factor = rhs._occupancy_factor;
  _temperature_factor = rhs._temperature_factor;

  return;
}

PDB_Stored_Atom_Information::PDB_Stored_Atom_Information (const PDB_Stored_Atom_Information & rhs)
{
  _do_copy(rhs);

  return;
}

PDB_Stored_Atom_Information &
PDB_Stored_Atom_Information::operator = (const PDB_Stored_Atom_Information & rhs)
{
  _do_copy (rhs);

  return *this;
}

template class resizable_array_p<PDB_Stored_Atom_Information>;
template class resizable_array_base<PDB_Stored_Atom_Information *>;

static resizable_array_p<PDB_Stored_Atom_Information> stored_pdb_atom_information;

const resizable_array_p<PDB_Stored_Atom_Information> &
stored_pdb_atom_information_last_molecule_read()
{
  return stored_pdb_atom_information;
}

resizable_array_p<PDB_Stored_Atom_Information> &
stored_pdb_atom_information_last_molecule_read(int notused)
{
  return stored_pdb_atom_information;
}

static int use_stored_atom_information_when_writing_pdb_files = 0;

void
set_use_stored_atom_information_when_writing_pdb_files (int s)
{
  use_stored_atom_information_when_writing_pdb_files = s;
}

static int write_pdb_files_in_fragment_order = 0;

void
set_write_pdb_files_in_fragment_order(int s)
{
  write_pdb_files_in_fragment_order = s;

  return;
}

static int number_by_element_count = 0;

void
set_pdb_number_by_element_count(int s)
{
  number_by_element_count = s;
}

static int number_within_sequence = 0;

void
set_pdb_number_within_sequence(int s)
{
  number_within_sequence = s;
}

static Atom *
atom_from_pdb_element (const const_IWSubstring & sym)
{
  int slen = sym.length();

  if (slen > 4)
  {
    cerr << "atom atom_from_pdb_element: cannot parse '" << sym << "'\n";
    return NULL;
  }
  
//  If there is just one character, assume that it is C N O F P S

  if (1 == slen)
  {
    const Element * e = get_element_from_symbol_no_case_conversion(sym);
    if (NULL == e)
    {
      cerr << "atom_from_pdb_element:unrecognised element '" << sym << "'\n";
      return NULL;

    }

    return new Atom (e);
  }

  const_IWSubstring mysym(sym);

// Strip off trailing numerics

  while (isdigit(mysym.last_item()))
  {
    mysym.chop();
    if (0 == mysym.length())
    {
      cerr << "atom_from_pdb_element:numeric element specifier invalid '" << sym << "'\n";
      return 0;
    }
  }

  const Element * e = get_element_from_symbol_no_case_conversion(mysym);

  if (NULL != e)
    return new Atom(e);

// Now things get ambiguous. We just focus on the first letter of the name. 

  if (mysym.starts_with("CL"))
    return new Atom(17);

  if (mysym.starts_with("BR"))
    return new Atom(35);

  if (mysym.starts_with("C"))
    return new Atom(6);

  if (mysym.starts_with("N"))
    return new Atom(7);

  if (mysym.starts_with("O"))
    return new Atom(8);

  if (mysym.starts_with("F"))
    return new Atom(7);

  if (mysym.starts_with("P"))
    return new Atom(15);

  if (mysym.starts_with("S"))
    return new Atom(16);

  cerr << "No element for '" << sym << "'\n";

  return NULL;
}

static int
parse_pdb_atom_record (IWString & buffer,     // not const
                       resizable_array<int> & anumber,
                       extending_resizable_array<int> & atom_number_to_atom_number,
                       resizable_array<Atom *> & atoms)
{
  const_IWSubstring token;
  int i = 6;            // note dangerous initialisation

  if (! buffer.nextword(token, i))   // sequence number
    return 0;

//cerr << "From " << buffer << " sequence number is " << token << endl;

  int atom_number;

  if (! token.numeric_value(atom_number))
    atom_number = -1;
  else if (atom_number > 0)
    ;
  else
  {
    cerr << "parse_pdb_atom_record:invalid sequence number '" << token << "'\n";
    return 0;
  }

  const_IWSubstring atom_name;
  if (! buffer.nextword(atom_name, i))
    return 0;

  const_IWSubstring residue_name;

  if (! buffer.nextword(residue_name, i))
    return 0;

  coord_t ix, iy, iz;
  if (3 != IW_SSCANF (buffer.chars () + 30, "%f %f %f", &ix, &iy, &iz))
  {
    cerr << "parse_pdb_atom_record:invalid coordinates\n";
    return 0;
  }

  Atom * a;
  if (buffer.length() >= 78 && isalpha(buffer[77]))
  {
    const_IWSubstring s;
    buffer.from_to(76, 77, s);     // columns 77 to 78
    s.strip_leading_blanks();
//  cerr << "From " << buffer << endl;
//  cerr << "Getting element from '" << s << "'\n";
    const Element * e = get_element_from_symbol_no_case_conversion(s);
    if (NULL == e)
    {
      cerr << "Molecule::read_molecule_pdb_ds:unrecognised element '" << s << "'\n";
      return 0;
    }
    a = new Atom(e);
  }
  else
  {
    a = atom_from_pdb_element (atom_name);

    if (NULL == a)
    {
      cerr << "read pdb: unrecognised element '" << atom_name << "'\n";
      return 0;
    }
  }

  a->setxyz (ix, iy, iz);

  if (atom_number > 0)
  {
    atom_number_to_atom_number[atom_number] = atoms.number_elements();
    anumber.add(atom_number);
  }

  atoms.add(a);

  if (store_pdb_atom_information)
  {
    PDB_Stored_Atom_Information * t = new PDB_Stored_Atom_Information();

    t->set_atom_number(atom_number + 1);
    t->set_atom_name(atom_name);
    t->set_residue_name(residue_name);

    IWString tmp(buffer);
    tmp.remove_leading_chars(55);
    int nw = tmp.nwords();

    if (nw)    // are there extra columns
    {
      IWString occupancy, temperature;

      tmp.word(0, occupancy);

      t ->set_occupancy(occupancy);
      if (nw > 1)
      {
        tmp.word(1, temperature);
        t->set_temperature(temperature);
      }
    }

    stored_pdb_atom_information.add (t);
  }

  return 1;
}

/*
  function to read a pdb file.
  We only process ATOM and CONECT records.
*/

int
Molecule::read_molecule_pdb_ds (iwstring_data_source & input)
{
  assert (ok ());
  assert (input.good ());

  input.set_strip_trailing_blanks(1);

  stored_pdb_atom_information.resize_keep_storage (0);

// Since atoms may come in in random order, we temporarily store the
// atoms in an array, and an array of corresponding atom numbers

  resizable_array<Atom *> atoms;
  resizable_array<atom_number_t> anumber;

  resizable_array_p<Bond> bonds;

/*
  A typical atom record looks like

ATOM      1  C   ACE     0     -37.000   8.810  17.821  1.00  0.00      3LDH    |216
*/
 
  extending_resizable_array<int> atom_number_to_atom_number;     // as read in to our index

  IWString buffer;
  while (input.next_record (buffer))
  {
    if (buffer.starts_with ("COMPND "))
    {
      set_name (buffer.substr (7));
      continue;
    }
    
    if (buffer.starts_with ("ATOM ") || buffer.starts_with ("HETATM"))
    {
      if (! parse_pdb_atom_record(buffer, anumber, atom_number_to_atom_number, atoms))
      {
        cerr << "Molecule::read_molecule_pdb_ds:invalid input '" << buffer << "', line " << input.lines_read() << endl;
        return 0;
      }
    }

    if (buffer.starts_with("TER"))
      continue;

/*
    A typical CONECT record looks like


CONECT 2556 2557 2558 2559 2578                                                 |3LDH2821

    At a minumum it must contain two atom numbers (specifying a bond)
    Being somewhat arbitrary, we impose a maximum of 8 connections.
    Illegal connections are simply ignored.
*/

#define PDB_MAX_CON 8
    
    if (buffer.starts_with ("CONECT "))
    {
      atom_number_t base;
      atom_number_t pdbc[PDB_MAX_CON];

      int ntokens = IW_SSCANF (buffer.chars (), "CONECT %d %d %d %d %d %d %d %d %d",
                    &base, &pdbc[0], &pdbc[1], &pdbc[2], &pdbc[3],
                           &pdbc[4], &pdbc[5], &pdbc[6], &pdbc[7]);
      if (ntokens < 2)
      {
        cerr << "read_pdb: bad record, line " << input.lines_read () <<
                ", '" << buffer << "'\n";
        return 0;
      }

      for (int i = 0; i < ntokens - 1; i ++)
      {
        if (pdbc[i] > base)
        {
          Bond * b = new Bond(base, pdbc[i], SINGLE_BOND);
          bonds.add(b);
        }
      }
    }
    else if (buffer.starts_with("ENDMDL"))
      ;
    else if (buffer.starts_with ("END"))
      break;
  }

  int na = atoms.number_elements();
  if (0 == na)
  {
    if (_molecule_name.length())
      cerr << "Molecule::read_molecule_pdb_ds:no atoms\n";

    return 0;
  }

  if (! resize(na))
    return 0;

// Make sure we have all the atoms

  const Atom ** tmp = new const Atom *[na]; std::unique_ptr<const Atom *[]> free_tmp(tmp);
  assert (NULL != tmp);

  for (int i = 0; i < na; i++)
  {
    tmp[i] = NULL;
  }

  for (int i = 0; i < na; i++)
  {
    atom_number_t j = anumber[i];

    int k = atom_number_to_atom_number[j];

    if (NULL != tmp[k])
    {
      cerr << "Molecule::read_molecule_pdb_ds:duplicate atom " << k << endl;
      return 0;
    }

    tmp[k] = atoms[i];
  }

// Did we get all the atoms

  for (int i = 0; i < na; i++)
  {
    if (NULL == tmp[i])
    {
      cerr << "Molecule::read_molecule_pdb_ds:no atom " << i << endl;
      return 0;
    }

    add(atoms[i]);
  }

  int nbonds = bonds.number_elements();

  for (int i = 0; i < nbonds; i++)
  {
    Bond * b = bonds[i];

    atom_number_t a1 = atom_number_to_atom_number[b->a1()];
    atom_number_t a2 = atom_number_to_atom_number[b->a2()];

    if (a1 < 0 || a2 < 0 || a1 >= na || a2 >= na)
    {
      cerr << "Molecule::read_molecule_pdb_ds:invalid bond " << b->a1() << " to " << b->a2() << endl;
      return 0;
    }

    if (are_bonded(a1, a2))
    {
      if (ignore_self_bonds())
        cerr << "Molecule::read_molecule_pdb_ds:ignoring self bond " << a1 << " to " << a2 << " in " << _molecule_name << endl;
      else
      {
        cerr << "Molecule::read_molecule_pdb_ds:atoms " << a1 << " and " << a2 << " in " << _molecule_name << " alread bonded\n";
        return 0;
      }
    }
    else
      add_bond(a1, a2, SINGLE_BOND, 1);
  }

  if (nbonds)
    check_bonding ();

  if (unconnect_covalently_bonded_non_organics_on_read ())
    _do_unconnect_covalently_bonded_non_organics ();

  if (store_pdb_atom_information && stored_pdb_atom_information.number_elements() != _number_elements)
  {
    cerr << "PDB stored name mismatch with atom count " << stored_pdb_atom_information.number_elements() << " and " << _number_elements << endl;
    abort();
  }

  if (discern_chirality_from_3d_coordinates() && 3 == highest_coordinate_dimensionality())
    (void) discern_chirality_from_3d_structure();

  return _number_elements;
}

int
Molecule::write_molecule_pdb (const char *fname, const IWString & comments)
{
  std::ofstream output (fname);
  if (! output.good ())
  {
    cerr << "Molecule::write_molecule_pdb: cannot open '" << fname << "'\n";
    return 0;
  }

  int rc = write_molecule_pdb (output, comments);

  if (stored_pdb_atom_information.number_elements ())
    stored_pdb_atom_information.resize_keep_storage (0);

  return rc;
}

static void
append_fixed_width_right_justified(const IWString & s, 
                   const int len, 
                   std::ostream & output)
{
  if (s.length() > len)
  {
    for (int i = 0; i < len; ++i)
    {
      output << s[i];
    }
  }
  else if (s.length() == len)
    output << s;
  else
  {
    for (int i = 0; i < (len - s.length()); ++i)
    {
      output << ' ';
    }
    output << s;
  }

  return;
}

static void
append_fixed_width_left_justified (const IWString & s,
                                   const int len,
                                   std::ostream & output)
{
  if (s.length() > len)
  {
    for (int i = 0; i < len; ++i)
    {
      output << s[i];
    }
  }
  else if (s.length() == len)
    output << s;
  else
  {
    output << s;
    for (int i = s.length(); i < len; ++i)
    {
      output << ' ';
    }
  }

  return;
}

static int
write_pdb_atom_common (const Atom * a,
                       const int atom_number,
                       const int fragment_number,
                       const PDB_Stored_Atom_Information * pdbsai,
                       std::ostream & output)
{
  output << "HETATM" << std::setw(5) << atom_number << ' ';

  if (NULL != pdbsai && pdbsai->atom_name().length() > 0)
    append_fixed_width_left_justified(pdbsai->atom_name(), 4, output);
  else
    append_fixed_width_left_justified(a->atomic_symbol(), 4, output);

  if (NULL != pdbsai && pdbsai->residue_name().length() > 0)
    append_fixed_width_right_justified(pdbsai->residue_name(), 4, output);
  else
    output << " UNK";

  output << " A";

  output << std::setw(4) << fragment_number;

  output << "    ";

//output_buffer << " UNK W   " << fragment_number << "    ";

  char buffer[32];

  IW_SPRINTF (buffer, "%8.3f", a->x ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->y ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->z ());
  output << buffer;

  if (NULL != pdbsai && pdbsai->occupancy_factor().length())
    append_fixed_width_right_justified(pdbsai->occupancy_factor(), 6, output);
  else
    output << "  1.00";

  if (NULL != pdbsai && pdbsai->temperature_factor().length())
    append_fixed_width_right_justified(pdbsai->temperature_factor(), 6, output);
  else
    output << " 20.00";

  output << "          ";

  append_fixed_width_right_justified(a->atomic_symbol(), 2, output);

  output << '\n';

  return output.good();
}

static int
write_pdb_atom(const Atom * a,
               int atom_number,
               const PDB_Stored_Atom_Information * pdbsai,
               int fragment_number,
               std::ostream & output)
{
  output << "HETATM" << std::setw(5) << (pdbsai->atom_number()+1) << ' ';

  output << ' ';           // first column of atom name seems not to be used
  if (pdbsai->atom_name().length() > 0)
    append_fixed_width_left_justified(pdbsai->atom_name(), 3, output);
  else
    append_fixed_width_left_justified(a->atomic_symbol(), 3, output);

  if (pdbsai->residue_name().length())
    append_fixed_width_right_justified(pdbsai->residue_name(), 4, output);
  else
    output << " UNK";

  output << " A";

  output << std::setw(4) << fragment_number;

  output << "    ";

//output_buffer << " UNK W   " << fragment_number << "    ";

  char buffer[32];

  IW_SPRINTF (buffer, "%8.3f", a->x ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->y ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->z ());
  output << buffer;

  if (pdbsai->occupancy_factor().length())
    append_fixed_width_right_justified(pdbsai->occupancy_factor(), 6, output);
  else
    output << "  1.00";

  if (pdbsai->temperature_factor().length())
    append_fixed_width_right_justified(pdbsai->temperature_factor(), 6, output);
  else
    output << " 20.00";

  output << "          ";

  append_fixed_width_right_justified(a->atomic_symbol(), 2, output);

  output << '\n';

  return output.good();
}

/*
  Make something like

HETATM    1  C1  UNK     1      -0.399   1.177  -2.204  1.00 20.00
HETATM    2 BR2  UNK     1       1.005   2.092  -3.080  1.00 20.00
HETATM    3  C3  UNK     1      -1.519   1.868  -1.734  1.00 20.00
HETATM   10  C10 UNK     1      -4.371   3.728  -2.158  1.00 20.00
HETATM   30 BR30 UNK     1      19.774  31.230  -0.019  1.00 20.00

*/

static int
write_pdb_atom (const Atom * a,
                int * ecount,
                int atom_number,
                int fragment_number,
                std::ostream & output)
{
  output << "HETATM" << std::setw(5) << atom_number << ' ';

  char buffer[32];

  IWString s = a->atomic_symbol ();

  if (2 == s.length ())
    s.to_uppercase ();

  IWString number_to_append;

  if (NULL == ecount)
    number_to_append << atom_number;
  else
  {
    atomic_number_t z = a->atomic_number();
    if (z >= 0)
    {
      ecount[z]++;
      number_to_append << ecount[z];
    }
    else
      number_to_append << '1';               // an arbitrary number
  }

// overall length must be 4, but by convention the first column seems to be not used

  output << ' ';
  if (s.length() + number_to_append.length() <= 3)
    append_fixed_width_left_justified(s, 3, output);
  else
  {
    s << '*';
    append_fixed_width_left_justified(s, 3, output);
  }

  output << " UNK A ";
  output << std::setw(3) << fragment_number;

  output << "    ";

  IW_SPRINTF (buffer, "%8.3f", a->x ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->y ());
  output << buffer;
  IW_SPRINTF (buffer, "%8.3f", a->z ());
  output << buffer;
  output << "  1.00 20.00";

  output << "          ";

  append_fixed_width_right_justified(a->atomic_symbol(), 2, output);

  output << '\n';

  return output.good ();
}

/*
  Some programmes insist on members of each fragment being
  written contiguously.
*/

int
Molecule::write_connection_table_pdb (std::ostream & os)
{
  assert (ok ());
  assert (os.good ());

  int matoms = natoms ();

  if (use_stored_atom_information_when_writing_pdb_files && stored_pdb_atom_information.number_elements() != matoms)
  {
    cerr << "Molecule::write_connection_table_pdb:stored name size mismatch, turned off\n";
    cerr << "Molecule has " << matoms << " atoms, stored names " << stored_pdb_atom_information.number_elements() << " items\n";
    use_stored_atom_information_when_writing_pdb_files = 0;
  }

//fmtflags ff = os.flags ();   // save state of os

  int * element_count = NULL;

  if (number_by_element_count)
    element_count = new_int(HIGHEST_ATOMIC_NUMBER + 1);

  std::unique_ptr<int[]> free_element_count(element_count);

  if (write_pdb_files_in_fragment_order)
  {
    int nf = number_fragments();
    for (int f = 0; f < nf; f++)
    {
      int atoms_this_fragment = 0;
      for (int j = 0; j < _number_elements; j++)
      {
        if (f != fragment_membership(j))
          continue;

        atoms_this_fragment++;

        if (use_stored_atom_information_when_writing_pdb_files)
          write_pdb_atom(_things[j], j + 1, stored_pdb_atom_information[j], f + 1, os);
        else if (NULL != element_count)
          write_pdb_atom(_things[j], element_count, (j + 1), f + 1, os);
        else if (number_within_sequence)
          write_pdb_atom(_things[j], NULL, atoms_this_fragment, f+ 1, os);
        else
          write_pdb_atom(_things[j], NULL, (j + 1), f+ 1, os);
      }
    }
  }
  else
  {
    for (int i = 0; i < matoms; i++)
    {
      int f;
      if (0 == _bond_list.number_elements())
        f = 1;
      else
        f = fragment_membership(i) + 1;

      if (use_stored_atom_information_when_writing_pdb_files)
        write_pdb_atom(_things[i], i + 1, stored_pdb_atom_information[i], f, os);
      else
        write_pdb_atom(_things[i], element_count, (i + 1), fragment_membership(i) + 1, os);
    }
  }

  for (int i = 0; i < matoms;i++)
  {
    const Atom * ai = _things[i];

    int ncon = ai->ncon();

    if (ncon > 0)
    {
      os << "CONECT" << std::setw (5) << i + 1;
      for (int j = 0; j < ncon; j++)
      {
        atom_number_t k = ai->other (i, j);
        os << std::setw (5) << k + 1;
      }
      os << '\n';
    }
  }

//os.setf (ff);

  return 1;
}

int
Molecule::write_molecule_pdb (std::ostream & os, const IWString & comments)
{
  assert (ok ());
  assert (os.good ());

  static const char * month [] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

  time_t t = time(NULL);
  struct tm * tm = gmtime(&t);

  os << "HEADER    UNK                                     " << tm->tm_mday << '-' << month[tm->tm_mon];
  os <<  '-' << std::setw(2) << (tm->tm_year - 100) << "  1UNK\n";
  os << "COMPND    " << name () << '\n';
  if (comments.length ())
    os << "REMARK    " << comments << endl;
  os << "REMARK   1 fileconv:" << __FILE__ << " compiled " << __DATE__ << " at " << __TIME__ << '\n';

  int rc = write_connection_table_pdb (os);

  os << "END\n";

  return rc;
}

// arch-tag: 40352ad4-d9bb-430c-adc6-52c128673138
