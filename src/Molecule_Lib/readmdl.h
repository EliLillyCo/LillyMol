#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <fstream>
#include <iomanip>
#include <memory>

// Be sure to define this symbol so all the private functions get defined

#define COMPILING_MDL_CC
#define COMPILING_CTB

#include "misc.h"
#include "iwcrex.h"

#include "mdl.h"
#include "atom_alias.h"
#include "molecule.h"
#include "misc2.h"
#include "chiral_centre.h"
#include "rwmolecule.h"
#include "aromatic.h"
#include "mdl_atom_record.h"

/*
  MDL files are somewhat strange in that there several varieties.
  Sometimes the connection table is terminated by 'M  END', sometimes by $$$$.
  Sometimes both are present.
  The way this works now is that by default, we insist on the $$$$.
  If return_on_m_end is set, then we return on 'M  END'
*/

template <typename T>
int
Molecule::_read_molecule_mdl_ds (T & input,
                                int return_on_m_end,
                                MDL_File_Supporting_Material & mdlfsm) 
{
  int nb = 0;
  int v30;
  if (! _read_mdl_atom_connection_table(input, nb, v30, mdlfsm))
  {
    skip_to_string(input, "$$$$", 1);

    return 0;
  }

  if(_bond_list.elements_allocated() < nb)
    _bond_list.resize(nb);

// Aromatic atoms are those found at the ends of aromatic bonds. Not ideal, but seems to work

  int * aromatic_atoms = NULL;
  int * aromatic_bonds = NULL;

  if (_number_elements > 0 && input_aromatic_structures())
  {
    aromatic_atoms = new_int(_number_elements);
    if (nb > 0)
      aromatic_bonds = new_int(nb);
  }

  int rc;
  if(v30)
    rc = _read_v30_bond_list(input, nb, aromatic_atoms, aromatic_bonds);
  else
    rc = _read_mdl_bond_list(input, nb, aromatic_atoms, aromatic_bonds, mdlfsm);

  if(v30)
  {
    if (! _process_v30_composite_records(input, mdlfsm))
    {
      cerr << "_read_molecule_mdl_ds:cannot process V30 composite records\n";
      return 0;
    }
    skip_to_string(input, "M  END", 0);    // 0 means do it quietly
  }

  if(unconnect_covalently_bonded_non_organics_on_read())
    _do_unconnect_covalently_bonded_non_organics();

  if (NULL != aromatic_atoms)
  {
    if (! _final_processing_of_aromatic_mdl_input(aromatic_atoms, aromatic_bonds))
      rc = 0;

    delete [] aromatic_atoms;
    delete [] aromatic_bonds;
  }

  if (v30 && return_on_m_end)    // possibly reading an rdf file
    return rc;

  if (! _read_molecule_mdl_trailing_records(input, return_on_m_end, mdlfsm))
    return rwmolecule_error("read_molecule_mdl_ds: bad stuff at end", input);

  return rc;
}

/*
  Main function for reading mdl files
*/

template <typename T>
int
Molecule::read_molecule_mdl_ds (T & input,
                                int return_on_m_end)
{
  assert(ok());

  resize(0);

  assert(input.good());
  if(input.eof())
    return 0;

  MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  mdlfos->reset_for_next_molecule();

  if (! _read_molecule_mdl_ds(input, return_on_m_end, *mdlfos))
    return 0;

// Clean up chirality and cis-trans stuff
// Aug 2001.  If there are no chiral centres, but wedge bonds present,
// use that data - files from Afferent are like that!  Note that this
// isn't robust.  We could have some centres that are not atom marked,
// but have wedge bonds

  if(ignore_all_chiral_information_on_input())
    _chiral_centres.resize(0);
  else if(mdlfos->discern_chirality_from_wedge_bonds())
   (void) discern_chirality_from_wedge_bonds();
  else if (0 == _chiral_centres.number_elements() && number_up_or_down_wedge_bonds())
     (void) discern_chirality_from_wedge_bonds();
  else if (discern_chirality_from_3d_coordinates() && 3 == highest_coordinate_dimensionality())
  {
    int d = discern_chirality_from_3d_coordinates();
    if (1 == d)
      (void) discern_chirality_from_3d_structure();
    else if (_number_elements <= d)
      (void) discern_chirality_from_3d_structure();
    else
      cerr << "Molecule::_read_molecule_mdl_ds:skipped d@3d for too many atoms '" << name() << "' " << _number_elements << '\n';
  }

  if (0 == _chiral_centres.number_elements())    // none to worry about
    ;
  else if(_complete_chiral_centres_from_mdl_files(*mdlfos))    // good
    ;
  else               // OOPS, bad chirality info
  {
    cerr << "Molecule::read_molecule_mdl_ds: erroneous chiral input\n";
    cerr << _molecule_name << endl;
    _chiral_centres.resize(0);
    if (! ignore_incorrect_chiral_input())
      return 0;
  }

// Do cis-trans bonds last because they depend on ring membership. Ran into a case, 583770, 
// where the ring perception forced a smiles ordering that was wrong because chirality hadn't
// been perceived

  if(discern_cis_trans_bonds())
   (void) discern_cis_trans_bonds_from_depiction();

  return 1;
}

template <typename T>
int
Molecule::_read_mdl_atom_connection_table (T & input, 
                                           int & nb,
                                           int & v30,
                                           MDL_File_Supporting_Material & mdlfos)
{
//  There are three header lines.
//  The line following will contain na and nb

// Note that the record with na and nb is being read here too!

  nb = 0;
  v30 = 0;

  const_IWSubstring buffer;
  for (int i = 0; i < 4; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");
    if (0 == i && ! mdlfos.discard_sdf_molecule_name())
      set_name(buffer);
  }

// buffer should now hold na and nb

  if (buffer.length() < 6)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table: the atoms/bond record must be at least 6 chars long, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  int na;
  if (2 != int3d(buffer, na, nb))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table: error from int3d '" << buffer << "'\n";
    cerr << "Line " << input.lines_read() << endl;
    return 0;
  }

//cerr << "Contains '" << na << " atoms and " << nb << " bonds\n";

  assert (na >= 0 && (nb >= 0));

  if (_elements_allocated < na)
    resize(na);

  if (buffer.contains(" V3000"))
  {
    v30 = 1;
    const auto isosave = mdlfos.read_isotopes_as_numbers_rather_than_differences_from_normal();
    mdlfos.set_read_isotopes_as_numbers_rather_than_differences_from_normal(1);
    const auto rc = _read_mdl_atom_connection_table_v30(input, nb, mdlfos);
    mdlfos.set_read_isotopes_as_numbers_rather_than_differences_from_normal(isosave);
    return rc;
  }

  MDL_Atom_Record mdl_atom_record;

  for (int i = 0; i < na; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

    if (! mdl_atom_record.build(buffer))
    {
      cerr << buffer << endl;
      return rwmolecule_error("read_molecule_mdl_ds:bad atom data", input);
    }

    Atom * a = mdl_atom_record.create_atom();

    if (NULL == a)
    {
      cerr << buffer << endl;
      cerr << rwmolecule_error("read_molecule_mdl_ds:bad element", input);
      return 0;
    }

    add(a);

    if (0 == mdl_atom_record.astere())     // the most common case, no chirality
      ;
    else if (! _mdl_atom_is_chiral_centre(_number_elements - 1, mdl_atom_record.astere(), mdlfos))
    {
      cerr << "Molecule::_read_mdl_atom_connection_table: invalid chirality on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

/*
  We always attempt to read NB records.
*/

template <typename T>
int
Molecule::_read_mdl_bond_list (T & input, int nb,
                               int * aromatic_atoms,
                               int * aromatic_bonds,
                               MDL_File_Supporting_Material & mdlfos)
{
  assert (nb >= 0);

  int na = _number_elements;

  int rc = 1;    // assume ok until proven otherwise

  const_IWSubstring buffer;
  for (int i = 0; i < nb; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

    if (0 == rc)     // just read the records if we have already failed
      continue;

    int a1, a2;
    int bond_type_read_in;
    int directionality;

    if (! mdlfos.parse_bond_record(buffer, na, a1, a2, bond_type_read_in, directionality))
    {
      cerr << "Molecule::_read_mdl_bond_list: bond record " << i << " is bad '" <<
              buffer << "'\n";
      rc = 0;
      continue;
    }

    bond_type_read_in = mdlfos.translate_input_bond_type(bond_type_read_in);

    if (mdlfos.ignore_self_bonds() && (a1 == a2))
    {
      cerr << "Molecule::_read_mdl_bond_list: ignoring self bonds " << a1 << " to " << a2 << endl;
      continue;
    }

    bond_type_t btype = INVALID_BOND_TYPE;

    if (4 == bond_type_read_in)
    {
      if (! input_aromatic_structures())
      {
        cerr << "Molecule::_read_mdl_bond_list:aromatic input not enabled, bond between atoms " << a1 << " and " << a2 << endl;
        rc = 0;
        continue;
      }

      btype = SINGLE_BOND;
      aromatic_atoms[a1] = 1;
      aromatic_atoms[a2] = 1;
      aromatic_bonds[i] = 1;
    }
    else if (! convert_from_mdl_number_to_bond_type(bond_type_read_in, btype))
    {
      cerr << "Molecule::_read_mdl_bond_list: bad bond type " << bond_type_read_in << endl;
      rc = 0;
      continue;
    }

    add_bond(a1, a2, btype, 1);     // 1 means partially built molecule
    if (directionality)
      _mdl_set_bond_directionality(a1, a2, directionality);
  }

  if (rc && nb > 0)
  {
    check_bonding();
  }

  return rc;
}

/*
  By convention, data following a tag may be multi-record, but will be
  terminated by a blank line
*/

extern int looks_like_sdf_tag (const const_IWSubstring & buffer);

template <typename T>
int
Molecule::_read_mdl_data_following_tag (T & input,
                                        const MDL_File_Supporting_Material & mdlfos)
{
  IWString buffer;

  for (int i = 0; input.next_record(buffer); i++)
  {
    if (looks_like_sdf_tag(buffer))
    {
      input.push_record();
      return 1;
    }

    if (' ' != mdlfos.gsub_mdl_file_data())
      buffer.gsub(' ', mdlfos.gsub_mdl_file_data());

    if (read_extra_text_info())   // even if buffer is empty
      _text_info.add(new IWString(buffer));

    if (0 == buffer.length() || (1 == buffer.length() && 13 == static_cast<int>(buffer[0])))    // 13 is ^M
      return 1;

    if (0 == i && mdlfos.prepend_sdfid())   // on first record, no space before extra info
      _molecule_name += buffer;
    else
      _molecule_name.append_with_spacer(buffer, mdlfos.insert_between_sdf_name_tokens());
  }

  cerr << "Molecule::_read_mdl_data_following_tag:premature eof\n";
  return 0;
}

extern int parse_m_sty_record (const const_IWSubstring & buffer, int & sgroups_present);

/*
  Function for reading all the "stuff" which comes between the bond
  list and the end of the molecule
*/

template <typename T>
int
Molecule::_read_molecule_mdl_trailing_records (T & input,
                                               int return_on_m_end,
                                               MDL_File_Supporting_Material & mdlfos)   // non const because of rx match
{
  if (read_extra_text_info() && 0 == _text_info.number_elements())
    _text_info.resize(10);

//Aprop atom_properties[MAX_PAIRS];

// We just read and ignore the S group information

  int sgroups_present = 0;
  int msal = 0;
  int msbl = 0;
  int msmt = 0;
  int msbv = 0;

  IWString buffer;
  int trailing_lines = 0;
  int extreg_found = 0;
  int got_dollar = 0;

  while (input.next_record(buffer))
  {
    trailing_lines++;

    buffer.strip_trailing_blanks();      // do this first, get rid of any garbage out the end

    if ("$$$$" == buffer)
    {
      got_dollar = 1;
      break;
    }

    int fatal;
    if (_common_parse_M_record(buffer, fatal))   // great, recognised and good
      continue;
    else if (fatal)
    {
      cerr << "Molecule::_read_molecule_mdl_trailing_records:invalid record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    if (buffer.starts_with("M  STY"))
    {
      if (! parse_m_sty_record(buffer, sgroups_present))
      {
        if(mdlfos.die_on_erroneous_m_input())
          return 0;
      }

      continue;
    }

    if (buffer.starts_with("M  SAL"))
    {
      msal--;
      continue;
    }

    if (buffer.starts_with("M SBL"))
    {
      msbl--;
      continue;
    }

    if (buffer.starts_with("M  SBV"))
    {
      msbv--;
      continue;
    }

    if (buffer.starts_with("M  SMT"))
    {
      msmt--;
      continue;
    }

    if (buffer.starts_with("M  S"))
      continue;

    if ("M  END" == buffer)
    {
      if (return_on_m_end)
        return 1;
      continue;
    }

    if (0 == extreg_found && mdlfos.name_in_m_tag().length() && buffer.nwords() > 2 &&
        buffer.starts_with("M ") && buffer.contains(mdlfos.name_in_m_tag()))
    {
      _molecule_name = buffer;
      _molecule_name.remove_leading_words(2);    // M REG
      extreg_found = 1;
      continue;
    }

    if (buffer.starts_with("M  "))
    {
      if (mdlfos.report_unrecognised_records() || ! mdlfos.ignore_unrecognised_m_records())
        cerr << "Unrecognised 'M  ' directive '" << buffer << "'\n";
      if (! mdlfos.ignore_unrecognised_m_records())
        return 0;
    }

    if (buffer.starts_with("A  "))
    {
      Atom_Alias * a = new Atom_Alias;
      if (! a->build(buffer, input))
      {
        cerr << "Invalid atom alias data, line " << input.lines_read() << endl;
        delete a;
        return 0;
      }

      mdlfos.add_alias(a);

     (void) input.next_record(buffer);

      if(buffer.length() > 0)
        input.push_record();

      continue;
    }

    if (buffer.starts_with("G  "))
    {
      mdlfos.extra_g_record_found();

      IWString g(buffer);   // needs two records, make a copy of first one
      (void) input.next_record(buffer);

      if (mdlfos.mdl_g_records_hold_atom_symbols())
      {
        if (! _process_mdl_g_record (g, buffer))
        {
          cerr << "Molecule::_read_molecule_mdl_trailing_records:cannot process G record '" << buffer << "'\n";
          return 0;
        }
      }

      continue;
    }

    if (read_extra_text_info())
      add_to_text_info(_text_info, buffer);

//  cerr << read_extra_text_info() << " now contains " << _text_info.number_elements() << endl;

    if (0 == buffer.length())
      continue;

//  Now all the various other identifiers possible in the file

    if (mdlfos.sdf_identifier_matches(buffer))
    {
      IWString id;
      extract_sdf_identifier(buffer, id);

      EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

      if(read_extra_text_info())
        add_to_text_info(_text_info, buffer);

      IWString tmp;
      if (mdlfos.replace_first_sdf_tag().length() > 0)
        tmp << mdlfos.replace_first_sdf_tag() << ':';
      else if(mdlfos.prepend_sdfid())
        tmp << id << ':';

      tmp += buffer;

      _molecule_name.append_with_spacer(tmp, mdlfos.insert_between_sdf_name_tokens());

      continue;
    }

//  Are we fetching all the SDF data from the file.
//  We run into lots of different formats.
//  >  <IDENTIFIER> 
//  data
//
// >  <IDENTIFIER> stuff
// data

// We need to fetch IDENTIFIER

    if (mdlfos.fetch_all_sdf_identifiers() && looks_like_sdf_tag(buffer))
    {
      IWString id;
      if (! extract_sdf_identifier(buffer, id))
        return rwmolecule_error("read_molecule_mdl_ds: cannor parse SDF identifier", input);

      if (mdlfos.replace_first_sdf_tag().length() > 0)
      {
        _molecule_name.append_with_spacer(mdlfos.replace_first_sdf_tag(), mdlfos.insert_between_sdf_name_tokens());
        _molecule_name.add(':');
      }
      else if (mdlfos.prepend_sdfid())
      {
        _molecule_name.append_with_spacer(id, mdlfos.insert_between_sdf_name_tokens());
        _molecule_name.add(':');
      }

      if (! _read_mdl_data_following_tag(input, mdlfos))
        return 0;

      continue;
    }

    if (0 == _molecule_name.length() && mdlfos.take_first_tag_as_name() &&
        looks_like_sdf_tag(buffer))
    {
      if (! _read_mdl_data_following_tag(input, mdlfos))
        return 0;

      continue;
    }

//  If we are storing the extra info, do so, otherwise silently ignore other
//  kinds of records.

    if (0 == extreg_found && mdlfos.extract_isis_extregno() &&
        buffer.starts_with(">  <") && 3 == buffer.nwords())
    {
      IWString extreg;
      buffer.word(2, extreg);
      if (extreg.starts_with('(') && extreg.ends_with(')'))
      {
        extreg.remove_leading_chars(1);
        extreg.chop(1);

        if (0 == _molecule_name.length())
          _molecule_name = extreg;
        else
        {
          extreg += ' ';     // allow a space before the current name
          _molecule_name.insert(extreg, 0);
        }

        extreg_found = 1;
        continue;
      }
    }

//  If we get to here, just ignore it.

    if(mdlfos.report_unrecognised_records())
      cerr << "Ignoring unrecognised form, line " << input.lines_read() << " '" << buffer << "'\n";
  }

  if (mdlfos.set_elements_based_on_atom_aliases() && mdlfos.number_aliases())
    _set_elements_based_on_atom_aliases(mdlfos.atom_aliases());

  if(got_dollar)
    return 1;

// If we come out here, it must be EOF. That's OK.

//if (trailing_lines > 6)     May 2005, this warning doesn't seem necessary
//  cerr << "mdl_read_ds: " << trailing_lines << " lines found between bonds and EOF\n";

  cerr << "mdl_read_ds returning at EOF without $$$$\n";

  return 1;
}

template <typename T>
int
skip_to_rdfile_start_of_record (T & input,
                                const IWString & rdfile_start_of_record)
{
  const_IWSubstring buffer;
  
  int records_read_here = 0;

  while (input.next_record(buffer))
  {
    records_read_here++;

    if (buffer.starts_with(rdfile_start_of_record))
    {
      input.push_record();
      return 1;;
    }
  }

  if (0 == records_read_here)
  {
    cerr << "read mol rdf eof\n";
    return 0;
  }

  cerr << "EOF reading RDFILE, cannot find start record '" << rdfile_start_of_record << "', tried " << records_read_here << "\n";
  return 0;
}

template <typename T>
int
Molecule::write_molecule_mdl (T & os,
                              const IWString & comments) const
{
  assert(ok());
  assert(os.good());

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  if (_number_elements > 999 || mdlfos->write_v30_mdl_files() || _contains_isotope_above (999))
    return write_molecule_mdl_v30(os, comments, 1);

  os << _molecule_name << newline_string();

  if(mdlfos->isis_standard_records())
  {
    int dim;
    if(highest_coordinate_dimensionality() > 2)
      dim = 3;
    else 
      dim = 2;

    os << "  -ISIS-  0516971354" << dim << "D 1   1.00000     0.00000     1" << newline_string();
    os << newline_string();
  }
  else
  {
    if(comments.length())
      os << comments << newline_string();
    else
      os << "Blank" << newline_string();

    os << "Blank" << newline_string();
  }

  int rc = write_connection_table_mdl(os);

//if(isis_standard_records)
//  os << "M  END\n";

  if (::write_extra_text_info())
    write_extra_text_info(os);

  if (mdlfos->write_mdl_dollars())
    os << "$$$$" << newline_string();

  if (flush_files_after_writing_each_molecule())
    os.flush();

  return rc;
}

template <typename T>
int
Molecule::write_extra_text_info (T & os) const
{
  int ne = _text_info.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const IWString * info = _text_info[i];
    os <<(*info) << newline_string();
  }

  return os.good();
}

/*
  Writing out the atoms and bonds record of an MDL file is complex
*/

template <typename T>
int
Molecule::_mdl_write_atoms_and_bonds_record (T & os,
                                             int nfc,
                                             int iat,
                                             MDL_File_Supporting_Material & mdlfos) const
{
  int write_stereo_info = mdlfos.write_mdl_chiral_flags();

  if(write_stereo_info)
    write_stereo_info = _chiral_centres.number_elements();

  int nb = _bond_list.number_elements();

  mdlfos.write_atoms_and_bonds (_number_elements, nb, os);

  if (mdlfos.isis_standard_records())
  {
    os << "  0";     // number of atoms lists
    os << "   ";     // obsolete
    if (write_stereo_info && _chiral_centres.number_elements())
      os << "  1";
    else
      os << "  0";
    os << "  0";     // number of stext entries
    os << "   ";     // reaction components + 1
    os << "   ";     // number of reactants
    os << "   ";     // number of products
    os << "   ";     // number of intermediates

//  Work out the number of M lines

    int mmm = 0;

    if(nfc)
    {
      if (0 == nfc % 8)
        mmm = nfc / 8;
      else
        mmm = nfc / 8 + 1;
    }

    if(iat)
    {
      if (0 == iat % 8)
        mmm += iat / 8;
      else
        mmm += iat / 8 + 1;
    }

    mmm++;    // don't forget the 'M  END' line

    os << mdlfos.digits3(mmm);

    os << " V2000";  // Ctab version
  }

  os << newline_string();

  return os.good();
}

/*
  This function just writes out the connection table.
*/

/*
  Variable write_stereo_info is really a counter of the number
  of remaining stereo centres to write.
*/

template <typename T>
int
Molecule::write_connection_table_mdl (T & os) const
{
  assert(ok());
  assert(os.good());

  MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  int nfc = number_formally_charged_atoms();
  int iat = number_isotopic_atoms();

  _mdl_write_atoms_and_bonds_record(os, nfc, iat, *mdlfos);

  int nb = _bond_list.number_elements();
  int width = os.width();
  IWString output_buffer;
  output_buffer.resize(80);

  int number_r_groups_present = 0;

  for (int i = 0; i < _number_elements && os.good(); i++)
  {
    const Atom * a = _things[i];

    a->write_coordinates(os);

    _write_mdl_atom_record_element_charge_and_chirality(i, output_buffer, *mdlfos);

    if(mdlfos->isis_standard_records())
      output_buffer += "  0  0  0           0  0  0";

    output_buffer += newline_string();

    os << output_buffer;

    if (! a->element()->is_in_periodic_table())    // is this an R# group
    {
      const IWString & s = a->element()->symbol();
      if ('R' == s[0] && 2 == s.length() && isdigit(s[1]))
        number_r_groups_present++;
    }
  }

  for (int i = 0; i < nb && os.good(); i++)
  {
    const Bond *b = bondi(i);
    os << mdlfos->digits3(b->a1() + 1);
    os << mdlfos->digits3(b->a2() + 1);

    if (mdlfos->mdl_write_aromatic_bonds() && (b->is_aromatic() || b->is_permanent_aromatic()))
      os << mdlfos->digits3(4);
    else if(b->is_single_bond())
      os << mdlfos->digits3(1);
    else if(b->is_double_bond())
      os << mdlfos->digits3(2);
    else if(b->is_triple_bond())
      os << mdlfos->digits3(3);
    else if(b->is_aromatic())
      os << mdlfos->digits3(4);
    else if (COORDINATION_BOND == b->btype())
      os << mdlfos->digits3(9);
    else
      os << mdlfos->digits3(1);       // defaults to single

    int directionality_written = 1;

    if(b->is_wedge_up())
      os << mdlfos->digits3(1);
    else if(b->is_wedge_down())
      os << mdlfos->digits3(6);
    else if(b->is_wedge_either())
      os << mdlfos->digits3(4);
    else if(b->is_cis_trans_either_double_bond())
      os << mdlfos->digits3(3);
    else
      directionality_written = 0;

    if(mdlfos->isis_standard_records())
    {
      if (! directionality_written)
        os << "  0";

      os << "     0  0";
    }

    os << newline_string();
  }

  os.width(width);

  int need_m_end = 0;

  if (nfc && mdlfos->write_mdl_charges_as_m_chg())
  {
    _write_m_chg_records(os, nfc);
    need_m_end = 1;
  }

  if(iat)
  {
    _write_m_iso_records(os, iat);
    need_m_end = 1;
  }

  if(number_r_groups_present)
    _write_M_RGP_records(*mdlfos, os);

  if (mdlfos->isis_standard_records() || (need_m_end && mdlfos->write_mdl_m_end_record()) || mdlfos->write_mdl_m_end_record() > 1)
    os << "M  END" << newline_string();

  if (! os.good())
  {
    cerr << "Molecule::write_connection_table_mdl: cannot write\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
Molecule::write_molecule_mdl_v30 (T & os, 
                                  const IWString & comments,
                                  int write_end_stuff) const
{
  MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  os << _molecule_name << '\n';

  os << "IWmolecule05289912212D 1   0.00000     0.00000      0\n";
  os << comments << '\n';
  os << "  0  0  0     0  0              0 V3000\n";
  os << "M  V30 BEGIN CTAB\n";
  os << "M  V30 COUNTS " << _number_elements << ' ' << _bond_list.number_elements() << " 0 0 " << (_chiral_centres.number_elements() > 0) << '\n';

  _write_molecule_atom_list_v30 (os, *mdlfos);
  _write_molecule_bond_list_v30 (os);

  if (write_end_stuff)
  {
    os << "M  V30 END CTAB\n";
    os << "M  END\n";
  }

  if (mdlfos->write_mdl_dollars())
    os << "$$$$\n";

  return os.good();
}

template <typename T>
int
Molecule::_write_molecule_bond_list_v30 (T & os) const
{
  os << "M  V30 BEGIN BOND\n";

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];
    os << "M  V30 " << (i + 1);

    if (b->is_single_bond())
      os << " 1 ";
    else if (b->is_double_bond())
      os << " 2 ";
    else if (b->is_triple_bond())
      os << " 3 ";
    else if (b->is_aromatic())
      os << " 4 ";
    else
    {
      cerr << "Molecule::_write_molecule_bond_list_v30: what kind of bond is this " << b->btype() << '\n';
      os << " ? ";
    }

    os << (b->a1() + 1) << ' ' << (b->a2() + 1);

    if (b->is_wedge_up())
      os << " CFG=1";
    else if (b->is_wedge_down())
      os << " CFG=3";
    else if (b->is_wedge_either())
      os << " CFG=2";

    os << '\n';
  }

  os << "M  V30 END BOND\n";

  return os.good();
}

template <typename T>
int
Molecule::_write_molecule_atom_list_v30 (T & os,
                                         MDL_File_Supporting_Material & mdlfos) const
{
  const Set_of_Atoms & unspecified = mdlfos.mdl_unspecified_chiral_atoms();

  os << "M  V30 BEGIN ATOM\n";

  int nc = _chiral_centres.number_elements();

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    os << "M  V30 " << (i + 1) << ' ' << a->atomic_symbol() << ' ' << a->x() << ' ' << a->y() << ' ' << a->z() << " 0";

    if (a->formal_charge())
    {
      os << " CHG=" << a->formal_charge();
    }

    if (a->is_isotope())
    {
      os << " MASS=" << a->isotope();
    }

    if (nc)    // still some chiral centres to process
    {
      for (int j = 0; j < _chiral_centres.number_elements(); j++)
      {
        const Chiral_Centre * cc = _chiral_centres[j];
        if (i == cc->a())
        {
          os << " CFG=" << cc->mdl_stereo_centre_value();
          nc--;     // one less chiral centre to process
          break;
        }
      }
    }

    if (unspecified.contains (i))
      os << " CFG=3";

    os << '\n';
  }

  os << "M  V30 END ATOM\n";

  return os.good();
}


/*
  Write M CHG record(s)

  There can be a maximum of 8 charge entries per record.
*/

#define M_CHG_PER_RECORD 8

template <typename T>
int
Molecule::_write_m_chg_records (T & os, int nc) const
{
  assert (nc > 0);

  int j = 0;    // index of atom being processed

  while(nc)
  {
    os << "M  CHG ";

    int items_this_record;
    if (nc > M_CHG_PER_RECORD)
      items_this_record = M_CHG_PER_RECORD;
    else
      items_this_record = nc;

    nc -= items_this_record;

    os << std::setw(2) << items_this_record;
    int items_written = 0;

    while (items_written < items_this_record)
    {
      formal_charge_t q = _things[j++]->formal_charge();
      if(q)
      {
        os << std::setw(4) << j << std::setw(4) << q;
        items_written++;
      }
    }

    os << newline_string();
  }

  return os.good();
}

template <typename T>
int
Molecule::_write_m_iso_records (T & os, int n) const
{
  assert (n > 0);

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
   
  int j = 0;        // index of atom being processed

  while(n)
  {
    os << "M  ISO ";

    int items_this_record;
    if (n > M_CHG_PER_RECORD)
      items_this_record = M_CHG_PER_RECORD;
    else
      items_this_record = n;

    n -= items_this_record;

    os << std::setw(2) << items_this_record;

    int items_written = 0;
    while (items_written < items_this_record)
    {
      const Atom * a = _things[j++];

      if (0 == a->isotope())
        continue;

      int to_be_written;
      if (mdlfos->write_M_isotopes_as_numbers_rather_than_differences_from_normal())
        to_be_written = a->isotope();
      else
        to_be_written = a->isotope() - a->element()->normal_isotope();

      os << std::setw(4) << j << std::setw(4) << to_be_written;

      items_written++;
    }

    os << newline_string();
  }

  return os.good();
}

template <typename T>
int
fetch_collection(T & input,
                 const IWString & collection_type,
                 resizable_array_p<IWString> & mgroup)
{
  const_IWSubstring line;

  IWString end_group = "END ";
  end_group << collection_type;

  IWString s;    // buffer in case of continuation lines

  int got_end_group = 0;

  while (input.next_record(line))
  {
    if (! line.starts_with("M  V30 "))
    {
      cerr << "fetch_collection:possibly invalid record prefix '" << line << "', continuing....\n";
      continue;
    }

    line.remove_leading_chars(7);

    if (line.ends_with(" -") || (79 == line.length() && line.ends_with('-')))
    {
      line.chop(1);
      s += line;
      continue;
    }

    if (line.starts_with(end_group))
    {
      got_end_group = 1;
      break;
    }

    s += line;      // S may be empty

    mgroup.add(new IWString(s));
    s.resize_keep_storage(0);
  }

  if (s.length())   // should not happen
    mgroup.add(new IWString(s));

  if (! got_end_group)
  {
    cerr << "fetch_collection:premature EOF. Did not find '" << end_group << "'\n";
  }

  return mgroup.number_elements();
}

template <typename T>
int
Molecule::_process_v30_composite_records(T & input,
                                         const MDL_File_Supporting_Material & mdlfsm)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with("M  V30 END CTAB"))
      break;

    if (buffer.starts_with("M  V30 BEGIN"))
    {
      resizable_array_p<IWString> s;
      IWString collection_type;

      collection_type = buffer;
      collection_type.remove_leading_chars(13);
      if (! fetch_collection(input, collection_type, s))
      {
        cerr << "Molecule::_process_v30_composite_records:cannot fetch '" << collection_type << " collection end\n";
        return 0;
      }

//    cerr << "Got group type " << collection_type << " with " << s.number_elements() << " lines\n";
      if ("SGROUP" == collection_type && mdlfsm.convert_single_atom_sgroup_to_element())
      {
        if (! _convert_sgroups_to_elements(s))
        {
          cerr << "Molecule::_process_v30_composite_records:cannot process SGROUP\n";
          for (int i = 0; i < s.number_elements(); ++i)
          {
            cerr << *s[i] << '\n';
          }
          return 0;
        }
      }
    }
  }

  return 1;
}
