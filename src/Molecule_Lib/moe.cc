// Reads MOE files

#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"

namespace moe {

using std::cerr;

// Moe files are unusual, in that they likely contain many molecules, probably a
// receptor and several ligands. The residue information will contain the
// per-atom name information, but we have no means of storing that, in scope,
// and returning it with the molecule.
// To make things easier, we do this horrible hacky thing and set up a file
// to which we can write the Moe fragments as they are read.

IWString_and_File_Descriptor stream_for_moe_fragments;

int
SetupMoeFragmentStream(const_IWSubstring fname) {
  IWString myfname(fname);
  if (! myfname.ends_with(".smi")) {
    myfname << ".smi";
  }
  if (! stream_for_moe_fragments.open(myfname.null_terminated_chars())) {
    cerr << "SetupMoeFragmentStream:cannot open '" << myfname << "'\n";
    return 0;
  }

  return 1;
}

}  // namespace moe

// many of these static things could be added to the moe namespace.
// do that if we ever create unit tests.

using std::cerr;

// Get the next token from `buffer`. If that fails, read the next
// record from `input` into `buffer` and extract the next word from
// there. `i` is the character pointer used by buffer.nextword.
static int
nextword(iwstring_data_source& input,
         const_IWSubstring& buffer,
         int& i,
         const_IWSubstring& result) {
  if (buffer.nextword(result, i)) {
    return 1;
  }
  if (! input.next_record(buffer)) {
    cerr << "moe:nextword:cannot get next record\n";
    return 0;
  }

  i = 0;
  return buffer.nextword(result, i);
}

struct Residue {
  int id;
  int natoms;
  IWString name;
};

int
GatherResidueData(iwstring_data_source& input,
                  const_IWSubstring& header,
                  std::vector<Residue>& results) {
  const_IWSubstring token;
  int i = 0;
  if (! header.nextword(token, i) || 
      ! header.nextword(token, i)) {
    cerr << "GatherResidueData::cannot extract tokens from header '" << header << "'\n";
    return 0;
  }
  int number_values;
  if (! token.numeric_value(number_values)) {
    cerr << "GatherResidueData::invalid number of residues '" << header << "'\n";
    return 0;
  }

  results.resize(number_values);

  const_IWSubstring buffer;
  i = 0;
  for (int j = 0; j < number_values; ++j) {
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueData:eof\n";
      return 0;
    }
    if (! token.numeric_value(results[j].id)) {
      cerr << "GatherResidueData:invalid id '" << token << "'\n";
      return 0;
    }

    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueData:eof\n";
      return 0;
    }
    if (! token.numeric_value(results[j].natoms)) {
      cerr << "GatherResidueData:invalid natoms '" << token << "'\n";
      return 0;
    }

    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueData:eof\n";
      return 0;
    }
    results[j].name = token;

    for (int notused = 0;  notused < 4; ++notused) {
      if (! nextword(input, buffer, i, token)) {
        cerr << "GatherResidueData:eof\n";
        return 0;
      }
    }
  }

  return number_values;
}

struct ResidueCount {
  int residue_count;
  IWString name;
  IWString header;
  IWString tag;
  // color_by
  // rgb
};

int
GatherResidueCount(iwstring_data_source& input,
                   const_IWSubstring& header,
                   std::vector<ResidueCount>& results) {
  const_IWSubstring token;
  int i = 0;
  if (! header.nextword(token, i) || 
      ! header.nextword(token, i)) {
    cerr << "GatherResidueCount::cannot extract tokens from header '" << header << "'\n";
    return 0;
  }
  int number_values;
  if (! token.numeric_value(number_values)) {
    cerr << "GatherResidueCount::invalid number of charges '" << header << "'\n";
    return 0;
  }

  results.resize(number_values);

  const_IWSubstring buffer;
  i = 0;
  for (int j = 0; j < number_values; ++j) {
    if (! nextword(input, buffer, i, token) ||
        ! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueCount:eof\n";
      return 0;
    }
    if (! token.numeric_value(results[j].residue_count)) {
      cerr << "GatherResidueCount::invalid residue_count '" << token << "'\n";
      return 0;
    }
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueCount:eof\n";
      return 0;
    }
    results[j].name = token;
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueCount:eof\n";
      return 0;
    }
    results[j].header = token;
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherResidueCount:eof\n";
      return 0;
    }
    results[j].tag = token;
    for (int notused = 0;  notused < 2; ++notused) {
      if (! nextword(input, buffer, i, token)) {
        cerr << "GatherResidueCount:eof\n";
        return 0;
      }
    }
  }

  return number_values;
}

int
GatherFormalCharges(iwstring_data_source& input,
                    const_IWSubstring& header,
                    std::vector<int> results) {

  const_IWSubstring token;
  int i = 0;
  if (! header.nextword(token, i) || 
      ! header.nextword(token, i)) {
    cerr << "GatherFormalCharges::cannot extract tokens from header '" << header << "'\n";
    return 0;
  }
  int number_values;
  if (! token.numeric_value(number_values)) {
    cerr << "GatherFormalCharges::invalid number of charges '" << header << "'\n";
    return 0;
  }

  const_IWSubstring buffer;
  i = 0;
  for (int j = 0; j < number_values; ++j) {
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherFormalCharges:eof\n";
      return 0;
    }
    int atom_number;
    if (! token.numeric_value(atom_number) || atom_number < 1) {
      cerr << "GatherFormalCharges:invalid atom " << token << '\n';
      return 0;
    }
    results.push_back(atom_number - 1);
  }

  return results.size();
}

int
GatherBonds(iwstring_data_source& input,
            const_IWSubstring& header,
            std::vector<int>& result) {

  const_IWSubstring token;
  int i = 0;
  if (! header.nextword(token, i) || 
      ! header.nextword(token, i)) {
    cerr << "GatherBonds::cannot extract tokens from header '" << header << "'\n";
    return 0;
  }
  int nbonds;
  if (! token.numeric_value(nbonds)) {
    cerr << "GatherBonds::invalid nbonds '" << header << "'\n";
    return 0;
  }

  const_IWSubstring buffer;

  if (! input.next_record(buffer)) {
    return 0;
  }

  i = 0;
  for (int j = 0; j < 2 * nbonds; ++j) {
    const_IWSubstring token;
    if (! nextword(input, buffer, i, token)) {
      cerr << "GatherBonds:cannot fetch next word\n";
      return 0;
    }
    int zatom;
    if (! token.numeric_value(zatom) || zatom < 1) {
      cerr << "GatherBonds:invalid atom number '" << token << "'\n";
      return 0;
    }
    result.push_back(zatom - 1);
  }

  if (result.size() % 2 != 0) {
    cerr << "GatherBonds:incorrect atoms for bonds " << result.size() << '\n';
    return 0;
  }

  return result.size();
}

int
Molecule::read_molecule_moe_ds(iwstring_data_source& input) {
  if (! input.skip_to("#system ")) {
    return 0;
  }
  if (! input.skip_to("#molecule ")) {
    return 0;
  }

  // Bonds are just stored as adjacent pairs in these arrays.
  std::vector<int> single_bonds, double_bonds, triple_bonds;
  std::vector<int> negative_charges, positive_charges;
  std::vector<Residue> residue_data;
  std::vector<ResidueCount> residue_count;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with("endmolecule")) {
      break;
    }
    int starting_line = input.lines_read();
    if (buffer.starts_with("#bond ") && buffer.ends_with("i=1")) {
      if (! GatherBonds(input, buffer, single_bonds)) {
        cerr << "Molecule::read_molecule_moe:bad bond list at " << starting_line << '\n';
        return 0;
      }
    }
    if (buffer.starts_with("#bond ") && buffer.ends_with("i=2")) {
      if (! GatherBonds(input, buffer, double_bonds)) {
        cerr << "Molecule::read_molecule_moe:bad bond list at " << starting_line << '\n';
        return 0;
      }
    }
    if (buffer.starts_with("#bond ") && buffer.ends_with("i=3")) {
      if (! GatherBonds(input, buffer, triple_bonds)) {
        cerr << "Molecule::read_molecule_moe:bad bond list at " << starting_line << '\n';
        return 0;
      }
    }
    if (buffer.starts_with("#attr ") &&
        buffer.contains(" ID i aName t aElement t aGeometry t aPosX r aPosY r aPosZ r")) {
      if (! MoeGatherAtoms(input, buffer)) {
        cerr << "Molecule::read_molecule_moe:bad atoms " << starting_line << '\n';
        return 0;
      }
    }

    if (buffer.starts_with("#attr ") &&
        buffer.contains(" ID i aIon i=-1")) {
      if (! GatherFormalCharges(input, buffer, negative_charges)) {
        cerr << "Molecule::read_molecule_moe:cannot gather -ve formal charges, line "
             << input.lines_read() << '\n';
        return 0;
      }
        }

    if (buffer.starts_with("#attr ") &&
        buffer.contains(" ID i aIon i=1")) {
      if (! GatherFormalCharges(input, buffer, positive_charges)) {
        cerr << "Molecule::read_molecule_moe:cannot gather +ve formal charges, line "
             << input.lines_read() << '\n';
        return 0;
      }
    }

    if (buffer.starts_with("#attr ") &&
        buffer.contains(" ID i rAtomCount i rName t rUID i rINS c rPos i rType t")) {
      if (! GatherResidueData(input, buffer, residue_data)) {
        cerr << "Molecule::read_molecule_moe:cannot gather residue data, line " <<
                input.lines_read() << '\n';
        return 0;
      }
    }

    if (buffer.starts_with("#attr ") &&
        buffer.contains(" ID i cResidueCount i cName t cHeader t cTag t cColorBy t cRGB ix")) {
      if (! GatherResidueCount(input, buffer, residue_count)) {
        cerr << "Molecule::read_molecule_moe:cannot gather residue count data, line " <<
                input.lines_read() << '\n';
        return 0;
      }
    }
  }

  // Place bonds and charges.
  MoePlaceBonds(single_bonds, SINGLE_BOND);
  MoePlaceBonds(double_bonds, DOUBLE_BOND);
  MoePlaceBonds(triple_bonds, TRIPLE_BOND);
  MoePlaceFormalCharges(negative_charges, -1);
  MoePlaceFormalCharges(positive_charges, +1);

  input.skip_past("#endsystem");

  // Try to assign a residue name to each atom.
  // Each molecule in the file will have a name, that is held in the residue_count
  // data. We could alternately assign names like ILE, LEU to the atoms, but that
  // is not the per fragment name.
  int atom_index = 0;
  int residue_index = 0;
  IWString mname;
  for (ResidueCount& rc : residue_count) {
#ifdef DEBUG_ASSIGN_ATOM_NAMES
    cerr << "Residue set contains " << rc.residue_count << " residues\n";
#endif
    for (int i = 0; i < rc.residue_count; ++i, ++residue_index) {
      const int atoms_in_residue = residue_data[residue_index].natoms;

#ifdef DEBUG_ASSIGN_ATOM_NAMES
      cerr << " residue " << residue_index << " contains " << atoms_in_residue << '\n';
#endif
      mname.append_with_spacer(rc.name, ':');

      for (int j = 0; j < atoms_in_residue; ++j, ++atom_index) {
#ifdef SETS_RESIDUE_NAME  // this works, but is not what I needed, ILE, LEU, CL, etc...
        //_things[atom_index]->set_user_specified_void_ptr(&residue_data[residue_index].name);
#endif
        _things[atom_index]->set_user_specified_void_ptr(&rc.name);
#ifdef DEBUG_ASSIGN_ATOM_NAMES
        cerr << " atom " << atom_index << " is " << rc.name << "\n";
#endif
      }
    }
  }
  
  set_name(mname);

  if (! moe::stream_for_moe_fragments.is_open()) {
    return 1;
  }

  atom_index = 0;
  cerr << "Will generate " << residue_data.size() << " residues\n";
  for (const Residue& residue : residue_data) {
    Molecule frag;
    cerr << "Creating fragment with " << residue.natoms << " atoms, atom_index " << atom_index << '\n';
    for (int i = 0; i < residue.natoms; ++i) {
      Atom* a = _things[atom_index + i];  // not const because we may set void_ptr
      frag.add(a->element());
      frag.setxyz(i, a->x(), a->y(), a->z());
      if (i == 0) {
        const IWString* s = reinterpret_cast<IWString*>(a->user_specified_void_ptr());
        if (s != nullptr) {
          frag.set_name(*s);
        }
      }
    }
    for (int i = 0; i < residue.natoms; ++i) {
      const Atom* a = _things[atom_index + i];
      for (const Bond* b : *a) {
        atom_number_t j = b->other(atom_index + i);
        if (j > atom_index + i) {
          continue;
        }
        atom_number_t a1 = b->a1() - atom_index;
        atom_number_t a2 = b->a2() - atom_index;
        if (! frag.ok_2_atoms(a1, a2)) {
          // inter residue bonds are fine - that is how a protein works.
          // cerr << "Molecule::read_molecule_moe:invalid bond atoms in frag with " << frag.natoms() << " atoms\n";
          // cerr << " atom index " << atom_index << ", i " << i << " bond " << *b << '\n';
          continue;
        }
        frag.add_bond(a1, a2, b->btype());
      }
    }
    atom_index += residue.natoms;
    moe::stream_for_moe_fragments << frag.smiles() << ' ' << frag.name() << '\n';
    moe::stream_for_moe_fragments.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// Add bonds of type `btype` to all adjacent atoms in `pairs`.
int
Molecule::MoePlaceBonds(const std::vector<int>& pairs,
                        bond_type_t btype) {
  const int n = pairs.size();
  for (int i = 0; i < n; i += 2) {
    atom_number_t a1 = pairs[i];
    atom_number_t a2 = pairs[i + 1];
    if (! add_bond(a1, a2, btype)) {
      cerr << "Molecule::PlaceBonds:cannot place bond between " << a1 << " and " << a2 << '\n';
      return 0;
    }
  }

  return 1;
}

int
Molecule::MoePlaceFormalCharges(const std::vector<int>& atoms,
                                int charge) {
  for (int a : atoms) {
    set_formal_charge(a, charge);
  }

  return 1;
}

struct MoeAtom {
  int atom_number;
  const Element* element;
  float x, y, z;
};

int
GetNextAtom(iwstring_data_source& input,
                const_IWSubstring& buffer,
                int& i,
                MoeAtom& moe_atom) {
  const_IWSubstring token;
  if (! nextword(input, buffer, i, token)) {
    cerr << "GetNextAtom:nextword failed\n";
    return 0;
  }

  if (! token.numeric_value(moe_atom.atom_number)) {
    cerr << "GetNextAtom:invalid atom number '" << token << "'\n";
    return 0;
  }
  // Convert to zero based numbering.
  --moe_atom.atom_number;

  if (! nextword(input, buffer, i, token) ||
      ! nextword(input, buffer, i, token)) {
    cerr << "GetNextAtom:cannot get atomic symbol\n";
    return 0;
  }

  const Element* e = get_element_from_symbol_no_case_conversion(token);
  if (e == nullptr) {
    cerr << "GetNextAtom:invalid atom '" << token << "'\n";
    return 0;
  }
  moe_atom.element = e;

  // Skip over hybridization state also.
  if (! nextword(input, buffer, i, token) ||
      ! nextword(input, buffer, i, token) ||
      ! token.numeric_value(moe_atom.x)) {
    cerr << "GetNextAtom:invalid or missing x coord '" << token << "'\n";
    return 0;
  }

  if (! nextword(input, buffer, i, token) || 
      ! token.numeric_value(moe_atom.y)) {
    cerr << "GetNextAtom:invalid or missing y coord '" << token << "'\n";
    return 0;
  }

  if (! nextword(input, buffer, i, token) || 
      ! token.numeric_value(moe_atom.z)) {
    cerr << "GetNextAtom:invalid or missing z coord '" << token << "'\n";
    return 0;
  }

  return 1;
}

int
Molecule::MoeGatherAtoms(iwstring_data_source& input,
                         const_IWSubstring& header) {
  int i = 0;
  const_IWSubstring token;
  if (! header.nextword(token, i) ||
      ! header.nextword(token, i)) {
    cerr << "Molecule::MoeGatherAtoms:cannot parse header '" << header << "'\n";
    return 0;
  }

  int number_atoms;
  if (! token.numeric_value(number_atoms)) {
    cerr << "Molecule::MoeGatherAtoms:invalid atom count " << header << '\n';
    return 0;
  }

  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    return 0;
  }

  resize(number_atoms);

  i = 0;

  MoeAtom moe_atom;
  for (int j = 0; j < number_atoms; ++j) {
    if (! GetNextAtom(input, buffer, i, moe_atom)) {
      cerr << "Molecule::MoeGatherAtoms:invalid atom '" << buffer << "' line " << input.lines_read() << '\n';
      return 0;
    }
    int atom_number = _number_elements;
    add(moe_atom.element);
    _things[atom_number]->setxyz(moe_atom.x, moe_atom.y, moe_atom.z);
  }

  return 1;
}
