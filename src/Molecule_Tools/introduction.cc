// An introduction to LillyMol.

// The most central object is the Molecule object. It has a large
// number of member functions, with many overloaded names.
// Most often a Molecule is introduced to a programme via
// a data_source_and_type<Molecule> producer object, which
// iteratively returns a newly allocated Molecule from an
// input stream.
// Alternatively, you can instantiate a Molecule, and either
// add atoms to it, or build it from a smiles string.

// A Molecule object contains a vector of Atoms, and these atoms
// can be iterated via a C++ range for loop. But almost all
// times an Atom is passed from the Molecule it will be const.
// The reason for this is that if certain atomic properties are
// altered, then certain molecular properties would need to be
// either changed or invalidated. For example, if an atom were
// to be assigned a new isotopic value, the smiles for the molecule
// must be recomputed.
// In reality, Molecule is lazy, and will only recompute the smiles
// if requested. Changing the isotope just invalidates the smiles.
// Similariy, removing a bond between atoms would need to invalidate
// fragment membership, ring membership, symmetry, aromaticity...
// But again, those quantities will only be recomputed if needed.

// One question that comes up is why is there so little use of
// standard data structures like std::string, std::vector etc.
// LillyMol began development in 1995, and at that time these
// concepts either did not exist or were poorly standardized.
// Over many years we fought with incompatible implementations
// of various things and in order to preserve sanity, implemented
// our own versions of things that offered cross-platform
// stability. Today the C++ world is in better shape.

// Nothing here raises exceptions. Generally LillyMol is
// thread safe. There are a great many optional settings that
// are stored in file scope static variables. Those are
// clearly not thread safe, but are usually set just once.
// With care, LillyMol has been successfully used in several
// multi-threaded applications.

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace lillymol_introcution {

using std::cerr;

// An empty molecule has no atoms or bonds, but can have a name.
void
DemoEmptyMolecule() {
  Molecule mol;
  mol.set_name("foo");
  // Should have zero atoms and molecular weight.
  cerr << "Empty molecule '" << mol.name() << "' has " << mol.natoms()
       << " atoms, amw " << mol.molecular_weight() << '\n';
  cerr << "Is the molecule empty " << mol.empty() << '\n';
}

// Using operator << with a string concatenates to the molecule name.
void
DemoStringConcat() {
  Molecule mol;
  mol << "hello";
  cerr << "Name updated to '" << mol.name() << "'\n";
  mol << " world";
  cerr << "Name updated to '" << mol.name() << "'\n";
}

// When presented with an invalid smiles, build_from_smiles will fail.
void
DemoCannotParseBadSmiles() {
  Molecule mol;
  if (mol.build_from_smiles("invalid")) {
    cerr << "Building from bad smiles succeeded, this should not happen.\n";
  } else {
    cerr << "Bad smiles cannot be parsed, good outcome.\n";
  }
}

void
DemoAtomicNumbers() {
  Molecule mol;
  if (!mol.build_from_smiles("C")) {
    cerr << "Unable to build methane!\n";
    return;
  }
  mol.set_name("methane");

  cerr << mol.name() << " has " << mol.natoms() << " atoms\n";

  // There are three ways to find the atomic number of atom atom.
  // Ask the molecule for the atomic number of a given atom.
  // We only have one atom in the molecule, so it is atom number 0.
  cerr << "Molecule's atomic number 0 " << mol.atomic_number(0) << '\n';

  // We can fetch an Atom and ask the atom for its atomic number.
  const Atom& a = mol.atom(0);
  cerr << "Atom's atomic number " << a.atomic_number() << '\n';

  // We can ask the atom for its Element and query the atomic number
  const Element* e = a.element();
  cerr << "Element's atomic number " << e->atomic_number() << '\n';
}

// Return a newly created molecule built from `smiles`.
// The molecule name is not set.
// If smiles interpretation fails, an empty molecule is returned.
// A better approach would be to use std::optional.
// Note that this may silently invoke a copy operation.
Molecule
MolFromSmiles(const char* smiles) {
  Molecule result;
  if (!result.build_from_smiles(smiles)) {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return Molecule();
  }

  return result;
}

// A very nice interface, but beware of copy operations happening.
std::optional<Molecule>
OptionalMolFromSmiles(const char* smiles) {
  Molecule result;
  if (!result.build_from_smiles(smiles)) {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return std::nullopt;
  }

  return result;
}

void
DemoOptionalMolFromSmiles() {
  std::optional<Molecule> m1 = OptionalMolFromSmiles("foo");
  if (m1) {  // should not happen.
    cerr << "OptionalMolFromSmiles suceeded for invalid smiles!!\n";
  } else {
    cerr << "Good, OptionalMolFromSmiles failed with invalid input\n";
  }

  std::optional<Molecule> m2 = OptionalMolFromSmiles("C1C2CC3CC1CC(C2)C3");
  if (! m2) {
    cerr << "OptionalMolFromSmiles did not parse adamantane\n";
    return;
  }

  cerr << "Adamantane has " << m2->nrings() << " SSSR rings\n";
}

// In LillyMol, there is a periodic table datastructure which contains
// an Element object for all elements known to the system. We can
// fetch pointers to those elements via various means.

// For historical reasons, to get an element from a two letter element
// we must use get_element_from_symbol_no_case_conversion. This should
// be fixed. We had some early mdl files where the symbols were all
// uppercase.
void
DemoElements() {
  const Element* carbon = get_element_from_atomic_number(6);
  const Element* uranium = get_element_from_symbol_no_case_conversion("U");
  const Element* chlorine = get_element_from_symbol_no_case_conversion("Cl");
  const Element* iron = get_element_from_atomic_number(26);

  cerr << "Carbon is " << carbon->symbol() << " atnum " << carbon->atomic_number() << '\n';
  cerr << "normal isotope " << carbon->normal_isotope() << '\n';
  cerr << "atomic_mass " << carbon->atomic_mass() << '\n';
  cerr << "exact_mass " << carbon->exact_mass() << '\n';
  cerr << "normal_valence " << carbon->normal_valence() << '\n';
  cerr << "number_alternate_valences " << carbon->number_alternate_valences() << '\n';
  cerr << "carbon organic " << carbon->organic() << '\n';
  cerr << "uranium organic " << uranium->organic() << '\n';
  cerr << "carbon needs_square_brackets " << carbon->needs_square_brackets() << '\n';
  cerr << "uranium needs_square_brackets " << uranium->needs_square_brackets() << '\n';
  cerr << "carbon is_halogen " << carbon->is_halogen() << '\n';
  cerr << "uranium is_halogen " << uranium->is_halogen() << '\n';
  cerr << "chlorine is_halogen " << chlorine->is_halogen() << '\n';
  cerr << "carbon is_in_periodic_table " << carbon->is_in_periodic_table() << '\n';
  cerr << "uranium is_in_periodic_table " << uranium->is_in_periodic_table() << '\n';
  cerr << "carbon is_metal " << carbon->is_metal() << '\n';
  cerr << "iron is_metal " << iron->is_metal() << '\n';

  // Computing the number of pi electrons and lone pairs from an element is
  // complex and you are better off asking a Molecule for the pi electron
  // or lone pair count on a given atom.
}

void
DemoAtomIteration() {
  Molecule mol = MolFromSmiles("CCC");

  // Loop through the atoms in `mol`, writing the atomic number.
  // Note that the iterator returns a pointer to an Atom.
  for (const Atom* atom : mol) {
    cerr << " atomic_number " << atom->atomic_number() << '\n';
  }

  // Or we can iterate based on atom numbers. This time we
  // can get a reference to an atom.
  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    cerr << " atomic number " << atom.atomic_number() << '\n';
  }
}

void
DemoBondIteration() {
  Molecule mol = MolFromSmiles("CC");
  for (const Bond* b : mol.bond_list()) {
    cerr << *b << '\n';
  }
}

// Molecules have an operator= overload and operator ==
void
DemoAssignment() {
  Molecule mol1 = MolFromSmiles("Oc1ccccc1");
  Molecule mol2 = mol1;
  cerr << "Smiles " << mol1.unique_smiles() << ' ' << mol2.unique_smiles() << '\n';

  // The operator == method is efficient, in that it only computes
  // the unique smiles if needed.
  cerr << "Same? " << (mol1 == mol2) << '\n';

  // Add a carbon atom to each molecule.
  mol1.add(get_element_from_atomic_number(6));
  cerr << "Same after adding " << (mol1 == mol2) << '\n';

  mol2.add(get_element_from_atomic_number(6));
  cerr << "Same after adding to both " << (mol1 == mol2) << '\n';
}

// The natoms method is overloaded so that if an atomic number, or
// Element, is given the number of atoms of that type is returned.
void
DemoNatoms() {
  Molecule mol = MolFromSmiles("OCN");

  cerr << "Molecule contains " << mol.natoms() << " atoms\n";

  cerr << "Oxygen atoms " << mol.natoms(8) << '\n';

  const Element* fluorine = get_element_from_atomic_number(9);

  cerr << "Molecule does not contain F " << mol.natoms(fluorine) << '\n';
}

// There are a wide variety of Molecule member functions for
// getting and setting isotopes.
void
DemoIsotopes() {
  Molecule mol = MolFromSmiles("CC");
  mol.set_isotope(0, 1);
  mol.set_isotope(1, 2);

  // Isotopes can be retrieved by querying either the molecule or the atom.
  cerr << "Isotope on atom 0 " << mol.isotope(0) << '\n';
  for (const Atom* atom : mol) {
    cerr << " isotope on atom " << atom->isotope() << '\n';
  }

  // Isotopes are arbitrary, unsigned, numbers.
  mol.set_isotope(0, 999);
  mol.set_isotope(1, std::numeric_limits<uint32_t>::max());
  cerr << "Ridiculous isotopes " << mol.smiles() << '\n';

  // We should be able to build a new molecule from `mol`
  // and their unique smiles should be the same.
  Molecule round_trip;
  round_trip.build_from_smiles(mol.smiles());
  cerr << "Same? " << (mol.unique_smiles() == round_trip.unique_smiles()) << '\n';

  // Remove all isotopes.
  // unset_isotopes() does the same thing.
  mol.transform_to_non_isotopic_form();
}

void
DemoIsotopesMolecularWeight() {
  Molecule mol = MolFromSmiles("OCN");
  float amw = mol.molecular_weight();
  cerr << "Non isotopic molecular weight " << amw << '\n';
  mol.set_isotope(0, 1);
  amw = mol.molecular_weight();
  cerr << "With isotopes, molecular weight is zero " << amw << '\n';

  amw = mol.molecular_weight_ignore_isotopes();
  cerr << "If isotopes are ignored, amw " << amw << '\n';

  // Isotope 1 gets counted as contributing 1.0 to the amw.
  amw = mol.molecular_weight_count_isotopes();
  cerr << "If isotopes included " << amw << '\n';
}

// Not all atoms need to have an atom map number. It is
// really just another arbitrary number that can be attached
// to an atom. They have particular meaning when dealing with
// reactions.
void
DemoAtomMapNumbers() {
  Molecule mol = MolFromSmiles("[CH3:5]-C");

  // Atom map numbers are not shown in a unique smiles. Atom map numbers
  // do not affect unique smiles formation.
  // Isotopes do and so they are included.
  cerr << "Smiles " << mol.smiles() << " usmi " << mol.unique_smiles() << '\n';

  // THe atom map number is a property of the Atom, so can be queried at
  // that level.

  for (const Atom* atom : mol) {
    cerr << "Atom map " << atom->atom_map() << '\n';
  }

  mol.reset_all_atom_map_numbers();
  // Note that removing the atom map numbers also removes the
  // implicit hydrogens known flag if possible - if the valence
  // is OK.
  cerr << "After removing atom map numbers " << mol.smiles() << '\n';

  // For many functions which return an atom number, returning INVALID_ATOM_NUMBER
  // is a common failure mode.
  atom_number_t twelve = mol.atom_with_atom_map_number(12);
  if (twelve != INVALID_ATOM_NUMBER) {
    cerr << "HUH, found atom map number 12\n";
  }
}

void
DemoImplicitAndExplicitHydrogens() {
  Molecule mol = MolFromSmiles("OCN");

  // Note that the ncon method, the number of connections to
  // an atom, does NOT reflect implicit Hydrogens.
  int matoms = mol.natoms();
  cerr << "Molecule begins with " << matoms << " atoms\n";
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " ncon " << mol.ncon(i) << " connections has " << mol.hcount(i) << " Hydrogens\n";
  }

  // Convert implicit Hydrogens to explicit.
  // Note that the hcount for the existing atoms, natoms, should
  // be unchaged, because the hcount() method includes both
  // explicit and implicit Hydrogens.
  // but now that the formerly implicit Hydrogens are now
  // explicit members of the connection table, the ncon values
  // for the existing atoms will change.

  cerr << "After implicit Hydrogens become explicit\n";
  mol.make_implicit_hydrogens_explicit();
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " ncon " << mol.ncon(i) << " connections has " << mol.hcount(i) << " Hydrogens\n";
  }

  cerr << "Molecule now contains " << mol.natoms() << " atoms\n";
}

// One of the most troublesome parts of LillyMol is the
// implicit hydrogens known attribute. If a smiles begins with
// C-[CH]-C where the centre Carbon is deficient one Hydrogen
// what happens when we add a Carbon atom to it, thereby
// satisfying the valence.
void 
DemoImplicitHydrogensKnown() {
  Molecule mol = MolFromSmiles("C-[CH]-C");
  cerr << "Valence? " << mol.valence_ok() << '\n';

  // Indeed the problem is with atom 1.
  for (int i = 0; i < mol.natoms(); ++i) {
    cerr << " atom " << i << " valence " << mol.valence_ok(i) << '\n';
  }

  // Add a carbon atom.
  // First fetch a pointer to the Element with atomic number 6.
  // Then call mol.add() with a pointer to the Element.
  // Internally, mol will instantiate an Atom which has
  // that element.
  const Element* carbon = get_element_from_atomic_number(6);
  mol.add(carbon);

  // The molecule now has two fragments.
  cerr << "Molecule has " << mol.number_fragments() << " fragments\n";

  // Add a bond from this last atom added to the middle atom to satisfy the valence.
  // The newly added atom will always have the highest atom number.
  mol.add_bond(1, 3, SINGLE_BOND);
  cerr << "Now " << mol.number_fragments() << " fragments\n";

  // Is the valence ok now?
  cerr << "After atom addition valence? " << mol.valence_ok() << '\n';

  // But the smiles shows that the centre atom still has its
  // implicit Hydrogens known attribute set. Even though it
  // is no longer necessary.
  cerr << "Smiles now " << mol.smiles() << '\n';
  cerr << "Is the value fixed " << mol.implicit_hydrogens_known(1) << '\n';

  // Remove the implicit Hydrogen known attribute.
  // Now the smiles does not have the square brackets.
  mol.set_implicit_hydrogens_known(1, 0);
  cerr << "Smiles now " << mol.smiles() << '\n';
  cerr << "Is the value fixed " << mol.implicit_hydrogens_known(1) << '\n';

  // This function can also be used, it examines all atoms and if possible
  // unsets the implicit hydrogen known flag.
  // In this case, the Carbon naturally has 3 Hydrogens, so the square
  // brackets are not necessary.
  Molecule mol2 = MolFromSmiles("C[CH3]");
  cerr << "Before removing unnecessary square brackets " << mol2.smiles() << '\n';
  mol2.unset_unnecessary_implicit_hydrogens_known_values();
  cerr << "After removing unnecessary square brackets " << mol2.smiles() << '\n';

  // If an implicit hydrogen flag is causing a valence error, this can
  // alleviate that problem.
  mol.build_from_smiles("C[C]C");
  cerr << "Before removing valence problem hydrogens known " << mol.smiles() << 
        " valence? " << mol.valence_ok() << '\n';
  mol.remove_hydrogens_known_flag_to_fix_valence_errors();
  cerr << "After removing valence problem hydrogens known " << mol.smiles() << 
        " valence? " << mol.valence_ok() << '\n';
}

// A data structure commonly encountered in LillyMol is a
// Set_of_Atoms. Despite the name, it is a vector of atom
// numbers, and there is no enforcement of uniqueness.
// It was a poor choice of name.
// Several Molecule member functions either consume or
// return Set_of_Atoms.
// Generally the ordering of the atoms in a Set_of_Atoms is
// unimportant. With an important exception.
// Rings inherit from a Set_of_Atoms, but they have the
// extra property that the atoms are ordered in a way
// that describes a circuit around the ring.
void
DemoSetOfAtoms() {
  Set_of_Atoms s{0, 5, 2};

  cerr << "Set contains " << s.size() << " atoms " << s << '\n';

  // Duplicate atom numbers are OK in a Set_of_Atoms, but
  // functions consuming them may be unhappy with duplicate
  // entries.
  s << 2;
  cerr << "See the dup " << s << '\n';

  // Remove the dup.
  s.chop();

  s.add_if_not_already_present(2);
  cerr << "Dup not added " << s << '\n';

  // These have a convenient `scatter` operation.

  std::unique_ptr<int[]> arbitrary = std::make_unique<int[]>(10);
  std::fill_n(arbitrary.get(), 10, 0);
  // For every atom number in `s`, set the corresponding array index to '2'.
  s.set_vector(arbitrary.get(), 2);
  // this also works with a std::vector.

  // There are many gather type operations associated with the
  // Set_of_Atoms class. See the header for info.
}

// Random smiles are useful for testing. Algorithms should not
// depend on the order of the smiles.
void
DemoRandomSmiles() {
  // Caffeine
  Molecule mol = MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");

  // Note that we make a copy of the unique smiles.
  // Generally we get a reference to the unique smiles,
  //   IWString& starting_unique_smiles = mol.unique_smiles()
  // But in this case the smiles of `mol` will be destroyed when the
  // random smiles are formed below, and the reference would
  // be invalid.

  const IWString starting_unique_smiles = mol.unique_smiles();

  for (int i = 0; i < 10; ++i) {
    cerr << " random smiles " << mol.random_smiles() << '\n';
    Molecule hopefully_the_same;
    hopefully_the_same.build_from_smiles(mol.random_smiles());
    cerr << "Same? " << (starting_unique_smiles == hopefully_the_same.unique_smiles()) << '\n';
  }
}

//
void
DemoRemoveAtom() {
  Molecule mol = MolFromSmiles("CCCCCCCCCCCCC");

  cerr << "Before atom removals, molecule contains " << mol.natoms() << ' ' << mol.smiles() << '\n';
  // Any atom can be removed with remove_atom.
  mol.remove_atom(0);
  cerr << "Molecule now contains " << mol.natoms() << " atoms " << mol.smiles() << "\n";

  // Atom removal may lead to multiple fragments
  mol.remove_atom(1);
  cerr << "Molecule now contains " << mol.natoms() << " atoms " << mol.smiles() << "\n";

  // If you wish to remove multiple atoms, that can be done one at a time.
  // It will be most efficient if you iterate from high to low atom numbers.
  // for (int i = mol.natoms() - 1; i >= 0; --i) {
  //   if (some condition) {
  //     mol.remove_atom(i);
  //   }
  // }

  // If you have an array with the atoms to be removed set, use that.
  const int matoms = mol.natoms();
  std::unique_ptr<int[]> to_remove = std::make_unique<int[]>(matoms);
  std::fill_n(to_remove.get(), matoms, 0);
  for (int i = 0; i < matoms; i += 2) {
    to_remove[i] = 1;
  }
  mol.remove_atoms(to_remove.get());
  cerr << "After removing every second atom " << mol.smiles() << '\n';
}

void
DemoFormalCharges() {
  Molecule mol = MolFromSmiles("CC[N+H3]");

  cerr << "Does the molecule have formal charges " << mol.has_formal_charges() <<
          " net " << mol.net_formal_charge() << '\n';
  cerr << "Charge " << mol.formal_charge(2) << '\n';
  cerr << "Charge on atom " << mol.atom(2).formal_charge() << '\n';

  mol.set_formal_charge(2, 0);
  cerr << "After reset " << mol.formal_charge(2) << '\n';
}

// Partial charges are a property of the molecule only.
// The molecule starts with the partial charge being null.
// The first time you set the partial charge on any atom,
// the partial charge array is allocated and set to zero.
void
DemoPartialCharges() {
  Molecule mol = MolFromSmiles("CC");
  cerr << "Should be no partial charges " << mol.has_partial_charges() << '\n';
  mol.set_partial_charge(1, 0.22);
  cerr << "After setting charge on atom 1 charge on 0 " << mol.partial_charge(0) << '\n';

  // Removing an atom will preserve partial charges.
  mol.remove_atom(0);
  cerr << "Charges present now? " << mol.has_partial_charges() << '\n';

  // Adding an atom will see that atom enter with zero partial charge.
  mol.add(get_element_from_atomic_number(6));
  cerr << "Partial charge on newly added atom " << mol.partial_charge(1) << '\n';

  // Partial charges can be dropped.
  mol.invalidate_partial_charges();

  // Some partial charge types are available.
  mol.compute_Abraham_partial_charges();
  cerr <<"Partial charge type " << mol.partial_charge_type() << '\n';
}

void
DemoNcon() {
  Molecule mol = MolFromSmiles("CC(C)(C)C");

  // The molecule is
  //              [C:2]
  //                |
  //                |
  //     [C:0] -- [C:1] -- [C:4]
  //                |
  //                |
  //              [C:3]
  // Atoms have different connectivity.

  // We can get the degree of any atom by asking either the
  // molecule, or the atom.
  // In addition get a smarts equivalent for that atom - this
  // is mostly useful for debugging.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    cerr << " atom " << i << " has " << mol.ncon(i) << " connections, atom "
         << atom.ncon() << " type " << mol.smarts_equivalent_for_atom(i) << '\n';
  }
}

void
DemoBondIterationAroundAtom() {
  Molecule mol = MolFromSmiles("CC(C)(C)C");

  // For each atom, enumerate the connected atoms.
  // Iteration is over Bond objects associated with each atom.
  // The Bond contains two atom numbers, we use the bond->other
  // method to determine what the other atom is.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    for (const Bond* bond : atom) {
      cerr << "atom " << i << " bonded to " << bond->other(i) <<
              " type " << bond->btype() << '\n';
    }
  }
}

// It is possible to iterate through the connections to an atom
// by iterating through the Molecule. It is more efficient and elegant
// to iterate via the Atom rather than via the Molecule.
void
DemoBondIterationMolecule() {
  Molecule mol = MolFromSmiles("CO");

  int ncon = mol.ncon(0);
  for (int i = 0; i < ncon; ++i) {
    atom_number_t j = mol.other(0, i);
    cerr << "Atom 0 is bonded to " << j << '\n';
  }
}

// While it is possible to get the connections to an
// atom, this is usually less efficient than iterating,
// since iterating does not require the creation of a
// Set_of_Atoms object.
void
DemoConnections() {
  Molecule mol = MolFromSmiles("C(C)O");

  const Set_of_Atoms connected = mol.connections(0);
  cerr  << " atoms connected to 0 " << connected << '\n';
}

// A common operation is the need to add all atoms and bonds
// from one molecule to another.
// There are variants of add_molecule that allow for
// excluding certain atoms from the added Molecule.
void
DemoAddMolecule() {
  Molecule mol1 = MolFromSmiles("CC");
  Molecule mol2 = MolFromSmiles("NN");

  cerr << "mol1 contains " << mol1.natoms() << " atoms\n";
  mol1.add_molecule(&mol2);
  cerr << "mol1 contains " << mol1.natoms() << " atoms and " <<
           mol1.number_fragments() << " fragments\n";

  mol1.add_bond(1, 2, SINGLE_BOND);
  cerr << "After bond addition " << mol1.number_fragments() << " fragments\n";
}

// An existing bond type can be changed.
// Note that the symbols SINGLE_BOND, DOUBLE_BOND, DOUBLE_BOND
// should not be interpreted as ints. The bits in these
// are used by LillyMol for various internal purposes.
void
DemoChangeBond() {
  Molecule mol = MolFromSmiles("CC");

  cerr << "Initial smiles " << mol.smiles() << '\n';

  const std::vector<int> bonds{SINGLE_BOND, DOUBLE_BOND, TRIPLE_BOND};
  for (auto btype : bonds) {
    mol.set_bond_type_between_atoms(0, 1, btype);
    cerr << "New bond " << mol.smiles() << '\n';
  }
}

// It is possible to change the element associated with an atom.
// Existing bonds are preserved.
// It is up to you to ensure that the resulting chemistry is valid.
void
DemoChangeElement() {
  Molecule mol = MolFromSmiles("CCC");

  // Change second atom to Fluorine!
  mol.set_atomic_number(1, 9);

  cerr << "Doubly bonded Fluorine! " << mol.smiles() << '\n';
  cerr << "ok valence? " << mol.valence_ok() << '\n';
}

// Start with ethane, verify one fragment, then remove the bond
// to generate two methanes.
void
DemoFragments() {
  Molecule mol = MolFromSmiles("CC");

  cerr << "Methane has " << mol.number_fragments() << " fragments\n";
  cerr << "Ethane contains only one edge " << mol.nedges() << '\n';

  // there are only two atoms, remove the bond between them.
  mol.remove_bond_between_atoms(0, 1);

  // We should have two fragments now.
  cerr << "After bond removal " << mol.number_fragments() << " fragments\n";

  // Fragment membership is only known by the Molecule.
  // Fragments are assigned sequential numbers, starting with zero.
  // We can ask the molecule, in which fragment is the i'th atom?
  // Atom's do not know anything about fragment membership.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << " in fragment " << mol.fragment_membership(i) << '\n';

    // const Atom& atom = mol.atom(i);
    // atom.fragment_membership() CANNOT WORK
  }

  // we can ask the molecule how many atoms are in each fragment.
  const int nf = mol.number_fragments();
  for (int i = 0; i < nf; ++i) {
    cerr << " Fragment " << i << " contains " << mol.atoms_in_fragment(i) << " atoms\n";
  }

  // To determine if two atoms are in the same fragment, compare 
  // mol.fragment_membership(i) == mol.fragment_membership(j);

  // Put the molecule back by adding a bond

  mol.add_bond(0, 1, SINGLE_BOND);
  cerr << "After recreating nfrag " << mol.number_fragments() << '\n';
  cerr << "Same frag " << (mol.fragment_membership(0) == mol.fragment_membership(1)) << '\n';
}

// when we have a multi-fragment molecule, we can create
// molecules from each fragment.
void
DemoCreateComponents() {
  Molecule mol = MolFromSmiles("C.C.C.C.C.C.C.C.C");
  resizable_array_p<Molecule> fragments;
  mol.create_components(fragments);
  cerr << "Created " << fragments.size() << " fragments from " << mol.smiles() << '\n';

  // The fragments should all have the same smiles.
  // Note that we cannot use 'const Molecule* m' here because generating
  // a smiles is a potentially non-const operation - see lazy evaluation.
  for (Molecule* m : fragments) {
    if (m->unique_smiles() != "C") {
      cerr << "Invalid fragment unique smiles\n";
    }
  }

  // We can remove a fragment by fragment number
  mol.delete_fragment(2);

  // Or remove whatever fragment contains atom 3
  mol.remove_fragment_containing_atom(3);
}

// Normally getting the largest fragment is straightforward and
// yields the expected result.
void
DemoReduceToLargestFragment() {
  Molecule mol = MolFromSmiles("CC.C");
  mol.reduce_to_largest_fragment();
  cerr << "after reduce_to_largest_fragment have " << mol.natoms() << " atoms " << mol.smiles() << '\n';

  // There are however cases where just choosing the largest
  // fragment is not the best choice. Prefer using
  // reduce_to_largest_fragment_carefully() which has
  // heuristics that help avoid some common problems.
}

// create_subset creates a new molecule that contains
// just some of the atoms in the starting molecule.
// Pass an int* 
void
DemoCreateSubset() {
  Molecule mol = MolFromSmiles("c1ccccc1OC");

  const int matoms = mol.natoms();

  std::unique_ptr<int[]> subset = std::make_unique<int[]>(matoms);
  std::fill_n(subset.get(), matoms, 0);

  // The first 6 atoms are the benzene ring. That is the subset.
  // Mark the first 6 entries with a non-zero integer.
  constexpr int kOne = 1;
  std::fill_n(subset.get(), 6, kOne);

  Molecule benzene;
  mol.create_subset(benzene, subset.get(), kOne);
  cerr << "Should be benzene " << benzene.unique_smiles() << '\n';
}

void
DemoDistanceMatrix() {
  Molecule mol = MolFromSmiles("CCCCCCCCCC");

  cerr << "Longest path " << mol.longest_path() << '\n';

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      int d = mol.bonds_between(i, j);
      if (d == 1) {
        cerr << " atoms " << i << " and " << j << " are bonded\n";
      } else {
        cerr << " atoms " << i << " and " << j << " are " << d << " bonds apart\n";
      }
    }
  }
}

void
DemoRing() {
  // Adjacent three and four membered rings.
  Molecule mol = MolFromSmiles("C1CCC1C1CC1");

  // We should get two rings.
  cerr << "Molecule contains " << mol.nrings() << " rings\n";

  // As listed here, the ring will have no status for aromaticity,
  // because aromaticity has not been determined.
  const int nr = mol.nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring* ring = mol.ringi(i);
    cerr << " ring " << i << " contains " << ring->size() << " atoms\n";
    cerr << " atoms " << *ring << '\n';
  }

  // We can check whether atoms are in the same ring
  cerr << "same ring 0,1 " << mol.in_same_ring(0, 1)  << '\n';
  cerr << "same ring 3,4 " << mol.in_same_ring(3, 4)  << '\n';
  cerr << "same ring 4,5 " << mol.in_same_ring(4, 5)  << '\n';

  // We can ask the molecule about the ring membership of an atom.
  // Note that Atom objects do not know anything about their ring status.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " Atom " << i << ' ' << mol.smarts_equivalent_for_atom(i) <<
         " in " << mol.nrings(i) << " rings, ring bond count " << mol.ring_bond_count(i) << '\n';
  }

  // Sometimes it is convenient to just get all the ring membership values in
  // an array - saves the overhead of querying the Molecule multiple times.
  // If `mol` is changed however, this pointer becomes invalid.
  const int * ring_membership = mol.ring_membership();

  // Querying the molecule for each atom is however safe if the molecule is
  // being altered.
  int ring_atoms = 0;
  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i] > 0) {
      ++ring_atoms;
    }
  }
  cerr << ring_atoms << " of " << matoms << " atoms in a ring\n";

  // If you have your own array, you can call mol.ring_membership(your_array) and
  // it will be filled with the ring membership of each atom. Here if the
  // molecule gets changed, you are on your own.
}

// There is a convenient Ring iterator.
void
DemoRingIteration () {
  // Cubane has 6 faces, but 5 sssr rings.
  Molecule mol = MolFromSmiles("C12C3C4C1C5C2C3C45");

  cerr << "Cubane has " << mol.nrings() << " SSSR rings\n";
  cerr << "Cubane has " << mol.non_sssr_rings() << " non SSSR rings\n";

  for (const Ring* r : mol.sssr_rings()) {
    cerr << " cubane ring " << *r << '\n';
  }

  int non_ssr_rings = mol.non_sssr_rings();
  for (int i = 0; i < non_ssr_rings; ++i) {
    cerr << " cubane NON SSSR " << *mol.non_sssr_ring(i) << '\n';
  }
}

void
DemoRingSystem() {
  // Fused 3 and 4 membered rings.

  //    C ---- C                         *
  //    |      |  \                      *
  //    |      |   C                     *
  //    |      |  /                      *
  //    C ---- C                         *

  Molecule mol = MolFromSmiles("C12CCC1C2");

  // Contains 2 rings. 
  cerr << "Molecule contains " << mol.nrings() << " rings\n";

  // The two rings must have the same fused system identifier
  // Generally rings are sorted by size so the smallest rings should be first.
  for (const Ring* r : mol.sssr_rings()) {
    cerr << *r << '\n';
    cerr << "fused rings:fused_system_identifier " << r->fused_system_identifier() << '\n';
    cerr << "is_fused " << r->is_fused() << '\n';
  }

  cerr << "Atoms in same ring system " << mol.in_same_ring_system(1, 2) << '\n';

  // only write a failure - cut down on uninteresting output.
  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (! mol.in_same_ring_system(i, j)) {
        cerr << "Atoms " << i << " and " << j << " not in same ring system\n";
      }
    }
  }

  // For each ring, we can fetch pointers to the Ring's attached.
  // In this case, where will be two atoms in common between the two rings.
  const int nr = mol.nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring* ri = mol.ringi(i);
    const int nfused = ri->fused_ring_neighbours();
    for (int j = 0; j < nfused; ++j) {
      const Ring* nbr = ri->fused_neighbour(j);
      cerr << "ring " << *ri << " fused to " << *nbr << '\n';
    }
  }

  // Two separated rings will have different fused_system_identifier

  mol.build_from_smiles("C1CC1C1CC1");
  for (const Ring* r : mol.sssr_rings()) {
    cerr << "separated rings: fused_system_identifier " <<r->fused_system_identifier() << '\n';
    cerr << "is_fused " << r->is_fused() << '\n';
  }
}

void
DemoAreBonded() {
  Molecule mol = MolFromSmiles("CC.C");

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (mol.are_bonded(i, j)) {
        cerr << "Atoms " << i << " and " << j << " bonded\n";
      } else {
        cerr << "Atoms " << i << " and " << j << " not bonded\n";
      }
    }
  }
}

// A LillyMol smiles extension is to allow coordinates in smiles.
// {{x,y,z}} is appended after each atom.
// This is more compact than a .sdf representation, and allows
// molecules to be processed one per line.

// Now shown here are the member functions
// set_bond_length, set_bond_angle and set_dihedral which are used
// to alter the geometry of a Molecule.

// rotate_atoms can rotate the whole molecule.
void
DemoCoordinatesInAtoms() {
  Molecule mol = MolFromSmiles("C{{-1,1,0}}C{{0,0,0}}C{{1,0,0}}C{{2,1,1}}");

  cerr << "Distance between 0 1 " << mol.distance_between_atoms(0, 1) << '\n';
  cerr << "Bond angle 0 1 2 " << (mol.bond_angle(0, 1, 2) * RAD2DEG) << '\n';
  cerr << "Dihedral angle 0 1 2 3 " << (mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG) << '\n';

  // Note that the distances above are geometric distances.
  // Topological distances are frequently of interest.
  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      cerr << "Bonds between " << i << " and " << j << ' ' << mol.bonds_between(i, j) << '\n';
    }
  }

  // We can ask the molecule or an individual atom for their coordinates.
  // Internally, the Molecule fetches the Atom.
  cerr << "Atom 0 at " << mol.x(0) << ',' << mol.y(0) << ',' << mol.z(0) << '\n';
  const Atom& atom = mol.atom(0);
  cerr << "Atom 0 at " << atom.x() << ',' << atom.y() << ',' << atom.z() << '\n';

  // We can get the spatial extremities of the molecule.
  coord_t xmin, xmax, ymin, ymax, zmin, zmax;
  mol.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);
  cerr << "X btw " << xmin << " and " << xmax << '\n';

  // Move the molecule to a quadrant of space.
  mol.translate_atoms(-xmin, -ymin, -zmin);
}

// The spinach is the complement of the scaffold.
void
DemoScaffold() {
  Molecule mol = MolFromSmiles("Fc1ccc(cc1)CC(=O)Nc1oncc1");

  const int matoms = mol.natoms();
  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);
  mol.identify_spinach(spinach.get());

  // Invert the spinach so we get the scaffold.
  for (int i = 0; i < matoms; ++i) {
    spinach[i] = ! spinach[i];
  }

  Molecule scaffold;
  mol.create_subset(scaffold, spinach.get(), 1);

  cerr << "Scaffold is " << scaffold.unique_smiles() << '\n';
}

// Performing a substructure search requires instantiation of
// a Substructure_Query object. That is then constructed
// from smarts, a query file, a proto, or other means.
// Once the Substructure_Query is built, use it to perform
// substructure searches on Molecule's.
void
DemoSubstuctureSearch() {
  Molecule mol = MolFromSmiles("OCN");

  Substructure_Query qry;
  qry.create_from_smarts("OCN");

  // For historical reasons, this method takes a pointer.
  const int nhits = qry.substructure_search(&mol);
  cerr << "Query makes " << nhits << " matches\n";
}

// The example above did not save the matched atoms. If we need
// the matched atoms, pass a Substructure_Results to the search.
void
DemoSubstuctureSearchSaveResults() {
  Molecule mol = MolFromSmiles("OCN");

  Substructure_Query qry;
  qry.create_from_smarts("OCN");

  Substructure_Results sresults;
  const int nhits = qry.substructure_search(mol, sresults);

  // Loop over the hits - one in this case. There is an iterator
  // for shown below.
  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);
    cerr << " atoms " << *e << '\n';
  }
}

// If all we care about is whether or not a query matches
// and we do not care about how many matches there might 
// be, and we do not want to save the matched atoms, we
// can limit the scope of the search.
// The initial molecule/smarts generates 559872 query
// matches. Every 10,000 matches found, there will be
// a message to stderr informing you of a potentially
// catastrophic matching.
// 559872 is 2^8 * 3^7
void
DemoParsimoniousSearch() {
  const char* smiles = "C(C)(C)(C)c1c(C(C)(C)(C))c(C(C)(C)(C))c(C(C)(C)(C))c(C(C)(C)(C))c1C(C)(C)(C)";
  Molecule mol = MolFromSmiles(smiles);

  Substructure_Query qry;
  qry.create_from_smarts(smiles);

  Substructure_Results sresults;
  int nhits = qry.substructure_search(mol, sresults);
  cerr << "Without constraints, find " << nhits << " substructure matches\n";

  qry.set_max_matches_to_find(1);
  qry.set_save_matched_atoms(0);

  nhits = qry.substructure_search(mol, sresults);
  cerr << "With constraints find " << nhits << " substructure matches\n";
}

// Using a Molecule_to_Match as the target for the query
// can be more efficient because computed properties are cached.
// Definitely recommended if the same molecule is being searched
// by multiple queries.
void
DemoSubstuctureSearchWithTarget() {
  Molecule mol;
  mol.build_from_smiles("Oc1ccccc1");

  Molecule_to_Match target(&mol);

  Substructure_Query qry;
  qry.create_from_smarts("O-c:c");

  Substructure_Results sresults;

  int nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (unconstrained)\n";

  // If only one embedding per start atom, just one match will be found.
  qry.set_find_one_embedding_per_atom(1);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (one embedding per start atom)\n";

  qry.set_find_one_embedding_per_atom(0);

  // Turn off perception of symmetry equivalent matches.
  qry.set_perceive_symmetry_equivalent_matches(0);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (no symmetry)\n";

  qry.set_perceive_symmetry_equivalent_matches(1);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (back to default)\n";
}

// Set up a vector of Substructure_Query objects.
// A resizable_array_p is a vector that contains pointers
// to objects. The resizable_array_p manages ownership of
// the objects.
void
DemoMultipleQueries() {
  resizable_array_p<Substructure_Query> queries;
  constexpr int kNqueries = 10;

  // Build up a set of smarts that are C, CC, CCC ...
  IWString smarts = 'C';
  for (int i = 0; i < kNqueries; ++i) {
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    qry->create_from_smarts(smarts);
    // The name of the query is the same as the smarts.
    qry->set_comment(smarts);
    queries << qry.release();

    smarts += 'C';
  }

  Molecule mol = MolFromSmiles("CCCCC");
  Molecule_to_Match target(&mol);

  // It is interesting to contemplate the number of matches found
  // to each of these queries.
  Substructure_Results sresults;
  for (Substructure_Query* qry : queries) {
    int nhits = qry->substructure_search(target, sresults);
    cerr << nhits << " matches to query " << qry->comment() << '\n';
  }
}

void
DemoQueryFromProto() {
  Molecule mol;
  mol.build_from_smiles("Oc1ccccc1");

  const std::string s = R"pb(
query {
    min_nrings: 1
    ring_specifier {
      aromatic: true
      fused: 0
      base {
        heteroatom_count: 0
      }
    }
    smarts: "O-c:c"
}
)pb";

  // Build a proto from the text.
  SubstructureSearch::SubstructureQuery proto;  
  if (! google::protobuf::TextFormat::ParseFromString(s, &proto)) {
    cerr << "Cannot parse proto\n";
    return;
  }

  // Now build a substructure query from the proto.

  Substructure_Query qry;
  qry.ConstructFromProto(proto);

  Substructure_Results sresults;
  qry.substructure_search(mol, sresults);
  for (const Set_of_Atoms* embedding : sresults.embeddings()) {
    cerr << "From proto matches atom " << *embedding << '\n';
  }
}

// If enabled, arbitrary two character combinations are valid elements.
void
DemoStrangeElements() {
  set_auto_create_new_elements(1);
  Molecule mol = MolFromSmiles("[Th][Eq][U]IC[K][Br]O[W]NFO[X][Ju][Mp]SO[Ve][Rt][La][Zy][D]O[G]");
  cerr << "Smiles " << mol.smiles() << '\n';

  cerr << "contains_non_periodic_table_elements " << mol.contains_non_periodic_table_elements() << '\n';
  cerr << "organic " << mol.organic_only() << '\n';

  Substructure_Query qry;
  qry.create_from_smarts("[Th][Eq][U]IC[K]");

  Substructure_Results sresults;
  int nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to '[Th][Eq][U]'\n";

  // But if there are element names that might collide with smarts
  // directives, we can enclose that in #{}
  qry.create_from_smarts("[#{Rt}][La]");
  nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to '[Rt][La]'\n";

  set_auto_create_new_elements(0);
}

// If enabled, any sequence of letters is an ok element
void
DemoAnyLengthElement() {
  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  Molecule mol = MolFromSmiles("[Ala]1[Arg][Asn][Asp]1");

  cerr << "Peptide " << mol.smiles() << '\n';

  // All regular smarts directives are available.
  Substructure_Query qry;
  qry.create_from_smarts("[#{Ala}D2R]1[#{Arg}D2x2][#{Asn}][#{Asp}D2]1");

  Substructure_Results sresults;
  const int nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to Ala-Arg-Asn-Asp\n";

  // This can be extended to searching peptides more generally.
  // [#{Ala}D1][#{Arg}]...{6-12}[#{Ala}]...{>4}[#{Tyr}D2]&&<4[#{Phe}]
  // which looks for a sequence that starts with Ala, then Arg,
  // between 6 and 12 other resides to Ala, and then at least
  // 4 resides to a non terminal Tyr (D2). AND there must be
  // fewer than 4 Phe residues in the sequence.

  // And of course logical expressions at the atomic level work
  // [#{Ala},#{Gly},#{Phe};D1][#{Arg}]...{6-12}[#{Ala}]...{>4}[#{Phe}]

  set_atomic_symbols_can_have_arbitrary_length(0);
  set_auto_create_new_elements(0);
}

// Aromaticity is only computed if requested. In the rings demo
// we saw that if the rings are requested, without aromaticity
// having been queried, it was not comptuted.
void
DemoAromaticity() {
  Molecule mol = MolFromSmiles("c1ccccc1");

  cerr << "Is the molecule aromatic " << mol.contains_aromatic_atoms() << '\n';
  cerr << "Count the aromatic atoms " << mol.aromatic_atom_count() << '\n';

  cerr << "First atom aromatic? " << mol.is_aromatic(0) << '\n';
  cerr << "Ring is now aromatic " << mol.ringi(0) << '\n';
  cerr << "We can ask if a ring is aromatic " << mol.ringi(0)->is_aromatic() << '\n';

  // If you are going to be querying rings for aromaticity,
  // call compute_aromaticity_if_needed which will be a no-op
  // if aromaticity is already computed.
  mol.compute_aromaticity_if_needed();

  // Beware: if you alter the molecule, all pointers
  // and references to internal structures, like rings,
  // fragments and smiles, are invalidated.
}

// During unique smiles determination, symmetry is perceived.
void
DemoSymmetry() {
  Molecule mol = MolFromSmiles("Oc1c(C)cc(O)cc1C");

  cerr << "Symmetric molecule contains " << mol.number_symmetry_classes() << " groups of symmetric atoms\n";

  // It can be convenient to see which atom number is which
  cerr << mol.isotopically_labelled_smiles() << '\n';

  // THe two methyls are the same
  cerr << "Methyl symmetry " << mol.symmetry_class(3) << " and " << mol.symmetry_class(9) << '\n';
  // and the two atoms in the ring adjacent to the methyls.
  cerr << "Methyl symmetry " << mol.symmetry_class(2) << " and " << mol.symmetry_class(8) << '\n';

  // But the two oxygens are not symmetric.
  cerr << "Oxygens not equivalent " << mol.symmetry_class(0) << " and " << mol.symmetry_class(6) << '\n';

  // We can get the symmetry equivalent atoms.
  // Note that the starting atom is not included - not sure if that was a good idea,
  // but that is how it works today.
  Set_of_Atoms same_as_3;
  mol.symmetry_equivalents(3, same_as_3);
  cerr << "Same as 3 " << same_as_3 << '\n';
}

// Generally we think of symmetry in terms of the atoms, but it can also
// be useful to examine symmetry as it relates to bonds.
void
DemoBondSymmetry() {
  Molecule mol = MolFromSmiles("FC(F)(F)Cl");

  const int nedges = mol.nedges();

  std::unique_ptr<int[]> bond_symm = std::make_unique<int[]>(nedges);

  // There are two variations on this function that make different
  // trade-offs
  mol.bond_symmetry_class_small_memory(bond_symm.get());

  // We should see that all the C-F bonds are one type
  // the the C-Cl bond is different.
  for (int i = 0; i < nedges; ++i) {
    const Bond* b = mol.bondi(i);
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    cerr << " bond btw " << mol.smarts_equivalent_for_atom(a1) << " and " << 
            mol.smarts_equivalent_for_atom(a2) << " symm " << bond_symm[i] << '\n';
  }
}

// A common operation is to ask how many heteroatoms are
// attached to an atom.
void
DemoAttachedHeteroatomCount() {
  Molecule mol = MolFromSmiles("CC(F)(F)F");

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << ' ' << mol.atomic_symbol(i) << " has " <<
             mol.attached_heteroatom_count(i) << " heteroatoms attached\n";
  }
}

// A common operation is to ask whether or not an atom
// is multiply bonded to a heteroatom.
// Note that it is a yes/no, so the Nitrogen that
// is part of the Nitro group reports 1.
// Note that there is a doubly_bonded_oxygen_count() member
// function.
void
DemoMultipleBondToHeteroatom() {
  Molecule mol = MolFromSmiles("OC(=O)CN(=O)=O");

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << ' ' << mol.smarts_equivalent_for_atom(i) << ' ' <<
            mol.multiple_bond_to_heteroatom(i) << " multiple bond to heteroatom\n";
  }
}

// Only tetrahedral chirality is supported.
void
DemoChirality() {
  Molecule mol = MolFromSmiles("F[C@H](C)N");

  const int nchiral = mol.chiral_centres();
  cerr << mol.smiles() << " contains " << nchiral << " chiral centres\n";

  // The chiral centre is described by up to 5 atoms.
  // The centre atom, top front and top back,
  // Then left down, and right down. Any of these atoms
  // may be an implicit hydrogen, or a pi electron.
  cerr << "initial smiles " << mol.unique_smiles() << '\n';
  for (const Chiral_Centre* c : mol.ChiralCentres()) {
    // Currently no operator<< for a Chiral_Centre object.
    c->debug_print(cerr);
  }

  mol.invert_chirality_on_atom(1);
  cerr << "After inversion " << mol.unique_smiles() << '\n';
  for (const Chiral_Centre* c : mol.ChiralCentres()) {
    c->debug_print(cerr);
  }

  // Get the chiral centre on a particular atom. If no chiral
  // centre on that atom, will be nullptr.
  // const Chiral_Centre* c = mol.chiral_centre_at_atom(1);

  // Chiral_Centre's can be added and removed.
}

// There are times when it is convenient to have an atom
// in a given position in the molecule.
void
DemoSwapAtoms() {
  Molecule mol = MolFromSmiles("CN");

  cerr << mol.unique_smiles() << " initial atomic numbers\n";
  for (const Atom* a : mol) {
    cerr << a->atomic_number() << '\n';
  }

  mol.swap_atoms(0, 1);
  cerr << mol.unique_smiles() << " after swap atomic numbers\n";
  for (const Atom* a : mol) {
    cerr << a->atomic_number() << '\n';
  }

  // It is possible to move an individual atom to the end of the
  // connection table,
  mol.move_atom_to_end_of_atom_list(0);
  // which is the same as mol.swap_atoms(i, mol.natoms() - 1);

  // It is also possible to move all atoms of a given type to the
  // end of the connection table. 
  mol.move_hydrogens_to_end_of_connection_table();
}

// A common task when breaking molecules is to need to know which atoms
// are on either side of a bond breakage.

void
DemoIdentifySideOfBond() {
  Molecule mol = MolFromSmiles("H[He][Li][Be]BCNOF[Ne]");

  const int matoms = mol.natoms();
  std::unique_ptr<int[]> sides = std::make_unique<int[]>(matoms);
  std::fill_n(sides.get(), matoms, 0);

  // We define two sides of the bond. In this case, the Boron and the
  // Carbon atoms.

  atom_number_t a1 = 4;
  atom_number_t a2 = 5;

  mol.identify_side_of_bond(sides.get(), a1, 1, a2);

  Molecule mol1, mol2;
  mol.create_subset(mol1, sides.get(), 0);
  mol.create_subset(mol2, sides.get(), 1);
  cerr << "Parent " << mol.smiles() << '\n';
  cerr << "Parts " << mol1.smiles() << ' ' << mol2.smiles() << '\n';
}

// An atom has a user specified pointer that can point to anything.
// Most of the time, it is not necessary to store information with atoms,
// but during complex processing, where atom numbers might be changing
// or where the provenance of atoms might be needed, storing something
// with the atom can be convenient. This keeps the Atom lightweight.

// The Molecule/Atom does NOT assume ownership of these pointers.
// Today this would likely be implemeted as a shared_ptr<void>.

// The molecule also has its own user_specified_void_ptr, not shown here.
void
DemoUserSpecifiedVoidPointer() {
  // Structures of arbitrary complexity can be used.
  // Or just a pointer to a numeric value that is retained in scope.
  struct Something {
    float number;
  };

  // Allocate a couple of these, and assign arbitrary numbers.
  std::unique_ptr<Something> s1 = std::make_unique<Something>();
  std::unique_ptr<Something> s2 = std::make_unique<Something>();
  s1->number = 1.0;
  s2->number = 2.0;

  Molecule mol = MolFromSmiles("CC");

  mol.set_user_specified_atom_void_ptr(0, reinterpret_cast<void*>(s1.get()));
  mol.set_user_specified_atom_void_ptr(1, reinterpret_cast<void*>(s2.get()));

  mol.swap_atoms(0, 1);

  // Since the atomic user
  for (const Atom* atom : mol) {
    const Something* s = reinterpret_cast<const Something*>(atom->user_specified_void_ptr());
    cerr << "Value is " << s->number << '\n';
  }

  // Note that during copy constructors, the user_specified void pointers are copied.

  Molecule mol2(mol);
  cerr << "In a copy\n";
  for (const Atom* atom : mol2) {
    const Something* s = reinterpret_cast<const Something*>(atom->user_specified_void_ptr());
    cerr << "Value is " << s->number << '\n';
  }
}

void
DemoOutput() {
  Molecule mol;
  mol.build_from_smiles("C");
  mol.set_name("methane");

  mol.write_molecule(std::cout, FILE_TYPE_SMI, "");
  // which in this case is the same as
  std::cout << mol.smiles() << ' ' << mol.name() << '\n';
}

// There are several methods that pass functions or lambda objects
// to a Molecule. Sometimes these are useful for the compactness
// they can enable, other times a function might work better.

// The each_atom method requires a struct that responds to operator()(const Atom&);
// This is fairly complex, and probably not where you should start.
// This may make sense if you had some complex processing to be
// done for each atom. But that could also be done in a function.
// Modern C++ does a lot with iterators like this, so....
// What you do get from this is a single line invocation that
// can be re-used.
void
DemoEachAtom() {
  Molecule mol = MolFromSmiles("OCN");

  // Sum the atomic_number() value for each atom.
  struct SumNcon {
    int sum = 0;
    void operator()(const Atom& a) {
      sum += a.atomic_number();
    }
  };

  SumNcon sum_ncon;
  mol.each_atom(sum_ncon);
  cerr << "Sum of ncon " << sum_ncon.sum << '\n';

  // Same as:
  int sum = 0;
  for (const Atom* a : mol) {
    sum += a->ncon();
  }
  cerr << "Via traditonal means " << sum << '\n';
}

// Find the highest atomic number in the molecule.
// Again, a complex case.
void
DemoEachAtomLambda() {
  Molecule mol = MolFromSmiles("OCN");

  int max_atomic_number = 0;
  mol.each_atom_lambda([&max_atomic_number] (const Atom& atom) {
    max_atomic_number = std::max(max_atomic_number, atom.atomic_number());
  });
  cerr << "highest atomic number " << max_atomic_number << '\n';
}

// Iterate through all bonds(edges) in the molecule.
void
DemoEachBond() {
  Molecule mol = MolFromSmiles("CC=CC");

  struct CountSingleBonds {
    int single_bonds = 0;
    void operator() (const Bond& b) {
      if (b.is_single_bond()) {
        ++single_bonds;
      }
    }
  };

  CountSingleBonds count_single_bonds;

  mol.each_bond(count_single_bonds);
  cerr << "Find " <<count_single_bonds.single_bonds << " single bonds\n";

  // This is the same as.
  int single_bonds = 0;
  for (const Bond* bond : mol.bond_list()) {
    if (bond->is_single_bond()) {
      ++single_bonds;
    }
  }
}

// Iterators through rings.
void
DemoEachRing() {
  Molecule mol = MolFromSmiles("C1CC1C1CCC1C1CCCC1");

  struct LargestRingSize {
    int largest_ring = 0;
    void operator() (const Ring& r) {
      largest_ring = std::max(r.number_elements(), largest_ring);
    }
  };

  LargestRingSize largest_ring;

  mol.each_ring(largest_ring);
  cerr << "Find " << largest_ring.largest_ring << " largest ring\n";

  // This is the same as.
  int big_ring = 0;
  for (const Ring* r : mol.sssr_rings()) {
    big_ring = std::max(r->number_elements(), big_ring);
  }
}

// Apply a Molecule member function to each atom.
// This particular demo is quite unusual, normally some
// kind of numeric accumulation would be done. THis
// works because IWString has a += operator.
void
DemoEachIndex() {
  Molecule mol = MolFromSmiles("OCN");

  Molecule::const_member_fn<const IWString&> fn = &Molecule::atomic_symbol;
  IWString result = mol.each_index<IWString>(fn);
  cerr << "Concatenated symbols '" << result << "'\n";
}

int
Main(int argc, char** argv) {
  DemoEmptyMolecule();
  DemoStringConcat();
  DemoCannotParseBadSmiles();
  DemoOptionalMolFromSmiles();
  DemoAtomicNumbers();
  DemoElements();
  DemoAtomIteration();
  DemoBondIteration();
  DemoAssignment();
  DemoNatoms();
  DemoIsotopes();
  DemoIsotopesMolecularWeight();
  DemoEachAtom();
  DemoEachAtomLambda();
  DemoEachBond();
  DemoEachRing();
  DemoEachIndex();
  DemoAtomMapNumbers();
  DemoImplicitAndExplicitHydrogens();
  DemoImplicitHydrogensKnown();
  DemoSetOfAtoms();
  DemoRandomSmiles();
  DemoRemoveAtom();
  DemoFormalCharges();
  DemoPartialCharges();
  DemoNcon();
  DemoBondIterationAroundAtom();
  DemoBondIterationMolecule();
  DemoConnections();
  DemoAddMolecule();
  DemoChangeBond();
  DemoChangeElement();
  DemoFragments();
  DemoCreateComponents();
  DemoReduceToLargestFragment();
  DemoCreateSubset();
  DemoDistanceMatrix();
  DemoRing();
  DemoRingIteration ();
  DemoRingSystem();
  DemoAreBonded();
  DemoCoordinatesInAtoms();
  DemoScaffold();
  DemoSubstuctureSearch();
  DemoSubstuctureSearchSaveResults();
  DemoParsimoniousSearch();
  DemoSubstuctureSearchWithTarget();
  DemoMultipleQueries();
  DemoQueryFromProto();
  DemoStrangeElements();
  DemoAnyLengthElement();
  DemoAromaticity();
  DemoSymmetry();
  DemoBondSymmetry();
  DemoAttachedHeteroatomCount();
  DemoMultipleBondToHeteroatom();
  DemoChirality();
  DemoSwapAtoms();
  DemoIdentifySideOfBond();
  DemoUserSpecifiedVoidPointer();
  DemoOutput();

  return 0;
}

}  // namespace lillymol_introcution

int
main(int argc, char** argv) {
  int rc = lillymol_introcution::Main(argc, argv);

  return rc;
}
