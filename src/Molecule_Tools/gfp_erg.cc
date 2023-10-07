// Generate fingerprints based on ErG reduced graphs.
// Mostly written by Nikolaus Stiefl, see
// See https://pubs.acs.org/doi/10.1021/ci050457y
// ErG:  2D Pharmacophore Descriptions for Scaffold Hopping
// Nikolaus Stiefl, Ian A. Watson, Knut Baumann, and Andrea Zaliani

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#define IWARAY_EACH_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

static isotope_t Ar_Isotope = 5;
static isotope_t Al_Isotope = 4;
static isotope_t PosCh_Isotope = 7;
static isotope_t NegCh_Isotope = 6;
#define NIK_ACCEPTOR_ISOTOPIC_LABEL 2
#define NIK_DUAL_ISOTOPIC_LABEL 1
static int Acc_Label_min_1 = NIK_ACCEPTOR_ISOTOPIC_LABEL - 1;
static int Don_Label_min_1 = DEFAULT_DONOR_ISOTOPIC_LABEL - 1;

static int NUMBER_OF_PROPERTIES = 6;
static int LFR = 6;
static int SFR = 4;
// static int fuse = 100;
static float FUZZYINCR = 0.3;
static float ON = 1;
static int LOW_DIST_FCUT = 1;
static int MAX_FF = 6;

static int truncated_flip_flops = 0;

// MAX_DISTANCE is always one larger than the real max_distance wanted because of
// the 0 start of arrays in C++
static int MAX_DISTANCE = 16;
static int MAX_STERIC_DISTANCE = 21;
static int MAX_DISTANCE_MINUS_ONE = MAX_DISTANCE - 1;
static int MAX_STERIC_DISTANCE_MINUS_ONE = MAX_STERIC_DISTANCE - 1;
static int VECTOR_SIZE =
    ((NUMBER_OF_PROPERTIES * (NUMBER_OF_PROPERTIES + 1)) / 2) * MAX_DISTANCE +
    MAX_STERIC_DISTANCE;
static int P_VECTOR_SIZE = VECTOR_SIZE - MAX_STERIC_DISTANCE;

static int ADSetOff = MAX_DISTANCE * 1;
static int DDSetOff = NUMBER_OF_PROPERTIES * MAX_DISTANCE;
static int StericSetOff =
    (((NUMBER_OF_PROPERTIES * (NUMBER_OF_PROPERTIES + 1)) / 2) * MAX_DISTANCE);

static int molecules_read = 0;

static int calc_set_theo = 0;

static Chemical_Standardisation chemical_standardisation;

static Charge_Assigner charge_assigner;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static resizable_array<int> property_offset;

// define them outside to reduce the overhead
static resizable_array<int> isotopes_to_handle_i;
static resizable_array<int> isotopes_to_handle_j;

static int work_as_filter = 0;

static std::ofstream unconvert_log;

const char* prog_name;

static int verbose = 0;

static int write_abstract_form = 0;

static int compress_consecutive_ch2 = 0;

static int one_fingerprint_per_molecule = 0;

static int flatten_fingerprints_to_01 = 0;

static Sparse_Fingerprint_Creator shared_sparse_fingerprint_creator;

static int replicates_in_shared_fingerprint_creator = 0;

static int normalise_fingerprint_counts = 0;

/*
  Handy when doing database loads
*/

static int write_unique_smiles_as_identifier = 0;

/*
  The smiles (unique or otherwise) of the stating molecule.
*/

static IWString origSmiles;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fingerprint_tag;

static IWString molecular_fingerprint_tag;

// we tend to get saturation with smaller bit vector sizes
static int molecular_fingerprint_nbits = 4096;

static int fingerprints_written_this_molecule = 0;

static extending_resizable_array<int> fingerprints_per_molecule;

static IWString descriptor_name_prefix = "erg_";

static Element_Transformations etrans;

static IWString convert_cf3_to_chloro = 0;
static IWString convert_tbutyl_to_chloro = 0;

static int min_ring_size_for_aliphatic_abstract_ring = 4;

static int convert_so2_to_sc = 0;

static void
preprocess(Molecule& m)
{
  m.remove_all_chiral_centres();

  m.transform_to_non_isotopic_form();

  m.reduce_to_largest_fragment();

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  m.revert_all_directional_bonds_to_non_directional();

  if (etrans.active()) {
    etrans.process(m);
  }

  return;
}

#ifdef NOT_BEING_USED
static int
do_convert_so2_to_sc(Molecule& m)
{
  int rc = 0;

  Set_of_Atoms to_be_deleted;

  for (auto i = m.natoms() - 1; i >= 0; i--) {
    const Atom* s = m.atomi(i);

    if (16 != s->atomic_number()) {
      continue;
    }

    if (4 != s->ncon()) {
      continue;
    }

    if (6 != s->nbonds()) {
      continue;
    }

    Set_of_Atoms oxygens;

    for (auto j = 0; j < 4; ++j) {
      const Bond* b = s->item(j);

      if (!b->is_double_bond()) {
        continue;
      }

      const atom_number_t o = b->other(i);

      if (8 != m.atomic_number(o)) {
        continue;
      }

      oxygens.add(o);
    }

    if (2 != oxygens.number_elements()) {
      continue;
    }

    m.remove_bonds_to_atom(oxygens[0]);
    m.remove_bonds_to_atom(oxygens[1]);
    to_be_deleted += oxygens;
    m.set_atomic_number(i, 21);
    rc++;
  }

  if (to_be_deleted.number_elements()) {
    m.remove_atoms(to_be_deleted);
  }

  return rc;
}
#endif

/*
  Works for CF3 and t_butyl
*/

static int
do_convert_cf3_to_chloro(Molecule& m, const atomic_number_t singly_attached_atomic_number)
{
  const auto matoms = m.natoms();

  Set_of_Atoms to_be_removed;

  for (auto i = 0; i < matoms; ++i) {
    const Atom* c = m.atomi(i);

    if (6 != c->atomic_number()) {
      continue;
    }

    if (4 != c->ncon()) {
      continue;
    }

    Set_of_Atoms attached;

    for (auto j = 0; j < 4; ++j) {
      const Bond* b = c->item(j);

      if (!b->is_single_bond()) {
        continue;
      }

      atom_number_t k = b->other(i);

      if (singly_attached_atomic_number != m.atomic_number(k)) {
        continue;
      }

      if (1 != m.ncon(k)) {
        continue;
      }

      attached.add(k);
    }

    if (3 != attached.number_elements()) {
      continue;
    }

    m.set_atomic_number(i, 17);
    m.set_isotope(i, 0);

    to_be_removed += attached;
  }

  if (to_be_removed.number_elements()) {
    m.remove_atoms(to_be_removed);
  }

  return 1;
}

static int
is_ch2(Molecule& m, const atom_number_t zatom, const int* ring_membership)
{
  if (ring_membership[zatom]) {
    return 0;
  }

  const Atom* a = m.atomi(zatom);

  if (6 != a->atomic_number()) {
    return 0;
  }

  if (2 != a->ncon()) {
    return 0;
  }

  if (2 != m.hcount(zatom)) {
    return 0;
  }

  return 1;
}

static int
is_ch2_or_ch3(Molecule& m, const atom_number_t zatom, const int* ring_membership)
{
  if (ring_membership[zatom]) {
    return 0;
  }

  const Atom* a = m.atomi(zatom);

  if (6 != a->atomic_number()) {
    return 0;
  }

  const int acon = a->ncon();

  if (acon > 2) {
    return 0;
  }

  if (2 == acon && 2 == m.hcount(zatom)) {
    ;
  } else if (1 == acon && 3 == m.hcount(zatom)) {
    ;
  } else {
    return 0;
  }

  return 1;
}

// #define DEBUG_COMPRESS_CONSECUTIVE

static int
do_compress_consecutive_ch2(Molecule& m)
{
#ifdef DEBUG_COMPRESS_CONSECUTIVE
  cerr << "Processing '" << m.name() << "'\n";
#endif

  const auto matoms = m.natoms();

  int* ring_membership = new int[matoms];
  std::unique_ptr<int[]> free_ring_membership(ring_membership);
  m.ring_membership(ring_membership);

  Set_of_Atoms to_be_removed;

  for (int i = m.natoms() - 1; i >= 0; i--) {
    if (!is_ch2(m, i, ring_membership)) {
      continue;
    }

    const Atom* a = m.atomi(i);

    //    lhs - I - ch2

    atom_number_t lhs = INVALID_ATOM_NUMBER;
    atom_number_t ch2 = INVALID_ATOM_NUMBER;

    for (auto j = 0; j < 2; ++j) {
      const Bond* b = a->item(j);

#ifdef DEBUG_COMPRESS_CONSECUTIVE
      cerr << "Bond from " << b->a1() << " to " << b->a2() << " single? "
           << b->is_single_bond() << endl;
#endif

      if (!b->is_single_bond()) {
        continue;
      }

      const atom_number_t k = b->other(i);

      if (INVALID_ATOM_NUMBER == ch2) {
        if (is_ch2_or_ch3(m, k, ring_membership)) {
          ch2 = k;
        } else {
          lhs = k;
        }
      } else {
        lhs = k;
      }
    }

#ifdef DEBUG_COMPRESS_CONSECUTIVE
    cerr << "Atom " << i << " " << m.smarts_equivalent_for_atom(i) << " is CH2, left "
         << lhs << " other ch2 " << ch2 << endl;
#endif

    if (INVALID_ATOM_NUMBER == ch2 || INVALID_ATOM_NUMBER == lhs) {
      continue;
    }

    to_be_removed.add(i);
    m.remove_bonds_to_atom(i);
    m.add_bond(lhs, ch2, SINGLE_BOND);
  }

  const auto nr = to_be_removed.number_elements();

  if (0 == nr) {
    return 0;
  }

  for (auto i = 0; i < nr; ++i) {
    m.remove_atom(to_be_removed[i]);
  }

  return nr;
}

static const Element* element_ar = nullptr;

// changed from IW's Al defintion to Hf - just in case there is an Aluminium
// in the DB. Hf is for Hydrof(!)obic

static const Element* element_al = nullptr;

// dummy atom type Dysprosium - hopefully there is never a Dysprosium containing
// compound in the data set ;-)

static const Element* element_dy = nullptr;

static int
add_aromatic_atom(Molecule& m)
{
  Atom* a = new Atom(element_ar);
  a->set_isotope(Ar_Isotope);
  m.add(a);

  return (1);
}

static int
add_aliphatic_atom(Molecule& m)
{
  Atom* a = new Atom(element_al);
  a->set_isotope(Al_Isotope);
  m.add(a);

  return (1);
}

static int
place_abstract_ring(const Molecule& m, const int* ring_membership, const Ring& r,
                    Molecule& m2, int* atoms_to_be_deleted, int* convert_to_dummy_atoms)
{
  atom_number_t attachment_point = m2.natoms();

  // if ring size is between 4 and 7 check if more than 50 % of the ring atoms
  // are sp2 hybridised - if so ad aromatic ... (e.g. Cox-2)
  int ring_size = r.number_elements();

  if (r.is_aromatic()) {
    add_aromatic_atom(m2);
  } else if ((ring_size >= min_ring_size_for_aliphatic_abstract_ring) &&
             (ring_size < 8)) {
    int flat_atoms = 0;
    for (int i = 0; i < ring_size; i++) {
      if (m.atomi(r[i])->ncon() < m.atomi(r[i])->nbonds()) {
        flat_atoms++;
      }
    }
    if (flat_atoms > ring_size / 2) {
      add_aromatic_atom(m2);
    }
    // added this to make sure that aliphatic rings fused to aromatic rings
    // are handled as one entity - here both as aromatic
    else if (r.is_fused()) {
      bool isAromatic(false);
      for (int i = 0; i < r.fused_ring_neighbours(); i++) {
        if ((r.fused_neighbour(i))->is_aromatic()) {
          isAromatic = true;
          break;
        }
      }
      if (isAromatic == true) {
        add_aromatic_atom(m2);
      } else {
        add_aliphatic_atom(m2);
      }
    } else {
      add_aliphatic_atom(m2);
    }
  } else {
    if (verbose > 1) {
      cerr << "Abstract ring not added, size " << ring_size << endl;
    }
    return 0;
  }

  for (int i = 0; i < ring_size; i++) {
    atom_number_t j = r[i];

    // nik - keep donor or acceptor atoms and attach to new Ar/Hf
    atoms_to_be_deleted[j] = 1;

    const Atom* aj = m.atomi(j);

    int jcon = aj->ncon();

    if ((m.isotope(j) > 0) &&
        (m.isotope(j) < static_cast<isotope_t>((NUMBER_OF_PROPERTIES + 2)))) {
      m2.add_bond(attachment_point, j, SINGLE_BOND);
      atoms_to_be_deleted[j] = 0;
    } else if (2 != jcon) {
      for (int k = 0; k < jcon; k++) {
        const Bond* b = aj->item(k);

        atom_number_t o = b->other(j);

        if (r.contains(o)) {
          continue;
        }

        // use attachment atom of ring as dummy atom for Ar atom
        if (b->is_aromatic()) {
          m2.add_bond(attachment_point, j,
                      DOUBLE_BOND);  // B is part of M, so it is still valid even if the
                                     // corresponding bond is removed from M2
        } else {
          m2.add_bond(attachment_point, j,
                      b->btype());  // B is part of M, so it is still valid even if the
                                    // corresponding bond is removed from M2
        }
        //      cerr << "After bond addition " << m2.smiles() << endl;
        // change attachment point to atom type dummy
        convert_to_dummy_atoms[j] = 1;
        atoms_to_be_deleted[j] = 0;
        // when reaching here, a bond was formed between this j and the attachment
        // point - get out of for loop for j's with 4 connections (e.g. sulphur)
        break;
      }
    }
  }  // end for i

  return 1;
}

// #define DEBUG_ERG

static int
aromatic_rings_to_ar_aliphatic_rings_to_al(Molecule& m, const int* ring_membership,
                                           Molecule& m2, int* atoms_to_be_deleted,
                                           int* ring_is_atom, int* convert_to_dummy_atoms)
{
  int nr = m.nrings();

  int* ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  m.compute_aromaticity_if_needed();

  int rings_processed = 0;

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    rings_processed++;

    ring_already_done[i] = 1;

    if (place_abstract_ring(m, ring_membership, *ri, m2, atoms_to_be_deleted,
                            convert_to_dummy_atoms)) {
      ring_is_atom[ri->ring_number()] = m2.natoms() - 1;
    }
  }

  const auto nb = m.nedges();

  m.compute_aromaticity_if_needed();

#ifdef DEBUG_ERG
  cerr << m.smiles() << endl;
#endif

  for (auto i = 0; i < nb; ++i) {
    const auto b = m.bondi(i);

    //  cerr << "bond " << *b << ", aromatic ? " << b->is_aromatic() << endl;

    if (!b->is_aromatic()) {
      continue;
    }

    m2.set_bond_type_between_atoms(b->a1(), b->a2(), DOUBLE_BOND);
  }

  if (nr == rings_processed) {  // finished!
    return nr;
  }

  cerr << m.name() << ": Processed " << rings_processed << " of " << nr
       << " rings: " << m.smiles() << endl;

  return 1;
}

static int
aromatic_rings_to_ar_aliphatic_rings_to_al(Molecule& m, const int* ring_membership,
                                           Molecule& m2, int* atoms_to_be_deleted,
                                           int* convert_to_dummy_atoms)
{
  const int matoms = m.natoms();

  const int nr = m.nrings();

  int rs;
  if (nr == 1) {
    rs = 2;
  } else {
    rs = nr * (nr - 1) / 2;
  }

  // For each ring, an indication of whether or not it is joined to some other ring entity

  int* joined_to_something = new_int(nr);
  std::unique_ptr<int[]> free_joined_to_something(joined_to_something);

  int* ring_to_delete = new_int(rs);
  std::unique_ptr<int[]> free_ring_to_delete(ring_to_delete);

  bond_type_t* bond_between_ring = new bond_type_t[nr * nr];
  std::unique_ptr<bond_type_t[]> free_bond_between_ring(bond_between_ring);

  set_vector(bond_between_ring, nr * nr,
             static_cast<bond_type_t>(
                 NOT_A_BOND));  // we'll use INVALID_BOND_TYPE for spiro fusions

  int* tmp = new int[matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  int bonds_between_rings_present = 0;
  int reduce_to_one_centroid = 0;
  int reduce_this_centroid = 0;

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    assert(i == ri->ring_number());

    set_vector(tmp, matoms, 0);

    ri->set_vector(tmp, 1);

    for (int j = i + 1; j < nr; j++) {
      const Ring* rj = m.ringi(j);

      int c = rj->count_members_set_in_array(tmp, 1);

      bond_type_t bt = NOT_A_BOND;

      if (1 == c) {  // spiro fused
        bt = INVALID_BOND_TYPE;
      } else if (2 == c) {  // single bond shared
        bt = DOUBLE_BOND;
      } else if (c == 0) {
        continue;
      } else  // strongly fused ring system
      {
        reduce_to_one_centroid = 1;
        reduce_this_centroid = 1;

        // check if we want to keep both centroids
        if (ri->number_elements() > rj->number_elements()) {
          if ((ri->number_elements() > LFR) && (rj->number_elements() > SFR)) {
            reduce_to_one_centroid = 0;
            reduce_this_centroid = 0;
          }
        } else {
          if ((rj->number_elements() > LFR) && (ri->number_elements() > SFR)) {
            reduce_to_one_centroid = 0;
            reduce_this_centroid = 0;
          }
        }
        bt = DOUBLE_BOND;
      }

      if (reduce_this_centroid) {
        ring_to_delete[(((2 * nr * i) - (i) - (i * i)) / 2) + j - 1 - i] = 1;
        reduce_this_centroid = 0;
      }

      joined_to_something[i] = 1;
      joined_to_something[j] = 1;  // that way no additional vector is needed

      bond_between_ring[i * nr + j] = bt;

      bonds_between_rings_present++;
    }
  }

  // Cross reference between ring numbers and abstract atoms placed

  int* ring_is_atom = new_int(nr, -1);
  std::unique_ptr<int[]> free_ring_is_atom(ring_is_atom);
  int* ring_deleted = new_int(nr);
  std::unique_ptr<int[]> free_ring_deleted(ring_deleted);

  int rc = aromatic_rings_to_ar_aliphatic_rings_to_al(
      m, ring_membership, m2, atoms_to_be_deleted, ring_is_atom, convert_to_dummy_atoms);

  if (0 == rc) {  // strange
    return 0;
  }

  if (0 == bonds_between_rings_present) {
    return rc;
  }

  // in case there are highly fused ring systems delete the unimportant rings and
  // reconnect atoms
  if (reduce_to_one_centroid) {
    // delete strongly fused rings which are too small
    for (int i = 0; i < nr; i++) {
      for (int j = i + 1; j < nr; j++) {
        if (ring_is_atom[j] < 0) {  // unconverted ring
          continue;
        }

        int* redirect_atoms = new_int(m.natoms());
        std::unique_ptr<int[]> free_redirect_atoms(redirect_atoms);

        if ((ring_to_delete[(((2 * nr * i) - (i) - (i * i)) / 2) + j - 1 - i]) &&
            (!ring_deleted[j])) {
          // find the atoms connected to j
          int matoms = m.natoms();
          for (int k = 0; k < matoms; k++) {
            if (m2.are_bonded(k, ring_is_atom[j])) {
              redirect_atoms[k] = 1;
            }
          }
          // check if deleted ring has Ar attribute. If so, set the ring kept to
          // aromatic
          if (m2.isotope(ring_is_atom[j]) == Ar_Isotope) {
            m2.set_element(ring_is_atom[i], element_ar);
            m2.set_isotope(ring_is_atom[i], Ar_Isotope);
          }
          // end check
          m2.remove_atom(ring_is_atom[j]);
          for (int k = j + 1; k < nr; k++) {
            ring_is_atom[k]--;
          }
          ring_is_atom[j] = ring_is_atom[i];

          ring_deleted[j] = i + 1;

          // redirect the atoms connected to j to i
          for (int k = 0; k < matoms; k++) {
            // in case we have a connection to two rings - check if not already
            // bonded
            if ((redirect_atoms[k]) && (!(m2.are_bonded(k, ring_is_atom[i])))) {
              m2.add_bond(ring_is_atom[i], k, SINGLE_BOND);
            }
          }
          joined_to_something[j] = 0;

          // connect separated rings
          bond_between_ring[i * nr + j] = NOT_A_BOND;
          for (int k = 0; k < nr; k++) {
            if (k != i) {
              if (k < j) {
                if (bond_between_ring[k * nr + j] != NOT_A_BOND) {
                  bond_between_ring[k * nr + j] = NOT_A_BOND;
                  if (k < i) {
                    if (ring_deleted[i] > 0) {
                      bond_between_ring[k * nr + (ring_deleted[i] - 1)] = DOUBLE_BOND;
                    } else {
                      bond_between_ring[k * nr + i] = DOUBLE_BOND;
                    }
                  } else if (ring_deleted[i] > 0) {
                    bond_between_ring[(ring_deleted[i] - 1) * nr + k] = DOUBLE_BOND;
                  } else {
                    bond_between_ring[i * nr + k] = DOUBLE_BOND;
                  }
                }
              } else {
                if (bond_between_ring[j * nr + k] != NOT_A_BOND) {
                  bond_between_ring[j * nr + k] = NOT_A_BOND;
                  if (k < i) {
                    if (ring_deleted[i] > 0) {
                      bond_between_ring[k * nr + (ring_deleted[i] - 1)] = DOUBLE_BOND;
                    } else {
                      bond_between_ring[k * nr + i] = DOUBLE_BOND;
                    }
                  } else if (ring_deleted[i] > 0) {
                    bond_between_ring[(ring_deleted[i] - 1) * nr + k] = DOUBLE_BOND;
                  } else {
                    bond_between_ring[i * nr + k] = DOUBLE_BOND;
                  }
                }
              }
            }
          }
        }
      }  // end for j
    }    // end for i
  }      // end if reduce_one_centroid

  return rc;
}

/*
  Interesting test molecule, OC1=C(CCC(C)C)C(O)=CC(O)=C1C(=O)CC(C)C PBCHM13411747
  Check Dy atoms sometime
*/

static int
aromatic_rings_to_ar_aliphatic_rings_to_al(Molecule& m, int* ring_membership,
                                           Molecule& m2)
{
  const int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int matoms = m.natoms();
  int matoms_nr = matoms + nr;

  // First task is to figure out which rings are joined to each other and how

  int* atoms_to_be_deleted = new_int(matoms_nr);
  std::unique_ptr<int[]> free_atoms_to_be_deleted(atoms_to_be_deleted);
  int* convert_to_dummy_atoms = new_int(matoms_nr);
  std::unique_ptr<int[]> free_convert_to_dummy_atoms(convert_to_dummy_atoms);

  int rc = aromatic_rings_to_ar_aliphatic_rings_to_al(
      m, ring_membership, m2, atoms_to_be_deleted, convert_to_dummy_atoms);
  // cerr << " m2 " << m2.smiles() << endl;

  if (rc) {
    // not nice but efficient
    for (int i = matoms - 1; i >= 0; i--) {
      if (convert_to_dummy_atoms[i]) {
        m2.set_element(i, element_dy);
      }
      // in strongly fused ring systems it can happen that an atom
      // is member of three rings! for such cases, it can occur that an atom is
      // dummy and to_be_deleted - for such cases we don't want the deletion, but
      // the ring closure - therefore: else if
      else if (atoms_to_be_deleted[i]) {
        m2.remove_atom(i);
      }
    }
  }
  return rc;
}

static int
do_produce_abstract_ring_forms(Molecule& m)
{
  int matoms = m.natoms();

  int* ring_membership = new int[matoms];
  std::unique_ptr<int[]> free_ring_membership(ring_membership);

  m.ring_membership(ring_membership);

  Molecule mcopy(m);
  mcopy.set_name(m.name());

  aromatic_rings_to_ar_aliphatic_rings_to_al(m, ring_membership, mcopy);
  // copying mcopy to m to process in next step
  m = mcopy;

  return 1;
}

static void
usage(int rc = 1)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "usage: " << prog_name << " -i <input type> file1 file2...\n";
  cerr << "The following options are recognised\n";
  cerr << "  -F <fIncr>     specify fuzzy increment (default:0.5).\n";
  cerr << "  -X <maxFF>     specify maximum number of increments for enumeration"
          "(default: 5).\n";
  cerr << "  -K             write header\n";
  cerr << "  -s             write set_theoretic form\n";
  cerr << "  -d <type>      specify maximum interproperty distance (default:15)\n";
  cerr << "  -U <fname>     write unconvertable molecules to <fname>\n";
  cerr << "  -L <fname>     write rejected molecules to <fname>\n";
  cerr << "  -m             write abstract molecular form instead of descriptor\n";
  cerr << "  -B ...         various abstraction options, enter '-B help' for info\n";
  cerr << "  -c             compress consecutive chain CH2 groups\n";
  cerr << "  -h             convert CF3 to Cl\n";
  cerr << "  -y             convert O=S=O to Sc\n";
  cerr << "  -q <query>     specify endcap query file\n";
  cerr << "  -q F:file      specify file of endcap queries\n";
  cerr << "  -f             work as a tdt filter\n";
  cerr << "  -J <tag>       generate fingerprints\n";
  cerr << "  -j             one (OR'd) fingerprint per molecule\n";
  cerr << "  -M <tag>       generate path based fingerprints of the reduced graph\n";
  cerr << "  -n             normalise fingerprint counts for flip-flop replicates (turns on -j)\n";
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  (void) display_standard_charge_assigner_options (cerr, 'N');
  cerr << "  -H <..>        donor acceptor assignment, enter '-H help' for info\n";
  cerr << "  -u             write the unique smiles as the identifier\n";
  cerr << "  -T El=El       standard element transformation options\n";
  cerr << "  -z <natoms>    min ring size for converting aliphatic ring to abstract form\n";
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
find_endcaps(Molecule& m, resizable_array_p<Substructure_Query>& endcap_queries)
{
  Molecule_to_Match target(&m);

  int nq = endcap_queries.number_elements();

  int nmatches = 0;
  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;
    nmatches = endcap_queries[i]->substructure_search(target, sresults);

    for (int j = 0; j < nmatches; j++) {
      const Set_of_Atoms* embedding = sresults.embedding(j);
      int na =
          embedding->number_elements();  // how many of the matched atoms do we process

      for (int k = 0; k < na; k++) {
        atom_number_t an = embedding->item(k);
        if (m.isotope(an) != 0) {
          cerr << "atom is already assigned isotope " << m.isotope(an) << ": "
               << m.smiles() << endl;
        }

        m.set_isotope(an, Al_Isotope);
      }
    }
  }

  return 1;
}

static int
change_charged_labels(Molecule& m)
{
  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++) {
    formal_charge_t qi = m.formal_charge(i);

    if (0 == qi) {
      ;
    } else if (qi > 0) {
      m.set_isotope(i, PosCh_Isotope);
    } else if (qi < 0) {
      m.set_isotope(i, NegCh_Isotope);
    }
  }

  return (1);
}

static int
generate_property_vector(Molecule& m, resizable_array<float>& mol_vector)
{
  int labeli(0);
  int labelj(0);
  int distance(0);
  int field(0);

  int natoms = m.natoms();

  for (int i = 0; i < natoms; i++) {
    labeli = m.isotope(i);

    if (0 == labeli) {
      continue;
    }

    // reduce label by 1 to handle isotopic labelling of FlipFlop
    labeli--;

    isotopes_to_handle_i.resize_keep_storage(0);

    if (labeli < 5) {
      isotopes_to_handle_i.add(labeli);
      // handling zero distance
      mol_vector[property_offset[labeli - 1] * MAX_DISTANCE]++;
    } else if (5 == labeli)  // acceptor and negative
    {
      // handling zero distance for both properties and mix
      mol_vector[property_offset[0] * MAX_DISTANCE]++;
      mol_vector[(property_offset[0] + 4) * MAX_DISTANCE]++;
      mol_vector[property_offset[4] * MAX_DISTANCE]++;

      isotopes_to_handle_i.add(Acc_Label_min_1);
      isotopes_to_handle_i.add(labeli);
    }
    // then it must be 6 -> donor and positive
    else {
      // handling zero distance for both properties and mix
      mol_vector[property_offset[1] * MAX_DISTANCE]++;
      mol_vector[(property_offset[1] + 4) * MAX_DISTANCE]++;
      mol_vector[property_offset[5] * MAX_DISTANCE]++;

      isotopes_to_handle_i.add(Don_Label_min_1);
      isotopes_to_handle_i.add(labeli);
    }

    // j starts at i+1 since zero distance is incremented already!
    for (int j = i + 1; j < natoms; j++) {
      labelj = m.isotope(j);
      if (0 == labelj) {
        continue;
      }

      labelj--;
      // minus one because of array offset (starting from 0)
      distance = m.bonds_between(i, j);

#ifdef DEBUG_ERG
      cerr << "Processing atoms " << i << " label " << labeli << " and " << j << " label "
           << labelj << " dist " << distance << endl;
#endif

      if (distance > MAX_DISTANCE_MINUS_ONE) {
        continue;
      }

      // changed the isotopic label for FlipFlop from Donor_acceptor.h!!
      // DUAL now 1 ACC 2
      //
      // now if FF then both labeli and labelj are 0
      // Acceptor labels are not included in calculatios since they are 0!!

      isotopes_to_handle_j.resize_keep_storage(0);

      if (labelj < 5) {
        isotopes_to_handle_j.add(labelj);
      } else if (labelj == 5)  // acceptor and negative
      {
        isotopes_to_handle_j.add(Acc_Label_min_1);
        isotopes_to_handle_j.add(labelj);
      }
      // then it must be 6 -> donor and positive
      else {
        isotopes_to_handle_j.add(Don_Label_min_1);
        isotopes_to_handle_j.add(labelj);
      }

      // generic approach
      if (calc_set_theo) {
        for (unsigned int iIso = 0; iIso < isotopes_to_handle_i.size(); iIso++) {
          for (unsigned int jIso = 0; jIso < isotopes_to_handle_j.size(); jIso++) {
            // either exchange labeli and labelj or do it this way ...
            if (isotopes_to_handle_i[iIso] < isotopes_to_handle_j[jIso]) {
              field = (property_offset[isotopes_to_handle_i[iIso] - 1] +
                       (isotopes_to_handle_j[jIso] - isotopes_to_handle_i[iIso])) *
                          MAX_DISTANCE +
                      distance;
            } else {
              field = (property_offset[isotopes_to_handle_j[jIso] - 1] +
                       (isotopes_to_handle_i[iIso] - isotopes_to_handle_j[jIso])) *
                          MAX_DISTANCE +
                      distance;
            }

            mol_vector[field] = ON;

            if (distance) {
              if ((distance > LOW_DIST_FCUT) && (mol_vector[field - 1] != ON)) {
                mol_vector[field - 1] = FUZZYINCR;
              }

              if ((distance != MAX_DISTANCE_MINUS_ONE) && (mol_vector[field + 1] != ON)) {
                mol_vector[field + 1] = FUZZYINCR;
              }
            }
          }
        }  // iIso
      } else {
        for (unsigned int iIso = 0; iIso < isotopes_to_handle_i.size(); iIso++) {
          for (unsigned int jIso = 0; jIso < isotopes_to_handle_j.size(); jIso++) {
            // either exchange labeli and labelj or do it this way ...
            if (isotopes_to_handle_i[iIso] < isotopes_to_handle_j[jIso]) {
              field = (property_offset[isotopes_to_handle_i[iIso] - 1] +
                       (isotopes_to_handle_j[jIso] - isotopes_to_handle_i[iIso])) *
                          MAX_DISTANCE +
                      distance;
            } else {
              field = (property_offset[isotopes_to_handle_j[jIso] - 1] +
                       (isotopes_to_handle_i[iIso] - isotopes_to_handle_j[jIso])) *
                          MAX_DISTANCE +
                      distance;
            }

            mol_vector[field]++;

#ifdef DEBUG_ERG
            cerr << " isotope " << isotopes_to_handle_i[iIso] << " and "
                 << isotopes_to_handle_i[jIso] << " field " << field
                 << " valud updated to " << mol_vector[field] << endl;
#endif

            if (distance) {
              if (distance > LOW_DIST_FCUT) {
                mol_vector[field - 1] += FUZZYINCR;
              }

              if (distance != MAX_DISTANCE_MINUS_ONE) {
                mol_vector[field + 1] += FUZZYINCR;
              }
            }
          }
        }
      }
    }
  }

#ifdef DEBUG_ERG
  for (auto i = 0; i < mol_vector.number_elements(); ++i) {
    cerr << " i " << i << " mol_vector " << mol_vector[i] << endl;
  }
#endif

  return 1;
}

static int
create_fingerprint(const resizable_array<float>& mol_vector,
                   Sparse_Fingerprint_Creator& sfp)
{
  for (int i = 0; i < P_VECTOR_SIZE; i++) {
    if (mol_vector[i] <= 0.0) {
      continue;
    }

    // don't want to include zero distance into descriptor
    if ((i % MAX_DISTANCE) != 0) {
      // make unsigned int with a max of 254 ...
      unsigned int b = (i - 1 - static_cast<unsigned int>(i / MAX_DISTANCE));

#ifdef DEBUG_ERG
      cerr << "Creating fingerprint, i = " << i << " mol_vector " << mol_vector[i]
           << " bit " << b << endl;
#endif

      if ((mol_vector[i] < 25.5)) {
        sfp.hit_bit(b, static_cast<int>(mol_vector[i] * 10.0f + 0.1f));
      } else {
        sfp.hit_bit(b, 254);
      }
    }
  }

#ifdef DEBUG_ERG
  sfp.debug_print(cerr);
#endif

  return 1;
}

static int
transfer_to_sparse_fp_creator(const IWMFingerprint& fp, Sparse_Fingerprint_Creator& sfc)
{
  const int* v = fp.vector();

  const int nbits = fp.nbits();

  for (int i = 0; i < nbits; ++i) {
    if (0 == v[i]) {
      continue;
    }

    sfc.hit_bit(i, v[i]);
  }

  return 1;
}

static int
create_fingerprint_and_write(const IWMFingerprint& fp, const IWString& tag,
                             IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator sfc;

  transfer_to_sparse_fp_creator(fp, sfc);

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tmp);

  output << tag << tmp << ">\n";

  return 1;
}

static int
write_molecular_fingerprint(Molecule& m, const int* atype, const IWString& tag,
                            IWString_and_File_Descriptor& output)
{
  IWMFingerprint fp;

  fp.construct_fingerprint(m, atype, nullptr);

  // cerr << "Fingerprint has " << fp.nset() << " bits set\n";

  if (tag.starts_with("FP")) {
    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << tag << tmp << ">\n";

    return 1;
  }

  return create_fingerprint_and_write(fp, tag, output);
}

static int
create_fingerprint_and_write(const resizable_array<float>& mol_vector,
                             IWString_and_File_Descriptor& output)
{
  if (one_fingerprint_per_molecule) {
    return create_fingerprint(mol_vector, shared_sparse_fingerprint_creator);
  }

  Sparse_Fingerprint_Creator fp;
  create_fingerprint(mol_vector, fp);

  if (fingerprints_written_this_molecule > 0) {
    output << ',';
  }

  IWString tmp;
  fp.daylight_ascii_form_with_counts_encoded(tmp);
  // cerr << "Length of fingerprint " << tmp.length() << '\n';

  output << tmp;

  fingerprints_written_this_molecule++;

  return 1;
}

static int
write_ascii_form(const IWString& mname, const float* mol_vector,
                 IWString_and_File_Descriptor& output)
{
  set_default_iwstring_float_concatenation_precision(5);

  if (write_unique_smiles_as_identifier) {
    output << origSmiles;
  } else {
    append_first_token_of_name(mname, output);
  }

  for (int i = 0; i < P_VECTOR_SIZE; i++) {
    if (0 == (i % MAX_DISTANCE)) {
      continue;
    }

    if (static_cast<float>(0.0) == mol_vector[i]) {
      output << " 0";
    } else {
      output << ' ' << mol_vector[i];
    }
  }

  output << '\n';

  return 1;
}

static void
update_molecular_fingerprint(Molecule& m, const int* atype,
                             Sparse_Fingerprint_Creator& sfc)
{
  IWMFingerprint fp;

  fp.construct_fingerprint(m, atype, nullptr);

  transfer_to_sparse_fp_creator(fp, sfc);

  return;
}

static int
write_fixed_width_fingerprint(const Sparse_Fingerprint_Creator& sfc, const IWString& tag,
                              IWString& output)
{
  IW_Bits_Base fp(molecular_fingerprint_nbits);

  const auto bits_found = sfc.bits_found();

  for (auto i = bits_found.begin(); i != bits_found.end(); ++i) {
    const auto b = (*i).first;

    fp.set(b % molecular_fingerprint_nbits);
  }

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);
  output << tag << tmp << ">\n";

  return 1;
}

static int
write_molecular_fingerprint(Molecule& m, int nff, Sparse_Fingerprint_Creator& sfc,
                            const IWString& tag, IWString_and_File_Descriptor& output)
{
  if (tag.starts_with("FP")) {
    return write_fixed_width_fingerprint(sfc, tag, output);
  }

  if (normalise_fingerprint_counts && nff > 0) {
    Sparse_Fingerprint_Creator::FPHash& b =
        const_cast<Sparse_Fingerprint_Creator::FPHash&>(sfc.bits_found());

    for (auto i = b.begin(); i != b.end(); ++i) {
      int c = (*i).second;
      if (c > nff) {
        c = c / nff;
        (*i).second = c;
      }
    }
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tmp);

  output << tag << tmp << ">\n";

  return 0;
}

static void
enumerate_DFS(Molecule& m, resizable_array<float>& mol_vector,
              Sparse_Fingerprint_Creator& molecular_fingerprint, int* atype,
              IWString_and_File_Descriptor& output, resizable_array<int>& fF, int level)
{
  // cerr << "DFS at level " << level << " of " << fF.size() << "
  // one_fingerprint_per_molecule " << one_fingerprint_per_molecule << endl; DFS
  for (int i = 0; i < NIK_ACCEPTOR_ISOTOPIC_LABEL; i++) {
    m.set_isotope(fF[level], i + NIK_ACCEPTOR_ISOTOPIC_LABEL);

    if (molecular_fingerprint_tag.length()) {
      atype[fF[level]] = i + NIK_ACCEPTOR_ISOTOPIC_LABEL;
    }

    if (level == static_cast<int>(fF.size() - 1)) {
#ifdef DEBUG_ERG
      cerr << "Generating fingerprint for '" << m.unique_smiles() << "'\n";
#endif
      if (molecular_fingerprint_tag.length() > 0) {
        update_molecular_fingerprint(m, atype, molecular_fingerprint);
        continue;
      }

      generate_property_vector(m, mol_vector);

      if (fingerprint_tag.length()) {
        create_fingerprint_and_write(mol_vector, output);
        mol_vector.set_all(0.0);
        //      mol_vector.each ([] (float & x) { x = 0.0;});
      } else if (one_fingerprint_per_molecule) {  // descriptors, but we are accumulating
                                                  // them
        ;
      } else {
        write_ascii_form(m.name(), mol_vector.rawdata(), output);
        mol_vector.set_all(0.0);
      }
    } else {
      enumerate_DFS(m, mol_vector, molecular_fingerprint, atype, output, fF, level + 1);
    }
  }

  return;
}

static int
canonical_handling_of_too_many_flipflops(Molecule& m, Set_of_Atoms& flipFlops)
{
  IWString s = m.unique_smiles();

  Molecule mcopy;

  set_input_aromatic_structures(1);

  if (!mcopy.build_from_smiles(s)) {
    cerr << "Yipes, cannot interpret smiles '" << m.name() << "'\n";
    return 0;
  }

  mcopy.set_name(m.name());

  int matoms = m.natoms();

  flipFlops.resize_keep_storage(0);

  for (int i = 0; i < matoms; i++) {
    if (mcopy.isotope(i) != NIK_DUAL_ISOTOPIC_LABEL) {
      continue;
    }

    if (flipFlops.number_elements() < MAX_FF) {
      flipFlops.add(i);
    } else {
      mcopy.set_isotope(i, NIK_ACCEPTOR_ISOTOPIC_LABEL);
    }
  }

  m = mcopy;

  truncated_flip_flops++;

  return 1;
}

static void
assign_atom_types(Molecule& m, int* atype)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);

    if (iso > 0) {
      atype[i] = 2000 + iso;
    } else {
      //    atype[i] = 16 * m.atomic_number(i) + 4 * m.ncon(i) + m.hcount(i);
      atype[i] = m.atomic_number(i);
    }
  }

  return;
}

static int
generate_property_vector_fF(Molecule& m, resizable_array<float>& mol_vector,
                            IWString_and_File_Descriptor& output)
{
  // find flip-flops
  int natoms = m.natoms();
  Set_of_Atoms flipFlops(natoms);

  for (int i = 0; i < natoms; i++) {
    if (m.isotope(i) == NIK_DUAL_ISOTOPIC_LABEL) {
      flipFlops.add(i);
    }
  }

  if (flipFlops.number_elements() > MAX_FF) {
    canonical_handling_of_too_many_flipflops(m, flipFlops);
  }
  /*{
      if (unconvert_log.rdbuf()->is_open() && unconvert_log.good())
            unconvert_log << "Molecule: " << m.name() << " could not be converted! "
                  << flipFlops.size() << " flip-flops!" << "  Smiles is: " << origSmiles
    << endl; else cerr << "Molecule: " << m.name() << " could not be converted! "
                  << flipFlops.size() << " flip-flops!\n";

      return (1);
    }*/

  if (verbose > 2) {
    cerr << m.name() << " has " << flipFlops.size() << " flip-flop atoms\n";
  }

  int* atype = nullptr;

  if (molecular_fingerprint_tag.length() > 0) {
    atype = new int[m.natoms()];
    assign_atom_types(m, atype);
  }

  if (flipFlops.size() > 0) {
    Sparse_Fingerprint_Creator molecular_fingerprint;
    enumerate_DFS(m, mol_vector, molecular_fingerprint, atype, output, flipFlops, 0);

    if (molecular_fingerprint_tag.length()) {
      write_molecular_fingerprint(m, flipFlops.size(), molecular_fingerprint,
                                  molecular_fingerprint_tag, output);
    } else if (0 == fingerprint_tag.length() && one_fingerprint_per_molecule) {
      if (normalise_fingerprint_counts) {
        const auto n = mol_vector.number_elements();

        const float nf = static_cast<float>(flipFlops.size());

        for (auto i = 0; i < n; ++i) {
          const auto mvi = mol_vector[i];
          if (mvi > 0.0f) {
            mol_vector[i] = 1.0f + mvi / nf;
          }
        }
      }

      write_ascii_form(m.name(), mol_vector.rawdata(), output);
    }
  } else {
    generate_property_vector(m, mol_vector);

    if (molecular_fingerprint_tag.length()) {
      write_molecular_fingerprint(m, atype, molecular_fingerprint_tag, output);
    } else if (fingerprint_tag.length()) {
      create_fingerprint_and_write(mol_vector, output);
    } else {
      write_ascii_form(m.name(), mol_vector.rawdata(), output);
    }
  }

  if (nullptr != atype) {
    delete[] atype;
  }

  return (1);
}

static int
generate_steric_vector(Molecule& m, resizable_array<float>& mol_vector, IWString& output)
{
  int distance(0);

  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++) {
    for (int j = i; j < natoms; j++) {
      if (i == j) {
        distance = 0;
      } else {
        // minus one because of array offset (starting from 0)
        distance = m.bonds_between(i, j);
      }

      if (distance > MAX_STERIC_DISTANCE_MINUS_ONE) {
        continue;
      }

      mol_vector[StericSetOff + distance]++;
    }
  }

  return (1);
}

static int
apply_donor_acceptor_and_switch_isotopes(Molecule& m)
{
  int rc = donor_acceptor_assigner.process(m);

  if (0 == rc) {
    return 0;
  }

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int iso = m.isotope(i);
    if (0 == iso) {
      ;
    } else if (DEFAULT_DUAL_ISOTOPIC_LABEL == iso) {
      m.set_isotope(i, NIK_DUAL_ISOTOPIC_LABEL);
    } else if (DEFAULT_ACCEPTOR_ISOTOPIC_LABEL == iso) {
      m.set_isotope(i, NIK_ACCEPTOR_ISOTOPIC_LABEL);
    }
  }

  return rc;
}

static void
do_normalise_fingerprint_counts(Sparse_Fingerprint_Creator& sfp, int replicates)
{
  auto h = sfp.bits_found();

  for (auto i = h.begin(); i != h.end(); ++i) {
    int c = (*i).second;

    c = (c + replicates) / replicates;

    (*i).second = c;
  }

  return;
}

static int
gfp_erg(Molecule& m, IWString_and_File_Descriptor& output,
        resizable_array_p<Substructure_Query>& endcap_queries)
{
  preprocess(m);

  if (one_fingerprint_per_molecule) {
    shared_sparse_fingerprint_creator.clear();
    replicates_in_shared_fingerprint_creator = 0;
  }

  if (write_unique_smiles_as_identifier) {
    origSmiles = m.unique_smiles();
  } else {
    origSmiles = m.smiles();
  }

  if (charge_assigner.active()) {
    charge_assigner.process(m);
  }

  // Nik swapped the default values

  if (donor_acceptor_assigner.active()) {
    apply_donor_acceptor_and_switch_isotopes(m);
  }

  // vector with steric properties attached at end

  resizable_array<float> mol_vector(VECTOR_SIZE, 0);

  // hydrophobe detection - better before acceptor/donor assignment - you never know

  generate_steric_vector(m, mol_vector, output);

  find_endcaps(m, endcap_queries);

  // change isotopic labels of charged atoms
  // positive charged -> label 6
  // negative charged -> label 7

  change_charged_labels(m);

#ifdef DEBUG_ERG
  cerr << m.smiles() << endl;
#endif

  if (convert_cf3_to_chloro) {
    do_convert_cf3_to_chloro(m, 9);
  }

  if (convert_tbutyl_to_chloro) {
    do_convert_cf3_to_chloro(m, 6);
  }

  if (m.nrings()) {
    do_produce_abstract_ring_forms(m);
  }

#ifdef DEBUG_ERG
  cerr << m.smiles() << endl;
#endif

  if (compress_consecutive_ch2) {
    do_compress_consecutive_ch2(m);
  }

  if (1 != m.number_fragments()) {
    cerr << "Fatal error processing '" << m.name() << "' multiple fragments '"
         << m.smiles() << "'\n";
    return 0;
  }

  if (fingerprint_tag.length() || molecular_fingerprint_tag.length()) {
    if (!work_as_filter) {
      output << smiles_tag << origSmiles << ">\n";
      output << identifier_tag << m.name() << ">\n";
    }

    fingerprints_written_this_molecule = 0;
  }

  if (write_abstract_form) {
    output << m.smiles() << ' ' << m.name() << '\n';
    return 1;
  }

  if (molecular_fingerprint_tag.length()) {
    output << "ERGSMI<" << m.smiles() << ">\n";
  }

  generate_property_vector_fF(m, mol_vector, output);

  if (fingerprint_tag.length()) {
    if (one_fingerprint_per_molecule) {
      if (flatten_fingerprints_to_01 > 1) {
        shared_sparse_fingerprint_creator.flatten_to_01();
      }

      if (normalise_fingerprint_counts && replicates_in_shared_fingerprint_creator > 1) {
        do_normalise_fingerprint_counts(shared_sparse_fingerprint_creator,
                                        replicates_in_shared_fingerprint_creator);
      }

      IWString tmp;
      shared_sparse_fingerprint_creator.daylight_ascii_form_with_counts_encoded(tmp);

      output << fingerprint_tag << tmp << ">\n";
    } else {
      fingerprints_per_molecule[fingerprints_written_this_molecule]++;
    }

    if (!work_as_filter) {
      output << "|\n";
    }
  } else if (molecular_fingerprint_tag.length()) {
    if (!work_as_filter) {
      output << "|\n";
    }
  }

  return 1;
}

static int
gfp_erg(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output,
        resizable_array_p<Substructure_Query>& endcap_queries)
{
  static int header_already_written = 0;

  // write the header
  if (fingerprint_tag.length() || molecular_fingerprint_tag.length()) {
    ;
  } else if (header_already_written) {
    ;
  } else if (write_abstract_form == 0) {
    output << "name ";
    for (int ii = 1; ii <= NUMBER_OF_PROPERTIES; ii++) {
      for (int jj = ii; jj <= NUMBER_OF_PROPERTIES; jj++) {
        for (int kk = 1; kk < MAX_DISTANCE; kk++) {
          output << descriptor_name_prefix << "P" << ii << "P" << jj << "D" << kk << " ";
        }
      }
    }
    output << '\n';
  }

  header_already_written = 1;

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (!gfp_erg(*m, output, endcap_queries)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
initialise_property_offset()
{
  property_offset.resize(NUMBER_OF_PROPERTIES);

  for (int i = 0; i < NUMBER_OF_PROPERTIES; i++) {
    if (i == 0) {
      property_offset.add(i);
    } else {
      property_offset.add(property_offset[i - 1] + NUMBER_OF_PROPERTIES - i + 1);
    }
  }

#ifdef ECHO_OFFSETS
  for (int i = 0; i < NUMBER_OF_PROPERTIES; i++) {
    cerr << "property_offset[" << i << "] = " << property_offset[i] << endl;
  }
#endif

  return 1;
}

static int
gfp_erg_filter(const const_IWSubstring& buffer,
               IWString_and_File_Descriptor& output_buffer,
               resizable_array_p<Substructure_Query>& endcap_queries)
{
  Molecule m;
  if (!m.build_from_smiles(buffer)) {
    cerr << "Cannot parse smiles '" << buffer << "'\n";
    return 0;
  }

  return gfp_erg(m, output_buffer, endcap_queries);
}

static int
gfp_erg_filter(iwstring_data_source& input, IWString_and_File_Descriptor& output,
               resizable_array_p<Substructure_Query>& endcap_queries)
{
  assert(fingerprint_tag.length() > 0);

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    if (!gfp_erg_filter(buffer, output, endcap_queries)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
gfp_erg_filter(const char* fname, IWString_and_File_Descriptor& output,
               resizable_array_p<Substructure_Query>& endcap_queries)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_erg_filter(input, output, endcap_queries);
}

int
gfp_erg(const char* fname, FileType input_type, IWString_and_File_Descriptor& output,
        resizable_array_p<Substructure_Query>& endcap_queries)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return gfp_erg(input, output, endcap_queries);
}

static void
display_abstraction_specifications(std::ostream& output)
{
  output << " -B cf3          convert CF3 to Chloro\n";
  output << " -B tb           convert t-butyl to Chloro\n";
  output << " -B so2          convert O=S=O groups to a single atom\n";
  output << " -B ch2          compress consecutive chain CH2 groups\n";
  output << " -B all          activate all other conversions\n";

  exit(1);
}

static int
gfp_erg(int argc, char** argv)
{
  Command_Line cl(argc, argv, "E:A:N:vg:d:q:H:mL:i:J:jU:sF:fKX:unM:T:z:B:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  if (argc == 1) {
    usage(4);
  }

  verbose = cl.option_count('v');

  if (verbose) {
    cerr << endl << argv[0] << " compiled " << __TIME__ << " " << __DATE__ << endl;
    cerr << endl << "Default settings used in compilation:" << endl;
    cerr << "\tFuzzy increment: 0.5" << endl;
    cerr << "\tMaximum interproperty distance: 15" << endl;
    cerr << "\tAll isotopic atoms will be converted to their non isotopic form" << endl;
    cerr << "\tMolecules containing elements not in periodic table "
            "will be discarded - problems with smiles!!"
         << endl;
    cerr << "\tThe largest fragment will be retained" << endl;
    cerr << endl << "Starting to work:" << endl;
  }

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot determine chemical standardisation from command line\n";
      usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!etrans.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Cannot initialise element transformations (-T)\n";
      return 2;
    }
  }

  if (cl.option_present('z')) {
    if (!cl.value('z', min_ring_size_for_aliphatic_abstract_ring) ||
        min_ring_size_for_aliphatic_abstract_ring < 3) {
      cerr << "The min aliphatic ring size for conversion to abstract form (-z) must be "
              "a valid ring size\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Aliphatic rings with " << min_ring_size_for_aliphatic_abstract_ring
           << " or more atoms converted to asbstract form\n";
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot determine charge assigner from command line\n";
      usage(77);
    }
  } else {
    cerr << "No charge assigner specifications supplied! Abort computation!" << endl;
    usage(6);
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose)) {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      usage(6);
    }
  } else {
    cerr << "No donor/acceptor assigner specifications supplied! Abort computation!"
         << endl;
    usage(6);
  }

  if (cl.option_present('f')) {
    work_as_filter = 1;
    if (verbose) {
      cerr << "Functioning as TDT filter\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (0 != input_type) {  // great, explicitly specified
    ;
  } else if (work_as_filter) {
    ;
  } else if (1 == cl.number_elements() &&
             0 == strcmp("-", cl[0]))  // reading from a pipe, assume smiles input
  {
    cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  } else if (all_files_recognised_by_suffix(cl)) {
    ;
  } else {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (cl.option_present('j')) {
    one_fingerprint_per_molecule = 1;
    if (verbose) {
      cerr << "Will generate one fingerprint/descriptor set per molecule\n";
    }
  }

  if (cl.option_present('n')) {
    normalise_fingerprint_counts = 1;
    one_fingerprint_per_molecule = 1;

    if (verbose) {
      cerr << "Will normalise fingerprint bit counts for flip-flop replicates\n";
    }
  }

  if (cl.option_present('J') && cl.option_present('M')) {
    cerr << "Sorry, cannot produce both ERG and reduced graph fingerprints\n";
    return 1;
  }

  if (cl.option_present('J')) {
    cl.value('J', fingerprint_tag);
    if (verbose) {
      cerr << "Fingerprints written with tag '" << fingerprint_tag << "\n";
    }

    if (!fingerprint_tag.ends_with('<')) {
      fingerprint_tag << '<';
    }
  }

  if (cl.option_present('M')) {
    cl.value('M', molecular_fingerprint_tag);
    if (verbose) {
      cerr << "Molecular fingerprints written with tag '" << molecular_fingerprint_tag
           << "\n";
    }

    if (!molecular_fingerprint_tag.ends_with('<')) {
      molecular_fingerprint_tag << '<';
    }

    one_fingerprint_per_molecule = 1;

    set_iwmfingerprint_nbits(molecular_fingerprint_nbits);
  }

  if (cl.empty()) {
    cerr << prog_name << ": insufficient arguments\n";
    usage(5);
  }

  if (cl.option_present('m')) {
    write_abstract_form = 1;
  }

  if (cl.option_present('B')) {
    IWString b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      b.to_lowercase();

      if ("ch2" == b) {
        compress_consecutive_ch2 = 1;
      } else if ("so2" == b) {
        convert_so2_to_sc = 1;
      } else if ("cf3" == b) {
        convert_cf3_to_chloro = 1;
      } else if ("tb" == b) {
        convert_tbutyl_to_chloro = 1;
      } else if ("all" == b) {
        compress_consecutive_ch2 = 1;
        convert_so2_to_sc = 1;
        convert_cf3_to_chloro = 1;
        convert_tbutyl_to_chloro = 1;
      } else if ("help" == b) {
        display_abstraction_specifications(cerr);
      } else {
        cerr << "Unrecognised abstraction qualifier '" << b << "'\n";
        display_abstraction_specifications(cerr);
      }
    }
  }

  if (cl.option_present('d')) {
    if (!cl.value('d', MAX_DISTANCE) || MAX_DISTANCE < 1) {
      cerr << "The max distance value (-d) option must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "MAX_DISTANCE set to " << MAX_DISTANCE << endl;
    }

    MAX_DISTANCE_MINUS_ONE = MAX_DISTANCE;
    MAX_DISTANCE += 1;
    VECTOR_SIZE =
        ((NUMBER_OF_PROPERTIES * (NUMBER_OF_PROPERTIES + 1)) / 2) * MAX_DISTANCE +
        MAX_STERIC_DISTANCE;
    P_VECTOR_SIZE = VECTOR_SIZE - MAX_STERIC_DISTANCE;

    ADSetOff = MAX_DISTANCE * 1;
    DDSetOff = NUMBER_OF_PROPERTIES * MAX_DISTANCE;
    StericSetOff =
        (((NUMBER_OF_PROPERTIES * (NUMBER_OF_PROPERTIES + 1)) / 2) * MAX_DISTANCE);
  }

  if (cl.option_present('F')) {
    if (!cl.value('F', FUZZYINCR) || FUZZYINCR <= 0.0) {
      cerr << "The -F option requires a valid +v floating point value\n";
      usage(4);
    }

    if (verbose) {
      cerr << "FUZZYINCR set to " << FUZZYINCR << endl;
    }
  }

  if (cl.option_present('s')) {
    calc_set_theo = 1;
  }

  if (cl.option_present('X')) {
    if (!cl.value('X', MAX_FF) || MAX_FF < 1) {
      cerr << "The max flip-flop option (-X) must be a whole +ve number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Max flip flop value " << MAX_FF << endl;
    }
  }

  if (cl.option_present('u')) {
    write_unique_smiles_as_identifier = 1;

    if (verbose) {
      cerr << "The unique smiles will be written as the identifier\n";
    }
  }

  if (cl.option_present('U')) {
    IWString unconvert_log_file_name;
    cl.value('U', unconvert_log_file_name);

    unconvert_log.open(unconvert_log_file_name.null_terminated_chars(), std::ios::out);
    if (!unconvert_log.good()) {
      cerr << "Cannot open unconvertable file '" << unconvert_log_file_name << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Unconverted structures will be written to '" << unconvert_log_file_name
           << "'\n";
    }
  }

  // endcap queries
  resizable_array_p<Substructure_Query> endcap_queries;

  endcap_queries.resize(256);  // hopefully large enough to avoid extra mallocs

  if (cl.option_present('q') && !process_queries(cl, endcap_queries, verbose)) {
    cerr << prog_name << ": cannot process queries from -q option(s)\n";
    return 6;
  }

  int nqueries = endcap_queries.number_elements();

  if (0 == nqueries) {
    cerr << prog_name << ": No endcap queries specified, use -q options\n";
    usage(12);
  }

  // assign global properties (endcaps) for queries and check queries
  for (int i = 0; i < nqueries; i++) {
    Substructure_Query* qi = endcap_queries[i];
    qi->set_find_unique_embeddings_only(1);

    assert(qi->ok());
  }

  element_ar = get_element_from_symbol_no_case_conversion("Ar");
  element_al = get_element_from_symbol_no_case_conversion("Hf");
  element_dy = get_element_from_symbol_no_case_conversion("Dy");

  initialise_property_offset();

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    const char* fname = cl[i];

    if (verbose) {
      cerr << "Processing '" << fname << "'\n";
    }

    int tmp;
    if (work_as_filter) {
      tmp = gfp_erg_filter(cl[i], output, endcap_queries);
    } else {
      tmp = gfp_erg(fname, input_type, output, endcap_queries);
    }

    if (0 == tmp) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << molecules_read << " molecules read" << endl;
    if (truncated_flip_flops) {
      cerr << truncated_flip_flops << " truncated flip flops " << MAX_FF << endl;
    }

    for (int i = 0; i < fingerprints_per_molecule.number_elements(); i++) {
      if (fingerprints_per_molecule[i]) {
        cerr << fingerprints_per_molecule[i] << " molecules generated " << i
             << " fingerprints\n";
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_erg(argc, argv);

  return rc;
}
