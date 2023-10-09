/*
  Computes ring size fingerprints
*/

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

const char* prog_name = nullptr;

class Ring_Size_Fingerprints
{
 private:
  int _verbose;
  int _molecules_read;
  Chemical_Standardisation _chemical_standardisation;
  int _reduce_to_largest_fragment;
  int _function_as_filter;

  int _must_be_macrocycle;
  int _discard_molecules_without_a_large_ring;
  int _non_macrocycle;

  int _examine_atom_types;
  int _examine_outside_ring_connections;

  IWString _smiles_tag;
  IWString _identifier_tag;
  IWString _tag;

  std::ofstream _stream_for_bits;

  // private functions

  void _default_values();
  void _usage(int rc);
  void _preprocess(Molecule& m);
  int _do_output(Molecule& m, const Sparse_Fingerprint_Creator&,
             IWString_and_File_Descriptor&);
  int _size_of_largest_ring_in_largest_ring_system(
      Molecule& m, const resizable_array<int>& largest_ring_system,
      const int* fsid) const;
  int _do_within_ring_unsaturation(Molecule& m, Sparse_Fingerprint_Creator&,
                               const resizable_array<int>& largest_ring_system,
                               const int* fsid) const;
  int _singly_fused_in_largest_system(Molecule& m, Sparse_Fingerprint_Creator& sfc,
                                  const resizable_array<int>& largest_ring_system,
                                  const int* fsid) const;
  int _outside_ring_connections(Molecule& m, Sparse_Fingerprint_Creator& sfc,
                            const resizable_array<int>& largest_ring_system,
                            const int* fsid) const;
  int _do_aromatic_in_largest_ring_system(Molecule& m, Sparse_Fingerprint_Creator& sfc,
                                      const resizable_array<int>& largest_ring_system,
                                      const int* fsid) const;
  int _within_system_adjacent_rings(Molecule& m, Sparse_Fingerprint_Creator& sfc,
                                const resizable_array<int>& largest_ring_system,
                                const int* fsid) const;
  int _two_connected_in_largest_ring(Molecule& m, Sparse_Fingerprint_Creator& sfc,
                                 const resizable_array<int>& largest_ring_system,
                                 const int* fsid) const;

  int _Ring_Size_Fingerprints_filter(iwstring_data_source& input,
                                 IWString_and_File_Descriptor& output);
  int _Ring_Size_Fingerprints_filter(const char* fname, IWString_and_File_Descriptor& output);

  int _Ring_Size_Fingerprints(Molecule& m, IWString_and_File_Descriptor& output);
  int _Ring_Size_Fingerprints(data_source_and_type<Molecule>& input,
                          IWString_and_File_Descriptor& output);
  int _Ring_Size_Fingerprints(const char* fname, FileType input_type,
                          IWString_and_File_Descriptor& output);

 public:
  Ring_Size_Fingerprints();

  int operator()(int argc, char** argv);
};

Ring_Size_Fingerprints::Ring_Size_Fingerprints()
{
  _default_values();

  return;
}

void
Ring_Size_Fingerprints::_default_values()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _function_as_filter = 0;

  _must_be_macrocycle = 0;
  _discard_molecules_without_a_large_ring = 0;
  _non_macrocycle = 0;

  _examine_atom_types = 0;
  _examine_outside_ring_connections = 0;

  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _tag = "NCRSZ<";

  return;
}

void
Ring_Size_Fingerprints::_usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Message about what the programme does\n";
  cerr << " -r <rsize>    only process macrocycles. Skip molecules without a ring at least <rsize>\n";
  cerr << " -t            include bits from ring atom types\n";
  cerr << " -e            include bits from outside ring connections\n";
  cerr << " -J <tag>      tag for fingerprints\n";
  cerr << " -f            work as a TDT filter\n";
  cerr << " -l            reduce to largest fragment\n";
  cerr << " -i <type>     input specification\n";
  cerr << " -g ...        chemical standardisation options\n";
  cerr << " -E ...        standard element specifications\n";
  cerr << " -A ...        standard aromaticity specifications\n";
  cerr << " -v            verbose output\n";
  cerr << " -x            exclude large rings, not compatible with -f\n";
  // clang-format on

  exit(rc);
}

void
Ring_Size_Fingerprints::_preprocess(Molecule& m)
{
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return;
}

int
Ring_Size_Fingerprints::_do_output(Molecule& m, const Sparse_Fingerprint_Creator& sfc,
                                   IWString_and_File_Descriptor& output)
{
  if (!_function_as_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_tag, tmp);
  output << tmp << '\n';

  if (!_function_as_filter) {
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  if (_stream_for_bits.is_open()) {
    //  _stream_for_bits << m.smiles() << ' ' << m.name();
    _stream_for_bits << m.name();
    sfc.to_svml(_stream_for_bits);
    _stream_for_bits << '\n';
  }

  return 1;
}

#define NRINGS_BIT 1
#define LARGEST_RING_SIZE 2
#define LARGEST_NON_SSSR_RING_SIZE 3
#define AROMATIC_RING_COUNT 4
#define ALIPHATIC_RING_COUNT 5
#define SMALLEST_RING 6
#define AVERAGE_SSSR_RING_SIZE 7
#define RINGS_IN_LARGEST_FUSED_SYSTEM 8
#define ISOLATED_RINGS 9
#define ISOLATED_AROMATIC_RINGS 10
#define ISOLATED_ALIPHATIC_RINGS 11
#define NUMBER_RING_SYSTEMS 12
#define FRACTION_RING_ATOMS 13
#define RING_AROM_DENSITY 14
#define NUMBER_NON_SSSR_RINGS 15
#define NUMBER_STRONGLY_FUSED_RINGS 16
#define AVERAGE_RING_SIZE_INC_NON_SSSR 17
#define NO_MACROCYCLE 18
#define AROMATIC_RING_ATOM_COUNT 19

#define MAX_ATOMS_IN_SYSTEM 20
#define MAX_RINGS_IN_SYSTEM 21

#define AVE_ATOMS_IN_SYSTEM 22
#define AVE_RINGS_IN_SYSTEM 23

#define FIVE_MEMBERED_AROMATIC_RINGS 25

#define DEGENERACY_LARGEST_RING_SYSTEM 27

// we keep track of rings in various bucket ranges. They overlap to avoid
// edge effects

#define RING_SIZE_8_12 100
#define RING_SIZE_13_17 101
#define RING_SIZE_18_22 102
#define RING_SIZE_23_27 103
#define RING_SIZE_28_32 104

#define RING_SIZE_10_14 104
#define RING_SIZE_15_19 105
#define RING_SIZE_20_24 106
#define RING_SIZE_25_29 107
#define RING_SIZE_30_34 108

#define ATOMS_IN_SYSTEM_8_12 100
#define ATOMS_IN_SYSTEM_13_17 101
#define ATOMS_IN_SYSTEM_18_22 102
#define ATOMS_IN_SYSTEM_23_27 103
#define ATOMS_IN_SYSTEM_28_32 104

#define ATOMS_IN_SYSTEM_10_14 106
#define ATOMS_IN_SYSTEM_15_19 107
#define ATOMS_IN_SYSTEM_20_24 108
#define ATOMS_IN_SYSTEM_25_29 109
#define ATOMS_IN_SYSTEM_30_34 110

#define RINGS_IN_SYSTEM_2_4 120
#define RINGS_IN_SYSTEM_5_7 121
#define RINGS_IN_SYSTEM_8_10 122
#define RINGS_IN_SYSTEM_11_13 123

#define RINGS_IN_SYSTEM_3_5 124
#define RINGS_IN_SYSTEM_6_8 125
#define RINGS_IN_SYSTEM_9_11 125
#define RINGS_IN_SYSTEM_12_14 125

#define TWO_CONNECTED_IN_LARGEST_SYSTEM_RING 130
#define TWO_CONNECTED_IN_LARGEST_SYSTEM_RING_FRACTION 131
#define TWO_CONNECTED_IN_LARGEST_SYSTEM 132
#define THREE_CONNECTED_IN_LARGEST_SYSTEM 133
#define FOUR_CONNECTED_IN_LARGEST_SYSTEM 134

#define SINGLY_FUSED_IN_LARGEST_SYSTEM 140

#define RING_FUSION_COUNT 210
#define RING_BIGGER_THAN_30 211

#define FULLY_SATURATED_ATOMS 220
#define UNSATURATED_ATOMS 221

#define LRS_FULLY_SATURATED_ATOMS 224
#define LRS_UNSATURATED_ATOMS 225

#define RING_CARBON_ATOMS 300
#define RING_HETERO_ATOMS 301
#define AROMATIC_CARBON 302
#define ALIPHATIC_CARBON 303

#define RING_NITROGEN_ATOMS 304
#define AROMATIC_NITROGEN 305
#define ALIPHATIC_NITROGEN 306

#define RING_OXYGEN_ATOMS 307
#define AROMATIC_OXYGEN 308
#define ALIPHATIC_OXYGEN 309

#define RING_SULPHUR_ATOMS 310
#define AROMATIC_SULPHUR 311
#define ALIPHATIC_SULPHUR 312

#define LRS_RING_CARBON_ATOMS 250
#define LRS_AROMATIC_CARBON 251
#define LRS_ALIPHATIC_CARBON 252
#define LRS_RING_NITROGEN_ATOMS 253
#define LRS_AROMATIC_NITROGEN 254
#define LRS_ALIPHATIC_NITROGEN 255
#define LRS_RING_OXYGEN_ATOMS 256
#define LRS_AROMATIC_OXYGEN 257
#define LRS_ALIPHATIC_OXYGEN 258
#define LRS_RING_SULPHUR_ATOMS 259
#define LRS_AROMATIC_SULPHUR 260
#define LRS_ALIPHATIC_SULPHUR 261
#define LRS_RING_HETERO_ATOMS 262

#define FSDR_RING_HETERO_ATOMS 270
#define FSDR_RING_CARBON_ATOMS 271
#define FSDR_AROMATIC_CARBON 272
#define FSDR_ALIPHATIC_CARBON 273
#define FSDR_RING_NITROGEN_ATOMS 274
#define FSDR_AROMATIC_NITROGEN 275
#define FSDR_ALIPHATIC_NITROGEN 276
#define FSDR_RING_OXYGEN_ATOMS 277
#define FSDR_AROMATIC_OXYGEN 278
#define FSDR_ALIPHATIC_OXYGEN 279
#define FSDR_RING_SULPHUR_ATOMS 280
#define FSDR_AROMATIC_SULPHUR 281
#define FSDR_ALIPHATIC_SULPHUR 282

#define WITHIN_RING_UNSATURATION 320
#define OUTSIDE_RING_UNSATURATION 321

#define RING_SYSTEM_NCON 330
#define RING_SYSTEM_SINGLY_CONNECTED 331
#define RING_SYSTEM_EXOCYCLIC_DOUBLE_BOND 332

#define GAP_FROM_LARGEST 340
#define FUSED_NBRS_LARGEST_RINGS 341
#define LARGEST_RING_IS_ISOLATED 342
#define LARGEST_RING_CH2 343
#define WITHIN_SYSTEM_ADJACENT_RINGS 343

#define AROMATIC_RINGS_IN_LARGEST_RING_SYSTEM 350

#define RING_BOND_COUNT_2 360
#define RING_BOND_COUNT_3 361

static void
do_hit_bit(Sparse_Fingerprint_Creator& sfc, const unsigned int b, const int c = 1)
{
  if (c > 0) {
    sfc.hit_bit(b, c);
  }

  return;
}

int
Ring_Size_Fingerprints::_do_aromatic_in_largest_ring_system(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  const int nr = m.nrings();

  int aromatic_rings_in_largest_ring_system = 0;

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (!ri->is_aromatic()) {
      continue;
    }

    if (!largest_ring_system.contains(fsid[ri->item(0)])) {
      continue;
    }

    aromatic_rings_in_largest_ring_system++;
  }

  do_hit_bit(sfc, AROMATIC_RINGS_IN_LARGEST_RING_SYSTEM,
             aromatic_rings_in_largest_ring_system);

  return 1;
}

int
Ring_Size_Fingerprints::_do_within_ring_unsaturation(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  const int matoms = m.natoms();

  int within_ring_unsaturation = 0;
  int outside_ring_unsaturation = 0;

  for (int i = 0; i < matoms; ++i) {
    if (!largest_ring_system.contains(fsid[i])) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const int acon = a->ncon();

    for (int j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      if (b->is_single_bond()) {
        continue;
      }

      if (b->is_aromatic()) {
        continue;
      }

      if (0 == b->nrings()) {
        within_ring_unsaturation++;
      } else {
        outside_ring_unsaturation++;
      }
    }
  }

  do_hit_bit(sfc, WITHIN_RING_UNSATURATION, within_ring_unsaturation);
  do_hit_bit(sfc, OUTSIDE_RING_UNSATURATION, outside_ring_unsaturation);

  return 1;
}

int
Ring_Size_Fingerprints::_outside_ring_connections(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
#ifdef ECHO_LARGEST_RING_SYSTEM_SIZE
  cerr << "Largest ring system size " << largest_ring_system.size();
  for (int i = 0; i < largest_ring_system.number_elements(); ++i) {
    cerr << ' ' << largest_ring_system[i];
  }
  cerr << '\n';
#endif

  const int matoms = m.natoms();

  int ncon = 0;
  int singly_connected = 0;
  int exocyclic_double_bond = 0;

  for (int i = 0; i < matoms; ++i) {
    if (!largest_ring_system.contains(fsid[i])) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const auto acon = a->ncon();

    if (2 == acon) {
      continue;
    }

    for (int j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      if (b->nrings()) {
        continue;
      }

      const auto k = b->other(i);

      ncon++;
      if (1 == m.ncon(k)) {
        singly_connected++;
        if (b->is_double_bond()) {
          exocyclic_double_bond++;
        }
      }
    }
  }

  do_hit_bit(sfc, RING_SYSTEM_NCON, ncon);
  do_hit_bit(sfc, RING_SYSTEM_SINGLY_CONNECTED, singly_connected);
  do_hit_bit(sfc, RING_SYSTEM_EXOCYCLIC_DOUBLE_BOND, exocyclic_double_bond);

  return 1;
}

int
Ring_Size_Fingerprints::_size_of_largest_ring_in_largest_ring_system(
    Molecule& m, const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  for (int i = m.nrings() - 1; i >= 0; --i) {
    const Ring* ri = m.ringi(i);

    if (largest_ring_system.contains(
            fsid[ri->item(0)])) {  // all atoms in the ring have the same fsid, so just
                                   // check the first one
      return ri->number_elements();
    }
  }

  cerr << "_size_of_largest_ring_in_largest_ring_system:should not happen\n";

  return 0;  // should not come here
}

int
Ring_Size_Fingerprints::_two_connected_in_largest_ring(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  int largest_ring_size_in_system =
      _size_of_largest_ring_in_largest_ring_system(m, largest_ring_system, fsid);

  int rc = 0;

  for (int i = m.nrings() - 1; i >= 0; --i) {
    const Ring* ri = m.ringi(i);

    const auto ring_size = ri->number_elements();

    if (ring_size > largest_ring_size_in_system) {
      continue;
    } else if (ring_size < largest_ring_size_in_system) {
      break;
    }

    if (!largest_ring_system.contains(fsid[ri->item(0)])) {
      continue;
    }

    if (2 == m.ncon(i)) {
      rc++;
    }
  }

  do_hit_bit(sfc, TWO_CONNECTED_IN_LARGEST_SYSTEM_RING, rc);

  float f = static_cast<float>(rc * 10) / static_cast<float>(largest_ring_size_in_system);

  do_hit_bit(sfc, TWO_CONNECTED_IN_LARGEST_SYSTEM_RING_FRACTION,
             static_cast<int>(f + 0.4999f));

  const int matoms = m.natoms();

  int rc2 = 0;
  int rc3 = 0;
  int rc4 = 0;

  for (int i = 0; i < matoms; ++i) {
    if (!largest_ring_system.contains(fsid[i])) {
      continue;
    }

    const auto icon = m.ncon(i);

    if (2 == icon) {
      rc2++;
    } else if (3 == icon) {
      rc3++;
    } else {
      rc4++;
    }
  }

  do_hit_bit(sfc, TWO_CONNECTED_IN_LARGEST_SYSTEM, rc2);
  do_hit_bit(sfc, THREE_CONNECTED_IN_LARGEST_SYSTEM, rc3);
  do_hit_bit(sfc, FOUR_CONNECTED_IN_LARGEST_SYSTEM, rc4);

  return 1;
}

int
Ring_Size_Fingerprints::_singly_fused_in_largest_system(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  int largest_ring_size_in_system =
      _size_of_largest_ring_in_largest_ring_system(m, largest_ring_system, fsid);

  int rc = 0;

  for (int i = m.nrings() - 1; i >= 0; --i) {
    const Ring* ri = m.ringi(i);

    const auto ring_size = ri->number_elements();

    if (ring_size > largest_ring_size_in_system) {
      continue;
    } else if (ring_size < largest_ring_size_in_system) {
      break;
    }

    if (!largest_ring_system.contains(fsid[ri->item(0)])) {
      continue;
    }

    if (1 == ri->fused_ring_neighbours()) {
      rc++;
    }
  }

  do_hit_bit(sfc, SINGLY_FUSED_IN_LARGEST_SYSTEM, rc);

  return 1;
}

/*
  Look for situations where we have a single bond joining two rings
  embedded within the largest ring system. For example a large macrocycle
  with adjacent pyridine rings within the large macrocycle.
  Such situtations are deliberately double counted.
*/

int
Ring_Size_Fingerprints::_within_system_adjacent_rings(
    Molecule& m, Sparse_Fingerprint_Creator& sfc,
    const resizable_array<int>& largest_ring_system, const int* fsid) const
{
  int largest_ring_size_in_system =
      _size_of_largest_ring_in_largest_ring_system(m, largest_ring_system, fsid);

  int rc = 0;

  for (int i = m.nrings() - 1; i >= 0; --i) {
    const Ring* ri = m.ringi(i);

    const int ring_size = ri->number_elements();

    if (ring_size > largest_ring_size_in_system) {
      continue;
    } else if (ring_size < largest_ring_size_in_system) {
      break;
    }

    if (!largest_ring_system.contains(fsid[ri->item(0)])) {
      continue;
    }

    for (int j = 0; j < ring_size; ++j) {
      const auto k = ri->item(j);

      const Atom* a = m.atomi(k);

      const auto acon = a->ncon();

      if (2 == acon) {
        continue;
      }

      if (1 == m.nrings(k)) {
        continue;
      }

      for (int l = 0; l < acon; ++l) {
        const Bond* b = a->item(l);

        if (1 != b->nrings()) {
          continue;
        }

        const auto o = b->other(k);

        assert(largest_ring_system.contains(fsid[o]));

        if (!ri->contains(o)) {
          continue;
        }

        if (m.nrings(o) > 1) {
          rc++;
        }
      }
    }
  }

  do_hit_bit(sfc, WITHIN_SYSTEM_ADJACENT_RINGS, rc);

  return rc;
}

static void
set_ring_size_bits(Sparse_Fingerprint_Creator& sfc, const int rsize)
{
  if (rsize > 30) {
    do_hit_bit(sfc, RING_BIGGER_THAN_30, 1);
    return;
  }

  if (rsize >= 8 && rsize <= 12) {
    do_hit_bit(sfc, RING_SIZE_8_12);
  } else if (rsize >= 13 && rsize <= 17) {
    do_hit_bit(sfc, RING_SIZE_13_17);
  } else if (rsize >= 18 && rsize <= 22) {
    do_hit_bit(sfc, RING_SIZE_18_22);
  } else if (rsize >= 23 && rsize <= 27) {
    do_hit_bit(sfc, RING_SIZE_23_27);
  } else if (rsize >= 28 && rsize <= 32) {
    do_hit_bit(sfc, RING_SIZE_28_32);
  }

  if (rsize >= 10 && rsize <= 14) {
    do_hit_bit(sfc, RING_SIZE_10_14);
  } else if (rsize >= 15 && rsize <= 19) {
    do_hit_bit(sfc, RING_SIZE_15_19);
  } else if (rsize >= 20 && rsize <= 24) {
    do_hit_bit(sfc, RING_SIZE_20_24);
  } else if (rsize >= 25 && rsize <= 29) {
    do_hit_bit(sfc, RING_SIZE_25_29);
  } else if (rsize >= 30 && rsize <= 34) {
    do_hit_bit(sfc, RING_SIZE_30_34);
  }

  return;
}

static void
set_ring_system_size_bits(Sparse_Fingerprint_Creator& sfc, const int rings_in_system,
                          const int atoms_in_system)
{
  if (atoms_in_system >= 8 && atoms_in_system <= 12) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_8_12);
  } else if (atoms_in_system >= 13 && atoms_in_system <= 17) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_13_17);
  } else if (atoms_in_system >= 18 && atoms_in_system <= 22) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_18_22);
  } else if (atoms_in_system >= 23 && atoms_in_system <= 27) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_23_27);
  } else if (atoms_in_system >= 28 && atoms_in_system <= 32) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_28_32);
  }

  if (atoms_in_system >= 10 && atoms_in_system <= 14) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_10_14);
  } else if (atoms_in_system >= 15 && atoms_in_system <= 19) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_15_19);
  } else if (atoms_in_system >= 20 && atoms_in_system <= 24) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_20_24);
  } else if (atoms_in_system >= 25 && atoms_in_system <= 29) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_25_29);
  } else if (atoms_in_system >= 30 && atoms_in_system <= 34) {
    do_hit_bit(sfc, ATOMS_IN_SYSTEM_30_34);
  }

  if (rings_in_system >= 2 && rings_in_system <= 4) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_2_4);
  } else if (rings_in_system >= 5 && rings_in_system <= 7) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_5_7);
  } else if (rings_in_system >= 8 && rings_in_system <= 10) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_8_10);
  } else if (rings_in_system >= 11 && rings_in_system <= 13) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_11_13);
  }

  if (rings_in_system >= 3 && rings_in_system <= 5) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_3_5);
  } else if (rings_in_system >= 6 && rings_in_system <= 8) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_6_8);
  } else if (rings_in_system >= 9 && rings_in_system <= 11) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_9_11);
  } else if (rings_in_system >= 12 && rings_in_system <= 14) {
    do_hit_bit(sfc, RINGS_IN_SYSTEM_12_14);
  }

  return;
}

static int
count_ch2(Molecule& m, const Ring& r)
{
  const int ring_size = r.number_elements();

  int rc = 0;

  for (int i = 0; i < ring_size; ++i) {
    const auto j = r[i];

    if (6 != m.atomic_number(j)) {
      continue;
    }

    if (2 != m.ncon(i)) {
      continue;
    }

    if (2 == m.hcount(i)) {
      rc++;
    }
  }

  return rc;
}

static int
only_in_isolated_small_ring(Molecule& m, const atom_number_t zatom,
                            const int max_ring_size)
{
  if (m.nrings(zatom) > 1) {
    return 0;
  }

  const auto nrings = m.nrings();

  for (int i = 0; nrings; ++i) {
    const Ring* ri = m.ringi(i);

    if (!ri->contains(zatom)) {
      continue;
    }

    if (ri->is_fused()) {
      return 0;
    }

    if (ri->number_elements() <= max_ring_size) {
      return 1;
    }

    return 0;
  }

  return 0;  // should never get here
}

class Element_Counter
{
 private:
  int _aromatic;
  int _aliphatic;

 public:
  Element_Counter();

  void
  extra(const int is_aromatic)
  {
    if (is_aromatic) {
      _aromatic++;
    } else {
      _aliphatic++;
    }
  }

  int
  natoms() const
  {
    return _aromatic + _aliphatic;
  }

  int
  aromatic() const
  {
    return _aromatic;
  }

  int
  aliphatic() const
  {
    return _aliphatic;
  }

  void
  reset()
  {
    _aromatic = _aliphatic = 0;
  }
};

Element_Counter::Element_Counter()
{
  _aromatic = 0;
  _aliphatic = 0;

  return;
}

int
Ring_Size_Fingerprints::_Ring_Size_Fingerprints(Molecule& m,
                                                IWString_and_File_Descriptor& output)
{
  const int nr = m.nrings();

  Sparse_Fingerprint_Creator sfc;

  do_hit_bit(sfc, NRINGS_BIT, nr + 1);  // add 1 so something with no rings gets something

  if (0 == nr) {
    return _do_output(m, sfc, output);
  }

  if (_must_be_macrocycle > 0 &&
      m.ringi(nr - 1)->number_elements() < _must_be_macrocycle) {
    _non_macrocycle++;
    if (_discard_molecules_without_a_large_ring) {
      return 1;
    }

    do_hit_bit(sfc, NO_MACROCYCLE);
    return _do_output(m, sfc, output);
  }

  m.compute_aromaticity_if_needed();

  Accumulator_Int<int> acc_rsize;

  int aromatic_ring_count = 0;
  int fused_rings = 0;
  int strongly_fused = 0;
  int isolated_rings = 0;
  int isolated_aliphatic_ring_count = 0;
  int largest_isolated_ring_size = 0;
  int isolated_aromatic_ring_count = 0;
  int five_membered_aromatic_rings = 0;

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    const int rsize = ri->number_elements();

    acc_rsize.extra(rsize);

    set_ring_size_bits(sfc, rsize);

    if (ri->is_aromatic()) {
      aromatic_ring_count++;
      if (5 == rsize) {
        five_membered_aromatic_rings++;
      }
    }

    if (ri->strongly_fused_ring_neighbours()) {
      strongly_fused++;
    }

    if (ri->is_fused()) {
      fused_rings++;
    } else {
      isolated_rings++;
      if (rsize > largest_isolated_ring_size) {
        largest_isolated_ring_size = rsize;
      }

      if (ri->is_aromatic()) {
        isolated_aromatic_ring_count++;
      } else {
        isolated_aliphatic_ring_count++;
      }
    }
  }

  do_hit_bit(sfc, SMALLEST_RING, acc_rsize.minval());
  do_hit_bit(sfc, LARGEST_RING_SIZE, acc_rsize.maxval());
  do_hit_bit(sfc, AVERAGE_SSSR_RING_SIZE, static_cast<int>(acc_rsize.average() + 0.4999));

  do_hit_bit(sfc, ISOLATED_AROMATIC_RINGS, isolated_aromatic_ring_count);
  do_hit_bit(sfc, ISOLATED_ALIPHATIC_RINGS, isolated_aliphatic_ring_count);

  do_hit_bit(sfc, RING_FUSION_COUNT, fused_rings);
  do_hit_bit(sfc, FIVE_MEMBERED_AROMATIC_RINGS, five_membered_aromatic_rings);

  do_hit_bit(sfc, AROMATIC_RING_COUNT, aromatic_ring_count);
  do_hit_bit(sfc, ALIPHATIC_RING_COUNT, nr - aromatic_ring_count);
  do_hit_bit(sfc, NUMBER_STRONGLY_FUSED_RINGS, strongly_fused);

  const int non_sssr_rings = m.non_sssr_rings();

  if (non_sssr_rings > 0) {
    do_hit_bit(sfc, NUMBER_NON_SSSR_RINGS, non_sssr_rings);

    int largest_non_sssr_ring_size = 0;

    for (int i = 0; i < non_sssr_rings; ++i) {
      const Ring* ri = m.non_sssr_ring(i);

      const int rsize = ri->number_elements();

      acc_rsize.extra(rsize);

      if (rsize > largest_non_sssr_ring_size) {
        largest_non_sssr_ring_size = rsize;
      }
    }

    do_hit_bit(sfc, AVERAGE_RING_SIZE_INC_NON_SSSR,
               static_cast<int>(acc_rsize.average() + 0.4999));
    do_hit_bit(sfc, LARGEST_NON_SSSR_RING_SIZE, largest_non_sssr_ring_size);
  }

  const int matoms = m.natoms();

  int aromatic_atom_count = 0;
  int non_aromatic_ring_atoms = 0;

  int saturated = 0;
  int unsaturated = 0;

  int ring_bond_count_2 = 0;
  int ring_bond_count_3 = 0;

  for (int i = 0; i < matoms; ++i) {
    if (!m.is_ring_atom(i)) {
      continue;
    }

    if (m.is_aromatic(i)) {
      aromatic_atom_count++;
    } else {
      non_aromatic_ring_atoms++;
    }

    if (m.ncon(i) == m.nbonds(i)) {
      saturated++;
    } else {
      unsaturated++;
    }

    const int rbc = m.ring_bond_count(i);
    if (2 == rbc) {
      ring_bond_count_2++;
    } else {
      ring_bond_count_3++;
    }
  }

  do_hit_bit(sfc, AROMATIC_RING_ATOM_COUNT, aromatic_atom_count);
  do_hit_bit(sfc, RING_AROM_DENSITY,
             static_cast<int>(
                 static_cast<float>(aromatic_atom_count * 10) /
                 static_cast<float>(aromatic_atom_count + non_aromatic_ring_atoms)));

  do_hit_bit(sfc, FRACTION_RING_ATOMS,
             static_cast<int>((aromatic_atom_count + non_aromatic_ring_atoms) * 10 /
                                  static_cast<float>(matoms) +
                              0.4999));

  do_hit_bit(sfc, FULLY_SATURATED_ATOMS, saturated);
  do_hit_bit(sfc, UNSATURATED_ATOMS, unsaturated);

  do_hit_bit(sfc, RING_BOND_COUNT_2, ring_bond_count_2);
  do_hit_bit(sfc, RING_BOND_COUNT_3, ring_bond_count_3);

  int* fsid = new int[m.natoms()];
  std::unique_ptr<int[]> free_fsid(fsid);

  m.label_atoms_by_ring_system_including_spiro_fused(fsid);

  resizable_array<int> already_done;

  Accumulator_Int<int> acc_atoms_in_ring_system, acc_rings_in_ring_system;

  int atoms_in_largest_ring_system = 0;
  resizable_array<int> largest_ring_system;

  int number_ring_systems = 0;

  for (int i = 0; i < matoms; ++i) {
    if (0 == fsid[i]) {
      continue;
    }

    if (already_done.contains(fsid[i])) {
      continue;
    }

    const int atoms_in_system = std::count(fsid, fsid + matoms, fsid[i]);

    acc_atoms_in_ring_system.extra(atoms_in_system);

    if (atoms_in_system > atoms_in_largest_ring_system) {
      atoms_in_largest_ring_system = atoms_in_system;
      largest_ring_system.resize_keep_storage(0);
      largest_ring_system.add(fsid[i]);
    } else if (atoms_in_system == atoms_in_largest_ring_system) {
      largest_ring_system.add(fsid[i]);
    }

    int rings_in_system = 0;

    for (int j = 0; j < nr; ++j) {
      const Ring* rj = m.ringi(j);

      if (!rj->any_members_set_in_array(fsid, fsid[i])) {
        continue;
      }

      rings_in_system++;
    }

    set_ring_system_size_bits(sfc, rings_in_system, atoms_in_system);

    acc_rings_in_ring_system.extra(rings_in_system);

    already_done.add(fsid[i]);

    number_ring_systems++;
  }

  _within_system_adjacent_rings(m, sfc, largest_ring_system, fsid);

  do_hit_bit(sfc, NUMBER_RING_SYSTEMS, number_ring_systems);
  do_hit_bit(sfc, MAX_ATOMS_IN_SYSTEM, acc_atoms_in_ring_system.maxval());
  do_hit_bit(sfc, MAX_RINGS_IN_SYSTEM, acc_rings_in_ring_system.maxval());

  do_hit_bit(sfc, AVE_ATOMS_IN_SYSTEM,
             static_cast<int>(acc_atoms_in_ring_system.average() + 0.4999));
  do_hit_bit(sfc, AVE_RINGS_IN_SYSTEM,
             static_cast<int>(acc_rings_in_ring_system.average() + 0.4999));
  do_hit_bit(sfc, DEGENERACY_LARGEST_RING_SYSTEM, largest_ring_system.number_elements());

  if (_examine_atom_types) {
    Element_Counter carbon, nitrogen, oxygen, sulphur;

    int ring_hetero_atoms = 0;

    for (int i = 0; i < matoms; ++i) {
      if (!m.is_ring_atom(i)) {
        continue;
      }

      const auto z = m.atomic_number(i);
      const auto arom = m.is_aromatic(i);

      if (6 == z) {
        carbon.extra(arom);
      } else {
        ring_hetero_atoms++;
        if (7 == z) {
          nitrogen.extra(arom);
        } else if (8 == z) {
          oxygen.extra(arom);
        } else if (16 == z) {
          sulphur.extra(arom);
        }
      }
    }

    do_hit_bit(sfc, RING_HETERO_ATOMS, ring_hetero_atoms);
    do_hit_bit(sfc, RING_CARBON_ATOMS, carbon.natoms());
    do_hit_bit(sfc, AROMATIC_CARBON, carbon.aromatic());
    do_hit_bit(sfc, ALIPHATIC_CARBON, carbon.aliphatic());
    do_hit_bit(sfc, RING_NITROGEN_ATOMS, nitrogen.natoms());
    do_hit_bit(sfc, AROMATIC_NITROGEN, nitrogen.aromatic());
    do_hit_bit(sfc, ALIPHATIC_NITROGEN, nitrogen.aliphatic());
    do_hit_bit(sfc, RING_OXYGEN_ATOMS, oxygen.natoms());
    do_hit_bit(sfc, AROMATIC_OXYGEN, oxygen.aromatic());
    do_hit_bit(sfc, ALIPHATIC_OXYGEN, oxygen.aliphatic());
    do_hit_bit(sfc, RING_SULPHUR_ATOMS, sulphur.natoms());
    do_hit_bit(sfc, AROMATIC_SULPHUR, sulphur.aromatic());
    do_hit_bit(sfc, ALIPHATIC_SULPHUR, sulphur.aliphatic());

    _do_within_ring_unsaturation(m, sfc, largest_ring_system, fsid);
    _do_aromatic_in_largest_ring_system(m, sfc, largest_ring_system, fsid);

    //  now do the same thing for just the largest ring system

    carbon.reset();
    nitrogen.reset();
    oxygen.reset();
    sulphur.reset();

    saturated = 0;
    unsaturated = 0;
    ring_hetero_atoms = 0;

    for (int i = 0; i < matoms; ++i) {
      if (!largest_ring_system.contains(fsid[i])) {
        continue;
      }

      if (m.ncon(i) == m.nbonds(i)) {
        saturated++;
      } else {
        unsaturated++;
      }

      const auto z = m.atomic_number(i);
      const auto arom = m.is_aromatic(i);

      if (6 == z) {
        carbon.extra(arom);
      } else {
        ring_hetero_atoms++;
        if (7 == z) {
          nitrogen.extra(arom);
        } else if (8 == z) {
          sulphur.extra(arom);
        } else if (16 == z) {
          sulphur.extra(arom);
        }
      }
    }

    do_hit_bit(sfc, LRS_RING_HETERO_ATOMS, ring_hetero_atoms);
    do_hit_bit(
        sfc, LRS_RING_CARBON_ATOMS,
        carbon.natoms() * 2);  // higher weight for being in the largest ring system
    do_hit_bit(sfc, LRS_AROMATIC_CARBON, carbon.aromatic() * 2);
    do_hit_bit(sfc, LRS_ALIPHATIC_CARBON, carbon.aliphatic() * 2);
    do_hit_bit(sfc, LRS_RING_NITROGEN_ATOMS, nitrogen.natoms() * 2);
    do_hit_bit(sfc, LRS_AROMATIC_NITROGEN, nitrogen.aromatic() * 2);
    do_hit_bit(sfc, LRS_ALIPHATIC_NITROGEN, nitrogen.aliphatic() * 2);
    do_hit_bit(sfc, LRS_RING_OXYGEN_ATOMS, oxygen.natoms() * 2);
    do_hit_bit(sfc, LRS_AROMATIC_OXYGEN, oxygen.aromatic() * 2);
    do_hit_bit(sfc, LRS_ALIPHATIC_OXYGEN, oxygen.aliphatic() * 2);
    do_hit_bit(sfc, LRS_RING_SULPHUR_ATOMS, sulphur.natoms() * 2);
    do_hit_bit(sfc, LRS_AROMATIC_SULPHUR, sulphur.aromatic() * 2);
    do_hit_bit(sfc, LRS_ALIPHATIC_SULPHUR, sulphur.aliphatic() * 2);

    do_hit_bit(sfc, LRS_FULLY_SATURATED_ATOMS, saturated);
    do_hit_bit(sfc, LRS_UNSATURATED_ATOMS, unsaturated);

    //  now omit isolated rings -

    carbon.reset();
    nitrogen.reset();
    oxygen.reset();
    sulphur.reset();
    ring_hetero_atoms = 0;

    for (int i = 0; i < matoms; ++i) {
      if (0 == fsid[i]) {
        continue;
      }

      if (largest_ring_system.contains(fsid[i])) {
        ;
      } else if (only_in_isolated_small_ring(m, i, 6)) {
        continue;
      }

      const auto z = m.atomic_number(i);
      const auto arom = m.is_aromatic(i);

      if (6 == z) {
        carbon.extra(arom);
      } else {
        ring_hetero_atoms++;
        if (7 == z) {
          nitrogen.extra(arom);
        } else if (8 == z) {
          oxygen.extra(arom);
        } else if (16 == z) {
          sulphur.extra(arom);
        }
      }
    }

    do_hit_bit(sfc, FSDR_RING_HETERO_ATOMS, ring_hetero_atoms);
    do_hit_bit(
        sfc, FSDR_RING_CARBON_ATOMS,
        carbon.natoms() * 2);  // higher weight for being in the largest ring system
    do_hit_bit(sfc, FSDR_AROMATIC_CARBON, carbon.aromatic() * 2);
    do_hit_bit(sfc, FSDR_ALIPHATIC_CARBON, carbon.aliphatic() * 2);
    do_hit_bit(sfc, FSDR_RING_NITROGEN_ATOMS, nitrogen.natoms() * 2);
    do_hit_bit(sfc, FSDR_AROMATIC_NITROGEN, nitrogen.aromatic() * 2);
    do_hit_bit(sfc, FSDR_ALIPHATIC_NITROGEN, nitrogen.aliphatic() * 2);
    do_hit_bit(sfc, FSDR_RING_OXYGEN_ATOMS, oxygen.natoms() * 2);
    do_hit_bit(sfc, FSDR_AROMATIC_OXYGEN, oxygen.aromatic() * 2);
    do_hit_bit(sfc, FSDR_ALIPHATIC_OXYGEN, oxygen.aliphatic() * 2);
    do_hit_bit(sfc, FSDR_RING_SULPHUR_ATOMS, sulphur.natoms() * 2);
    do_hit_bit(sfc, FSDR_AROMATIC_SULPHUR, sulphur.aromatic() * 2);
    do_hit_bit(sfc, FSDR_ALIPHATIC_SULPHUR, sulphur.aliphatic() * 2);
  }

  for (int i = 0; i < matoms; ++i) {
    if (!m.is_ring_atom(i)) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const auto acon = a->ncon();

    if (2 == acon) {
      continue;
    }
  }

  const int largest_ring_size = m.ringi(nr - 1)->number_elements();

  Accumulator_Int<int> nbrs_largest_ring;
  int gap_to_next_largest = 0;
  int ch2_largest_ring = 0;

  int largest_ring_is_isolated = 1;

  for (int i = nr - 1; i >= 0; --i) {
    const Ring* ri = m.ringi(i);

    if (ri->number_elements() != largest_ring_size) {
      gap_to_next_largest = largest_ring_size - ri->number_elements();
      break;
    }

    nbrs_largest_ring.extra(ri->fused_ring_neighbours());

    if (ri->is_fused()) {
      largest_ring_is_isolated = 0;
    }

    ch2_largest_ring += count_ch2(m, *ri);
  }

  do_hit_bit(sfc, GAP_FROM_LARGEST, gap_to_next_largest);
  do_hit_bit(sfc, FUSED_NBRS_LARGEST_RINGS,
             static_cast<int>(nbrs_largest_ring.average() + 1.4999f));
  do_hit_bit(sfc, LARGEST_RING_IS_ISOLATED, largest_ring_is_isolated);
  do_hit_bit(sfc, LARGEST_RING_CH2, ch2_largest_ring);

  if (_examine_outside_ring_connections) {
    _outside_ring_connections(m, sfc, largest_ring_system, fsid);
  }

  _two_connected_in_largest_ring(m, sfc, largest_ring_system, fsid);
  _singly_fused_in_largest_system(m, sfc, largest_ring_system, fsid);

  return _do_output(m, sfc, output);
}

int
Ring_Size_Fingerprints::_Ring_Size_Fingerprints(data_source_and_type<Molecule>& input,
                                                IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    _molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    if (!_Ring_Size_Fingerprints(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
Ring_Size_Fingerprints::_Ring_Size_Fingerprints_filter(
    iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if (!buffer.starts_with(_smiles_tag)) {
      continue;
    }

    buffer.remove_leading_chars(_smiles_tag.length());
    buffer.chop();

    Molecule m;

    if (!m.build_from_smiles(buffer)) {
      cerr << "Ring_Size_Fingerprints::filter:bad smiles '" << buffer << "'\n";
      return 0;
    }

    if (!_Ring_Size_Fingerprints(m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Ring_Size_Fingerprints::_Ring_Size_Fingerprints_filter(
    const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Ring_Size_Fingerprints:filter:cannot open '" << fname << "'\n";
    return 0;
  }

  return _Ring_Size_Fingerprints_filter(input, output);
}

int
Ring_Size_Fingerprints::_Ring_Size_Fingerprints(const char* fname, FileType input_type,
                                                IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 1) {
    input.set_verbose(1);
  }

  return _Ring_Size_Fingerprints(input, output);
}

int
Ring_Size_Fingerprints::operator()(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lfr:J:tD:ex");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    _usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, _verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      _usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, _verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;

    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('f')) {
    _function_as_filter = 1;

    if (_verbose) {
      cerr << "Will work as a TDT filter\n";
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', _must_be_macrocycle) || _must_be_macrocycle < 6) {
      cerr << "The min ring size for a macrocycle (-r) must be a valid ring size\n";
      _usage(1);
    }

    if (_verbose) {
      cerr << "Will only process rings having a ring of size " << _must_be_macrocycle
           << " or larger\n";
    }
  }

  if (cl.option_present('x')) {
    if (cl.option_present('f')) {
      cerr << "Sorry, the discard molecules without a suitable ring size (-x) and -f "
              "options are incompatible\n";
      _usage(1);
    }

    _discard_molecules_without_a_large_ring = 1;

    if (_verbose) {
      cerr << "Will discard molecules that do not have a ring large enough\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);

    if (!_tag.ends_with('<')) {
      _tag << '<';
    }

    if (_verbose) {
      cerr << "Fingerprints written with tag '" << _tag << "'\n";
    }
  }

  if (cl.option_present('t')) {
    _examine_atom_types = 1;

    if (_verbose) {
      cerr << "Will include bits based on ring atom types\n";
    }
  }

  if (cl.option_present('e')) {
    _examine_outside_ring_connections = 1;
    if (_verbose) {
      cerr << "Will include bits based on connections outside the largest ring system\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (_function_as_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      _usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  if (cl.option_present('D')) {
    const char* d = cl.option_value('D');

    _stream_for_bits.open(d);

    if (!_stream_for_bits.good()) {
      cerr << "_Ring_Size_Fingerprints:cannot open -D file '" << d << "'\n";
      return 1;
    }

    if (_verbose) {
      cerr << "Bit contributions written to '" << d << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    int rc;
    if (_function_as_filter) {
      rc = _Ring_Size_Fingerprints_filter(cl[i], output);
    } else {
      rc = _Ring_Size_Fingerprints(cl[i], input_type, output);
    }

    if (0 == rc) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (_verbose) {
    cerr << "Read " << _molecules_read << " molecules, " << _non_macrocycle
         << " were not macrocycles\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  Ring_Size_Fingerprints Ring_Size_Fingerprints;

  return Ring_Size_Fingerprints(argc, argv);
}
