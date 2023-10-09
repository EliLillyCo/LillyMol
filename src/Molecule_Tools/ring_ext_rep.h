#ifndef MOLECULE_TOOLS_RING_EXT_REP_H
#define MOLECULE_TOOLS_RING_EXT_REP_H

#include <cstdint>
#include <fstream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"

namespace ring_replacement {

#define RING_4A 0
#define RING_5a 1
#define RING_5A 2
#define RING_6a 3
#define RING_6A 4
#define RING_7a 5
#define RING_7A 6
#define RING_8A 7

#define FUSED_4A4A 8

#define FUSED_5a5a 9
#define FUSED_5a5A 10
#define FUSED_5A5A 11

#define FUSED_5a6a 12
#define FUSED_5a6A 13
#define FUSED_5A6a 14
#define FUSED_5A6A 15

#define FUSED_6a6a 16
#define FUSED_6a6A 17
#define FUSED_6A6A 18

#define FUSED_4A5A 19
#define FUSED_4A6a 20
#define FUSED_4A6A 21
#define FUSED_4a6a 22

#define FUSED_6a7A 23
#define FUSED_5A7A 24
#define FUSED_5A7a 25
#define FUSED_6a7a 26
#define FUSED_5a7a 27
#define FUSED_5a7A 28
#define FUSED_4A7A 29
#define FUSED_6A7A 30
#define FUSED_6A7a 31

#define FUSED_5a8A 32
#define FUSED_6a8A 33
#define FUSED_4A8A 34
#define FUSED_5A8A 35
#define FUSED_6A8A 36
#define FUSED_7A8A 37
#define FUSED_7a8A 38
#define FUSED_8A8A 39

#define FUSED_7A7A 40

#define FUSED_4A5a 42
#define FUSED_4A7a 45
#define FUSED_7a7a 46
#define FUSED_7a7A 47

#define FUSED_4a4a 48
#define FUSED_4a6A 49
#define FUSED_4a5a 50
#define FUSED_4a7a 51
#define FUSED_4a4A 52
#define FUSED_4a7A 53

// 3 membere fused systems.
#define FUSED_4A4A4A 54
#define FUSED_4A4A5A 55
#define FUSED_4A4A6A 56
#define FUSED_4A4A7A 57

#define FUSED_4A4A4a 58
#define FUSED_4A4A5a 59
#define FUSED_4A4A6a 60
#define FUSED_4A4A7a 61

#define FUSED_4A4a4A 62
#define FUSED_4A4a5A 63
#define FUSED_4A4a6A 64
#define FUSED_4A4a7A 65

#define FUSED_4A4a4a 66
#define FUSED_4A4a5a 67
#define FUSED_4A4a6a 68
#define FUSED_4A4a7a 69

#define FUSED_4A5A5A 70
#define FUSED_4A5A6A 71
#define FUSED_4A5A7A 72

#define FUSED_4A5A5a 73
#define FUSED_4A5A6a 74
#define FUSED_4A5A7a 75

#define FUSED_4A5a5A 76
#define FUSED_4A5a6A 77
#define FUSED_4A5a7A 78

#define FUSED_4A5a5a 79
#define FUSED_4A5a6a 80
#define FUSED_4A5a7a 81

#define FUSED_4A6A6A 82
#define FUSED_4A6A7A 83

#define FUSED_4A6A6a 84
#define FUSED_4A6A7a 85

#define FUSED_4A7A7A 86
#define FUSED_4A7A7a 87

#define FUSED_4A7a7A 88
#define FUSED_4A7a7a 89

#define FUSED_4a4A4a 90
#define FUSED_4a4A5a 91
#define FUSED_4a4A6a 92
#define FUSED_4a4A7a 93

#define FUSED_4a4a4A 94
#define FUSED_4a4a5A 95
#define FUSED_4a4a6A 96
#define FUSED_4a4a7A 97

#define FUSED_4a4a4a 98
#define FUSED_4a4a5a 99
#define FUSED_4a4a6a 100
#define FUSED_4a4a7a 101

#define FUSED_4a5A5A 102
#define FUSED_4a5A6A 103
#define FUSED_4a5A7A 104

#define FUSED_4a5A5a 105
#define FUSED_4a5A6a 106
#define FUSED_4a5A7a 107

#define FUSED_4a5a5A 108
#define FUSED_4a5a6A 109
#define FUSED_4a5a7A 110

#define FUSED_4a5a5a 111
#define FUSED_4a5a6a 112
#define FUSED_4a5a7a 113

#define FUSED_4a6A6A 114
#define FUSED_4a6A7A 115

#define FUSED_4a6A6a 116
#define FUSED_4a6A7a 117

#define FUSED_4a7A7A 118
#define FUSED_4a7A7a 119

#define FUSED_4a7a7A 120
#define FUSED_4a7a7a 121

// Must be the number of ring types above - one larger than the largest enum...

#define RING_ARRAY_SIZE 122

class Ring_Extraction_Replacement_Conditions
{
  protected:
    int _ring_size_needed;

    int _ring_aromaticity_needed;

    isotope_t _isotope_for_ring_fusion;

    const Element * _ring_fusion_element;

    isotope_t _isotope_for_substitution_points;

    int _include_substituents;

    int _remove_chirality;

    int _only_process_fused_rings;

    int _only_process_unfused_rings;

    int _fused_neighbours_allowed;

//  private functions
   
  private:

    void _add_substituents (Molecule & m, 
                  atom_number_t ring_atom,
                  atom_number_t outside_ring,
                  int * include_atom) const;

  public:
    Ring_Extraction_Replacement_Conditions ();
    ~Ring_Extraction_Replacement_Conditions ();

    int initialise (Command_Line &, int);

    int remove_chirality() const { return _remove_chirality;}

    int can_be_processed (Molecule & m, const Ring &) const;

    int isotope_for_substitution_points() const { return _isotope_for_substitution_points;}
    int isotope_for_ring_fusion() const { return _isotope_for_ring_fusion;}
    
    int identify_atoms_associated_with_ring (Molecule & m,
                                              const Ring & r,
                                              const int * in_same_ring,
                                              int * include_atom) const;
    int identify_atoms_associated_with_ring_system (Molecule & m,
                                              int * include_atom) const;
    int append_connectivity_smarts (Molecule & m, atom_number_t zatom,
                                    const int * include_atom,
                                    int aromatic,
                                    IWString & smarts) const;

    void append_connectivity_smarts (Molecule & m, atom_number_t zatom, int aromatic, IWString & smarts) const;
};

// When producing partial smiles or smarts, we have lots of 
// info to be passed, so we build a class

int display_standard_ring_ext_rep_options (std::ostream & os);

int initialise_in_same_ring_array (Molecule & m, int * in_same_ring);

uint32_t RingHash(const extending_resizable_array<int>& aliph, const extending_resizable_array<int>& arom);

// Generate a name like '4A5a' from lists of aliphatic and aromatic ring sizes.
IWString RingHashName(const extending_resizable_array<int>& aliph, const extending_resizable_array<int>& arom);
// same as RingHashName, except now we can specify the 'a' and 'A' strings.
IWString RingHashName(const extending_resizable_array<int>& aliph, 
                      const IWString& aliph_suffix,
                      const extending_resizable_array<int>& arom,
                      const IWString& arom_suffix);

}  // namespace ring_replacement
#endif  // MOLECULE_TOOLS_RING_EXT_REP_H
