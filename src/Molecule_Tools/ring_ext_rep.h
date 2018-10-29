#ifndef RING_EXT_REP_H
#define RING_EXT_REP_H

#include <fstream>

#include "cmdline.h"

#include "molecule.h"

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

// Must be the number of ring types above - one larger than the largest enum...

#define RING_ARRAY_SIZE 54

class Ring_Extraction_Replacement_Conditions
{
  protected:
    int _ring_size_needed;

    int _ring_aromaticity_needed;

    int _isotope_for_ring_fusion;

    const Element * _ring_fusion_element;

    int _isotope_for_substitution_points;

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

/*
  When producing partial smiles or smarts, we have lots of 
  info to be passed, so we build a class
*/

extern int display_standard_ring_ext_rep_options (std::ostream & os);

extern int
initialise_in_same_ring_array (Molecule & m, int * in_same_ring);
#endif
