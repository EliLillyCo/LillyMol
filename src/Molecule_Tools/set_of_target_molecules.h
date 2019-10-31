#ifndef SET_OF_TARGET_MOL_H
#define SET_OF_TARGET_MOL_H

#include "accumulator.h"

#include "atom_typing.h"
#include "mpr.h"
#include "maccskeys_fn5.h"
#include "istream_and_type.h"
#include "linear_path_fingerprint.h"


/*
  The various random molecule creation tools might be trying to get as
  close as possible to one or more target molecules
*/

template <typename F>
class Set_of_Target_Molecules
{
  private:
    resizable_array_p<Molecule> _m;

    F * _fp;

    Molecular_Properties_Generator _mpr;
    MACCSKeys _mk;

//  for the linear fingerprints

    Atom_Typing_Specification _atom_typing_specification;
    LFP::Linear_Fingerprint_Defaults _lfpd;
    LFP::Linear_Fingerprint_Creator _lfp;

    int _tmp[NUMBER_MACCS_KEYS];    // working storage, no thread safety here

//  Once we have found a match to a given molecule, we might want to "turn it off"

    float _turn_off_if_within;    // if we achieve a match within this distance, turn that molecule off

    int * _match_already_found;

//  private functions

    void _compute_fingerprint(Molecule & m, F & fp);

  public:
    Set_of_Target_Molecules();
    ~Set_of_Target_Molecules();

    int build (const char * fname);
    int build (data_source_and_type<Molecule> & input);

    int set_turn_off_molecules_matched_to_within(const float dist);   // matched 

    int number_molecules() const { return _m.number_elements();}

    float closest_distance(F &) const;
    float closest_distance(Molecule &);
    void closest_distance(Molecule & m, Accumulator<float> &);
};

#endif
