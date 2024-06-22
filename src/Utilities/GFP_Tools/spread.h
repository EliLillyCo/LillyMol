#ifndef SPREADV2_H
#define SPREADV2_H

/*
  All the variants of spread need a class like this
*/

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "gfp.h"
#include "smiles_id_dist.h"

class Tversky;

class Spread_Object : public IW_General_Fingerprint
{
  private:
    int _selected;

    Smiles_ID_Dist _nearest_selected_neighbour;

//  Our own smiles

    IWString _smiles;

    float _scale;     // the demerit value

    float _scale_of_nearest_selected_nbr;

    float _scaled_distance;

//  private functions

    int _determine_scale_from_hash (const IW_STL_Hash_Map_float & id_to_scale);
    int _determine_scale_from_tag (const IW_TDT & tdt, const IWString & scale_tag);
    int _determine_scale_from_column (int c);

    void _update_nsn_stuff (const Spread_Object & fpsel, similarity_type_t d);

  public:
    Spread_Object ();

    Spread_Object & operator= (const Spread_Object & rhs);

    void set_scale (float s) { assert (s > 0.0); _scale = s;}
    float scale () const { return _scale;}

//  similarity_type_t distance () const { return _scale * _nearest_selected_neighbour.distance (); }
    similarity_type_t distance () const { return _scaled_distance;}
    void set_distance (similarity_type_t d) { _nearest_selected_neighbour.set_distance (d); _scaled_distance = _scale * d;}

    int construct_from_tdt (IW_TDT &, int &);

    const IWString & smiles () const { return _smiles;}

    int  selected () const { return _selected;}
    void set_selected () { assert (! _selected);  _selected = 1;}
    void set_selected (int s) { _selected = s;}

    const Smiles_ID_Dist & nsn () const { return _nearest_selected_neighbour;}
    int has_a_nearest_selected_neighbour () const { return _nearest_selected_neighbour.id ().length ();}

    void set_nearest_previously_selected_neighbour (const IWString &, const IWString &, similarity_type_t);

    int set_distance_to_previously_selected_from_column (int);

    int object_has_been_selected (Spread_Object &);
    int object_has_been_selected (Spread_Object &, float blurr_distances);
    int object_has_been_selected (Spread_Object &, const Tversky & tv);
    int object_has_been_selected (Spread_Object &, const Tversky & tv, float blurr_distances);
    int object_has_been_selected_max_distance (Spread_Object &, float max_dist);
};

extern int set_scale_tag (const const_IWSubstring &);
extern int read_scaling_data (const char *, int verbose);
extern void set_every_object_must_have_a_scale_factor (int s);
extern void set_scaling_factor_column (int c);

extern const Accumulator<float> & scale_factor_statistics ();

extern similarity_type_t do_blurring (similarity_type_t d, float);

extern void set_include_scale_of_nsn_with_scale (int);

#endif
