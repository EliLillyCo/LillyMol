#ifndef LEADER_RESULTS_H
#define LEADER_RESULTS_H

#include <iostream>

#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::ostream;
using std::cerr;
using std::endl;
using std::pair;

class IWString_STL_Hash_Set;

/*
  Several programmes work on the output from gfp_nearneighbours
*/


/*
  The template'd object N, will typically include smiles, id and distance.

  The items that own the neighbours
*/

template <typename N>
class NN_Item_Base : public resizable_array_p<N>
{
  protected:
    IWString _smiles;
    IWString _id;

#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array_p<N>::_number_elements;
    using resizable_array_p<N>::_things;
#endif

  public:

    const IWString & smiles () const { return _smiles;}

    const IWString & id () const { return _id;}
    void  set_id (const IWString & s) { _id = s;}

    int debug_print (ostream &, int) const;

    int build (iwstring_data_source & input);

//  Version for when we have either a max or min distance specified

    int build (iwstring_data_source & input, similarity_type_t, similarity_type_t);

    int build_from_neighbour_list (const pair<int, similarity_type_t> *, int);

    template <typename F> void change_distances (F &);
};

template <typename N>
class NN_Results_Base
{
  protected:
    N * _results;
    int _number_items;

    similarity_type_t _min_distance_to_keep;
    similarity_type_t _max_distance_to_keep;

//  private functions

    int _resize ();

    int _count_items_in_input_and_resize (iwstring_data_source & input);

    int _build_from_distance_matrix (IWDistanceMatrixFloat & dm,
                                     const DM_to_NN_Conditions<float> & dmc,
                                     const int * only_use,
                                     pair<int, float> * p,
                                     similarity_type_t *);

  public:
    NN_Results_Base ();
    ~NN_Results_Base ();

    int number_results () const { return _number_items;}
    int number_elements () const { return _number_items;}

    void set_min_distance_to_keep (similarity_type_t m) {_min_distance_to_keep = m;}
    void set_max_distance_to_keep (similarity_type_t m) {_max_distance_to_keep = m;}

    int debug_print (ostream & os, int verbose) const;

    const N & operator [] (int i) const { return _results[i];}
    N & operator [] (int i) { return _results[i];}

    int build (const char * fname);
    int build (iwstring_data_source & input);

    int build_from_distance_matrix (const char * fname, const DM_to_NN_Conditions<float> &, const IWString_STL_Hash_Set & only_use);
    int build_from_distance_matrix (iwstring_data_source & input, const DM_to_NN_Conditions<float> &, const IWString_STL_Hash_Set & only_use);

    int build_from_gfp_nearneighbours_dash_o (const char * fname);
    int build_from_gfp_nearneighbours_dash_o (iwstring_data_source & input);

    template <typename F> void change_distances (F & f);
};

#endif
