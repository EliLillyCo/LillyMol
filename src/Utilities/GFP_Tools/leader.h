#ifndef GFP_LEADER_H
#define GFP_LEADER_H

#include "gfp.h"
#include "set_or_unset.h"

typedef float score_t;

class GFP_L : public IW_GFP_D
{
  private:
    int _selected;

    Set_or_Unset<similarity_type_t> _threshold;

//  Each item can limit the number of items in its cluster

    Set_or_Unset<int> _max_cluster_size;

    score_t _score;

//  shortest distance to a previously selected cluster centre

    similarity_type_t _shortest_distance_to_cluster_centre;

  public:
    GFP_L ();

    int selected () const { return _selected;}
    int & selected () { return _selected;}
    void set_selected (int s) { _selected = s;}
//  void set_selected (similarity_type_t s) { _selected = 1; _distance = s;}

    score_t score () const { return _score;}

    similarity_type_t shortest_distance_to_cluster_centre () const { return _shortest_distance_to_cluster_centre;}
    void set_shortest_distance_to_cluster_centre (similarity_type_t d) { _shortest_distance_to_cluster_centre = d;}

    int construct_from_tdt (IW_TDT &, int &);

    int threshold (similarity_type_t & t) const { return _threshold.value (t);}
    int set_threshold (similarity_type_t t) { return _threshold.set (t);}

    int max_cluster_size (int & m) const { return _max_cluster_size.value (m);}
};

#endif
