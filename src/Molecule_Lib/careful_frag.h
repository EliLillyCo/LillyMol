#ifndef MOLECULE_LIB_CAREFUL_FRAG_H_
#define MOLECULE_LIB_CAREFUL_FRAG_H_
/*
  Some extra private functions needed by careful_frag
*/


    int _reduce_to_largest_fragment_carefully (Fragment_Data * fc, int * already_counted);
    int _is_nitro (atom_number_t, int *) const;
    int _is_sulphate_like (atom_number_t, int *) const;
    int _identify_fragment_undesirable_groups (int * exclude) const;

#endif  // MOLECULE_LIB_CAREFUL_FRAG_H_
