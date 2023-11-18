#ifndef IW_GFP_CLUSTER_J
#define IW_GFP_CLUSTER_J

#include "Utilities/GFP_Tools/gfp.h"

class GFP_C : public IW_General_Fingerprint
{
  private:
    int _cluster;

  public:

    GFP_C ();

    int cluster () const { return _cluster;}

    void set_cluster (int c) { _cluster = c;}

    int construct_from_tdt (IW_TDT &, int &);
};

extern int set_cluster_id_tag (const const_IWSubstring &);

#endif
