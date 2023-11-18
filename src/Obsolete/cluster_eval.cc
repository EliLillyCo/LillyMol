#include <stdlib.h>
#include <iostream>

#include "Foundational/iw_tdt/iw_tdt.h"

#include "cluster_eval.h"

using std::cerr;
using std::endl;

static IWString cluster_id_tag ("CLUSTER");

int
set_cluster_id_tag (const const_IWSubstring & t)
{
  if (0 == t.length ())
  {
    cluster_id_tag.resize (0);
    return 1;
  }

  cluster_id_tag.resize_keep_storage (0);

  cluster_id_tag = t;
  if (! cluster_id_tag.ends_with ('<'))
    cluster_id_tag += '<';

  return 1;
}

GFP_C::GFP_C ()
{
  _cluster = -1;
}

int
GFP_C::construct_from_tdt (IW_TDT & tdt, int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt (tdt, fatal))
  {
    return 0;
  }

  if (0 == cluster_id_tag.length ())
    return 1;

  if (! tdt.dataitem_value (cluster_id_tag, _cluster) || _cluster < 0)
  {
    cerr << "GFP_C::construct_from_tdt: missing or invalid '" << cluster_id_tag << "' value\n";
    cerr << tdt;
    return 0;
  }

  return 1;
}
