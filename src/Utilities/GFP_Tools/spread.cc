#include <iostream>

#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"

#include "spread.h"

using std::cerr;
using std::endl;

static IWString scale_tag;

static IWString smiles_tag("$SMI<");

int
set_scale_tag(const const_IWSubstring & tag)
{
  scale_tag = tag;

  if (0 == scale_tag.length())   // has been turned off
    ;
  else if (! scale_tag.ends_with('<'))
    scale_tag += '<';

  return 1;
}

static Numeric_Data_From_File<float> id_to_scale;

static int scaling_factor_column = -1;

void
set_scaling_factor_column(int c)
{
  if (c >= 1)
    scaling_factor_column = c - 1;   // convert to word index
  else if (0 == c)
    cerr << "Column zero is invalid for scaling factor column\n";
  else
    scaling_factor_column = c;

  return;
}

static int every_object_must_have_a_scale_factor = 1;

void
set_every_object_must_have_a_scale_factor(int s)
{
  every_object_must_have_a_scale_factor = s;
}

static int include_scale_of_nsn_with_scale = 0;

void
set_include_scale_of_nsn_with_scale(int s)
{
  include_scale_of_nsn_with_scale = s;
}

int
read_scaling_data(const char * fname,
                  int verbose)
{
  if (! id_to_scale.read_data(fname))
  {
    cerr << "Cannot read scaling factor data from '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Read " << id_to_scale.size() << " id->scale values\n";

  return 1;
}

/*
  We keep track of statistics on scale factors
*/

static Accumulator<float> scale_stats;

const Accumulator<float> &
scale_factor_statistics()
{
  return scale_stats;
}

Spread_Object::Spread_Object()
{
  _selected = 0;

  _scale = 1.0;

  _scale_of_nearest_selected_nbr = 1.0;

  _scaled_distance = 0.0;

  return;
}

Spread_Object &
Spread_Object::operator=(const Spread_Object & rhs)
{
  IW_General_Fingerprint::operator=(rhs);

  _nearest_selected_neighbour = rhs._nearest_selected_neighbour;

  _smiles = rhs._smiles;

  _selected = rhs._selected;
  _scale = rhs._scale;
  _scale_of_nearest_selected_nbr = rhs._scale_of_nearest_selected_nbr;
  _scaled_distance = rhs._scaled_distance;
  return *this;
}

void
Spread_Object::_update_nsn_stuff(const Spread_Object & fpsel,
                                 similarity_type_t d)
{
  _nearest_selected_neighbour.set_distance(d);
  _nearest_selected_neighbour.set_smiles(fpsel.smiles());
  _nearest_selected_neighbour.set_id(fpsel.id());

  if (include_scale_of_nsn_with_scale)
  {
    _scale_of_nearest_selected_nbr = fpsel.scale();

    _scaled_distance = _scale * _scale_of_nearest_selected_nbr * d;
  }
  else
    _scaled_distance = _scale * d;

  return;
}

int
Spread_Object::object_has_been_selected(Spread_Object & fpsel)
{
  similarity_type_t new_distance = IW_General_Fingerprint::distance(fpsel);

  if (new_distance >= _nearest_selected_neighbour.distance())
    return 0;

  _update_nsn_stuff(fpsel, new_distance);

  return 1;
}

int
Spread_Object::object_has_been_selected(Spread_Object & fpsel,
                                        float blurr_distances)
{
  similarity_type_t new_distance = IW_General_Fingerprint::distance(fpsel);

  new_distance = do_blurring(new_distance, blurr_distances);

  if (new_distance >= _nearest_selected_neighbour.distance())
    return 0;

  _update_nsn_stuff(fpsel, new_distance);

  return 1;
}

int
Spread_Object::object_has_been_selected(Spread_Object & fpsel,
                                        const Tversky & tversky)
{
  similarity_type_t new_distance = static_cast<similarity_type_t>(1.0) - IW_General_Fingerprint::tversky(fpsel, tversky);

  if (new_distance >= _nearest_selected_neighbour.distance())
    return 0;

  _update_nsn_stuff(fpsel, new_distance);

  return 1;
}

int
Spread_Object::object_has_been_selected(Spread_Object & fpsel,
                                        const Tversky & tversky,
                                        float blurr_distances)
{
  similarity_type_t new_distance = static_cast<similarity_type_t>(1.0) - IW_General_Fingerprint::tversky(fpsel, tversky);

  new_distance = do_blurring(new_distance, blurr_distances);

  if (new_distance >= _nearest_selected_neighbour.distance())
    return 0;

  _update_nsn_stuff(fpsel, new_distance);

  return 1;
}


int
Spread_Object::object_has_been_selected_max_distance(Spread_Object & fpsel,
                                         similarity_type_t max_dist)
{
  similarity_type_t new_distance = IW_General_Fingerprint::distance(fpsel);

  if (new_distance > max_dist)
    new_distance = max_dist;

  if (new_distance >= _nearest_selected_neighbour.distance())
    return 0;

  _update_nsn_stuff(fpsel, new_distance);

  return 1;
}

int
Spread_Object::_determine_scale_from_column(int c)
{
  if (_id.nwords() > c)   // great, our token is there
    ;
  else if (every_object_must_have_a_scale_factor)
  {
    cerr << "Spread_Object::_determine_scale_from_column:no column data '" << _id << "'\n";
    cerr << "Column " << c << endl;
    return 0;
  }
  else
  {
    _scale = static_cast<float>(1.0);
    return 1;
  }

  const_IWSubstring token = _id.word(c);

  if (! token.numeric_value(_scale) || _scale < 0.0)   // even though a zero value makes no sense.
  {
    cerr << "Spread_Object::_determine_scale_from_column:invalid scaling value '" << _id << "'\n";
    return 0;
  }

//cerr << "Converted '" << token << "' to " << _scale << endl;

  scale_stats.extra(_scale);

  return 1;
}

int
Spread_Object::_determine_scale_from_hash(const IW_STL_Hash_Map_float & id_to_scale)
{
  IW_STL_Hash_Map_float::const_iterator f;
  
  if (_id.nwords() > 1)
  {
    IWString tmp(_id);
    tmp.truncate_at_first(' ');
    f = id_to_scale.find(tmp);
  }
  else
    f = id_to_scale.find(_id);

  if (f != id_to_scale.end())
  {
    _scale = (*f).second;

    scale_stats.extra(_scale);

    return 1;
  }

  if (every_object_must_have_a_scale_factor)
  {
    cerr << "No scaling data for '" << _id << "'\n";
    return 0;
  }

  return 1;
}

int
Spread_Object::_determine_scale_from_tag(const IW_TDT & tdt,
                                         const IWString & scale_tag)
{
  const_IWSubstring s;
  if (! tdt.dataitem_value(scale_tag, s))   // no scale tag present
    return 1;

  if (! s.numeric_value(_scale) || _scale <= 0.0)
  {
    cerr << "Spread_Object::create_from_tdt: invalid scale/demerit value\n";
    return 0;
  }

  scale_stats.extra(_scale);

  return 1;
}

int
Spread_Object::construct_from_tdt(IW_TDT & tdt, int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal))
    return 0;

  _scale = 1.0;

  if (! tdt.dataitem_value(smiles_tag, _smiles))
  {
    cerr << "Spread_Object::construct_from_tdt: no smiles tag '" << smiles_tag << "' in TDT\n";
    return 0;
  }

  if (scale_tag.length())
  {
    if (! _determine_scale_from_tag(tdt, scale_tag))
    {
      fatal = 1;
      return 0;
    }
  }
  else if (id_to_scale.size())
  {
    if (! _determine_scale_from_hash(id_to_scale))
    {
      fatal = 1;
      return 0;
    }
  }
  else if (scaling_factor_column > 0)
  {
    if (! _determine_scale_from_column(scaling_factor_column))
    {
      fatal = 1;
      return 0;
    }
  }

//cerr << "Object " << _id << " with scale " << _scale << endl;

  return 1;
}

void
Spread_Object::set_nearest_previously_selected_neighbour(const IWString & nnsmiles,
                            const IWString & nnid,
                            similarity_type_t d)
{
  _nearest_selected_neighbour.set_smiles(nnsmiles);
  _nearest_selected_neighbour.set_id(nnid);
  _nearest_selected_neighbour.set_distance(d);

  _scaled_distance = d;
  _scale_of_nearest_selected_nbr = d;

  return;
}

similarity_type_t
do_blurring(similarity_type_t d,
            float blurr_distances)
{
  float tmp = d * blurr_distances + static_cast<float>(0.499999);

  tmp = static_cast<float>(static_cast<int>(tmp));

  return tmp / blurr_distances;
}

int
Spread_Object::set_distance_to_previously_selected_from_column(int col)
{
  const_IWSubstring token;

  if (! _id.word(col, token))
  {
    cerr << "Spread_Object::set_distance_to_previously_selected_from_column:cannot extract column " << (col+1) << endl;
    return 0;
  }

  similarity_type_t d;
  if (! token.numeric_value(d) || d < 0.0)
  {
    cerr << "Spread_Object::set_distance_to_previously_selected_from_column:invalid distance " << token << "'\n";
    return 0;
  }

  _nearest_selected_neighbour.set_distance(d);

  return 1;
}
