/*
  A spread variant where items must be selected in groups.
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Utilities/GFP_Tools/smiles_id_dist.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = nullptr;

static int verbose = 0;

static int groups_to_select = 0;

static int groups_selected = 0;

static int items_selected = 0;    // each group may contain multiple items

static similarity_type_t too_close = static_cast<similarity_type_t>(0.0);

static similarity_type_t chop_long_distances_at = static_cast<similarity_type_t> (1.0);

static int previously_selected_file_present = 0;

static int previously_selected_items_processed = 0;

static int report_processing_previously_selected = 0;

static IWString precomputed_distance_to_previously_selected_tag;

static int dist_to_prev_column = -1;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Spread variant where each item has a number of neighbours that also get selected\n";
  cerr << "Input file has sets of identifiers what are grouped together\n";
  cerr << " -A <fname>     fingerprint file of already selected objects\n";
  cerr << " -r <number>    report progress of the previously selected items\n";
  cerr << " -I <fname>     fingerprint file for items to be selected\n";
  cerr << " -R <tag>       distance to previously selected items in <tag> in -I file\n";
  cerr << " -R col=<col>   distance to previously selected items column <col> of name field\n";
  cerr << " -t <dist>      define 'zero' distance. These are avoided\n";
  cerr << " -T <dist>      truncate all long distances to <dist>\n";
  cerr << " -n <number>    number of groups to select\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

class Default_Desirability_Measure
{
  private:
    int _nzero;    // number of zero distance instances

    Accumulator<similarity_type_t> _acc;

    double _scale;

  public:
    Default_Desirability_Measure();

    int debug_print(ostream &) const;

    void extra (similarity_type_t);

    int nzero() const { return _nzero;}

    void set_scale(double s) { _scale = s;}

    int is_set() const { return (_nzero > 0 || _acc.n() > 0);}

    int compare (const Default_Desirability_Measure &) const;

    void reset();
};

Default_Desirability_Measure::Default_Desirability_Measure()
{
  _nzero = 0;

  _scale = 1.0;

  return;
}

void
Default_Desirability_Measure::reset()
{
  _nzero = 0;
  _acc.reset();

  return;
}

int
Default_Desirability_Measure::debug_print(ostream & os) const
{
  os << "Desirability: scale " << _scale << " nzero " << _nzero;
  if (0 == _acc.n())
    ;
  else if (1 == _acc.n())
    os << " dist " << _acc.minval();
  else
  {
    os << " dist btw " << _acc.minval() << " and " << _acc.maxval() << " ave " << _acc.average_if_available_minval_if_not();
  }

  os << '\n';

  return os.good();
}

void
Default_Desirability_Measure::extra(similarity_type_t d)
{
  if (d <= too_close)
    _nzero++;
  else if (d > chop_long_distances_at)
    _acc.extra(chop_long_distances_at);
  else
    _acc.extra(d);

  return;
}

/*
  Return -1 if THIS is better, 1 if RHS is better
*/

int
Default_Desirability_Measure::compare(const Default_Desirability_Measure & rhs) const
{

  if (_nzero == rhs._nzero)
    ;
  else if (_nzero < rhs._nzero)
    return -1;
  else 
    return 1;

// The two items have the same number of zero distance items

  if (0 == _acc.n())
  {
    if (0 == rhs._acc.n())
      return 0;

    return 1;
  }
  else if (0 == rhs._acc.n())
    return -1;

// We want the longest possible average distance to previously selected items

  assert (_acc.n() > 0);
  assert (rhs._acc.n() > 0);

  double d1 = _acc.average_if_available_minval_if_not() * _scale;
  double d2 = rhs._acc.average_if_available_minval_if_not() * rhs._scale;

  if (d1 > d2)
    return -11;
  else if (d1 < d2)
    return 1;

  return 0;
}

/*
  The previously selected items will be simple FP_and_Smiles objects
*/

/* Ingvar: Class definition  moved elsewhere
class FP_and_Smiles : public IW_General_Fingerprint
{
  protected:
    IWString _smiles;

  public:
    int construct_from_tdt (IW_TDT &, int &);

    const IWString & smiles() const { return _smiles;}
};
*/

int
FP_and_Smiles::construct_from_tdt (IW_TDT & tdt,
                                   int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt (tdt, fatal))
    return 0;

  if (! tdt.dataitem_value(smiles_tag, _smiles))
  {
    cerr << "FP_and_Smiles::construct_from_tdt:cannot find smiles in\n";
    cerr << tdt;
    return 0;
  }

  return 1;
}

class Gspread_Object : public FP_and_Smiles
{
  private:
    int _selected;

    Smiles_ID_Dist _nearest_selected_neighbour;

  public:
    Gspread_Object();

    int selected() const { return _selected;}
    void set_selected () { _selected = 1;}

    int object_has_been_selected (FP_and_Smiles &);

    similarity_type_t distance() const { return _nearest_selected_neighbour.distance();}

    void set_distance_to_previously_selected(similarity_type_t d);

    int write_smiles_and_id(ostream &) const;
    int write_nearest_previously_selected(ostream &) const;
};

static Gspread_Object * pool = nullptr;

static int pool_size = 0;

Gspread_Object::Gspread_Object()
{
  _selected = 0;

  return;
}

int
Gspread_Object::object_has_been_selected(FP_and_Smiles & fp)
{
  similarity_type_t new_distance = IW_General_Fingerprint::distance (fp);

  if (new_distance >= _nearest_selected_neighbour.distance ())
    return 0;

//cerr << "Updating near neighbour distance\n";

  _nearest_selected_neighbour.set_distance (new_distance);
  _nearest_selected_neighbour.set_smiles (fp.smiles ());
  _nearest_selected_neighbour.set_id (fp.id ());

  return 1;
}

int
Gspread_Object::write_smiles_and_id(ostream & os) const
{
  os << smiles_tag << _smiles << ">\n";
  os << identifier_tag << id() << ">\n";

  return os.good();
}

int
Gspread_Object::write_nearest_previously_selected(ostream & os) const
{
  os << smiles_tag << _nearest_selected_neighbour.smiles() << ">\n";
  os << identifier_tag << _nearest_selected_neighbour.id() << ">\n";
  os << distance_tag << _nearest_selected_neighbour.distance() << ">\n";
  
  return os.good();
}

void
Gspread_Object::set_distance_to_previously_selected(similarity_type_t d)
{
  _nearest_selected_neighbour.set_distance(d);

  return;
}

/*
  Really should be a template for the desirability measure, but this is quicker...
  Besides, we'll probably never have different desirability measures
*/

class Spread_Group
{
  private:
    Gspread_Object ** _fp;
    int _nfp;

    int _selected;

    Default_Desirability_Measure _desirability_measure;

  public:
    Spread_Group();

    int debug_print (ostream &) const;

    int build(const const_IWSubstring & buffer, const IW_STL_Hash_Map_int & id_to_ndx);

    int number_fingerprints() const { return _nfp;}

    int unselected_items() const;

    int items_in_group() const { return _nfp;}

    int selected() const { return _selected;}
    void set_selected();

    int group_has_been_selected(const Spread_Group &);

    int compare (const Spread_Group & rhs) const;

//  After we process the -A file, we need to establish the desirability measure

    int establish_desirability_measure();

    int write(ostream &) const;

    int output_for_nplotnn(ostream &) const;
};

static Spread_Group * group = nullptr;
static int ngroup = 0;

Spread_Group::Spread_Group()
{
  _fp = nullptr;
  _nfp = 0;

  _selected = 0;

  return;
}

int
Spread_Group::debug_print(ostream & os) const
{
  os << "Group with " << _nfp << " fingerprints, " << unselected_items() << " unselected\n";

  return _desirability_measure.debug_print (os);
}

int
Spread_Group::unselected_items() const
{
  int rc = 0;

  for (int i = 0; i < _nfp; i++)
  {
    if (! _fp[i]->selected())
      rc++;
  }

  return rc;
}

int
Spread_Group::compare(const Spread_Group & rhs) const
{
  int u1 = unselected_items();
  int u2 = rhs.unselected_items();

  if (0 == u1 && 0 == u2)
    return 0;

  if (0 == u1)
    return 1;
  
  if (0 == u2)
    return -1;

// the number of items that have finite distances

  int f1 = _nfp - _desirability_measure.nzero();
  int f2 = rhs._nfp - rhs._desirability_measure.nzero();

  if (0 == f1 && 0 == f2)
    return 0;

  if (0 == f1)
    return 1;

  if (0 == f2)
    return -1;

  return _desirability_measure.compare(rhs._desirability_measure);
}

int
Spread_Group::build(const const_IWSubstring & buffer,
                    const IW_STL_Hash_Map_int & id_to_ndx)
{
  _nfp = buffer.nwords();
  if (_nfp < 1)
  {
    cerr << "FP_and_NBRS::build:empty input\n";
    return 0;
  }

  _fp = new Gspread_Object *[_nfp];

  if (NULL == _fp)
  {
    cerr << "Spread_Group::build:cannot allocate " << _nfp << " fingerprint pointers\n";
    return 0;
  }

  _nfp = 0;

  IWString token;
  int i = 0;

  while (buffer.nextword(token, i))
  {
    if (token.starts_with("SCALE="))
    {
      token.remove_leading_chars(6);
      double scale;
      if (! token.numeric_value(scale) || scale <= 0.0)
      {
        cerr << "Spread_Group::build:invalid scale specification '" << token << "'\n";
        return 0;
      }

      _desirability_measure.set_scale(scale);

      continue;
    }

    IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find(token);

    if (f == id_to_ndx.end())
    {
      cerr << "Spread_Group::build:no fingerprint for '" << token << "'\n";
      return 0;
    }

    int k = (*f).second;

    _fp[_nfp] = &(pool[k]);
    _nfp++;
  }

  return _nfp;
}

void
Spread_Group::set_selected()
{
  for (int i = 0; i < _nfp; i++)
  {
    if (! _fp[i]->selected())    // may be owned by some other group
      _fp[i]->set_selected();
  }

  _selected = 1;

  return;
}

int
Spread_Group::group_has_been_selected(const Spread_Group & rhs)
{
  for (int i = 0; i < _nfp; i++)
  {
    Gspread_Object * fpi = _fp[i];
    if (fpi->selected())
      continue;

//  cerr << "i = " << i << " initial distance " << fp->distance() << endl;

    for (int j = 0; j < rhs._nfp; j++)
    {
      const Gspread_Object * fpj = rhs._fp[j];

      if (fpj->selected())
        continue;

      if (fpi == fpj)   // otherwise we get everything with 0 distance
        continue;

#ifdef DEBUG_GROUP_HAS_BEEN_SELECTED
      cerr << "D = " << fp->IW_General_Fingerprint::distance(*(rhs._fp[j])) << endl;
#endif
      fpi->object_has_been_selected(*(rhs._fp[j]));
    }
  }

#ifdef DEBUG_GROUP_HAS_BEEN_SELECTED
  cerr << "After scanning " << distances_updated << " of " << _nfp << " distances updated, measure " << _desirability_measure.is_set() << endl;
#endif

  establish_desirability_measure();

  return 1;
}

int
Spread_Group::write(ostream & os) const
{
  os << "'";
  for (int i = 0; i < _nfp; i++)
  {
    if (i > 0)
      os << ' ';

    os << _fp[i]->id();
  }

  os << "'\n";
  _desirability_measure.debug_print(os);

  return os.good();
}

/*
  nplotnn has a hard time with a target that is multiple molecules, so
  we create a multi-fragment smiles in those cases
*/

int
Spread_Group::output_for_nplotnn(ostream & os) const
{
  if (1 == _nfp)
    _fp[0]->write_smiles_and_id(os);
  else
  {
    IWString smiles, id;        // combined form
    for (int i = 0; i < _nfp; i++)
    {
      const Gspread_Object * f = _fp[i];

      if (0 == i)
      {
        smiles = f->smiles();
        id = f->id();
      }
      else
      {
        smiles.append_with_spacer(f->smiles(), '.');
        id.append_with_spacer(f->id(), '+');
      }
    }

    os << smiles_tag << smiles << ">\n";
    os << identifier_tag << id << ">\n";
  }

  for (int i = 0; i < _nfp; i++)
  {
    const Gspread_Object * f = _fp[i];
    f->write_nearest_previously_selected(os);
  }

  os << "|\n";

  return os.good();
}

int
Spread_Group::establish_desirability_measure()
{
  _desirability_measure.reset();

  for (int i = 0; i < _nfp; i++)
  {
    if (! _fp[i]->selected())
      _desirability_measure.extra(_fp[i]->distance());
  }

  return 1;
}

static int
process_previously_selected_item (FP_and_Smiles & fp)
{
  for (int i = 0; i < pool_size; i++)
  {
    pool[i].object_has_been_selected(fp);
  }

  return 1;
}

static int
process_previously_selected_items(iwstring_data_source & input)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    FP_and_Smiles fp;
    int fatal;
    if(! fp.construct_from_tdt(tdt, fatal))
    {
      if (! fatal)
        continue;

      cerr << "Fatal error reading previously selected fingerprint, line " << input.lines_read() << endl;
      return 0;
    }

    process_previously_selected_item (fp);

    previously_selected_items_processed++;

    if (0 == report_processing_previously_selected)
      ;
    else if (0 == previously_selected_items_processed % report_processing_previously_selected)
      cerr << "Processed " << previously_selected_items_processed << " previously selected items\n";
  }

  return 1;
}

static int
process_previously_selected_items(const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open previously selected file '" << fname << "' for input\n";
    return 0;
  }

  return process_previously_selected_items(input);
}

static int
determine_distance_to_previously_selected_dataitem(Gspread_Object & gs,
                                                   const IW_TDT & tdt)
{
  similarity_type_t d;
  if (! tdt.dataitem_value(precomputed_distance_to_previously_selected_tag, d))
  {
    cerr << "Yipes, cannot extract '" << precomputed_distance_to_previously_selected_tag << "' value from tdt\n";
    return 0;
  }

  gs.set_distance_to_previously_selected(d);

  return 1;
}

static int
determine_distance_to_previously_selected_name_col(Gspread_Object & gs)
{
  const IWString & id = gs.id();

  const_IWSubstring token;
  if (! id.word(dist_to_prev_column, token))
  {
    cerr << "Cannot extract column " << (dist_to_prev_column + 1) << " from '" << id << "'\n";
    return 0;
  }

  float d;
  if (! token.numeric_value(d) || d < static_cast<float> (0.0))
  {
    cerr << "Invalid distance '" << id << "'\n";
    return 0;
  }

  gs.set_distance_to_previously_selected(d);

  return 1;
}

static int
determine_distance_to_previously_selected(Gspread_Object & gs,
                                          const IW_TDT & tdt)
{
  if (precomputed_distance_to_previously_selected_tag.length())
    return determine_distance_to_previously_selected_dataitem(gs, tdt);
  else
    return determine_distance_to_previously_selected_name_col(gs);
}

static int
build_pool (iwstring_data_source & input)
{
  assert (pool_size > 0);

  int i = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    if (! pool[i].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    if (dist_to_prev_column > 0 || precomputed_distance_to_previously_selected_tag.length())
      determine_distance_to_previously_selected(pool[i], tdt);

    if (pool[i].id().contains(' '))
    {
      const_IWSubstring tmp(pool[i].id());
      tmp.truncate_at_first(' ');
      pool[i].set_id(tmp);
    }

    i++;

    if (i >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }
  }

  if (verbose)
    cerr << i << " fingerprint objects added to pool\n";

  pool_size = i;


  return 1;
}

static int
build_pool (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)     // must grep the file to find out how many
  {
    assert (NULL == pool);

    IWString tmp = "^|";

    pool_size = input.grep (tmp);
    if (0 == pool_size)
    {
      cerr << "Yipes, cannot find any '" << tmp << "' in the input\n";
      return 0;
    }

    if (verbose)
      cerr << "Input contains " << pool_size << " fingerprints\n";

    pool = new Gspread_Object[pool_size];

    if (NULL == pool)
    {
      cerr << "Yipes, cannot allocate space for " << pool_size << " fingerprints\n";
      return 0;
    }
  }

  return build_pool (input);
}

int
read_groups(iwstring_data_source & input,
            const IW_STL_Hash_Map_int & id_to_ndx)
{
  int ndx = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! group[ndx].build(buffer, id_to_ndx))
    {
      cerr << "Invalid input '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    ndx++;
  }

  return ndx;
}

int
read_groups (const char * fname,
             const IW_STL_Hash_Map_int & id_to_ndx)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == group)
  {
    ngroup = input.records_remaining();
    if (ngroup < 2)
    {
      cerr << "Group file must contain more than 1 items\n";
      return 0;
    }

    group = new Spread_Group[ngroup];
    if (NULL == group)
    {
      cerr << "Cannot allocate " << ngroup << " group items\n";
      return 0;
    }

    if (verbose)
      cerr << "Identifier file contains data on " << ngroup << " groups\n";
  }

  return read_groups (input, id_to_ndx);
}

static int
choose_next_group()
{
  int best = -1;

  for (int i = 0; i < ngroup; i++)
  {
    Spread_Group & g = group[i];

    if (g.selected())
      continue;

    if (best < 0)
    {
      if (g.unselected_items())
        best = i;
    }
    else if (group[i].compare(group[best]) < 0)
      best = i;
  }

  return best;
}

static int
group_spread(ostream & os)
{
  if (0 == groups_to_select)
    groups_to_select = ngroup;

  int next_group_selected;
  if (previously_selected_file_present)
    next_group_selected = choose_next_group();
  else
    next_group_selected = 0;

  while (groups_selected < groups_to_select &&
         next_group_selected >= 0)
  {
    if(verbose > 1)
      cerr << "Next item selected is " << next_group_selected << ", " << group[next_group_selected].number_fingerprints() << " items\n";

    Spread_Group & sg = group[next_group_selected];

    items_selected += sg.unselected_items();

    for (int i = 0; i < ngroup; i++)
    {
      if (group[i].selected())
        continue;

      if (i != next_group_selected)
        group[i].group_has_been_selected(sg);
    }

    if (verbose > 2)
      sg.debug_print(cerr);

    sg.set_selected();    // very important to do this after group_has_been_selected()
    groups_selected++;
    cerr << groups_selected << " Selected ";
    sg.write(cerr);
    sg.output_for_nplotnn(os);

    next_group_selected = choose_next_group();
  }

  return os.good();
}

static int
group_spread (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vI:F:P:W:Q:G:V:n:A:r:T:t:R:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage (23);
    }
  }

  if (! cl.option_present('I'))
  {
    cerr << "Must specify fingerprint file via the -I option\n";
    usage(4);
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring fname = cl.string_value('I');
    if (! build_pool(fname))
    {
      cerr << "Cannot read fingerprint file '" << fname << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Read " << pool_size << " fingerprints from '" << fname << "'\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', too_close) || too_close < 0.0)
    {
      cerr << "The too close value (-t) must be a valid distance\n";
      usage(4);
    }

    if (verbose)
      cerr << "Distances less than " << too_close << " will be avoided as much as possible\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', chop_long_distances_at) || chop_long_distances_at <= too_close)
    {
      cerr << "The chop long distances value (-T) must be a valid distance\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will truncate distances greater than " << chop_long_distances_at << endl;
  }

  if (cl.option_present('R'))
  {
    const_IWSubstring r = cl.string_value('R');

    if (r.starts_with("col="))
    {
      r.remove_leading_chars(4);
      if (! r.numeric_value(dist_to_prev_column) || dist_to_prev_column < 2)
      {
        cerr << "The distance to previously selected column directive must be a valid column number\n";
        usage(3);
      }

      if (verbose)
        cerr << "Distance to nearest previously selected is column " << dist_to_prev_column << " of the name\n";

      dist_to_prev_column--;
    }
    else
    {
      precomputed_distance_to_previously_selected_tag = r;
      if (verbose)
        cerr << "Precomputed distances in '" << precomputed_distance_to_previously_selected_tag << "' tag of fingerprint\n";

      if (! precomputed_distance_to_previously_selected_tag.ends_with('<'))
        precomputed_distance_to_previously_selected_tag << '<';
    }
  }

  IW_STL_Hash_Map_int id_to_ndx;

  for (int i = 0; i < pool_size; i++)
  {
    const IWString & id = pool[i].id();

    id_to_ndx[id] = i;
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, can only process one file at a time\n";
    usage(4);
  }

  if (! read_groups(cl[0], id_to_ndx))
  {
    cerr << "Cannot read group assignments from '" << cl[0] << "'\n";
    return 3;
  }

  if (verbose)
  {
    extending_resizable_array<int> gsize;
    for (int i = 0; i < ngroup; i++)
    {
      int items_in_group = group[i].items_in_group();

      gsize[items_in_group]++;
    }

    for (int i = 0; i < gsize.number_elements(); i++)
    {
      if (gsize[i])
        cerr << gsize[i] << " groups had " << i << " members\n";
    }
  }

  if (cl.option_present('A'))
  {
    if (cl.option_present('r'))
    {
      if (! cl.value('r', report_processing_previously_selected) || report_processing_previously_selected < 1)
      {
        cerr << "The report progress of previously selected (-r) option must be a whole +ve number\n";
        usage(4);
      }

      if (verbose)
        cerr << "Will report progress of the previously selected molecules every " << report_processing_previously_selected << " items\n";
    }

    int i = 0;
    const_IWSubstring a;
    while (cl.value('A', a, i++))
    {
      if (! process_previously_selected_items(a))
      {
        cerr << "Cannot process previously selected file '" << a << "'\n";
        return i + 1;
      }
    }
  }

  for (int i = 0; i < ngroup; i++)
  {
    group[i].establish_desirability_measure();
  }

  if (verbose)
    cerr << "Start spread selections\n";

  group_spread(cout);

  if (verbose)
  {
    cerr << "Selected " << groups_selected << " groups, containing " << items_selected << " items\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = group_spread (argc, argv);

  return rc;
}
