/*
  An implementation of nearest neighbour predictions.
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <utility>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"
#include "nn_specification.h"
#include "prediction_weight_function.h"

using std::cerr;
using std::endl;
using std::pair;

static int verbose = 0;

static Tversky tversky;

/*
  It can often be easier to match things up if we strip leading zeros
*/

static int remove_leading_zeros = 0;

/*
  If we are reading experimental values from a file, it may be a descriptor file
*/

static int experimental_data_file_is_a_descriptor_file = 0;

static IWString missing_value ('.');

static int write_molecules_with_no_predictions = 0;

/*
  We can decide that we will only consider neighbours within a given distance
*/

static similarity_type_t maximum_distance = 1.0;

/*
  When doing leave-one-out scoring it is useful to drop exact matches
*/

static similarity_type_t minimum_distance = -1.0;

/*
  We can also decide that in order to get a valid prediction, we must
  have at least this number of valid neighbours - probably within maximum_distance
*/

static int min_neighbours_needed = 1;

/*
  If we want to do a manual leave-one-out computation, we need to be able to exclude a molecule
  from its neighbour list
*/

static int do_not_compare_molecules_with_themselves = 0;

static int classification_model = 0;

/*
  When reading model(s) from file(s), we need to make sure all the model types are consistent
*/

static int model_type_specified = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Does GFP nearneighbour predictions and uses a neighbour definition to make predictions\n";
  cerr << " -p <file>        training set fingerprint file\n";
  cerr << " -N <type>        what type of estimate to do\n";
  cerr << "                  2N     use the 2 nearest neighbours\n";
  cerr << "                  0.2    use all neighbours within a radius of 0.2\n";
  cerr << "                  0.1-3N all neighbours within 0.1, but\n";
  cerr << "                    with a minumum of 3 neighbours\n";
  cerr << " -T <dist>        discard distances greater than cutoff\n";
  cerr << " -t <dist>        discard distances of <dist> or shorter - useful for\n";
  cerr << "                  doing leave-one-out estimates of a training set\n";
  cerr << " -m <number>      require at least <number> valid neighbours for a prediction\n";
  cerr << " -M <string>      missing value specifier (default " << missing_value << ")\n";
  cerr << " -E file=<fname>  file with training set experimental data\n";
  cerr << " -E col=<col>     experimental data is always column <col> in the training set names\n";
  cerr << " -E YYY=          experimental data is a token in the name starting with 'YYY='\n";
  cerr << " -C ...           classification model\n";
  cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -j               skip first line in experimental data file\n";
  cerr << " -V ...           Tversky options, enter '-T help' for details\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

/*
  We impute activity to each target based on:
    the activity of each neighbour
    the distance from that neighbour
    a weighting function
*/

typedef float experimental_value_t;

/*
  Each member of the training set is a fingerprint and associated activity or category
  The experimental result will be either a continuous number or a category
*/

typedef int category_t;
typedef int category_type_t;
typedef double probability_type_t;

template <typename E>
class Training_Set_Member : public IW_General_Fingerprint
{
  protected:
    E _experimental_result;

  public:
    Training_Set_Member();

    int set_experimental_value_from_id_token (int token);
    int set_experimental_value_from_token (const const_IWSubstring & zprefix);

    void set_observed_activity (E e) { _experimental_result = e;}
    E observed_activity() const { return _experimental_result;}
};

class Classification_Training_Set_Member : public Training_Set_Member<category_type_t>
{
  private:
  public:
    Classification_Training_Set_Member();

};

static int highest_category_number = 0;

class Activity_Training_Set_Member : public Training_Set_Member<experimental_value_t>
{
  private:
  public:
    Activity_Training_Set_Member();
};

template <typename E>
Training_Set_Member<E>::Training_Set_Member()
{
  return;
}

Activity_Training_Set_Member::Activity_Training_Set_Member()
{
  _experimental_result = 0.0;

  return;
}

Classification_Training_Set_Member::Classification_Training_Set_Member()
{
  _experimental_result = -1;
}

class NN_Computation;

/*
  P is what the pool is composed of, and E will be either experimental_value_t or category_type_t
*/

template <typename P, typename E, typename R>
class Training_Set_Pool
{
  protected:
    int _pool_size;

    P * _pool;

//  When identifying neighbours, we need the list of activities and their distance in an array
//  that can be sorted by distance
    
    int _neighbours_found;
    pair<E, similarity_type_t> * _nbr;

//  Each pool needs to keep track of its prediction(s) for the current molecule

    R * _prediction;

//  We can add a constant value to the predicted values. Classification problems may need this

    E _offset;

//  private functions

    int _associate_experimental_values_with_training_set_members_column (const const_IWSubstring &);
    int _associate_experimental_values_with_training_set_members_file (const const_IWSubstring &);
    int _associate_experimental_values_with_training_set_members_prefix (const const_IWSubstring &);
    int _read_experimental_data (iwstring_data_source & input, IW_STL_Hash_Map<IWString, E> &);

    int _make_predictions (const IWString &, NN_Computation * prediction, int number_models, IWString_and_File_Descriptor &);
    int _make_prediction (const NN_Computation & prediction, int ndx);

    int _write_no_predictions (const IWString &, int number_computations, IWString_and_File_Descriptor & output);
    int _write_predicted_value (const IWString & id, int number_models, IWString_and_File_Descriptor & output);

  public:
    Training_Set_Pool();
    ~Training_Set_Pool();

    int build_pool (iwstring_data_source &);
    int build_pool (const const_IWSubstring &);

    int allocate_prediction_arrays (int);

    int associate_experimental_values_with_training_set_members (const const_IWSubstring & e);
    int associate_experimental_values_with_training_set_members_column (int);
    int associate_experimental_values_with_training_set_members_file (const const_IWSubstring & fname);
    int associate_experimental_values_with_training_set_members_file (iwstring_data_source & input);

    int make_predictions (IW_General_Fingerprint & fp, NN_Computation * prediction, int number_models, IWString_and_File_Descriptor & output);
    int make_predictions (IW_TDT & tdt, NN_Computation * prediction, int number_models, IWString_and_File_Descriptor & output);
    int make_predictions (iwstring_data_source & input, NN_Computation * prediction, int number_models, IWString_and_File_Descriptor & output);
    int make_predictions (const char * fname, NN_Computation * prediction, int number_models, IWString_and_File_Descriptor & output);

    int report_prediction_statistics (int, std::ostream &) const;
    int report_prediction_statistics (int number_computations, NN_Computation * computation, std::ostream & output) const;
};

template <typename P, typename E, typename R>
Training_Set_Pool<P, E, R>::Training_Set_Pool()
{
  _pool_size = 0;
  _pool = nullptr;

  _nbr = nullptr;

  _prediction = nullptr;

  _offset = static_cast<E>(0);

  return;
}

template <typename P, typename E, typename R>
Training_Set_Pool<P, E, R>::~Training_Set_Pool()
{
  if (NULL != _pool)
    delete [] _pool;

  if (NULL != _nbr)
    delete [] _nbr;

  if (NULL != _prediction)
    delete [] _prediction;

  return;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::allocate_prediction_arrays (int number_computations)
{
  assert (NULL == _prediction);

  _prediction = new R[number_computations];

  return NULL != _prediction;
}

static IWString identifier_tag ("PCN<");

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::build_pool (iwstring_data_source & input)
{
  if (0 == _pool_size)     // must grep the file to find out how many
  {
    assert (NULL == _pool);

    IWString tmp = '^';
    tmp << identifier_tag;

    _pool_size = input.grep (tmp);
    if (0 == _pool_size)
    {
      cerr << "Yipes, cannot find any '" << tmp << "' in the input\n";
      return 0;
    }

    if (verbose)
      cerr << "Input contains " << _pool_size << " fingerprints\n";

    _pool = new P[_pool_size];

    if (NULL == _pool)
    {
      cerr << "Yipes, cannot allocate space for " << _pool_size << " fingerprints\n";
      return 0;
    }
  }

  int i = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    if (! _pool[i].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    i++;

    if (i >= _pool_size)
    {
      cerr << "Pool is full, max " << _pool_size << endl;
      break;
    }
  }

  if (verbose)
  {
    cerr << i << " fingerprint objects added to pool\n";
  }

  _pool_size = i;

  _nbr = new pair<E, similarity_type_t> [_pool_size];

  if (NULL == _nbr)
    cerr << "Training_Set_Pool::build_pool:cannot allocate array of size " << _pool_size << endl;

  return 1;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::build_pool (const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  return build_pool (input);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::associate_experimental_values_with_training_set_members (const const_IWSubstring & e)
{
  if (e.starts_with ("col="))
    return _associate_experimental_values_with_training_set_members_column (e);
  if (e.starts_with ("file="))
    return _associate_experimental_values_with_training_set_members_file (e);

  return _associate_experimental_values_with_training_set_members_prefix (e);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_associate_experimental_values_with_training_set_members_column (const const_IWSubstring & e)
{
  const_IWSubstring ecopy (e);
  ecopy.remove_up_to_first ('=');

  int col;
  if (! ecopy.numeric_value (col) || col < 1)
  {
    cerr << "Training_Set_Pool::_associate_experimental_values_with_training_set_members_column: invalid column specifier '" << e << "'\n";
    return 0;
  }

  col--;

  return associate_experimental_values_with_training_set_members_column (col);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::associate_experimental_values_with_training_set_members_column (int col)
{
  int rc = 1;
  for (int i = 0; i < _pool_size; i++)
  {
    if (! _pool[i].set_experimental_value_from_id_token (col))
      rc = 0;
  }

  return rc;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_associate_experimental_values_with_training_set_members_prefix (const const_IWSubstring & e)
{
  int rc = 1;
  for (int i = 0; i < _pool_size; i++)
  {
    if (! _pool[i].set_experimental_value_from_token (e))
      rc = 0;
  }

  return rc;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_associate_experimental_values_with_training_set_members_file (const const_IWSubstring & e)
{
  const_IWSubstring ecopy (e);
  ecopy.remove_up_to_first ('=');

  return associate_experimental_values_with_training_set_members_file (ecopy);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::associate_experimental_values_with_training_set_members_file (const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Training_Set_Pool::associate_experimental_values_with_training_set_members_file: cannot open '" << fname << "'\n";
    return 0;
  }

  return associate_experimental_values_with_training_set_members_file (input);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::associate_experimental_values_with_training_set_members_file (iwstring_data_source & input)
{
  IW_STL_Hash_Map<IWString, E> zhash;
  if (! _read_experimental_data (input, zhash))
  {
    cerr << "Training_Set_Pool::associate_experimental_values_with_training_set_members_file: cannot read experimental data\n";
    return 0;
  }

  if (verbose)
    cerr << "Read " << zhash.size() << " experimental data values\n";

  int rc = 1;
  for (int i = 0; i < _pool_size; i++)
  {
    const IWString & id = _pool[i].id();

    if (zhash.contains (id))
    {
      _pool[i].set_observed_activity (zhash[id]);
      continue;
    }

    if (remove_leading_zeros || id.contains (' '))
    {
      IWString tmp (id);
      if (remove_leading_zeros)
        tmp.remove_leading_chars ('0');
      tmp.truncate_at_first (' ');

      if (zhash.contains (tmp))
      {
        _pool[i].set_observed_activity (zhash[tmp]);
        continue;
      }
    }

    cerr << "Training_Set_Pool::associate_experimental_values_with_training_set_members_file: no data for '" << id << "'\n";
    rc = 0;
  }

  if (0 == rc)
  {
    cerr << "Failed to match ids to experimental data\n";
    for (typename IW_STL_Hash_Map<IWString, E>::const_iterator i = zhash.begin(); i != zhash.end(); ++i)
    {
      cerr << (*i).first << " value " << (*i).second << endl;
    }
  }

  return rc;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_read_experimental_data (iwstring_data_source & input,
                                                  IW_STL_Hash_Map<IWString, E> & zhash)
{
  const_IWSubstring buffer;

  if (experimental_data_file_is_a_descriptor_file)
    input.next_record (buffer);

  while (input.next_record (buffer))
  {
    const_IWSubstring id, activity;
    if (! buffer.split (id, ' ', activity))
    {
      cerr << "Training_Set_Pool::_read_experimental_data: cannot parse '" << buffer << "'\n";
      return 0;
    }

    if (remove_leading_zeros)
      id.remove_leading_chars ('0');

    activity.strip_trailing_blanks();
    activity.strip_leading_blanks();

    E obs;
    if (! activity.numeric_value (obs))
    {
      if (1 == input.lines_read())
        continue;

      cerr << "Training_Set_Pool::_read_experimental_data: invalid observed value '" << buffer << "'\n";
      cerr << " token is '" << activity << "'\n";
      return 0;
    }

    zhash[id] = obs;
  }

  return zhash.size();
}

/*
  If we are going to do different types of computation, we need an object to hold
  the different types
*/

class NN_Computation
{
  private:
    int _number_neighbours;

    double _radius;

    IWString _type;

    NN_Weight_Function _weight_function;

//  For classification models, we may have a special class and a threshold probability

    category_type_t _special_class;
    probability_type_t _special_class_probability;

//  private functions

    int _is_neighbours (const const_IWSubstring &);
    int _is_distance   (const const_IWSubstring &);
    int _build_weight_function (const const_IWSubstring &);
    int _build_neighbour_specification (const const_IWSubstring & c);
    int _build         (const const_IWSubstring &);
    int _parse_special_class_specification (const const_IWSubstring &);

  public:
    NN_Computation();

    int build (const const_IWSubstring &);

    int  number_neighbours() const { return _number_neighbours;}

    double radius() const { return _radius;}

    category_type_t special_category() const { return _special_class;}
    probability_type_t special_class_probability() const { return _special_class_probability;}

    const NN_Weight_Function & weight_function() const { return _weight_function;}

    const IWString & estimate_type() const { return _type;}
};

NN_Computation::NN_Computation()
{
  _number_neighbours = -1;
  _radius = -1.0;

  _special_class = -1;
  _special_class_probability = 0.0;

  return;
}

static NN_Computation * computation = nullptr;
static int number_computations = 0;

static IWString no_predictions;

/*
  With predictions, we output confidence values. Our initial implementation is
  to produce the sum of the reciprocal distances used in the estimate. In order
  to avoid large numbers, we have a cutoff distance
*/

typedef float confidence_t;

static confidence_t min_distance_in_confidence_computations = 0.1;

/*
  A prediction consists of the actual prediction (which may be unset) and
  the distance of the nearest neighbour in the training set. We also keep
  track of the statistics on all molecules predicted

  We have two fundamentally different types of prediction, continuous response
  and classification
*/

template <typename E>
class Prediction_Base : public Set_or_Unset<E>
{
  protected:

    int _number_missing_values;

//  This distance is for the current molecule

    similarity_type_t _distance_to_nearest_training_set_member;

    int _number_neighbours_used_for_estimate;

//  We record various confidence values. This one is the sum of the
//  reciprocal distances used for the activity estimation

    confidence_t _confidence;

//  private functions

    void _extra_data (confidence_t confidence, int neighbours_used, similarity_type_t mindist);

  public:
    Prediction_Base();

    int report (std::ostream &) const;

    int write_predicted_value (IWString_and_File_Descriptor &, E = 0) const;

    int extra_missing_value();

    similarity_type_t distance_to_nearest_training_set_member() const { return _distance_to_nearest_training_set_member;}
    void set_distance_to_nearest_training_set_member ( similarity_type_t s) {_distance_to_nearest_training_set_member = s;}
};

template <typename E>
Prediction_Base<E>::Prediction_Base()
{
  _number_missing_values = 0;
}

template <typename E>
int
Prediction_Base<E>::extra_missing_value()
{
  Set_or_Unset<E>::unset();

  _number_missing_values++;

  return 1;
}

template <typename E>
int
Prediction_Base<E>::write_predicted_value (IWString_and_File_Descriptor & output, E offset) const
{
  output << ' ';

  E p;
  if (Set_or_Unset<E>::value (p))
  {
    output << (p + offset) << ' ' << _confidence << ' ' << (_confidence / _number_neighbours_used_for_estimate);
  }
  else
    output << missing_value << ' ' << missing_value << ' ' << missing_value;

  return output.good();
}

template <typename E>
void
Prediction_Base<E>::_extra_data (confidence_t confidence, int neighbours_used, similarity_type_t mindist)
{
  _confidence = confidence;

  _number_neighbours_used_for_estimate = neighbours_used;

  _distance_to_nearest_training_set_member = mindist;

  return;
}

class Prediction_Continuous_Response : public Prediction_Base<experimental_value_t>
{
  private:

//  We gather statistics on the predictions for every molecule

    Accumulator<double> _stats;

  public:

    int report (std::ostream &) const;

    int compute_predicted_values (int neighbours_found, pair<experimental_value_t, similarity_type_t> * adp, int, similarity_type_t, const NN_Computation &);
    int compute_predicted_values (int neighbours_to_use, pair<category_type_t, similarity_type_t> * cdp, const NN_Computation &);
};

int
Prediction_Continuous_Response::report (std::ostream & os) const
{
  if (0 == _stats.n())
  {
    os << "no predictions\n";
    return os.good();
  }

  os << "values between " << _stats.minval() << " and " << _stats.maxval();
  if (_stats.n() > 1)
    os << " average " << _stats.average() << " variance " << _stats.variance();
  if (_number_missing_values)
    os << '\n' << _number_missing_values << " missing values";

  os << '\n';

  return os.good();
}

class Prediction_Categorical_Response : public Prediction_Base<category_t>
{
  private:
    int * _put_in_class;

    weight_t * _category_sumw;

    confidence_t * _category_confidence;

    float * _probability_of_category;

    int _compute_predicted_category (weight_t & w, category_type_t &, const NN_Computation &);
    int _compute_predicted_category (weight_t & w, category_type_t &);

  public:
    Prediction_Categorical_Response();
    ~Prediction_Categorical_Response();

    int initialise_category_arrays();

    int report (std::ostream &) const;

    int compute_predicted_values (int neighbours_found, pair<category_type_t, similarity_type_t> * cdp, int, similarity_type_t, const NN_Computation &);
    int compute_predicted_values (int neighbours_to_use, pair<category_type_t, similarity_type_t> * cdp, const NN_Computation &);

    int write_predicted_value (IWString_and_File_Descriptor &, int = 0) const;
};

Prediction_Categorical_Response::Prediction_Categorical_Response()
{
  cerr << "In constructoer,  " << highest_category_number << endl;
  if (highest_category_number > 0)
    initialise_category_arrays();
  else
  {
    _put_in_class = nullptr;
    _category_sumw = nullptr;
    _category_confidence = nullptr;
    _probability_of_category = nullptr;
  }

  return;
}

int
Prediction_Categorical_Response::initialise_category_arrays()
{
  assert (highest_category_number > 0);

  _put_in_class = new_int (highest_category_number + 1);
  _category_sumw = new weight_t[highest_category_number + 1];
  _category_confidence = new confidence_t[highest_category_number + 1];
  _probability_of_category = new float[highest_category_number + 1];

  return 1;
}

Prediction_Categorical_Response::~Prediction_Categorical_Response()
{
  if (NULL != _put_in_class)
    delete _put_in_class;
  if (NULL != _category_sumw)
    delete _category_sumw;
  if (NULL != _category_confidence)
    delete _category_confidence;
  if (NULL != _probability_of_category)
    delete _probability_of_category;

  return;
}

int
Prediction_Categorical_Response::report (std::ostream & os) const
{
  if (NULL == _put_in_class)
  {
    cerr << "Prediction_Categorical_Response::report:no values\n";
    return 0;
  }

  int n = 0;
  for (int i = 0; i <= highest_category_number; i++)
  {
    n += _put_in_class[i];
  }

  if (0 == n)
  {
    cerr << "No Predictions\n";
    return os.good();
  }

  os << n << " predictions\n";
  for (int i = 0; i <= highest_category_number; i++)
  {
    if (_put_in_class[i])
      os << " placed in class " << i << " " << _put_in_class[i] << " times\n";
  }

  if (_number_missing_values)
    os << _number_missing_values << " missing values\n";

  return os.good();
}

int
Prediction_Categorical_Response::write_predicted_value (IWString_and_File_Descriptor & os,
                        int offset) const
{
  if (! Prediction_Base<category_t>::write_predicted_value (os, offset))
    return 0;

  if (Set_or_Unset<category_t>::is_set())
  {
    for (int i = 0; i <= highest_category_number; i++)
    {
      os << ' ' << _probability_of_category[i];
    }
  }
  else
  {
    for (int i = 0; i <= highest_category_number; i++)
    {
      os << ' ' << missing_value;
    }
  }

  return os.good();
}

class Classification_Training_Set_Pool : public Training_Set_Pool<Classification_Training_Set_Member, category_type_t, Prediction_Categorical_Response>
{
  private:

  public:

    int initialise_category_data();
};

class Activity_Training_Set_Pool : public Training_Set_Pool<Activity_Training_Set_Member, experimental_value_t, Prediction_Continuous_Response>
{
  private:
  public:
};

int
Classification_Training_Set_Pool::initialise_category_data()
{
  assert (0 == highest_category_number);

  int lowest_category_number = _pool[0].observed_activity();

  highest_category_number = _pool[0].observed_activity();

  for (int i = 1; i < _pool_size; i++)
  {
    int c = _pool[i].observed_activity();

    if (c < lowest_category_number)
      lowest_category_number = c;
    else if (c > highest_category_number)
      highest_category_number = c;
  }

  if (lowest_category_number < 0)
  {
    _offset = lowest_category_number;

    highest_category_number -= lowest_category_number;

    for (int i = 0; i < _pool_size; i++)
    {
      int c = _pool[i].observed_activity();
      _pool[i].set_observed_activity(c - lowest_category_number);
    }
  }

  return 1;
}

/*
  If we don't satisfy the number of neighbours condition, we may get
  all missing values
*/

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_write_no_predictions (const IWString & id,
                                                   int number_computations,
                                                   IWString_and_File_Descriptor & output)
{
  append_first_token_of_name (id, output);

  output << no_predictions;

  if (verbose)
  {
    for (int i = 0; i < number_computations; i++)
    {
      _prediction[i].extra_missing_value();
    }
  }

  return output.good();
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_write_predicted_value (const IWString & id,
                                           int number_models,
                                           IWString_and_File_Descriptor & output)
{
  append_first_token_of_name (id, output);

  for (int i = 0; i < number_computations; i++)
  {
    _prediction[i].write_predicted_value (output, _offset);
  }

  output << '\n';

  return output.good();
}

/*
  When doing continuous response things we need to keep track of all
  the activity and distance of all the neighbours
*/

template <typename E>
class Observed_Distance_Pair : public pair <E, similarity_type_t>
{
  private:

#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using pair<E, similarity_type_t>::first;
    using pair<E, similarity_type_t>::second;
#endif

  public:

    E observed_activity() const { return first;}
    similarity_type_t    distance() const { return second;}

    void set_observed_activity (E e) { first = e;}
    void set_distance (similarity_type_t e) { second = e;}
};

class Activity_Distance_Pair : public Observed_Distance_Pair<experimental_value_t>
{
  private:
  public:
};

class Category_Distance_Pair : public Observed_Distance_Pair<category_type_t>
{
  private:
  public:
};


int
Prediction_Categorical_Response::_compute_predicted_category (weight_t & sumw,
                                     category_type_t & pred)
{
  pred = 0;
  sumw = _category_sumw[0];

  weight_t highest_weight = _category_sumw[0];
  for (int i = 1; i <= highest_category_number; i++)
  {
    sumw += _category_sumw[i];

    if (_category_sumw[i] > highest_weight)
    {
      highest_weight = _category_sumw[i];
      pred = i;
    }
  }

  return 1;
}

/*
  We are doing a model with error balancing
*/

int
Prediction_Categorical_Response::_compute_predicted_category (weight_t & sumw,
                                     category_type_t & pred,
                                     const NN_Computation & computation)
{
  category_type_t special_class = computation.special_category();
  assert (special_class >= 0 && special_class <= highest_category_number);

  sumw = sum_vector (_category_sumw, highest_category_number + 1);

  probability_type_t p = _category_sumw[special_class] / sumw;

  if (p >= computation.special_class_probability())
  {
    pred = special_class;
    return 1;
  }

// We didn't exceed the threshold for our special class, just do it the regular way

  return _compute_predicted_category (sumw, pred);    // sumw will be recomputed, oh well..
}

int
Prediction_Categorical_Response::compute_predicted_values (int number_neighbours,
                          pair<category_type_t, similarity_type_t> * cdp,
                          int min_neighbours_needed,
                          similarity_type_t max_dist_this_estimate,
                          const NN_Computation & computation)
{
  assert (NULL != cdp);

  const NN_Weight_Function & weight_function = computation.weight_function();

  if (NULL == _category_sumw)
    initialise_category_arrays();

  int neighbours_used_in_estimate = 0;

  set_vector (_category_sumw, highest_category_number + 1, static_cast<weight_t> (0.0));
  set_vector (_category_confidence, highest_category_number + 1, static_cast<confidence_t> (0.0));

  for (int i = 0; i < number_neighbours; i++)
  {
    const pair<category_type_t, similarity_type_t> & cdpi = cdp[i];

    similarity_type_t d = cdpi.second;

    if (d > max_dist_this_estimate)   // max_dist is the radius for this particular estimation
    {
      if (neighbours_used_in_estimate >= min_neighbours_needed)
        break;
    }

    neighbours_used_in_estimate++;

    confidence_t conf;

    if (d <= min_distance_in_confidence_computations)
      conf = 1.0 / min_distance_in_confidence_computations;
    else
      conf = 1.0 / d;

    weight_t w = weight_function.weight (d);

    category_type_t c = cdpi.first;

    _category_sumw[c] += w;
    _category_confidence[c] += conf;
  }

  if (neighbours_used_in_estimate < min_neighbours_needed)   // min_neighbours_needed must be initialised to a non-zero value
    cerr << "ONly found " << neighbours_used_in_estimate << " neighbours\n";
  if (neighbours_used_in_estimate < min_neighbours_needed)   // min_neighbours_needed must be initialised to a non-zero value
    return extra_missing_value();

  weight_t sumw;
  category_t pred = 0;

  if (computation.special_category() >= 0)
    _compute_predicted_category (sumw, pred, computation);
  else
    _compute_predicted_category (sumw, pred);

  Set_or_Unset<category_t>::set (pred);

  Prediction_Base<category_t>::_extra_data (_category_confidence[pred], neighbours_used_in_estimate, cdp[0].second);

  _put_in_class[pred]++;

  for (int i = 0; i <= highest_category_number; i++)
  {
    _probability_of_category[i] = _category_sumw[i] / sumw;
  }

  return 1;
}

int
Prediction_Continuous_Response::compute_predicted_values (int number_neighbours,
                                    pair<experimental_value_t, similarity_type_t> * adp,
                                    int min_neighbours_needed,
                                    similarity_type_t max_dist_this_estimate,
                                    const NN_Computation & computation)
{
  assert (NULL != adp);

  const NN_Weight_Function & weight_function = computation.weight_function();

  experimental_value_t pred = 0.0;

  confidence_t confidence = 0.0;

  int neighbours_used_in_estimate = 0;

  weight_t sumw = 0.0;

  for (int i = 0; i < number_neighbours; i++)
  {
    const pair<experimental_value_t, similarity_type_t> & adpi = adp[i];

    similarity_type_t d = adpi.second;

    if (d > max_dist_this_estimate)   // max_dist is the radius for this particular estimation
    {
      if (neighbours_used_in_estimate >= min_neighbours_needed)
        break;
    }

    neighbours_used_in_estimate++;

    if (d <= min_distance_in_confidence_computations)
      confidence += 1.0 / min_distance_in_confidence_computations;
    else
      confidence += 1.0 / d;

    weight_t w = weight_function.weight (d);

//  cerr << "Distance " << d << " weight " << w << endl;
    assert (w >= 0.0 && w <= 1.0);

    experimental_value_t e = adpi.first;

    pred += e * w;
    sumw += w;
  }

  if (neighbours_used_in_estimate < min_neighbours_needed)   // min_neighbours_needed must be initialised to a non-zero value
    return extra_missing_value();

  assert (0.0 != sumw);

  pred = pred / sumw;

  Set_or_Unset<experimental_value_t>::set (pred);

  _stats.extra (static_cast<double> (pred));

  Prediction_Base<experimental_value_t>::_extra_data (confidence, neighbours_used_in_estimate, adp[0].second);

  return 1;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_make_prediction (const NN_Computation & prediction,
                                              int ndx)
{
  int min_neighbours_needed = prediction.number_neighbours();

  if (min_neighbours_needed > _neighbours_found)
    min_neighbours_needed = _neighbours_found;

  similarity_type_t radius = prediction.radius();

//cerr << "Needs " << min_neighbours_needed << " and radius " << radius << endl;

  if (radius <= static_cast<similarity_type_t> (0.0))    // just a number of neighbours
    return _prediction[ndx].compute_predicted_values (min_neighbours_needed, _nbr, min_neighbours_needed, 1.0, prediction);

// This must be the case of a distance and a minimum number of neighbours

  return _prediction[ndx].compute_predicted_values (_neighbours_found, _nbr, min_neighbours_needed, radius, prediction);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::_make_predictions (const IWString & id,
                                           NN_Computation * prediction,
                                           int number_models,
                                           IWString_and_File_Descriptor & output)
{
  assert (_neighbours_found > 0);

  for (int i = 0; i < number_models; i++)
  {
    _make_prediction (prediction[i], i);
  }

  return _write_predicted_value (id, number_models, output);
}

/*
  Somewhat dangerous because we use the same comparitor for both category_type_t and experimental_value_t.
  Make sure those things remain the same size
*/

static int
OD_comparitor (const void * p1, const void * p2)
{
  typedef pair<experimental_value_t, similarity_type_t> OD;

  const OD * odp1 = (const OD *) p1;
  const OD * odp2 = (const OD *) p2;

  similarity_type_t d1 = odp1->second;
  similarity_type_t d2 = odp2->second;

//cerr << "Comparing " << d1 << " with " << d2 << endl;

  if (d1 < d2)
    return -1;
  if (d1 > d2)
    return 1;

  return 0;
}

#define CHECK_SORT
#ifdef CHECK_SORT

template <typename T>
void
check_ordering (const T * xdp, int n)
{
  int out_of_order = 0;

  similarity_type_t previous_distance = xdp[0].second;
  for (int i = 1; i < n; i++)
  {
    similarity_type_t d = xdp[i].second;

    if (d < previous_distance)
    {
      cerr << "Distances out of order, i = " << i << " prev " << previous_distance << " d = " << d << endl;
      out_of_order = 1;
    }
  }

  if (out_of_order)
  {
    cerr << "Distances out of order, found " << n << " neighbours\n";
    for (int i = 0; i < n; i++)
    {
      cerr << " i = " << i << " d = " << xdp[i].second << endl;
    }

    abort();
  }

  return;
}

#endif

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::make_predictions (IW_General_Fingerprint & fp,
                NN_Computation * prediction,
                int number_models,
                IWString_and_File_Descriptor & output)
{
  _neighbours_found = 0;

  for (int i = 0; i < _pool_size; i++)
  {
    P & pi = _pool[i];

    if (! can_be_compared (pi, fp))
      continue;

    similarity_type_t d;
    if (tversky.active())
      d = static_cast<similarity_type_t> (1.0) - pi.tversky (fp, tversky);
    else
      d = fp.distance (pi);

//  cerr << "pool item " << i << " '" << pi.id() << "' at distance " << d << " from " << fp.id() << endl;

    if (d <= minimum_distance)
      continue;
    if (d > maximum_distance)    // these two molecules too far apart to influence each other
      continue;

    if (do_not_compare_molecules_with_themselves && static_cast<similarity_type_t> (0.0) == d && pi.id() == fp.id())
      continue;

    _nbr[_neighbours_found].first = pi.observed_activity();
    _nbr[_neighbours_found].second = d;

    _neighbours_found++;
  }

//cerr << "Found " << _neighbours_found << " neighbours, need " << min_neighbours_needed << endl;

  if (_neighbours_found < min_neighbours_needed)
  {
    if (write_molecules_with_no_predictions)
      return _write_no_predictions (fp.id(), number_models, output);

    return 1;
  }

  if (_neighbours_found > 1)
  {
    assert (sizeof (similarity_type_t) == sizeof (category_type_t));

    qsort (_nbr, _neighbours_found, sizeof (pair<E, similarity_type_t>), OD_comparitor);
  }

#ifdef CHECK_SORT
  check_ordering (_nbr, _neighbours_found);
#endif

  return _make_predictions (fp.id(), prediction, number_models, output);
}

static int tdts_read = 0;

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::make_predictions (IW_TDT & tdt,
                NN_Computation * prediction,
                int number_models,
                IWString_and_File_Descriptor & output)
{
  IW_General_Fingerprint fp;
  int fatal;
  if (! fp.construct_from_tdt (tdt, fatal))
  {
    if (fatal)
    {
      cerr << "Training_Set_Pool::make_predictions: cannot parse TDT\n";
      cerr << tdt;
      return 0;
    }
  }

  return make_predictions (fp, prediction, number_models, output);
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::make_predictions (iwstring_data_source & input,
                NN_Computation * prediction,
                int number_models,
                IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    if (! make_predictions (tdt, prediction, number_models, output))
    {
      cerr << "Training_Set_Pool::make_predictions: prediction failure near line " << input.lines_read() << endl;
      return 0;
    }
  }

  return output.good();
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::make_predictions (const char * fname,
                NN_Computation * prediction,
                int number_models,
                IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Training_Set_Pool::make_predictions: cannot open '" << fname << "'\n";
    return 0;
  }

  return make_predictions (input, prediction, number_models, output);
}

/*
  Common code for extracting the experimental value from a token in a string
*/

template <typename B, typename A>
int
get_experimental_value (const B & buffer,
                        int data_column,
                        A & expt)
{
  const_IWSubstring sexpt;
  if (! buffer.word (data_column, sexpt))
  {
    cerr << "Cannot extract experimental value, column " << data_column << " in '" << buffer << "'\n";
    return 0;
  }

  if (! sexpt.numeric_value (expt))
  {
    cerr << "Invalid numeric value for experimental activity '" << sexpt << "'\n";
    return 0;
  }

  return 1;
}

template <typename E>
int
Training_Set_Member<E>::set_experimental_value_from_id_token (int c)
{
  if (! get_experimental_value (_id, c, _experimental_result))
  {
    cerr << "Training_Set_Member::set_experimental_value_from_id_token: cannot set activity from id token\n";
    return 0;
  }

  return 1;
}

template <typename E>
int
Training_Set_Member<E>::set_experimental_value_from_token (const const_IWSubstring & zprefix)
{
  assert (zprefix.length());

  int i = 0;
  const_IWSubstring token;

  while (_id.nextword (token, i))
  {
    if (! token.starts_with (zprefix))
      continue;

    token.remove_leading_chars (zprefix.length());

    if (! token.numeric_value (_experimental_result))
    {
      cerr << "Invalid '" << zprefix << "' activity '" << _id << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "Training_Set_Member::set_experimental_value_from_token: no identifier tokens start with '" << zprefix << "'\n";
  return 0;
}

template int get_experimental_value (const const_IWSubstring &, int, experimental_value_t &);
template int get_experimental_value (const IWString &, int, experimental_value_t &);

int
NN_Computation::_is_distance (const const_IWSubstring & d)
{
  if (d.numeric_value (_radius) && _radius >= 0.0 && _radius <= 1.0)
    return 1;

  return 0;
}

int
NN_Computation::_is_neighbours (const const_IWSubstring & n)
{
  if (! n.ends_with ('N'))
    return 0;

  const_IWSubstring myn = n;
  myn.chop (1);

  if (myn.numeric_value (_number_neighbours) && _number_neighbours > 0)
    return 1;

  return 0;
}

int
NN_Computation::build (const const_IWSubstring & c)
{
  if (! _build (c))
    return 0;

  _type = c;

  return 1;
}

/*
  The threshold for the special class will be of the form
    S=0.nn
*/

int
NN_Computation::_parse_special_class_specification (const const_IWSubstring & token)
{
  assert (token.starts_with ('S'));

  const_IWSubstring mytoken (token);
  mytoken.remove_leading_chars (1);    // get rid of the 'S'

  const_IWSubstring sclass, sthreshold;
  if (! mytoken.split (sclass, '=', sthreshold))
  {
    cerr << "NN_Computation::_parse_special_class_specification: invalid special class '" << token << "'\n";
    return 0;
  }

  category_type_t c;     // store in temporary variable first

  if (! sclass.numeric_value (c) || c < 0)
  {
    cerr << "NN_Computation::_parse_special_class_specification: invalid special class value '" << token << "'\n";
    return 0;
  }

  if (! sthreshold.numeric_value (_special_class_probability) || _special_class_probability <= 0.0 || _special_class_probability > 1.0)
  {
    cerr << "NN_Computation::_parse_special_class_specification: invalid probability '" << token << "'\n";
    return 0;
  }

  _special_class = c;

  return 1;
}

int
NN_Computation::_build_weight_function (const const_IWSubstring & w)
{
  assert (w.starts_with ("W:"));
  const_IWSubstring myw = w;
  myw.remove_leading_chars (2);

  return _weight_function.construct_from_command_line_token (myw, ' ');
}

/*
  A computation specification looks like

  0.22
  0.18-3N
  0.22,W=linear
  0.18-3N,W=step0.5
  3N,W=def,S1=0.4
*/

int
NN_Computation::_build (const const_IWSubstring & c)
{
  const_IWSubstring weight_function_specification;
  if (! parse_model_specifications (c,
                            weight_function_specification,
                            _special_class,
                            _special_class_probability,
                            _radius,
                            _number_neighbours))

  {
    cerr << "NN_Computation:_build: cannot parse specification '" << c << "'\n";
    return 0;
  }

  if (0 == weight_function_specification.length())
    _weight_function.initialise_default_weight_function();
  else if (! _weight_function.construct_from_command_line_token (weight_function_specification, 'x'))
  {
    cerr << "Cannot parse weight function specification '" << c << "' to '" << weight_function_specification << "'\n";
    return 0;
  }

  if (verbose)
  {
    cerr << "Model with";
    if (_number_neighbours > 0)
      cerr << ' ' << _number_neighbours << " neighbours";
    if (_radius > 0.0)
      cerr << ' ' << _radius << " radius";
    if (_special_class >= 0)
      cerr << " errbal = " << _special_class_probability << " class " << _special_class;

    cerr << " W:" << _weight_function.weight_function_type();

    cerr << endl;
  }

  return 1;
}

int
NN_Computation::_build_neighbour_specification (const const_IWSubstring & c)
{
  const_IWSubstring c1, c2;

  if (c.split (c1, '-', c2))
  {
    if (0 == c1.length() || 0 == c2.length())
    {
      cerr << "NN_Computation::_build_neighbour_specification: invalid specification '" << c << "'\n";
      return 0;
    }

    if (_is_distance (c1) && _is_neighbours (c2))
      return 1;

    cerr << "NN_Computation::_build_neighbour_specification: cannot process both '" << c1 << "' and '" << c2 << "'\n";
    return 0;
  }

  if (_is_distance (c))
    ;
  else if (_is_neighbours (c))
    ;
  else
  {
    cerr << "NN_Computation::_build_neighbour_specification: cannot process '" << c << "'\n";
    return 0;
  }

  return 1;
}

template <typename P, typename E, typename R>
int
Training_Set_Pool<P, E, R>::report_prediction_statistics (int number_computations,
                                          NN_Computation * computation,
                                          std::ostream & output) const
{
  for (int i = 0; i < number_computations; i++)
  {
    const NN_Computation & ci = computation[i];

    output << "Prediction '" << ci.estimate_type() << "' ";

    const R & pi = _prediction[i];

    pi.report (output);
  }

  return output.good();
}

template <typename N>
int
nn_predictions (Command_Line & cl,
                NN_Computation * computation,
                int number_computations,
                N & pool,
                IWString_and_File_Descriptor & output)
{

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! pool.make_predictions (cl[i], computation, number_computations, output))
    {
      cerr << "Could not make predictions for something in '" << cl[i] << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  if (verbose)
    pool.report_prediction_statistics (number_computations, computation, cerr);

  return output.good();
}

static int
write_header (IWString_and_File_Descriptor & os, int number_computations)
{
  os << "name";
  for (int i = 0; i < number_computations; i++)
  {
    os << " pred" << i << " conf" << i << " nconf" << i;
  }

  return os.good();
}

template <typename N>
int
build_pool_and_initialise_expt_data (const IWString & trgfp,
                                     const IWString & dash_e,
                                     N & pool,
                                     int number_computations)
{
  assert (trgfp.length());

  if (! pool.build_pool (trgfp))
  {
    cerr << "Cannot build pool from '" << trgfp << "'\n";
    return 0;
  }

  pool.allocate_prediction_arrays (number_computations);

  if (! pool.associate_experimental_values_with_training_set_members (dash_e))
  {
    cerr << "Cannot associate experimental values with taining set members via '" << dash_e << "' specification\n";
    return 0;
  }

  return 1;
}

static int
nn_predictions_classification (Command_Line & cl,
                               const IWString & trgfp,
                               const IWString & dash_e,
                               NN_Computation * computation,
                               int number_computations, 
                               IWString_and_File_Descriptor & output)
{
  Classification_Training_Set_Pool pool;

  if (! build_pool_and_initialise_expt_data (trgfp, dash_e, pool, number_computations))
  {
    cerr << "Cannot initialise pool\n";
    return 0;
  }

  pool.initialise_category_data();

  if (verbose)
    cerr << "Highest category is " << highest_category_number << endl;

  write_header (output, number_computations);

  for (int i = 0; i <= highest_category_number; i++)
  {
    output << " prob" << i;
  }

  output << '\n';

  return nn_predictions (cl, computation, number_computations, pool, output);
}

static int
examine_and_check_model_type (const const_IWSubstring & mtype)
{
  if (KNN_CONTINUOUS == mtype)
  {
    if (0 == model_type_specified)
      model_type_specified = 1;
    else if (1 != model_type_specified)
    {
      cerr << "Inconsistent model types encountered\n";
      return 0;
    }
  }
  else if (KNN_CLASSIFICATION == mtype)
  {
    if (0 == model_type_specified)
      model_type_specified = 2;
    else if (2 != model_type_specified)
    {
      cerr << "Inconsistent model types encountered\n";
      return 0;
    }
  }
  else
  {
    cerr << "Unrecognised model type '" << mtype << "'\n";
    return 0;
  }

  return 1;
}

static int
count_models_in_file (iwstring_data_source & input,
                      IWString & trgfp,
                      IWString & dash_e,
                      resizable_array_p<IWString> & string_models)
{
  int rc = 0;

  const_IWSubstring buffer;

  if (! input.next_record (buffer))
  {
    cerr << "Empty model specification file\n";
    return 0;
  }

  if (buffer != KNN_MODEL_FILE_VERSION_1)
  {
    cerr << "Invalid KNN model file '" << buffer << "'\n";
    return 0;
  }

  while (input.next_record (buffer))
  {
    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    const_IWSubstring directive, qualifiers;

    if (! buffer.split (directive, ' ', qualifiers))
      continue;

    if (KNN_MODEL_TYPE == directive)
    {
      if (! examine_and_check_model_type (qualifiers))
        return 0;
    }
    else if (TR_GFP_DIRECTIVE == directive)
    {
      if (0 == trgfp.length())
        trgfp = qualifiers;
      else if (trgfp != qualifiers)
      {
        cerr << "Training set gfp mismatch '" << trgfp << "' vs '" << qualifiers << "'\n";
        return 0;
      }
    }
    else if (KNN_MODEL_PREFIX == directive)
    {
      IWString * tmp = new IWString (qualifiers);
      string_models.add (tmp);

      rc++;
    }
    else if (EXPT_SOURCE == directive)
    {
      dash_e = qualifiers;
    }
    else if ("gfpm:" == directive)
      ;
    else
    {
      cerr << "Unrecognised directive '" << buffer << "'\n";
//    return 0;
    }
  }

  return rc;
}

static int
count_models_in_file (const const_IWSubstring & fname,
                      IWString & trgfp,
                      IWString & dash_e,
                      resizable_array_p<IWString> & string_models)
{
  iwstring_data_source input (fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return count_models_in_file (input, trgfp, dash_e, string_models);
}

/*
  Models can be specified with the -N option, or 
  -N FILE:
  which is a file of models. To determine the number of models, we
  must read those files. We save the string representation of each
  model in the array
*/

static int
count_models_to_evaluate (const Command_Line & cl,
                          char cflag,
                          IWString & trgfp,
                          IWString & dash_e,
                          resizable_array_p<IWString> & string_models)
{
  int rc = 0;

  int i = 0;
  const_IWSubstring token;
  while (cl.value (cflag, token, i++))
  {
    if (token.starts_with ("FILE:"))
    {
      token.remove_leading_chars (5);
      int tmp = count_models_in_file (token, trgfp, dash_e, string_models);
      if (0 == tmp)
      {
        cerr << "Cannot determine number of models in file '" << token << "'\n";
        return 0;
      }

      rc += tmp;
    }
    else
    {
      IWString * tmp = new IWString (token);
      string_models.add (tmp);
      rc++;
    }
  }

  return rc;
}

static int
nn_predictions_continuous (Command_Line & cl,
                           const IWString & trgfp,
                           const IWString & dash_e,
                           NN_Computation * computation,
                           int number_computations, 
                           IWString_and_File_Descriptor & output)
{
  Activity_Training_Set_Pool pool;

  if (! build_pool_and_initialise_expt_data (trgfp, dash_e, pool, number_computations))
  {
    cerr << "Cannot initialise pool\n";
    return 0;
  }

  write_header (output, number_computations);

  output << '\n';

  return nn_predictions (cl, computation, number_computations, pool, output);
}

static int
nn_predictions (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "F:P:W:vjN:M:E:p:t:T:m:zV:C:wh");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised options encountered\n";
    usage (4);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('C'))
  {
    classification_model = 1;

    const_IWSubstring c = cl.string_value ('C');

    if (verbose)
      cerr << "Will treat as a classification model\n";
  }

  if (cl.option_present ('F') || cl.option_present ('P') || cl.option_present ('W'))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage (23);
    }
  }

  if (cl.option_present ('V'))
  {
    if (! tversky.parse_command_line (cl, 'V', verbose))
    {
      cerr << "Cannot initialise Tversky conditions\n";
      usage (7);
    }
  }

  if (cl.option_present ('M'))
  {
    missing_value = cl.string_value ('M');

    if (verbose)
      cerr << "Missing values written as '" << missing_value << "'\n";

    write_molecules_with_no_predictions = 1;
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', minimum_distance) || minimum_distance < 0.0 || minimum_distance >= 1.0)
    {
      cerr << "The minimum neighbour distance (-t) option must have a valid distance\n";
      usage (13);
    }

    if (verbose)
      cerr << "Neighbours closer than " << minimum_distance << " or shorter will be discarded\n";
  }

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', maximum_distance) || maximum_distance <= 0.0 || maximum_distance > 1.0)
    {
      cerr << "The maximum neighbour distance (-T) option must have a valid distance\n";
      usage (13);
    }

    if (verbose)
      cerr << "Only neighbours closer than " << maximum_distance << " will be considered for making predictions\n";

    if (cl.option_present ('m'))
    {
      if (! cl.value ('m', min_neighbours_needed) || min_neighbours_needed < 1)
      {
        cerr << "The minimum valid neighbours (-m) option must be a whole positive number\n";
        usage (18);
      }

      if (verbose)
        cerr << "Unless there are " << min_neighbours_needed << " valid neighbours, no prediction will be made\n";
    }
  }

  if (cl.option_present ('w'))
  {
    write_molecules_with_no_predictions = 1;

    if (verbose)
      cerr << "Will write molecules with no predictions\n";
  }

  if (cl.option_present ('h'))
  {
    do_not_compare_molecules_with_themselves = 1;

    if (verbose)
      cerr << "Will discard neighbours with zero distance and the same id as the target\n";
  }

  if (cl.option_present ('z'))
  {
    remove_leading_zeros = 1;

    if (verbose)
      cerr << "For an identifier match, leading zero's will be removed\n";
  }

  if (cl.option_present ('j'))
  {
    experimental_data_file_is_a_descriptor_file = 1;

    if (verbose)
      cerr << "Experimental data file is a descriptor file\n";
  }

  if (! cl.option_present ('N'))
  {
    cerr << "Must specify the models to evaluate via the -N option\n";
    usage (4);
  }

  resizable_array_p<IWString> string_models;

  IWString trgfp, dash_e;

  number_computations = count_models_to_evaluate (cl, 'N', trgfp, dash_e, string_models);

  if (0 == number_computations)
  {
    cerr << "Fatal error determining number of computations\n";
    usage (8);
  }

  if (trgfp.length())
    ;
  else if (cl.option_present ('p'))
    cl.value ('p', trgfp);
  else
  {
    cerr << "Must specify training set fingerprints via the -p option\n";
    usage (6);
  }

// We must have a source of experimental data

  if (dash_e.length())
    ;
  else if (cl.option_present ('E'))
    cl.value ('E', dash_e);
  else
  {
    cerr << "Must specify source of experimental data via the -E option\n";
    usage (4);
  }

  for (int i = 0; i < number_computations; i++)
  {
    no_predictions << ' ' << missing_value << ' ' << missing_value << ' ' << missing_value;
  }
  no_predictions += '\n';

  computation = new NN_Computation[number_computations];

  for (int i = 0; i < number_computations; i++)
  {
    const IWString & n = *(string_models[i]);

    assert (n.length());

    if (! computation[i].build (n))
    {
      cerr << "Invalid computation specification '" << n << "'\n";
      usage (8);
    }

    if (verbose)
      cerr << "number " << i << " is " << computation[i].estimate_type() << endl;
  }

  string_models.resize (0);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (4);
  }

  IWString_and_File_Descriptor output(1);

  int rc;
  if (classification_model)
    rc = nn_predictions_classification (cl, trgfp, dash_e, computation, number_computations, output);
  else
    rc = nn_predictions_continuous (cl, trgfp, dash_e, computation, number_computations, output);

  output.flush();

  if (verbose)
  {
    cerr << "Read " << tdts_read << " fingerprints for prediction\n";
  }

  delete [] computation;

  return (! rc);
}

int
main (int argc, char ** argv)
{
  int rc = nn_predictions (argc, argv);

  return rc;
}
