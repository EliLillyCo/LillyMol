/*
  Confusion matrix for classification problems

  Complicated by the need to have the data input be flexible.

  it can either be a single file, with real values and predicted values in that file,
  or it can be two separate files
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
using std::cerr;
using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Utilities/GFP_Tools/Metric.h"

const char * prog_name = NULL;

static int verbose = 0;

static int display_matrix = 1;

//static int convert_experimental_to_class = 0;
//static float convert_experimental_to_class_threshold = 0.0f;

static int prediction_file_has_header = 1;

static char input_separator = ' ';

static IWString positive_class;

static int enrichment_specified = 0;
static double bedroc_alpha = 0.20;
static double enrichment_cutoff = 0.0;
static double enrichment_factor_fraction = 0.5;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes classification accuracies\n";
  cerr << " -E <fname>     experimental values file\n";
  cerr << " -e <col>       experimental values in column <col>\n";
  cerr << " -p <col>       predicted values column <col>\n";
  cerr << " -np <col>      numeric predicted values in column <col>\n";
  cerr << " -n <c>         convert floating point data to class labels, < c and > c\n";
  cerr << " -n <l,h>       convert floating point data to class labels, < l and > h\n";
  cerr << " -pos <class>   compute PPV and NPV measures assuming 'class' is the positive class\n";
  cerr << " -nm            suppress display of the matrix\n";
  cerr << " -o <float>     cutoff for active/inactive (for enrichment metrics BEDROC EF ...)\n";
  cerr << " -a <float>     BEDROC alpha value (default " << bedroc_alpha << ")\n";
  cerr << " -f <float>     Enrichment Factor default fraction (default " << enrichment_factor_fraction << ")\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class AClass
{
  private:
    const IWString _name;
    int _n;
    IW_STL_Hash_Map_int _classified;

    float _sensitivity;
    float _ppv;
    float _specificity;
    float _npv;

  public:
    AClass(const IWString &);

    int debug_print(std::ostream & output) const;

    const IWString & name () const { return _name;}

    int times_predicted_in_class(const IWString & cname) const;

    void found_in(const IWString & c);

    void another_instance() { _n++;}

    int  n () const { return _n;}

    void compute_ppv(const AClass & negative_class, std::ostream & output);
    void compute_npv(const AClass & positive_class, std::ostream & output);

    double sum_measures_positive_class() const;
    double sum_measures_negative_class() const;

    double report(const IW_STL_Hash_Set & classes, std::ostream & output) const;

    void compute_tn_fn(int & tp, int & tn) const;
    void compute_tp_fp(int & tp, int & fp) const;
};

typedef IW_STL_Hash_Map<IWString, AClass *> IW_STL_Hash_Map_AClass;

AClass::AClass(const IWString & s) : _name(s)
{
  _n = 1;

  return;
}

int
AClass::debug_print(std::ostream & output) const
{
  output << "Class " << _name << " found " << _n << " times\n";

  return 1;
}

void
AClass::found_in(const IWString & c)
{
  auto f = _classified.find(c);
  if (f == _classified.end())
    _classified[c] = 1;
  else
    f->second = f->second + 1;

  return;
}


double
AClass::report(const IW_STL_Hash_Set & class_names,
               std::ostream & output) const
{
  int longest_class_name = 0;
  for (auto c : class_names)
  {
    if (c.length() > longest_class_name)
      longest_class_name = c.length();
  }

  output << _n << " instances of class " << _name << ", ";

  double fraction_class_correct = 0.0 ;   // our return code

  auto f = _classified.find(_name);
  if (f != _classified.end())
  {
    fraction_class_correct = static_cast<float>(f->second) / static_cast<float>(_n);
    output << "correct " << f->second << " " << fraction_class_correct << '\n';
  }
  else
    output << "correct 0 0.0\n";

  if (! display_matrix)
    return fraction_class_correct;

  for (auto c: class_names)
  {
    auto f = _classified.find(c);
    if (f != _classified.end())
    {
      const double ratio = static_cast<float>(f->second) / static_cast<float>(_n);
      output << std::setw(longest_class_name) << c << ' ' << f->second << ' ' << ratio << '\n';
    }
    else
      output << std::setw(longest_class_name) << c << " 0 0.0\n";
  }

  return fraction_class_correct;
}

double
AClass::sum_measures_positive_class() const
{
  return _sensitivity + _ppv;
}

double
AClass::sum_measures_negative_class() const
{
  return _specificity + _npv;
}

int
AClass::times_predicted_in_class(const IWString & cname) const
{
  const auto f = _classified.find(cname);

  if (f == _classified.end())
    return 0;

  return f->second;
}

void
AClass::compute_ppv(const AClass & negative_class, std::ostream & output)
{
  if (_classified.size() > 2)
  {
    cerr << "AClass::compute_ppv:not a two class problem " << _classified.size() << endl;
    return;
  }

  int true_positives = 0;
  const auto f = _classified.find(_name);
  if (f != _classified.end())
    true_positives = f->second;

  _sensitivity = static_cast<float>(true_positives)/ static_cast<float>(_n);
  output << "Sensitivity " << _sensitivity;

//const int false_positives = _n - true_positives;

  const double g = true_positives + negative_class.times_predicted_in_class(_name);

  if (g > 0)
  {
    _ppv = static_cast<float>(true_positives) / g;
    output << " PPV " << _ppv << "\n";
  }
  else
    output << " PPV 0\n";

  return;
}

void
AClass::compute_npv(const AClass & positive_class, std::ostream & output)
{
  if (_classified.size() > 2)
  {
    cerr << "AClass:;compute_npv:not a two class problem\n";
    return;
  }

  int true_negatives = 0;
  auto f = _classified.find(_name);
  if (f != _classified.end())
    true_negatives = f->second;

//const int false_negatives = _n - true_negatives;

  _specificity = static_cast<float>(true_negatives) / static_cast<float>(_n);
  output << "Specificity " << _specificity;

  const int h = positive_class.times_predicted_in_class(_name) + true_negatives;

  if (h > 0)
  {
    _npv = static_cast<float>(true_negatives) / h;
    output << " NPV " << _npv << "\n";
  }
  else
    output << " NPV undefined\n";

  return;
}

void
AClass::compute_tn_fn(int & tn, int & fn) const
{
  tn = 0;
  fn = 0;

  const auto f = _classified.find(_name);
  if (f != _classified.end())
    tn = f->second;

  fn = _n - tn;

  return;
}

void
AClass::compute_tp_fp(int & tp, int & fp) const
{
  tp = 0;
  fp = 0;
  const auto f = _classified.find(_name);
  if (f != _classified.end())
    tp = f->second;

  fp = _n - tp;

  return;
}

template <typename T>
int
get_next_token(const const_IWSubstring & buffer,
               T & token,
               int & i,
               char token_separator)
{
  if (' ' == token_separator)
    return buffer.nextword(token, i);

  return buffer.nextword_single_delimiter(token, i, token_separator);
}

class Class_Cutoff
{
  private:
    float _cutoff;
    IWString _low, _high;
  public:
    Class_Cutoff();

    int active() const { return _low.length();}

    float cutoff() const { return _cutoff;}
    const IWString & low_label() const { return _low;}
    const IWString & high_label() const { return _high;}

    int build(const Command_Line_v2 & cl, const char * flag);
};

Class_Cutoff::Class_Cutoff()
{
  _cutoff = 0.0f;
}

int
Class_Cutoff::build(const Command_Line_v2 & cl,
                    const char * flag)
{
  int got_cutoff = 0;

  const_IWSubstring s;
  for (int i = 0; cl.value(flag, s, i); ++i)
  {
    if (s.contains(','))
    {
      if (! s.split(_low, ',', _high) || 0 == _low.length() || 0 == _high.length())
      {
        cerr << "Class_Cutoff::build:invalid class labels '" << s << "'\n";
        return 0;
      }
    }
    else if (! s.numeric_value(_cutoff))
    {
      cerr << "Class_Cutoff::build:invalid cutoff " << s << "'\n";
      return 0;
    }
    else
      got_cutoff = 1;
  }

  if (! got_cutoff)
  {
    cerr << "Class_Cutoff::build:no cutoff specified\n";
    return 0;
  }

  if (0 == _low.length() || 0 == _high.length())
  {
    cerr << "Class_Cutoff::build:no labels specified\n";
    return 0;
  }

  return 1;
}

class DataPoint
{
  private:
    IWString _id;
    IWString _expt;
    IWString _pred;

    float _numeric_expt;    // don't really use this
    float _numeric_pred;

  public:
    DataPoint();

    const IWString & id() const { return _id;}

    float numeric_pred() const { return _numeric_pred;}

    void set_predicted_class (const IWString & s) { _pred = s;}
    void set_expt_class (const IWString & s) { _expt = s;}

    const IWString & expt() const { return _expt;}
    const IWString & pred() const { return _pred;}

    int build(const const_IWSubstring & buffer,
              const int numeric_expt_col, const int expt_col,
              const int numeric_pred_col, const int pred_col);

    int debug_print(std::ostream & output) const;
};

DataPoint::DataPoint()
{
  _numeric_expt = 0.0f;
  _numeric_pred = 0.0f;
}

int
DataPoint::debug_print(std::ostream & output) const
{
  output << _id << " expt " << _expt << " pred " << _pred << " numeric " << _numeric_pred << endl;
  return 1;
}

int
DataPoint::build(const const_IWSubstring & buffer,
                 const int numeric_expt_col, const int expt_col,
                 const int numeric_pred_col, const int pred_col)
{
  int i = 0;
  const_IWSubstring token;

  for (int c = 0; get_next_token(buffer, token, i, input_separator); ++c)
  {
    if (0 == c)
      _id = token;
    else if (expt_col == c)
      _expt = token;
    else if (pred_col == c)
      _pred = token;
    else if (numeric_expt_col == c)
    {
      if (! token.numeric_value(_numeric_expt))
      {
        cerr << "DataPoint::build:invalid numeric experimental value '" << token << "'\n";
        return 0;
      }
    }
    else if (numeric_pred_col == c)
    {
      if (! token.numeric_value(_numeric_pred))
      {
        cerr << "DataPoint::build:invalid numeric predicted value '" << token << "'\n";
        return 0;
      }
    }
  }

  return 1;
}

#ifdef NOT_USED_JJJJJJJJJJJJJJJJJJJJJJJJ
class Prediction
{
  private:
    IWString _id;
    float _v;
    IWString _c;

  public:
    Prediction();

    int build (const const_IWSubstring & buffer, const int value_col, const int class_col);

};

Prediction::Prediction()
{
  _v = 0.0f;

  return;
}

int
Prediction::build(const const_IWSubstring & buffer,
                  const int value_col,
                  const int class_col)
{
  int i = 0;
  const_IWSubstring token;

  for (int c = 0; get_next_token(buffer, token, i, input_separator); ++c)
  {
    if (0 == c)
      _id = token;
    else if (value_col == c)
    {
      if (! token.numeric_value(_v))
      {
        cerr << "Prediction:;build:invalid numeric '" << token << "'\n";
        return 0;
      }
    }
    else if (class_col == c)
    {
      _c = token;
    }
  }
  return 1;
}
#endif

/*
  Predicted data has
*/

#ifdef NOT_USED_JJJJJJJJJJJJJJJJJJJJJJJJ
class Predicted_Data
{
  private:
    IWString _response;

    resizable_array_p<Prediction> _pred;

  public:
    int build (const char * fname, const int value_col, const int class_col);
    int build (iwstring_data_source & input, const int value_col, const int class_col);
};

int
Predicted_Data::build(const char * fname,
                      const int value_col, 
                      const int class_col)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Predicted_Data::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input, value_col, class_col);
}
int
Predicted_Data::build(iwstring_data_source & input,
                      const int value_col, 
                      const int class_col)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Predicted_Data::build:cannot read header record\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    Prediction * p = new Prediction();

    if (! p->build(buffer, value_col, class_col))
    {
      cerr << "Predicted_Data::build:invalid data '" << buffer << "'\n";
      delete p;
      return 0;
    }

    _pred.add(p);
  }

  return _pred.number_elements();
}
#endif

static void
compute_matthews(const int tp,
                 const int tn,
                 const int fp,
                 const int fn,
                 std::ostream & output)
{
//cerr << "Matthews tp " << tp << " tn " << tn << " fp " << fp << " fn " << fn << endl;
  if (0 == (tp+fn) ||
      0 == (tn+fp) ||
      0 == (tp+fp) ||
      0 == (tn+fn))
    output << "MatthewsCoefficient 0\n";
  else
  {
    const float m = static_cast<float>(tp*tn - fp*fn) /
                        sqrt(static_cast<double>(tp+fn)*static_cast<double>(tn+fp)*static_cast<double>(tp+fp)*static_cast<double>(tn+fn));
    output << "MatthewsCoefficient " << m << "\n";
  }
}

static void
compute_accuracy(const int tp,
                 const int tn,
                 const int fp,
                 const int fn,
                 std::ostream & output)
{
  output << "Accuracy " << static_cast<float>(tp + tn) / static_cast<float>(tp + tn + fp + fn) << "\n";
}

#ifdef NOT_USED_QQQQQQQQQQQQQA
static int
read_predicted_data(const const_IWSubstring & buffer,
                    const int expt_col,
                    const int class_col,
                    IW_STL_Hash_Map_String & expt)
{
  IWString id;
  int i = 0;
  if (! buffer.nextword(id, i))
    return 0;

  const_IWSubstring token;

  for (int j = 0; j < expt_col; ++j)
  {
    if (! buffer.nextword(token, i))
      return 0;
  }

  if (! convert_experimental_to_class)
  {
    expt[id] = token;
    return 1;
  }

  float v;
  if (! token.numeric_value(v))
  {
    cerr << "Invalid numeric '" << token << "'\n";
    return 0;
  }

  if (v <= convert_experimental_to_class_threshold)
    expt[id] = "-1";
  else
    expt[id] = "1";

  return 1;
}
#endif

static int
read_predicted_data(iwstring_data_source & input,
                    const int expt_class_col,
                    const int pred_value_col,
                    const int pred_class_col,
                    resizable_array_p<DataPoint> & mydata,
                    const IW_STL_Hash_Map_String & expt)
{
  const_IWSubstring buffer;

  if (prediction_file_has_header && ! input.next_record(buffer))
  {
    cerr << "Cannot read experimental data header\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    DataPoint * d = new DataPoint;

    if (! d->build(buffer, -1, expt_class_col, pred_value_col, pred_class_col))
    {
      cerr << "Invalid experimental data record '" << buffer << "'\n";
      return 0;
    }

    mydata.add(d);

    if (0 == expt.size())
      continue;

    const auto f = expt.find(d->id());

    if (f == expt.end())
    {
      cerr << "read_predicted_data:no experimental data for '" << d->id() << "'\n";
      return 0;
    }

    d->set_expt_class(f->second);
  }

  return mydata.number_elements();
}

static int
read_predicted_data(const char * fname,
                    const int expt_class_col,
                    const int pred_value_col,
                    const int pred_class_col,
                    resizable_array_p<DataPoint> & mydata,
                    const IW_STL_Hash_Map_String & expt)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_predicted_data:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_predicted_data(input, expt_class_col, pred_value_col, pred_class_col, mydata, expt);
}

#ifdef NOT_USED_PPPPPPPPPPPPPPPPP
static DataPoint *
get_data_point(const const_IWSubstring & buffer,
               const IW_STL_Hash_Map_int & id_to_ndx,
               const resizable_array_p<DataPoint> & mydata)
{
  IWString id;
  int i = 0;
  if (! get_next_token(buffer, id, i, input_separator))
  {
    cerr << "get_next_token:cannot extract identifier '" << buffer << "'\n";
    return nullptr;
  }

  const auto f = id_to_ndx.find(id);

  if (f == id_to_ndx.end())
  {
    cerr << "get_next_token:no experimental data for '" << id << "'\n";
    return nullptr;
  }

  return mydata[f->second];
}
#endif

static int
read_experimental_data_record(const const_IWSubstring & buffer,
                              const int class_col,
                              IW_STL_Hash_Map_String & expt)
{
  int i = 0;
  const_IWSubstring token;

  IWString id;

  for (int c = 0; get_next_token(buffer, token, i, input_separator); ++c)
  {
    if (0 == c)
      id = token;
    else if (c == class_col)
    {
      expt[id] = token;
      return 1;
    }
  }

  return 0;
}

static int
read_experimental_data(iwstring_data_source & input,
                       const int class_col,
                       IW_STL_Hash_Map_String & expt)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    cerr << "read_experimental_data:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    if (! read_experimental_data_record(buffer, class_col, expt))
    {
      cerr << "Invalid experimental input " << buffer << "'\n";
      return 0;
    }
  }

  return expt.size();
}

static int
read_experimental_data(const char * fname,
                       const int class_col,
                       IW_STL_Hash_Map_String & expt)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_experimental_data(input, class_col, expt);
}

static int
gather_classes(const resizable_array_p<DataPoint> & mydata, 
               IW_STL_Hash_Map_AClass & classes)
{
  for (auto f : mydata)
  {
//  f->debug_print(cerr);
    const IWString & p = f->expt();
//  cerr << " predicted value " << p << endl;

    auto j = classes.find(p);

    if (j == classes.end())
      classes[p] = new AClass(p);
    else
      j->second->another_instance();
  }

  if (verbose)
  {
    for (auto i : classes)
    {
      cerr << "Found class '" << i.first << "' in " << i.second->n() << " items\n";
    }
  }

  return 1;
}

class DataPoint_Comparator_Numeric
{
  private:
  public:
    int operator() (const DataPoint * d1, const DataPoint * d2) const;
};

int
DataPoint_Comparator_Numeric::operator() (const DataPoint * d1, const DataPoint * d2) const
{
  const auto v1 = d1->numeric_pred();
  const auto v2 = d2->numeric_pred();

  if (v1 < v2)
    return 1;

  if (v1 > v2)
    return -1;

  return 0;
}

class DataPoint_Comparator
{
  private:
  public:
    int operator() (const DataPoint * d1, const DataPoint * d2) const;
};

int
DataPoint_Comparator::operator() (const DataPoint * d1, const DataPoint * d2) const
{
  const auto & v1 = d1->pred();
  const auto & v2 = d2->pred();

  return v1.strcmp(v2);
}

/*
  The tricky part here is preparing the F array for buildFromSortedArray

  It needs a floating point array of predicted values and a cutoff. 

  
*/

static int
do_enrichment_measures(resizable_array_p<DataPoint> & mydata,
                       std::ostream & output)
{
  DataPoint_Comparator_Numeric dpc;

  mydata.iwqsort(dpc);

  const int n = mydata.number_elements();

  float * f = new float[n]; std::unique_ptr<float[]> free_f(f);

  for (int i = 0; i < n; ++i)
  {
    const auto * x = mydata[i];

    const auto v = x->numeric_pred();

    if (v < enrichment_cutoff)     // would be predicted negative
    {
      if (positive_class == x->expt())     // incorrect prediction
        f[i] = -1.0f;
      else
        f[i] = 1.0f;
    }
    else                                   // would be predicted positive
    {
      if (positive_class == x->expt())     // good prediction
        f[i] = 1.0f;
      else
        f[i] = -1.0f;
    }
  }

  int above = 0;
  int below = 0;
  for (int i = 0; i < n; ++i)
  {
    if (f[i] < 0.0f)
      below++;
    else
      above++;
  }

  if (0 == below || 0 == above)
  {
    cerr << "No values around enrichment cutoff " << enrichment_cutoff << ", " << below << " points below and " << above << " above. No enrichment measures\n";
    return 0;
  }

  if (verbose)
    cerr << "Given enrichment cutoff " << enrichment_cutoff << ", " << below << " values below and " << above << " values above\n";

  Enrichment enrichment;

//enrichment.buildFromSortedArray(f, n, enrichment_cutoff);
  enrichment.buildFromSortedArray(f, n, 0.0f);

  calculateBEDROC bedroc(bedroc_alpha);

  const auto bedroc_value = bedroc(enrichment);

  output << "BEDROC " << setprecision(4) << static_cast<float>(bedroc_value) << endl;

  calculateROC roc;

  const double roc_value = roc(enrichment);

  output << "ROC " << setprecision(4) << static_cast<float>(roc_value) << endl;

  calculateEF ef(enrichment_factor_fraction);
  const double ef_value = ef(enrichment);

  output << "EF " << static_cast<float>(ef_value) << endl;

  return 1;
}

static int
confusion_matrix(resizable_array_p<DataPoint> & mydata, 
                 IW_STL_Hash_Map_AClass & classes,
                 std::ostream & output)
{
  int number_correct = 0;

  IW_STL_Hash_Set class_names;

  for (auto i : classes)
  {
    class_names.insert(i.first);
  }

  for (auto i : mydata)
  {
    const IWString & e = i->expt();

    const auto f = classes.find(e);

    const IWString & p = i->pred();

    f->second->found_in(p);

    if (e == p)
      number_correct++;
  }

  if (verbose)
    cerr << "Found info on " << mydata.number_elements() << " items in " << classes.size() << " classes\n";

  resizable_array<AClass *> classes_sorted;
  for (auto i : classes)
  {
    classes_sorted.add(i.second);
  }

  classes_sorted.iwqsort_lambda([] (const AClass * c1, const AClass * c2) { 
     if (c1->n() < c2->n())
       return -1;
     if (c1->n() > c2->n())
       return 1;

     return c1->name().strcmp(c2->name());
  });

  double total_f = 0.0;
  double min_fraction_correct = 1.0;

  for (auto i : classes_sorted)
  {
    const double f = i->report(class_names, output);

    total_f += f;
    if (f < min_fraction_correct)
      min_fraction_correct = f;
  }

  output << "Fraction correct " << setw(3) << setprecision(3) << static_cast<float>(number_correct)/ static_cast<float>(mydata.number_elements()) << endl;
  output << "Average fraction correct " << setw(3) << setprecision(3) <<  static_cast<float>(total_f) / static_cast<float>(classes.size()) << endl;
  output << "Min fraction correct " << setw(3) << setprecision(3) <<  min_fraction_correct << endl;

  if (enrichment_specified)
    do_enrichment_measures(mydata, output);

  if (0 == positive_class.length())
    return 1;

  if (! classes.contains(positive_class))
  {
    cerr << "Positive class is '" << positive_class << "' but no instances it\n";
    output << "Sensitivity 0.0 PPV 0.0\n";
    output << "Specificity 0.0 NPV 0.0\n";
    return 0;
  }

  if (2 != classes.size())
  {
    cerr << "Positive class for sensitivity/specificity computation, but " << classes.size() << " classes\n";
    return 0;
  }

  AClass * p = nullptr;
  AClass * n = nullptr;

  for (auto i : classes)
  {
    if (positive_class == i.first)
      p = i.second;
    else
      n = i.second;
  }

  if (nullptr == p || nullptr == n)
  {
    cerr << "No positive class '" << positive_class << "'\n";
    return 0;
  }

  p->compute_ppv(*n, output);
  n->compute_npv(*p, output);

  const double s = p->sum_measures_positive_class() + n->sum_measures_negative_class();
  output << "Sum all measures " << setw(4) << s << "\n";

  int tp, fp, tn, fn;
  p->compute_tp_fp(tp, fp);
  n->compute_tn_fn(tn, fn);
  compute_matthews(tp, tn, fp, fn, output);
  compute_accuracy(tp, tn, fp, fn, output);

  return 1;
}

static void
reset_class_assignments(resizable_array_p<DataPoint> & mydata, 
                        const Class_Cutoff & cutoff)
{
  int set_high = 0;

  const int n = mydata.number_elements();

  for (int i = 0; i < n; ++i)
  {
    DataPoint * p = mydata[i];

    if (p->numeric_pred() > cutoff.cutoff())
    {
      p->set_predicted_class(cutoff.high_label());
      set_high++;
    }
    else
    {
      p->set_predicted_class(cutoff.low_label());
    }
  }

  if (verbose)
    cerr << "At cutoff " << cutoff.cutoff() << " set " << set_high << " to class " << cutoff.high_label() << " " << (static_cast<float>(set_high) / static_cast<float>(n)) << endl;

  return;
}

static int
confusion_matrix (int argc, char ** argv)
{
  Command_Line_v2 cl(argc, argv, "-v-E=sfile-e=ipos-p=ipos-np=ipos-pos=s-nm-i=s-n=s-o=f-a=fraction-f=fraction");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  int expt_col = -1;

  if (cl.option_present('e'))
  {
    cl.value('e', expt_col);

    if (verbose)
      cerr << "Experimental data in column " << expt_col << endl;

    expt_col--;
  }

  IW_STL_Hash_Map_String expt;

  if (cl.option_present('E'))
  {
    int myexptcol = expt_col;
    if (myexptcol < 0)
      myexptcol = 1;

    const char * fname = cl.option_value('E');

    if (! read_experimental_data(fname, myexptcol, expt))
    {
      cerr << "Cannot read experimental data from '" << fname << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << expt.size() << " experimental data values from '" << fname << "'\n";

    expt_col = -1;    // so nobody goes looking for it in the main file
  }

// Typical input will be ID NUMERIC CLASS

  int predicted_value_col = 1;
  int predicted_class_col = 2;

  if (cl.option_present('p'))
  {
    cl.value('p', predicted_class_col);

    if (verbose)
      cerr << "Predicted class in col " << predicted_class_col << endl;

    predicted_class_col--;
  }

  if (cl.option_present("np"))
  {
    cl.value("np", predicted_value_col);

    if (verbose)
      cerr << "Predicted numeric values in col " << predicted_value_col << endl;

    predicted_value_col--;
  }

  if (predicted_class_col == predicted_value_col)
  {
    cerr << "Predicted class column and predicted class numeric value column identical\n";
    return 1;
  }

  if (0 == expt.size() && predicted_class_col == expt_col)
  {
    cerr << "Reading predictions and experimental from same file, but columns identical, impossible\n";
    return 1;
  }

  if (cl.option_present("nm"))
  {
    display_matrix = 0;
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "Unrecognised input separator specification '" << i << "'\n";
      return 1;
    }

    input_separator = i[0];
  }

  if (cl.option_present("pos"))
  {
    cl.value("pos", positive_class);

    if (verbose)
      cerr << "Positive class " << positive_class << endl;
  }

  Class_Cutoff cutoff;

  if (cl.option_present('n'))
  {
    if (! cutoff.build(cl, "n"))
    {
      cerr << "Cannot initialise class cutoff specification (-n)\n";
      return 1;
    }
  }

  if (cl.option_present('a'))
  {
    cl.value('a', bedroc_alpha);

    if (verbose)
      cerr << "BEDROC alpha set to " << bedroc_alpha << endl;
  }

  if (cl.option_present('f'))
  {
    if (! cl.value('f', enrichment_factor_fraction) || enrichment_factor_fraction <= 0.0 || enrichment_factor_fraction >= 1.0)
    {
      cerr << "The enrichment factor fraction (-f) must be a valid fraction\n";
      usage(2);
    }

    if (verbose)
      cerr << "Enrichment factor " << enrichment_factor_fraction << endl;

    enrichment_specified = 1;
  }

  if (cl.option_present('o'))
  {
    cl.value('o', enrichment_cutoff);

    if (verbose)
      cerr << "Enrichment cutoff " << enrichment_cutoff << endl;

    enrichment_specified = 1;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    resizable_array_p<DataPoint> mydata;

    if (! read_predicted_data(cl[0], expt_col, predicted_value_col, predicted_class_col, mydata, expt))
    {
      cerr << "Cannot read predicted values from '" << cl[i] << "'\n";
      return 1;
    }

    if (cutoff.active())
      reset_class_assignments(mydata, cutoff);

    IW_STL_Hash_Map_AClass classes;

    if (! gather_classes(mydata, classes))
    {
      cerr << "Cannot determine classes in '" << cl[0] << "'\n";
      return 1;
    }

    if (! confusion_matrix(mydata, classes, std::cout))
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return 1;
    }

    for (auto c : classes)
    {
      delete c.second;
    }
  }

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = confusion_matrix(argc, argv);

  return rc;
}
