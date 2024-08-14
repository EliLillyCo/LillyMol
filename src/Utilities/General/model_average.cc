/*
  Performs model averaging for regression models
*/

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static IWString missing_value(".");

static Accumulator<double> score_acc;

#define MODEL_AVERAGE_TYPE_AVERAGE 1
#define MODEL_AVERAGE_TYPE_MAXVAL 2
#define MODEL_AVERAGE_TYPE_VOTING 4
#define MODEL_AVERAGE_TYPE_SUM 8

#define MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_ONE_VS_ALL 5
#define MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_PAIRS 6

static int classification_model = 0;

static float class_cutoff = 0.0;

static int write_composite_prediction_first = 0;

static int just_write_composite_prediction = 0;

static IWString not_model_string("QJNOTZ_");

static int discern_class_labels_from_input_file = 0;
static IWString class_label[2] = {" -1", " 1"};

static int include_range_average_std = 0;

static char input_separator = ' ';
static char output_separator = ' ';

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(
Performs model averaging for regression models. All values in one file
 -c <col>       model in column <col>. c1,c2,c3 and c1-c4 also recognised.
 -c <col>,w=0.2 model in column <col>, specify relative weight
 -a ...         type of averaging to do
 -a ave         average of all values (default)
 -a max         final result is max of all values encountered
 -a sum         final result is sum of all values encountered
 -a pairs       many-class model, pair-wise models built
 -a ovsall      many-class model, individual built
 -C <cutoff>    classification model, convert result to class
 -d             discern class labels from input
 -f             write composite prediction first rather than last
 -n             just write the composite prediction
 -r             include range/average/standard deviation of component models
 -v             verbose output
)";
  // clang-format on

  exit(rc);
}

class Model
{
 private:
  int _column;
  double _weight;

  int _from_arbor;

  //  private functions

  int
  _parse_numeric_column_specification(const const_IWSubstring& c);

 public:
  Model();

  int
  build(const const_IWSubstring&);

  int
  from_arbor() const
  {
    return _from_arbor;
  }

  void
  set_weight(double w)
  {
    _weight = w;
  }

  double
  weight() const
  {
    return _weight;
  }

  int
  column() const
  {
    return _column;
  }

  void
  set_column_and_weight(int c, double w)
  {
    _column = c;
    _weight = w;
  }
};

Model::Model()
{
  _column = -1;
  _weight = static_cast<double>(0.0);

  _from_arbor = 0;

  return;
}

int
Model::_parse_numeric_column_specification(const const_IWSubstring& c)
{
  if (!c.numeric_value(_column) || _column < 1) {
    cerr << "Model::build:invalid column specification '" << c << "'\n";
    return 0;
  }

  if (1 == _column) {
    cerr << "Model::column 1 is assumed to be the identifier column\n";
    return 0;
  }

  _column--;

  return 1;
}

#ifdef NOT_USED_ANY_MORE
int
Model::build(const const_IWSubstring& buffer)
{
  const_IWSubstring mybuffer(buffer);

  if (mybuffer.starts_with("A:")) {
    mybuffer.remove_leading_chars(2);
    _from_arbor = 1;
  }

  int comma = mybuffer.index(',');

  if (comma < 0) {
    return _parse_numeric_column_specification(mybuffer);
  }

  const_IWSubstring c, w;
  if (!mybuffer.split(c, ',', w) || 0 == c.length() || w.length() < 3) {
    cerr << "Model::build:invalid column,weight specification '" << mybuffer << "'\n";
    return 0;
  }

  if (!_parse_numeric_column_specification(c)) {
    return 0;
  }

  if (!w.starts_with("w=")) {
    cerr << "Model::build:the weight specification must start with 'w=', '" << mybuffer
         << "' is invalid\n";
    return 0;
  }

  w.remove_leading_chars(2);

  if (!w.numeric_value(_weight) || _weight <= 0.0 || _weight >= 1.0) {
    cerr << "Model::build:the weight must be a valid fraction\n";
    return 0;
  }

  return 1;
}
#endif

static int nmodels = 0;

static Model* model = NULL;

/*
  We need a quick way of knowing which column corresponds to which model
*/

static extending_resizable_array<int> col_to_model(-1);

class Average
{
 private:
  int _n;
  double _sum;
  double _wsum;

  Accumulator<double> _acc;

  //  private functions

  void
  _default_values();

 public:
  Average();

  void
  reset();

  const char*
  text_description() const
  {
    return "avescore";
  }

  int
  extra(double v, double w, int);

  double
  final_result();

  template <typename T>
  int
  write_range_ave_std(T&);
};

void
Average::_default_values()
{
  _n = 0;
  _sum = 0.0;
  _wsum = 0.0;
}

Average::Average()
{
  _default_values();
}

void
Average::reset()
{
  _acc.reset();

  _default_values();
}

int
Average::extra(double v, double w, int from_arbor)
{
  _acc.extra(v);

  _n += 1;

  if (from_arbor) {
    v = 2.0 * v - 1.0;
  }

  _sum += v * w;
  _wsum += w;

  return 1;
}

double
Average::final_result()
{
  if (0 == _n) {
    cerr << "Average::final_result:no data\n";
    return 0.0;
  }

  if (0.0 == _wsum) {
    cerr << "Average::final_result:no weight\n";
    return 0.0;
  }

  return _sum / _wsum;
}

template <typename T>
int
Average::write_range_ave_std(T& output)
{
  output << output_separator << _acc.minval() << output_separator << _acc.maxval()
         << output_separator << static_cast<float>(_acc.maxval() - _acc.minval())
         << output_separator << static_cast<float>(_acc.average()) << output_separator
         << static_cast<float>(sqrt(_acc.variance()));

  return 1;
}

template int
Average::write_range_ave_std(IWString_and_File_Descriptor&);

class MaxVal
{
 private:
  int _n;
  double _max;

  Accumulator<double> _acc;

  //  private functions

  void
  _default_values();

 public:
  MaxVal();

  void
  reset();

  int
  extra(double v, double w, int);

  const char*
  text_description() const
  {
    return "maxscore";
  }

  double
  final_result() const;

  template <typename T>
  int
  write_range_ave_std(T&);
};

void
MaxVal::_default_values()
{
  _n = 0;

  return;
}

MaxVal::MaxVal()
{
  _default_values();
}

void
MaxVal::reset()
{
  _acc.reset();

  _default_values();
}

int
MaxVal::extra(double v, double w, int from_arbor)
{
  if (from_arbor) {
    v = 2.0 * v - 1.0;
  }

  _acc.extra(v);

  if (0 == _n) {
    _max = v;
  } else if (v > _max) {
    _max = v;
  }

  _n++;

  return 1;
}

double
MaxVal::final_result() const
{
  if (0 == _n) {
    cerr << "MaxVal::final_result:no values\n";
    return 0.0;
  }

  return _max;
}

template <typename T>
int
MaxVal::write_range_ave_std(T& output)
{
  output << output_separator << _acc.minval() << output_separator << _acc.maxval()
         << output_separator << static_cast<float>(_acc.maxval() - _acc.minval())
         << output_separator << static_cast<float>(_acc.average()) << output_separator
         << static_cast<float>(sqrt(_acc.variance()));

  return 1;
}

template int
MaxVal::write_range_ave_std(IWString_and_File_Descriptor&);

class Voting
{
 private:
  resizable_array<float> _v;
  resizable_array<int> _c;

 public:
  void
  reset();

  int
  extra(double v, double w, int);

  const char*
  text_description() const
  {
    return "vote";
  }

  double
  final_result() const;

  template <typename T>
  int
  write_range_ave_std(T&)
  {
    return 0;
  }  // not implemented
};

void
Voting::reset()
{
  _v.resize_keep_storage(0);
  _c.resize_keep_storage(0);

  return;
}

int
Voting::extra(double v, double w, int from_arbor)
{
  for (int i = 0; i < _v.number_elements(); i++) {
    if (fabs(v - _v[i]) < 1.0e-05) {
      _c[i]++;
      return 1;
    }
  }

  _v.add(v);
  _c.add(1);

  return 1;
}

double
Voting::final_result() const
{
  int n = _v.number_elements();

  if (0 == n) {
    cerr << "Voting::final_result:no data\n";
    return 0.0;
  }

  int max_count = _c[0];
  double item_with_most_votes = _v[0];

  for (int i = 1; i < n; i++) {
    if (_c[i] > max_count) {
      max_count = _c[i];
      item_with_most_votes = _v[i];
    }
  }

  return item_with_most_votes;
}

class Sum
{
 private:
  double _sum;

  Accumulator<double> _acc;

 public:
  Sum();

  void
  reset();

  int
  extra(double v, double w, int);

  const char*
  text_description() const
  {
    return "sum";
  }

  double
  final_result() const;

  template <typename T>
  int
  write_range_ave_std(T&);
};

Sum::Sum()
{
  _sum = 0.0;
}

void
Sum::reset()
{
  _acc.reset();
  _sum = 0.0;
}

int
Sum::extra(double v, double w, int)
{
  _sum += v;

  _acc.extra(v);

  return 1;
}

double
Sum::final_result() const
{
  return _sum;
}

template <typename T>
int
Sum::write_range_ave_std(T& output)
{
  output << output_separator << _acc.minval() << output_separator << _acc.maxval()
         << output_separator << static_cast<float>(_acc.maxval() - _acc.minval())
         << output_separator << static_cast<float>(_acc.average()) << output_separator
         << static_cast<float>(sqrt(_acc.variance()));

  return 1;
}

template int
Sum::write_range_ave_std(IWString_and_File_Descriptor&);

template <typename T>
int
model_average2(const const_IWSubstring& buffer, T& c,
               IWString_and_File_Descriptor& output)
{
  c.reset();

  const_IWSubstring token;
  int i = 0;

  int n = 0;

  for (int col = 0; buffer.nextword(token, i, input_separator); col++) {
    int ndx = col_to_model[col];

    if (ndx < 0) {
      continue;
    }

    if (missing_value == token) {
      continue;
    }

    double v;
    if (!token.numeric_value(v)) {
      cerr << "Invalid numeric '" << token << "'\n";
      return 0;
    }

    const Model& m = model[ndx];

    c.extra(v, m.weight(), m.from_arbor());

    n++;
  }

  if (0 == n) {
    output << ' ' << missing_value;

    return 1;
  }

  float result = c.final_result();

  score_acc.extra(result);

  output << ' ' << static_cast<float>(result);

  if (include_range_average_std) {
    c.write_range_ave_std(output);
  }

  if (!classification_model) {  // no class label to write
    ;
  } else if (result >= class_cutoff) {
    output << class_label[1];
  } else {
    output << class_label[0];
  }

  return 1;
}

template int
model_average2(const const_IWSubstring&, Average&, IWString_and_File_Descriptor&);
template int
model_average2(const const_IWSubstring&, MaxVal&, IWString_and_File_Descriptor&);
template int
model_average2<Voting>(const_IWSubstring const&, Voting&, IWString_and_File_Descriptor&);
template int
model_average2<Sum>(const_IWSubstring const&, Sum&, IWString_and_File_Descriptor&);

/*
  We have a classification problem. We need to discern that the class labels are
  We assume that every column being processed is followed by a class label assignment
*/

static int
do_discern_class_labels_from_input_file2(iwstring_data_source& input,
                                         IWString* class_label)
{
  class_label[0] = "";
  class_label[1] = "";

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    float previous_numeric = 0.0;

    const_IWSubstring token;
    int i = 0;

    for (int col = 0; buffer.nextword(token, i, input_separator); col++) {
      if (0 == col) {
        continue;
      }

      int ndx = col_to_model[col];

      if (ndx >= 0)  // has numeric values, not processed
      {
        if (!token.numeric_value(previous_numeric)) {
          cerr << "discern_class_labels_from_input_file:non numeric '" << buffer << "'\n";
          return 0;
        }

        continue;
      }

      if (class_label[0] == token) {
        ;
      } else if (class_label[1] == token) {
        ;
      } else if (previous_numeric >= class_cutoff && 0 == class_label[1].length()) {
        class_label[1] = token;
      } else if (previous_numeric < class_cutoff && 0 == class_label[0].length()) {
        class_label[0] = token;
      } else {
        cerr << "More than two class labels, found '" << class_label[0] << "' and '"
             << class_label[1] << "'\n";
        return 0;
      }
    }
  }

  if (verbose) {
    cerr << "Class label translation, values above " << class_cutoff << " '"
         << class_label[1] << "' below '" << class_label[0] << "'\n";
  }

  for (auto i = 0; i < 2; ++i) {
    class_label[i].insert(' ', 0);
  }

  return 1;
}

static int
do_discern_class_labels_from_input_file(iwstring_data_source& input,
                                        IWString* class_label)
{
  off_t o = input.tellg();

  int rc = do_discern_class_labels_from_input_file2(input, class_label);

  if (!input.seekg(o)) {
    cerr << "discern_class_labels_from_input_file2:yipes, cannot seek back to start of "
            "file\n";
    return 0;
  }

  return rc;
}

template <typename C>
int
model_average(iwstring_data_source& input, C& c, IWString_and_File_Descriptor& output)
{
  IWString buffer;  // must save state!

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (discern_class_labels_from_input_file) {
    do_discern_class_labels_from_input_file(input, class_label);
  }

  if (write_composite_prediction_first) {
    const_IWSubstring token;
    int i = 0;
    buffer.nextword(token, i, input_separator);

    output << token;
    output << ' ' << c.text_description();
    if (include_range_average_std) {
      output << " MA.min MA.max MA.range MA.ave MA.std";
    }
    if (classification_model) {
      output << " C.CLASS";
    }

    if (just_write_composite_prediction) {
      ;
    } else {
      while (buffer.nextword(token, i, input_separator)) {
        output << ' ' << token;
      }
    }
  } else {
    output << buffer << ' ' << c.text_description();
    if (include_range_average_std) {
      output << " MA.min MA.max MA.range MA.std";
    }
    if (classification_model) {
      output << " C.CLASS";
    }
  }

  output << '\n';

  if (write_composite_prediction_first) {
    while (input.next_record(buffer)) {
      int i = buffer.index(' ');

      output.strncat(buffer, i);

      if (!model_average2(buffer, c, output)) {
        return 0;
      }

      buffer.remove_leading_chars(i + 1);

      if (just_write_composite_prediction) {
        output << '\n';
      } else {
        output << ' ' << buffer << '\n';
      }

      output.write_if_buffer_holds_more_than(32768);
    }
  } else {
    while (input.next_record(buffer)) {
      output << buffer;

      if (!model_average2(buffer, c, output)) {
        return 0;
      }

      output << '\n';

      output.write_if_buffer_holds_more_than(32768);
    }
  }

  return 1;
}

template int
model_average(iwstring_data_source&, Average&, IWString_and_File_Descriptor&);
template int
model_average(iwstring_data_source&, MaxVal&, IWString_and_File_Descriptor&);
template int
model_average<Voting>(iwstring_data_source&, Voting&, IWString_and_File_Descriptor&);
template int
model_average<Sum>(iwstring_data_source&, Sum&, IWString_and_File_Descriptor&);

static int
parse_column_weight(const_IWSubstring token,  // our own copy
                    resizable_array<int>& col, resizable_array<double>& w)
{
  int comma = token.index(',');
  if (comma < 0) {
    return 0;
  }

  const_IWSubstring sc;
  token.from_to(0, comma - 1);

  int c;
  if (!sc.numeric_value(c) || c < 1) {
    cerr << "Invalid column specification '" << token << "'\n";
    return 0;
  }

  token.remove_leading_chars(comma);

  if (!token.starts_with("w=")) {
    cerr << "Invalid weight\n";
    return 0;
  }

  token.remove_leading_chars(2);  // remove 'w='

  double wt;

  if (!token.numeric_value(wt) || wt <= 0.0 || wt >= 1.0) {
    cerr << "invalid weight\n";
    return 0;
  }

  col.add(c);
  w.add(wt);

  return 1;
}

/*
  Gets kind of complicated

  -c 2,3

  means two models, in columns 2 and 3, equally weighted.

  -c 2,w=0.4 -c 3,w=0.6

  means models in columns 2 and 3, but this time, weights are different

  We also handle

  -c 2-8
*/

static int
gather_columns_and_weights(const Command_Line& cl, resizable_array<int>& col,
                           resizable_array<double>& w)
{
  const_IWSubstring s;

  for (int i = 0; cl.value('c', s, i); ++i) {
    int c;
    if (s.numeric_value(c) && c > 0) {
      col.add(c - 1);
      w.add(-1.0);
      continue;
    }

    if (s.contains(',')) {
      if (s.contains(",w=")) {
        int j = 0;
        const_IWSubstring token;

        while (s.nextword(token, j)) {
          if (!parse_column_weight(token, col, w)) {
            cerr << "Invalid column,weight specification '" << token << "'\n";
            return 0;
          }
        }
      } else {
        const_IWSubstring token;
        int j = 0;

        while (s.nextword(token, j, ',')) {
          int c;
          if (!token.numeric_value(c) || c <= 0) {
            cerr << "Invalid column specification '" << token << "'\n";
            return 0;
          }

          col.add(c);
          w.add(-1.0);
        }
      }
    } else if (s.contains('-')) {
      const_IWSubstring sc1, sc2;
      int c1, c2;
      if (s.split(sc1, '-', sc2) && sc1.numeric_value(c1) > 0 &&
          sc2.numeric_value(c2) > 0 && c1 < c2) {
        for (auto j = c1; j <= c2; ++j) {
          col.add(j - 1);
          w.add(-1.0);
        }
      } else {
        cerr << "Invalid column range specification '" << s << "'\n";
        return 0;
      }
    } else {
      cerr << "Do not know how to process model column specification '" << s << "'\n";
      return 0;
    }
  }

  const int nc = col.number_elements();
  const int nw = w.number_elements();

  if (nc != nw)  // gack, I have messed up
  {
    cerr << "Mismatch between columns entered " << nc << " and weights " << nw << '\n';
    return 0;
  }

  if (1 == nc) {
    cerr << "Only one column specified\n";
    return 1;
  }

  // Default behaviour is to assign a default weight to all unweighted models

  double weights_specified = 0.0;
  int items_without_weights = 0;

  for (int i = 0; i < nw; i++) {
    if (w[i] < 0.0) {
      items_without_weights++;
    } else {
      weights_specified += w[i];
    }
  }

  if (weights_specified > 1.0) {
    cerr << "Sorry, weights need to sum to 1.0, " << weights_specified << '\n';
    return 0;
  }

  if (items_without_weights) {
    double weight_for_unassigned_columns =
        (1.0 - weights_specified) / static_cast<double>(items_without_weights);

    for (int i = 0; i < nw; i++) {
      if (w[i] < 0.0) {
        w[i] = weight_for_unassigned_columns;
      }
    }
  }

  return nw;
}

template <typename T>
int
model_average(const char* fname, T& c, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  return model_average(input, c, output);
}

template int
model_average(const char*, Average&, IWString_and_File_Descriptor&);
template int
model_average(const char*, MaxVal&, IWString_and_File_Descriptor&);
template int
model_average<Voting>(char const*, Voting&, IWString_and_File_Descriptor&);
template int
model_average<Sum>(char const*, Sum&, IWString_and_File_Descriptor&);

static int
model_average_many_class_one_vs_all2(const const_IWSubstring& buffer,
                                     IWString_and_File_Descriptor& output)
{
  const_IWSubstring token;
  int i = 0;

  if (!buffer.nextword(token, i, input_separator)) {
    cerr << "Cannot read identifier\n";
    return 0;
  }

  output << token << ' ';

  float max_pred = -std::numeric_limits<float>::max();
  IWString pred;

  while (buffer.nextword(token, i, input_separator)) {
    float p;
    if (!token.numeric_value(p)) {
      cerr << "Invalid floating point input '" << token << "'\n";
      return 0;
    }

    if (!buffer.nextword(token, i, input_separator)) {
      cerr << "Huh, cannot extract class after numeric prediction\n";
      return 0;
    }

    if (p < class_cutoff) {
      p = class_cutoff - p;
    } else {
      p = p - class_cutoff;
    }

    if (p < max_pred) {
      continue;
    }

    if (token.contains(not_model_string)) {  // predicted to be in the NOT class
      continue;
    }

    max_pred = p;
    pred = token;
  }

  if (0 == pred.length()) {
    output << class_cutoff << " NO_PREDICTION\n";
  } else {
    output << (class_cutoff + max_pred) << ' ' << pred << '\n';
  }

  return 1;
}

static int
model_average_many_class_pairs2(const const_IWSubstring& buffer,
                                IWString_and_File_Descriptor& output)
{
  const_IWSubstring token;
  int i = 0;

  if (!buffer.nextword(token, i, input_separator)) {
    cerr << "Cannot read identifier\n";
    return 0;
  }

  output << token << ' ';

  float max_pred = -std::numeric_limits<float>::max();
  IWString pred;

  while (buffer.nextword(token, i, input_separator)) {
    float p;
    if (!token.numeric_value(p)) {
      cerr << "Invalid floating point input '" << token << "'\n";
      return 0;
    }

    if (!buffer.nextword(token, i, input_separator)) {
      cerr << "Huh, cannot extract class after numeric prediction\n";
      return 0;
    }

    if (p < class_cutoff) {
      p = class_cutoff - p;
    } else {
      p = p - class_cutoff;
    }

    if (p < max_pred) {
      continue;
    }

    max_pred = p;

    pred = token;
  }

  output << (class_cutoff + max_pred) << ' ' << pred << '\n';

  return 1;
}

static int
model_average_many_class_pairs(iwstring_data_source& input,
                               IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record of input\n";
    return 0;
  }

  output << "ID SVMMC PRED\n";

  while (input.next_record(buffer)) {
    if (!model_average_many_class_pairs2(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
model_average_many_class_one_vs_all(iwstring_data_source& input,
                                    IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record of input\n";
    return 0;
  }

  output << "ID SVMMC PRED\n";

  while (input.next_record(buffer)) {
    if (!model_average_many_class_one_vs_all2(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
model_average_many_class_pairs(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  return model_average_many_class_pairs(input, output);
}

static int
model_average_many_class_one_vs_all(const char* fname,
                                    IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  return model_average_many_class_one_vs_all(input, output);
}

static int
model_average(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vc:a:C:fndr");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('C')) {
    if (!cl.value('C', class_cutoff)) {
      cerr << "The class cutoff value (-C) must be a valid floating point number\n";
      usage(3);
    }

    classification_model = 1;

    if (verbose) {
      cerr << "Class cutoff at " << class_cutoff << '\n';
    }

    if (cl.option_present('d')) {
      discern_class_labels_from_input_file = 1;

      if (verbose) {
        cerr << "Will discern class labels from input file\n";
      }
    }
  }

  if (cl.option_present('f')) {
    write_composite_prediction_first = 1;

    if (verbose) {
      cerr << "Will write the composite prediction first\n";
    }
  }

  if (cl.option_present('n')) {
    just_write_composite_prediction = 1;
    write_composite_prediction_first = 1;

    if (verbose) {
      cerr << "Will only write the composite prediction\n";
    }
  }

  if (cl.option_present('r')) {
    include_range_average_std = 1;

    if (verbose) {
      cerr << "Will include the range/average and standard deviation of the component "
              "model values\n";
    }
  }

  // I have set this up to do multiple types of combinations,
  // but not so easy to implement below, so for now, this really
  // does not work as written here, just one combination type is
  // possible.

  int many_class_model = 0;

  int model_average_type = MODEL_AVERAGE_TYPE_AVERAGE;  // default

  if (cl.option_present('a')) {
    model_average_type = 0;
    int i = 0;
    const_IWSubstring a;
    while (cl.value('a', a, i++)) {
      if ("ave" == a) {
        model_average_type |= MODEL_AVERAGE_TYPE_AVERAGE;
        if (verbose) {
          cerr << "Result will be model average\n";
        }
      } else if ("max" == a) {
        model_average_type |= MODEL_AVERAGE_TYPE_MAXVAL;
        if (verbose) {
          cerr << "Result will be max value\n";
        }
      } else if ("vote" == a) {
        model_average_type |= MODEL_AVERAGE_TYPE_VOTING;
        if (verbose) {
          cerr << "Result will be class that occurs most often\n";
        }
      } else if ("sum" == a) {
        model_average_type |= MODEL_AVERAGE_TYPE_SUM;
        if (verbose) {
          cerr << "Result will be sum of all model values\n";
        }
      } else if ("pairs" == a) {
        model_average_type = MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_PAIRS;
        many_class_model = 1;
        if (verbose) {
          cerr << "Many class model, pair-wise models built\n";
        }
      } else if ("ovsall" == a) {
        model_average_type = MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_ONE_VS_ALL;
        many_class_model = 1;
        if (verbose) {
          cerr << "Many class model, one vs all models built\n";
        }
      } else {
        cerr << "Unrecognised -a qualifier '" << a << "'\n";
        usage(3);
      }
    }
  }

  if (!many_class_model) {
    resizable_array<int> col;
    resizable_array<double> w;

    nmodels = gather_columns_and_weights(cl, col, w);

    if (nmodels < 2) {
      cerr << "Must specify two or model results via the -c option\n";
      usage(4);
    }

    model = new Model[nmodels];

    for (int i = 0; i < nmodels; i++) {
      model[i].set_column_and_weight(col[i], w[i]);

      col_to_model[model[i].column()] = i;

      if (verbose > 1) {
        cerr << "Data in column " << (model[i].column() + 1) << " is model " << i << '\n';
      }
    }

    if (verbose > 1) {
      for (int i = 0; i < nmodels; i++) {
        const Model& mi = model[i];

        cerr << "Model in column " << (mi.column() + 1) << " weight " << mi.weight()
             << '\n';
      }
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, don't know how to handle multiple files\n";
    return 8;
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    int rc;

    //  cerr << "mat " << model_average_type << '\n';

    if (MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_PAIRS == model_average_type) {
      rc = model_average_many_class_pairs(cl[i], output);
    } else if (MODEL_AVERAGE_MANY_CLASS_CLASSIFICATION_ONE_VS_ALL == model_average_type) {
      rc = model_average_many_class_one_vs_all(cl[i], output);
    } else if (MODEL_AVERAGE_TYPE_AVERAGE == model_average_type) {
      Average a;
      rc = model_average(cl[i], a, output);
    } else if (MODEL_AVERAGE_TYPE_MAXVAL == model_average_type) {
      MaxVal m;
      rc = model_average(cl[i], m, output);
    } else if (MODEL_AVERAGE_TYPE_VOTING == model_average_type) {
      Voting v;
      rc = model_average(cl[i], v, output);
    } else if (MODEL_AVERAGE_TYPE_SUM == model_average_type) {
      Sum s;
      rc = model_average(cl[i], s, output);
    } else {
      cerr << "Not sure what to do with " << model_average_type << '\n';
      return 8;
    }

    if (0 == rc) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Made predictions for " << score_acc.n() << " items\n";
    if (score_acc.n() > 1) {
      cerr << "Values between " << score_acc.minval() << " and " << score_acc.maxval()
           << " ave " << static_cast<float>(score_acc.average()) << '\n';
    }
  }

  if (NULL != model) {
    delete[] model;
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = model_average(argc, argv);

  return rc;
}
