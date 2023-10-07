/*
  We have a file with arbitrary class labels. We want to translate
  these into known class labels for svm-lite

  Forward translation will be to translate from user input class labels
  to +1 and -1

  Backward translation will be from the numeric range back to initial class
  labels

  Feb 2008. Want several different kinds of translations
*/

#include <stdlib.h>

#include <iostream>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int column_for_class_label = 1;

static resizable_array_p<IWString> label_to_assign;

static int any_number_of_class_labels = 0;

static int labels_assigned = 0;

static int header_records_to_skip = 0;

static float upper_class_cutoff = 0.0;
static float lower_class_cutoff = 0.0;

/*
  There is ambiguity about how cutoffs should be treated when there
  is rescaling via the -R option. In the current implementation, if
  the -R option is used, then the class cutoff is what should be
  used directly. What a mess....
*/

static float class_cutoff_from_command_line = 0.0;

static int rarest_class_becomes_class1 = 0;

/*
  Various kinds of output
*/

static int write_value_then_class = 1;
static int write_class_and_score = 0;
static int write_class_and_two_scores = 0;
static int rescale_to_01 = 0;

static int rescale_via_user_specified_range = 0;
static float user_specified_min;
static float user_specified_max;

static IWString negative_class("-1");
static IWString positive_class("1");

static IWString ambiguous_class("*");

static int assigned_positive_class = 0;
static int assigned_negative_class = 0;
static int ambiguous_assignments = 0;
static int assigned_user_specified_class = 0;

static int write_ambiguous_predictions = 1;

static int truncated_to_lower_end = 0;
static int truncated_to_upper_end = 0;

static int flush_after_every_molecule = 0;

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
  cerr << "Handles class label translations to and from svm-lite\n";
  cerr << "Forward translation options marked with >, backward with <\n";
  cerr << " -L lab1,lab2,... >  preferred class labels (default '-1' and '1')\n";
  cerr << " -a                  allow any number of classes\n";
  cerr << " -C <fname>       >  write the class translation table to <fname>\n";
  cerr << " -s <number>      <> header records to skip\n";
  cerr << " -c <col>         <  lables to translate are in column <col>\n";
  cerr << " -U <fname>       <  reverse translation using file created by -C option\n";
  cerr << " -n <float>       <  cutoff value for converting svm-lite numbers to classes\n";
  cerr << " -n <float,float> <  values within range not assigned\n";
  cerr << " -x                  suppress writing ambiguous predictions\n";
  cerr << " -r 01               scale scores to the range 0-1\n";
  cerr << " -O ...              various kinds of output\n";
  cerr << " -O SC               write score then class (default)\n";
  cerr << " -O CS               write class then score\n";
  cerr << " -O CSS              write class then score for each class\n";
  cerr << " -R min-max          specify min and max for scaling\n";
  cerr << " -S <fname>       <  file name for range to string specification(s)\n";
  cerr << " -f               >  make sure smallest class is class '1' \n";
  cerr << " -Y ...              more options, enter '-Y help' for info\n";
  cerr << " -v               <> verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  June 2010. Folks want to be able to call certain ranges different things.
  So, values in the range 0.4 to 0.6 might be called "intedeterminate"
*/

class Range_to_String
{
 private:
  float _min;
  float _max;
  IWString _label;

 public:
  Range_to_String();

  int build(const const_IWSubstring&);

  int do_any_rescaling_needed();

  int matches(float, IWString&) const;
};

Range_to_String::Range_to_String()
{
  _min = 1.0f;  //
  _max = 0.0f;

  return;
}

int
Range_to_String::build(const const_IWSubstring& buffer)
{
  int i = 0;
  const_IWSubstring token;

  if (buffer.nwords() < 3) {
    cerr << "Range_to_String::build:input must have at least 3 tokens\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (!token.numeric_value(_min)) {
    cerr << "Range_to_String::build:invalid min specification '" << buffer << "'\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (!token.numeric_value(_max)) {
    cerr << "Range_to_String::build:invalid max specification '" << buffer << "'\n";
    return 0;
  }

  if (_min > _max) {
    cerr << "Range_to_String::build:inconsistent specification '" << buffer << "'\n";
    return 0;
  }

  buffer.nextword(_label, i);

  return 1;
}

/*
  The values read from the file will be in the range people are expecting the output,
  so we need to scale back to the -1,1 range
*/

int
Range_to_String::do_any_rescaling_needed()
{
  if (rescale_to_01) {
    _min = 2.0f * _min - 1.0f;
    _max = 2.0f * _max - 1.0f;
  } else if (rescale_via_user_specified_range) {
    cerr << "Range_to_String::do_any_rescaling_needed:not implemented for user specified "
            "range, see Ian\n";
    return 0;
  }

  return 1;
}

int
Range_to_String::matches(float v, IWString& s) const
{
  if (v < _min) {
    return 0;
  }

  if (v > _max) {
    return 0;
  }

  s = _label;

  return 1;
}

class Class_Label_Translation_Table
{
 private:
  IW_STL_Hash_Map_String _xref;
  IW_STL_Hash_Map_int _class_count;

  resizable_array_p<Range_to_String> _range_to_string;

  //  private functions

  int _read_buffer(const const_IWSubstring& buffer);

 public:
  int store_label(const IWString&);
  int
  do_translate(const IWString&, IWString&) const;

  int items_in_xref() const
  {
    return _xref.size();
  }

  int report_counts(std::ostream&) const;

  //  If things come in as -1 and +1, then we don't need to do anything
  //  But we may need to swap the mapping, so non-const

  int translation_needed();

  int skip_and_echo_header_records_reverse(iwstring_data_source& input,
                                       const Class_Label_Translation_Table& ct,
                                       IWString_and_File_Descriptor& output) const;
  int echo_header_record_reverse(const const_IWSubstring& buffer,
                             IWString_and_File_Descriptor& output) const;

  int do_read(const char*);
  int do_read(iwstring_data_source&);
  int do_write(const char*) const;
  int do_write(IWString_and_File_Descriptor&) const;

  int read_range_to_string_specifications(const char*);
  int read_range_to_string_specifications(iwstring_data_source&);

  int range_to_string_specifications() const {
    return _range_to_string.number_elements();
  }

  void do_any_range_to_string_rescaling_needed();

  int matches_a_range_to_string_specifier(float v, IWString& translated_form) const;

  int number_classses_identified() const {
    return _xref.size();
  }

  int sort_class_labels(const resizable_array_p<IWString>&);
  int rarest_class_becomes_class1();
};

int
Class_Label_Translation_Table::do_translate(const IWString& sfrom, IWString& sto) const
{
  IW_STL_Hash_Map_String::const_iterator f = _xref.find(sfrom);

  if (f != _xref.end()) {
    sto = (*f).second;
    return 1;
  }

  return 0;
}

int
Class_Label_Translation_Table::store_label(const IWString& s)
{
  IW_STL_Hash_Map_String::const_iterator f = _xref.find(s);

  if (f != _xref.end())  // seen this before
  {
    _class_count[s] += 1;
    return 1;
  }

  if (labels_assigned < label_to_assign.number_elements()) {
    _xref[s] = *(label_to_assign[labels_assigned]);
  } else if (any_number_of_class_labels) {
    IWString* tmp = new IWString;
    tmp->append_number(any_number_of_class_labels);
    _xref[s] = *tmp;
    label_to_assign.add(tmp);
    any_number_of_class_labels++;
  } else {
    cerr << "Class_Label_Translation_Table::store_label:all labels assigned, cannot "
            "classify '"
         << s << "'\n";
    return 0;
  }

  labels_assigned++;

  _class_count[s] = 1;

  return 2;
}

int
Class_Label_Translation_Table::_read_buffer(const const_IWSubstring& buffer)
{
  if (2 != buffer.nwords()) {
    cerr << "Class_Label_Translation_Table::_read_buffer:must be two tokens '" << buffer
         << "'\n";
    return 0;
  }

  IWString t1, t2;

  if (!buffer.split(t1, ' ', t2))  // how could this happen?
  {
    cerr << "Class_Label_Translation_Table::_read_buffer:cannot split\n";
    return 0;
  }

  _xref[t2] = t1;

  return 1;
}

int
Class_Label_Translation_Table::do_read(iwstring_data_source& input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (!_read_buffer(buffer)) {
      cerr << "Class_Label_Translation_Table::do_read:invalid input '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Class_Label_Translation_Table::do_read(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Class_Label_Translation_Table::do_read:cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return do_read(input);
}

int
Class_Label_Translation_Table::read_range_to_string_specifications(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Class_Label_Translation_Table::read_range_to_string_specifications:cannot "
            "open '"
         << fname << "'\n";
    return 0;
  }

  return read_range_to_string_specifications(input);
}

int
Class_Label_Translation_Table::read_range_to_string_specifications(
    iwstring_data_source& input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (0 == buffer.length() || buffer.starts_with('#')) {
      continue;
    }

    Range_to_String* r = new Range_to_String;

    if (!r->build(buffer)) {
      cerr << "Class_Label_Translation_Table::read_range_to_string_specifications:"
              "invalid range to string '"
           << buffer << "'\n";
      return 0;
    }

    _range_to_string.add(r);
  }

  // cerr << "Build " << _range_to_string.number_elements() << " range to string
  // specifications\n";

  return 1;
}

void
Class_Label_Translation_Table::do_any_range_to_string_rescaling_needed()
{
  for (int i = 0; i < _range_to_string.number_elements(); i++) {
    _range_to_string[i]->do_any_rescaling_needed();
  }

  return;
}

int
Class_Label_Translation_Table::matches_a_range_to_string_specifier(
    float v, IWString& translated_form) const
{
  for (int i = 0; i < _range_to_string.number_elements(); i++) {
    if (_range_to_string[i]->matches(v, translated_form)) {
      return 1;
    }
  }

  return 0;
}

int
Class_Label_Translation_Table::do_write(const char* fname) const
{
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Class_Label_Translation_Table::do_write:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(output);
}

int
Class_Label_Translation_Table::do_write(IWString_and_File_Descriptor& output) const
{
  output << "# Class_Label_Translation_Table\n";
  output << "#\n";

  for (IW_STL_Hash_Map_String::const_iterator i = _xref.begin(); i != _xref.end(); i++) {
    output << (*i).first << ' ' << (*i).second << '\n';
  }

  return 1;
}

int
Class_Label_Translation_Table::report_counts(std::ostream& os) const
{
  for (IW_STL_Hash_Map_int::const_iterator i = _class_count.begin();
       i != _class_count.end(); i++) {
    os << "Found " << (*i).second << " instances of class '" << (*i).first << "'\n";
  }

  return 1;
}

class String_Comparator
{
 private:
 public:
  int
  operator()(const IWString* s1, const IWString* s2) const;
};

int
String_Comparator::operator()(const IWString* s1, const IWString* s2) const
{
  int l1 = s1->length();
  int l2 = s2->length();

  if (l1 < l2) {
    return -1;
  }

  if (l1 > l2) {
    return 1;
  }

  return s1->strcmp(*s2);
}

int
Class_Label_Translation_Table::sort_class_labels(
    const resizable_array_p<IWString>& label_to_assign)
{
  if (_xref.size() < 2) {
    cerr << "Class_Label_Translation_Table::sort_class_labels:only " << _xref.size()
         << " classes\n";
    return 0;
  }

  if (_xref.size() != static_cast<unsigned int>(label_to_assign.number_elements())) {
    cerr << "Class_Label_Translation_Table::sort_class_labels:found " << _xref.size()
         << " tokens, but " << label_to_assign.number_elements()
         << " classes to assign\n";
    return 0;
  }

  resizable_array_p<IWString> tmp;

  for (IW_STL_Hash_Map_String::const_iterator i = _xref.begin(); i != _xref.end(); i++) {
    tmp.add(new IWString((*i).first));
  }

  String_Comparator sc;

  tmp.iwqsort(sc);

  _xref.clear();

  for (int i = 0; i < tmp.number_elements(); i++) {
    const IWString* s = tmp[i];

    _xref[*s] = *(label_to_assign[i]);

    if (verbose > 1) {
      cerr << "Class '" << (*s) << "' translated to '" << *(label_to_assign[i]) << "'\n";
    }
  }

  return 1;
}

class Label_and_Count
{
 private:
  const IWString _label;
  const int _count;

 public:
  Label_and_Count(const IWString& s, int c) : _label(s), _count(c){};

  const IWString&
  label() const
  {
    return _label;
  }

  int
  count() const
  {
    return _count;
  }
};

class Label_and_Count_Comparator
{
 private:
 public:
  int
  operator()(const Label_and_Count*, const Label_and_Count*) const;
};

int
Label_and_Count_Comparator::operator()(const Label_and_Count* lc1,
                                       const Label_and_Count* lc2) const
{
  int c1 = lc1->count();
  int c2 = lc2->count();

  if (c1 > c2) {
    return -1;
  }

  if (c1 < c2) {
    return 1;
  }

  return 0;
}

int
Class_Label_Translation_Table::rarest_class_becomes_class1()
{
  assert(2 == _xref.size());

  resizable_array_p<Label_and_Count> lc;

  for (IW_STL_Hash_Map_int::const_iterator i = _class_count.begin();
       i != _class_count.end(); ++i) {
    const IWString& l = (*i).first;
    int c = (*i).second;

    Label_and_Count* t = new Label_and_Count(l, c);
    lc.add(t);
  }

  Label_and_Count_Comparator lcc;

  lc.iwqsort(lcc);

  _xref.clear();

  for (int i = 0; i < lc.number_elements(); i++) {
    const Label_and_Count* lci = lc[i];

    cerr << lci->count() << " occurrences of '" << lci->label() << "'\n";

    _xref[lci->label()] = *(label_to_assign[i]);
  }

  return 1;
}

/*
  Do we need to translate the labels from something else to get to '1' and '-1'

  Tried to do something to account for a label of '+1', but it doesn't
  work because we end up writing '1' to the translation file. Let's not
  worry about that now...
*/

int
Class_Label_Translation_Table::translation_needed()
{
  IW_STL_Hash_Map_String::const_iterator fplus = _xref.find("1");
  if (fplus == _xref.end()) {
    fplus = _xref.find("+1");

    if (fplus == _xref.end()) {
      return 1;  // needs translating
    }

    IWString tmp = (*fplus).second;
    _xref.erase("+1");
    _xref["1"] = tmp;

    fplus = _xref.find("1");
  }

  IW_STL_Hash_Map_String::const_iterator fminus = _xref.find("-1");
  if (fminus == _xref.end()) {
    return 1;
  }

  if ((*fplus).second == "1" && (*fminus).second == "-1") {
    return 0;
  }

  // Both '1' and '-1' are present, just need to change the translation

  _xref["1"] = "1";
  _xref["-1"] = "-1";

  return 0;
}

#ifdef NOT_USED_ASDASDASD
static float
rescaled_if_necessary(float v, float class_cutoff)
{
  if (!rescale_to_01) {
    return v;
  }

  if (v < class_cutoff) {
    return 1.0 - v;
  } else {
    return v;
  }
}
#endif

/*
  Put back into the -1 to +1 range
*/

static void
do_rescale_via_user_specified_range(float& v, float user_specified_min,
                                    float user_specified_max)
{
  if (v <= user_specified_min) {
    v = static_cast<float>(-1.0);
    return;
  }

  if (v >= user_specified_max) {
    v = static_cast<float>(1.0);
    return;
  }

  // cerr << "Rescaling " << v;
  v = -1.0 + 2.0 * (v - user_specified_min) / (user_specified_max - user_specified_min);
  // cerr << " to " << v << endl;

  return;
}

static int
do_reverse_translation_record(const const_IWSubstring& input_buffer,
                              const Class_Label_Translation_Table& ct,
                              IWString_and_File_Descriptor& output)
{
  IWString output_buffer;  // if an ambiguous prediction, we may not output anything
  output_buffer.resize(input_buffer.length() + 10);

  int i = 0;
  IWString token;

  for (int col = 0; input_buffer.nextword(token, i); col++) {
    if (col > 0) {
      output_buffer << ' ';
    }

    if (col != column_for_class_label) {
      output_buffer << token;
      continue;
    }

    float v;

    if (!token.numeric_value(v)) {
      cerr << "do_reverse_translation_record:invalid numeric '" << token << "'\n";
      return 0;
    }

    if (rescale_via_user_specified_range) {
      do_rescale_via_user_specified_range(v, user_specified_min, user_specified_max);
    }

    //  the only output format that tolerates out of range values is the default

    if (write_value_then_class) {
      ;
    } else if (v < static_cast<float>(-1.0)) {
      v = static_cast<float>(-1.0);
      truncated_to_lower_end++;
    } else if (v > static_cast<float>(1.0)) {
      v = static_cast<float>(1.0);
      truncated_to_upper_end++;
    }

    int class_this_item_assigned = 0;

    IWString translated_form;
    if (ct.items_in_xref() > 2) {
    } else if (ct.matches_a_range_to_string_specifier(v, translated_form)) {
      assigned_user_specified_class++;
    } else if (v < lower_class_cutoff) {
      assigned_negative_class++;
      ct.do_translate(negative_class, translated_form);
      class_this_item_assigned = -1;
    } else if (v >= upper_class_cutoff)  // not equality sign, arbitrary choice
    {
      assigned_positive_class++;
      ct.do_translate(positive_class, translated_form);
      class_this_item_assigned = 1;
    } else {
      ambiguous_assignments++;
      cerr << "AMbiguous " << v << endl;
      if (!write_ambiguous_predictions) {
        return 1;
      }
      translated_form = ambiguous_class;
      //    ct.do_translate(ambiguous_class, translated_form);
    }

    //  Now do any rescaling for output

#ifdef DEBUG_UNSCALING
    cerr << "Before rescaling v = " << v << endl;
#endif

    if (rescale_to_01) {
      //    cerr << "Translating " << v;
      v = (v + 1.0) / 2.0;
      if (v < 0.0) {
        v = 0.0;
        token = '0';
      } else if (v > 1.0) {
        v = 1.0;
        token = '1';
      } else {
        token.resize_keep_storage(0);
        token << v;
      }
      //    cerr << " to " << v << " cf cutoff " << class_cutoff << endl;
    }

#ifdef DEBUG_UNSCALING
    if (rescale_to_01) {
      cerr << "Rescaled to " << v << endl;
    }
#endif

    if (write_value_then_class) {
      if (rescale_via_user_specified_range) {  // must write new numeric form
        output_buffer << v << ' ' << translated_form;
      } else {
        output_buffer << token << ' ' << translated_form;
      }
    } else if (write_class_and_score) {
      output_buffer << translated_form << ' ';
      if (!rescale_to_01) {
        output_buffer << v;
      } else if (class_this_item_assigned < 0) {
        output_buffer << static_cast<float>(1.0 - v);
      } else {
        output_buffer << v;
      }
    } else if (write_class_and_two_scores) {
#ifdef DEBUG_UNSCALING
      cerr << "Value is " << v << " cutoff " << lower_class_cutoff << endl;
#endif
      output_buffer << translated_form << ' ';
      output_buffer << static_cast<float>(1.0 - v) << ' ' << v;
    } else {
      output_buffer << translated_form << ' ' << token;
    }
  }

  output << output_buffer << '\n';

  return 1;
}

static int
echo_input_to_output(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  off_t s = input.file_size();

  return input.echo(output, s);
}

static int
echo_input_to_output(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return echo_input_to_output(input, output);
}

// After each molecule is processed, we might need to flush
// `output`, or write if it gets too full.
void
MaybeFlush(IWString_and_File_Descriptor& output)
{
  if (flush_after_every_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(32768);
  }
}

int
Class_Label_Translation_Table::echo_header_record_reverse(
    const const_IWSubstring& buffer, IWString_and_File_Descriptor& output) const
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); col++) {
    if (col > 0) {
      output << ' ';
    }

    if (col != column_for_class_label) {
      output << token;
      continue;
    }

    if (write_value_then_class) {
      output << token << ".SCORE " << token;
    } else if (write_class_and_score) {
      output << token << ' ' << token << ".SCORE";
    } else if (write_class_and_two_scores) {
      //    output << token << ' ' << token << ".SCORE1 " << token << ".SCORE2";
      IWString tmp;
      do_translate(negative_class, tmp);
      output << token << ' ' << token << '.' << tmp << ' ' << token << '.';
      do_translate(positive_class, tmp);
      output << tmp;
    } else {
      output << token << ' ' << token << ".SCORE ";
    }
  }

  output << '\n';

  MaybeFlush(output);

  return 1;
}

/*
  This is actually wrong, it doesn't work properly unless the number of
  header records to skip is 1
*/

int
Class_Label_Translation_Table::skip_and_echo_header_records_reverse(
    iwstring_data_source& input, const Class_Label_Translation_Table& ct,
    IWString_and_File_Descriptor& output) const
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++) {
    if (!input.next_record(buffer)) {
      cerr << "Premature end of file reading header records\n";
      return 0;
    }

    echo_header_record_reverse(buffer, output);
  }

  return 1;
}

static int
do_reverse_translation(iwstring_data_source& input,
                       const Class_Label_Translation_Table& ct,
                       IWString_and_File_Descriptor& output)
{
  if (!ct.skip_and_echo_header_records_reverse(input, ct, output)) {
    return 0;
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!do_reverse_translation_record(buffer, ct, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    MaybeFlush(output);
  }

  return 1;
}

static int
do_forward_translation_record(const const_IWSubstring& buffer,
                              const Class_Label_Translation_Table& ct,
                              IWString_and_File_Descriptor& output)
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); col++) {
    if (col > 0) {
      output << ' ';
    }

    if (col != column_for_class_label) {
      output << token;
    } else {
      IWString translated_form;
      if (ct.do_translate(token, translated_form)) {
        output << translated_form;
      } else {
        cerr << "No class translation available for '" << token << "'\n";
        return 0;
      }
    }
  }

  output << '\n';

  return 1;
}

static int
echo_header_record_forward(const const_IWSubstring& buffer,
                           IWString_and_File_Descriptor& output)
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); col++) {
    if (col > 0) {
      output << ' ';
    }

    if (col == column_for_class_label) {
      output << token;
    } else {
      output << token;
    }
  }

  output << '\n';

  return 1;
}

static int
skip_and_echo_header_records_forward(iwstring_data_source& input,
                                     IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++) {
    if (!input.next_record(buffer)) {
      cerr << "Premature end of file reading header records\n";
      return 0;
    }

    echo_header_record_forward(buffer, output);
  }

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
do_forward_translation(iwstring_data_source& input,
                       const Class_Label_Translation_Table& ct,
                       IWString_and_File_Descriptor& output)
{
  if (!skip_and_echo_header_records_forward(input, output)) {
    return 0;
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!do_forward_translation_record(buffer, ct, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    MaybeFlush(output);
  }

  return 1;
}

static int
do_reverse_translation(const char* fname, const Class_Label_Translation_Table& ct,
                       IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return do_reverse_translation(input, ct, output);
}

static int
do_forward_translation(const char* fname, const Class_Label_Translation_Table& ct,
                       IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return do_forward_translation(input, ct, output);
}

static int
identify_labels_present_in_file_record(const const_IWSubstring& buffer,
                                       Class_Label_Translation_Table& ct)
{
  IWString id;

  if (!buffer.word(column_for_class_label, id)) {
    return 0;
  }

  return ct.store_label(id);
}

static int
identify_labels_present_in_file(iwstring_data_source& input,
                                Class_Label_Translation_Table& ct)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++) {
    if (!input.next_record(buffer)) {
      cerr << "Premature end of file reading header records\n";
      return 0;
    }
  }

  while (input.next_record(buffer)) {
    if (!identify_labels_present_in_file_record(buffer, ct)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return ct.number_classses_identified();
}

static int
identify_labels_present_in_file(const char* fname, Class_Label_Translation_Table& ct)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return identify_labels_present_in_file(input, ct);
}

static int
do_forward_translation(Command_Line& cl, IWString_and_File_Descriptor& output)
{
  Class_Label_Translation_Table ct;

  if (cl.option_present('L')) {
    int i = 0;
    const_IWSubstring l;
    while (cl.value('L', l, i++)) {
      int j = 0;
      const_IWSubstring token;
      while (l.nextword(token, j, ',')) {
        label_to_assign.add(new IWString(token));
      }
    }

    if (verbose) {
      cerr << "Assigned " << label_to_assign.number_elements() << " labels\n";
    }
  } else if (any_number_of_class_labels) {
    ;
  } else {
    label_to_assign.add(new IWString("-1"));
    label_to_assign.add(new IWString("1"));
  }

  if (cl.option_present('f')) {
    rarest_class_becomes_class1 = 1;

    if (verbose) {
      cerr << "The rarest class will become class '1'\n";
    }
  }

  if (!identify_labels_present_in_file(cl[0], ct)) {
    cerr << "Cannot determine labels present in file '" << cl[0] << "'\n";
    return 4;
  }

  if (0 == ct.number_classses_identified()) {
    cerr << "No classes\n";
    return 0;
  }

  if (verbose) {
    ct.report_counts(cerr);
  }

  const char* c = cl.option_value('C');  // the output file we will create

  // cerr << "Found " << ct.number_classses_identified() << " Rare " <<
  // rarest_class_becomes_class1 << endl;

  if (1 == ct.number_classses_identified()) {
    cerr << "Warning, only one class in the input\n";
  } else if (!ct.translation_needed())  // already -1 and 1
  {
    if (verbose) {
      cerr << "No translation needed\n";
    }

    ct.do_write(c);

    return echo_input_to_output(cl[0], output);
  } else if (rarest_class_becomes_class1) {
    ct.rarest_class_becomes_class1();
  } else {
    ct.sort_class_labels(label_to_assign);
  }

  ct.do_write(c);

  return do_forward_translation(cl[0], ct, output);
}

static int
do_reverse_translation(Command_Line& cl, IWString_and_File_Descriptor& output)
{
  Class_Label_Translation_Table ct;

  const char* u = cl.option_value('U');
  if (!ct.do_read(u)) {
    cerr << "Cannot read cross reference file '" << u << "'\n";
    return 4;
  }

  if (verbose) {
    cerr << "Read " << ct.number_classses_identified() << " class labels from '" << u
         << "'\n";
  }

  if (ct.number_classses_identified() > 2) {
    cerr << "Sorry, don't know how to handle " << ct.number_classses_identified()
         << " classes\n";
    return 4;
  }

  if (0 == ct.number_classses_identified()) {
    cerr << "No classes\n";
    return 8;
  }

  IWString tmp;
  if (ct.do_translate(positive_class, tmp) && ct.do_translate(negative_class, tmp)) {
    ;
  } else {
    cerr << "Cannot use translation table, no translation for '" << positive_class
         << "' and/or '" << negative_class << "'\n";
    return 0;
  }

  if (cl.option_present('S')) {
    const char* s = cl.option_value('S');

    if (!ct.read_range_to_string_specifications(s)) {
      cerr << "Cannot read range to string specifications from '" << s << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Read " << ct.range_to_string_specifications()
           << " range to label specifications\n";
    }

    ct.do_any_range_to_string_rescaling_needed();  // after we have initialised all the
                                                   // various scaling things
  }

  int rc = do_reverse_translation(cl[0], ct, output);

  if (verbose) {
    IWString tmp;
    ct.do_translate(positive_class, tmp);
    cerr << "assigned " << assigned_positive_class << " to the positive class '" << tmp
         << "'\n";
    ct.do_translate(negative_class, tmp);
    cerr << "assigned " << assigned_negative_class << " to the negative class '" << tmp
         << "'\n";
    if (ambiguous_assignments) {
      cerr << ambiguous_assignments << " ambiguous not assigned\n";
    }
    if (assigned_user_specified_class) {
      cerr << assigned_user_specified_class << " assigned to user specified classes\n";
    }
  }

  return rc;
}

static void
DisplayDashYOptions(std::ostream& output)
{
  output << " -Y flush         flush after each molecule processed\n";

  ::exit(0);
}

static int
class_label_translation(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vc:C:U:L:s:n:abr:O:R:xfS:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    if (!cl.value('c', column_for_class_label) || column_for_class_label < 1) {
      cerr << "The column number option (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Class labels found in column " << column_for_class_label << '\n';
    }

    column_for_class_label--;
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', header_records_to_skip) || header_records_to_skip < 0) {
      cerr << "The header records to skip option (-s) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will skip " << header_records_to_skip << " header records\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, don't know how to handle multiple input files\n";
    usage(4);
  }

  if (cl.option_present('L') && cl.option_present('U')) {
    cerr << "The use existing cross reference (-U) and -L options don't make sense "
            "together\n";
    usage(4);
  }

  if (cl.option_present('L') && cl.option_present('a')) {
    cerr << "The -L and -a options are mutually inconsistent\n";
    usage(4);
  }

  if (cl.option_present('a')) {
    any_number_of_class_labels = 1;

    if (verbose) {
      cerr << "Will assign an arbitrary number of class labels\n";
    }
  }

  if (cl.option_present('r')) {
    const_IWSubstring r = cl.string_value('r');

    if ("01" == r) {
      rescale_to_01 = 1;
      if (verbose) {
        cerr << "Will rescale to 01 range\n";
      }
    } else if ("11" == r) {  // default behaviour
      ;
    } else {
      cerr << "The only rescaling I can do is '01', don't recognise '" << r << "'\n";
      usage(3);
    }
  }

  if (cl.option_present('R')) {
    const_IWSubstring r = cl.string_value('R');

    const_IWSubstring mi, ma;
    if (!r.split(mi, ':', ma) || 0 == mi.length() || 0 == ma.length()) {
      cerr << "The response scaling option (-R) must have a valid range like "
              "'float:float'\n";
      usage(4);
    }

    if (!mi.numeric_value(user_specified_min)) {
      cerr << "The minimum for the user specified raw score scaling (-R) must be a valid "
              "numeric '"
           << r << "'\n";
      return 6;
    }

    if (!ma.numeric_value(user_specified_max)) {
      cerr << "The maximum for the user specified raw score scaling (-R) must be a valid "
              "numeric '"
           << r << "'\n";
      return 6;
    }

    if (user_specified_min < user_specified_max) {
      ;
    } else {
      cerr << "The user specified raw score scaling min must be less than the max '" << r
           << "'\n";
      usage(5);
    }

    rescale_via_user_specified_range = 1;

    if (verbose) {
      cerr << "Raw scores scaled to the range " << user_specified_min << " to "
           << user_specified_max << endl;
    }
  }

  if (cl.option_present('O')) {
    int i = 0;
    const_IWSubstring o;
    while (cl.value('O', o, i++)) {
      if ("SC" == o) {
        write_value_then_class = 1;
        write_class_and_score = 0;
        write_class_and_two_scores = 0;
        if (verbose) {
          cerr << "Will write the numeric value then the class\n";
        }
      } else if ("CS" == o) {
        write_value_then_class = 0;
        write_class_and_score = 1;
        write_class_and_two_scores = 0;
        if (verbose) {
          cerr << "Will write class and then score\n";
        }

        //      rescale_to_01 = 1;
      } else if ("CSS" == o) {
        write_value_then_class = 0;
        write_class_and_score = 0;
        write_class_and_two_scores = 1;
        if (verbose) {
          cerr << "Will write class followed by two scores\n";
        }

        rescale_to_01 = 1;  // always
      } else {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        return 6;
      }
    }
  }

  if (cl.option_present('n')) {
    const_IWSubstring n = cl.string_value('n');

    if (n.contains(',')) {
      const_IWSubstring n1, n2;
      if (!n.split(n1, ',', n2) || 0 == n1.length() || 0 == n2.length()) {
        cerr << "The class cutoff range (-n f1,f2) must contain two valid floating point "
                "numbers '"
             << n << "' invalid\n";
        usage(2);
      }

      if (!n1.numeric_value(lower_class_cutoff) ||
          !n2.numeric_value(upper_class_cutoff) ||
          lower_class_cutoff > upper_class_cutoff) {
        cerr << "Invalid class cutoff specification '" << n << "'\n";
        usage(2);
      }

      if (rescale_to_01) {
        upper_class_cutoff = (upper_class_cutoff * 2.0) - 1.0;
        lower_class_cutoff = (lower_class_cutoff * 2.0) - 1.0;
        if (verbose) {
          cerr << "Rescaling to 01, class cutoffs adjusted to " << lower_class_cutoff
               << ',' << upper_class_cutoff << endl;
        }
      }
    } else {
      if (!cl.value('n', upper_class_cutoff)) {
        cerr << "The class cutoff value (-n) must be a valid number\n";
        usage(4);
      }

      if (verbose) {
        cerr << "Class cutoff set to " << upper_class_cutoff << endl;
      }

      class_cutoff_from_command_line = upper_class_cutoff;

      //    If they are operating in the 01 range, change the cutoff to the -1 to +1 range

      if (rescale_to_01) {
        upper_class_cutoff = (upper_class_cutoff * 2.0) - 1.0;
        lower_class_cutoff = upper_class_cutoff;
        if (verbose) {
          cerr << "Rescaling to 01, class cutoff adjusted to " << upper_class_cutoff
               << endl;
        }
      }
    }
  }

  if (cl.option_present('x')) {
    write_ambiguous_predictions = 0;

    if (verbose) {
      cerr << "Will not write ambiguous predictions\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "flush") {
        flush_after_every_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after every molecule\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  // There really should be two separate programmes to do this, but instead we have two
  // functions

  IWString_and_File_Descriptor output(1);

  int rc;
  if (cl.option_present('C')) {
    rc = do_forward_translation(cl, output);
  } else if (cl.option_present('U')) {
    rc = do_reverse_translation(cl, output);
  } else {
    cerr << "Must either create (-C) or use (-U) a class label transation table\n";
    usage(4);
  }

  output.flush();

  if (verbose) {
    if (truncated_to_upper_end) {
      cerr << truncated_to_upper_end << " predictions truncated to upper limit\n";
    }
    if (truncated_to_lower_end) {
      cerr << truncated_to_lower_end << " predictions truncated to lower limit\n";
    }
  }

  return !rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = class_label_translation(argc, argv);

  return rc;
}
