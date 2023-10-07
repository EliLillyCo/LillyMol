/*
  Splits a classifcation data set into train/test/splits
*/

#include <stdlib.h>

#include <iostream>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "ttitem.h"

using std::cerr;

const char *prog_name = NULL;

static int verbose = 0;

static int number_splits = 1;

static int *items_in_training_set = NULL;

static float fraction_in_training_set = 0.5;

static IWString output_stem("TTC");

static IWString suffix(".dat");

static IWString header_record;

static int write_class_memberships = 1;

static void usage(int rc) {
  // clang-format off
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  // clang-format on
  // clang-format off
  cerr << "Produce train/test splits for a classification problem\n";
  cerr << "Input is an activity file, ID class\n";
  cerr << " -p <pct>       percent of molecules in the training set (default 50)\n";
  cerr << " -n <number>    number of splits to create\n";
  cerr << " -S <stem>      write results to files starting with <stem> (default '" << output_stem << ")\n";
  cerr << " -i             just write ID's and not class membership to output file(s)\n";
  cerr << " -U <suffix>    create files with suffix <suffix> (default '" << suffix << "')\n";
  // cerr << " -s <seed>      random number seed\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

class Class_Members : public resizable_array_p<TTItem> {
private:
  IWString _class_label;

  int _splits_generated;

public:
  Class_Members();

  void set_class_label(const const_IWSubstring &s) { _class_label = s; }
  const IWString &class_label() const { return _class_label; }

  int extra(const const_IWSubstring &);

  int generate_split(double f);

  int write_split(IWString_and_File_Descriptor &train_output,
                  IWString_and_File_Descriptor &test_output) const;

  int items_in_training_set() const;

  int report(std::ostream &) const;
};

Class_Members::Class_Members() {
  _splits_generated = 0;

  return;
}

int Class_Members::extra(const const_IWSubstring &s) {
  TTItem *t = new TTItem(s);

  resizable_array_p<TTItem>::add(t);

  return _number_elements;
}

/*
  We want something to impose a largely random ordering of the items
*/

class Random_Number_Comparitor {
private:
public:
  int operator()(const TTItem *, const TTItem *) const;
};

int Random_Number_Comparitor::operator()(const TTItem *t1,
                                         const TTItem *t2) const {
  double h1 = t1->rnum();
  double h2 = t2->rnum();

  h1 += t1->times_in_training_set();
  h2 += t2->times_in_training_set();

  if (h1 < h2)
    return -1;

  if (h1 > h2)
    return 1;

  return 0;
}

int Class_Members::generate_split(double fraction_in_training_set) {
  if (_splits_generated > 0) {
    Random_Number_Comparitor rnc; // first randomise it
    resizable_array_p<TTItem>::iwqsort(rnc);
  }

  Times_in_Training_Set_Comparitor ttsc;
  resizable_array_p<TTItem>::iwqsort(ttsc);

  int to_be_assigned_training_set =
      static_cast<int>(_number_elements * fraction_in_training_set + 0.499);

  if (0 == to_be_assigned_training_set)
    to_be_assigned_training_set = 1;

  int put_into_training_set = 0;

  for (int i = 0; i < _number_elements; i++) {
    if (put_into_training_set > to_be_assigned_training_set)
      _things[i]->set_in_training_set(0);
    else {
      _things[i]->set_in_training_set(1);
      put_into_training_set++;
    }

    _things[i]->reset_random_number();
  }

  _splits_generated++;

  return 1;
}

int Class_Members::write_split(
    IWString_and_File_Descriptor &train_output,
    IWString_and_File_Descriptor &test_output) const {
  for (int i = 0; i < _number_elements; i++) {
    const TTItem *it = _things[i];

    if (it->in_training_set()) {
      train_output << it->id();
      if (write_class_memberships)
        train_output << ' ' << _class_label;
      train_output << '\n';
    } else {
      test_output << it->id();
      if (write_class_memberships)
        test_output << ' ' << _class_label;
      test_output << '\n';
    }

    train_output.write_if_buffer_holds_more_than(8192);
    test_output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int Class_Members::report(std::ostream &os) const {
  os << "Report for class '" << _class_label << "' " << _number_elements
     << " items\n";

  int *times_in_training_set = new_int(number_splits + 1);
  std::unique_ptr<int[]> free_times_in_training_set(times_in_training_set);

  for (int i = 0; i < _number_elements; i++) {
    times_in_training_set[_things[i]->times_in_training_set()]++;
  }

  for (int i = 0; i < number_splits + 1; i++) {
    os << times_in_training_set[i] << " items in training set " << i
       << " times\n";
  }

  return 1;
}

int Class_Members::items_in_training_set() const {
  int rc = 0;

  for (int i = 0; i < _number_elements; i++) {
    if (_things[i]->in_training_set())
      rc++;
  }

  return rc;
}

static int open_split_file(IWString &fname,
                           IWString_and_File_Descriptor &output) {
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

static void form_file_names(IWString &train_fname, IWString &test_fname,
                            int ndx) {
  train_fname = output_stem;
  train_fname << ndx;
  train_fname << 'R';

  if (suffix.length())
    train_fname << suffix;

  test_fname = output_stem;
  test_fname << ndx;
  test_fname << 'E';

  if (suffix.length())
    test_fname << suffix;

  return;
}

static int test_train_split_classification(Class_Members &class1,
                                           Class_Members &class2) {
  for (int i = 0; i < number_splits; i++) {
    class1.generate_split(fraction_in_training_set);
    class2.generate_split(fraction_in_training_set);

    IWString train_fname, test_fname;
    form_file_names(train_fname, test_fname, i);
    IWString_and_File_Descriptor train_output, test_output;
    if (!open_split_file(train_fname, train_output) ||
        !open_split_file(test_fname, test_output)) {
      cerr << "Cannot open one or more of the split files '" << train_fname
           << "' and/or '" << test_fname << "'\n";
      return 0;
    }

    train_output << header_record << '\n';
    test_output << header_record << '\n';

    class1.write_split(train_output, test_output);
    class2.write_split(train_output, test_output);

    if (verbose) {
      int c1 = class1.items_in_training_set();
      int c2 = class2.items_in_training_set();
      cerr << "Created split " << i << " " << (c1 + c2)
           << " items in training set, " << c1 << " from '"
           << class1.class_label() << "' " << c2 << " from '"
           << class2.class_label() << "'\n";
      items_in_training_set[i] = c1 + c2;
    }
  }

  return 1;
}

static int read_identifier(const const_IWSubstring &buffer,
                           Class_Members &class1, Class_Members &class2) {
  if (buffer.nwords() < 2) {
    cerr << "Input must have at least ID and class label\n";
    return 0;
  }

  const_IWSubstring id;
  int i = 0;

  buffer.nextword(id, i);

  const_IWSubstring cm; // class membership
  buffer.nextword(cm, i);

  if (cm == class1.class_label())
    class1.extra(id);
  else if (cm == class2.class_label())
    class2.extra(id);
  else if (0 == class1.size()) {
    class1.set_class_label(cm);
    class1.extra(id);
  } else if (0 == class2.size()) {
    class2.set_class_label(cm);
    class2.extra(id);
  } else {
    cerr << "More than two classes detected, cannot process '" << cm << "'\n";
    return 0;
  }

  return 1;
}

static int read_identifiers(iwstring_data_source &input, Class_Members &class1,
                            Class_Members &class2) {
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (0 == header_record.size())
    header_record = buffer;

  while (input.next_record(buffer)) {
    if (!read_identifier(buffer, class1, class2)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return class1.size() + class2.size();
}

static int read_identifiers(const char *fname, Class_Members &class1,
                            Class_Members &class2) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_identifiers(input, class1, class2);
}

static int test_train_split_classification(int argc, char **argv) {
  Command_Line cl(argc, argv, "vn:p:s:S:iU:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n')) {
    if (!cl.value('n', number_splits) || number_splits < 1) {
      cerr << "The number of splits to create option (-n) must be a whole "
              "positive number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will create " << number_splits << " random splits\n";
  } else
    cerr << "Only creating one split by default\n";

  items_in_training_set = new_int(number_splits);
  std::unique_ptr<int[]> free_items_in_training_set(items_in_training_set);

  if (cl.option_present('p')) {
    int p;
    if (!cl.value('p', p) || p < 1 || p > 99) {
      cerr << "The percent in training set option (-p) must be a valid "
              "percentage\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will put " << p << " percent of the items in the training set\n";

    fraction_in_training_set =
        static_cast<float>(p) / static_cast<float>(100.0);
  }

#ifdef SEED_NO_LONGER_SETTABLE
  if (cl.option_present('s')) {
    unsigned int s;
    if (!cl.value('s', s)) {
      cerr << "The random number seed option (-s) must be a valid random "
              "number seed\n";
      usage(5);
    }

    if (verbose)
      cerr << "Random number seed " << s << '\n';

    iw_set_rnum_seed(s);
  } else
    iw_random_seed();
#endif

  if (cl.option_present('S')) {
    output_stem = cl.string_value('S');

    if (verbose)
      cerr << "Will create file(s) will stem '" << output_stem << "'\n";
  }

  if (cl.option_present('i')) {
    write_class_memberships = 0;

    if (verbose)
      cerr << "Will not write class membership data to output file(s)\n";
  }

  if (cl.option_present('U')) {
    cl.value('U', suffix);

    if (verbose)
      cerr << "Files created with suffix '" << suffix << "'\n";
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  Class_Members class1, class2;

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!read_identifiers(cl[i], class1, class2)) {
      cerr << "Cannot read data from '" << cl[i] << "'\n";
      return 0;
    }
  }

  if (verbose)
    cerr << "Read " << class1.size() << " items from class '"
         << class1.class_label() << "' and " << class2.size()
         << " items from class '" << class2.class_label() << "'\n";

  if (class1.empty() || class2.empty()) {
    cerr << "Zero or only one class present, cannot process\n";
    return 3;
  }

  test_train_split_classification(class1, class2);

  if (verbose) {
    class1.report(cerr);
    class2.report(cerr);
  }

  return 0;
}

int main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = test_train_split_classification(argc, argv);

  return rc;
}
