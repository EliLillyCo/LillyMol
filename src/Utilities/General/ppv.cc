/*
  Computes ROC scores
  Today this would be better done with R/Python/Julia
*/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/accumulator/accumulator.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/data_source/iwstring_data_source.h"
#define ACTIVITY_DATA_IMPLEMENATION_H
#include "Foundational/iwmisc/activity_data_from_file.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static IWString missing_value('.');

static IWString_and_File_Descriptor stream_for_xmgrace;

static int write_npv = 1;
static int write_sum_sens_ppv = 0;
static int write_sum_sens_ppv_xmgrace = 0;

static int display_npv_curves_in_xmgrace_output = 1;
static int display_specificity_curves_in_xmgrace_output = 0;

/*
  We will make either 3 or 4 plots for each file
*/

static int ngraphs = 3;

/*
  If we are doing output for xmgrace, we need to accumulate all
  the values computed and then write them out at the end
*/

static resizable_array<float> all_sensitivity;
static resizable_array<float> all_ppv;
static resizable_array<float> all_npv;
static resizable_array<float> all_specificity;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes ROC scores\n";
  cerr << " -E <fname>     file with experimental data\n";
  cerr << " -e <col>       experimental values in column <col> (default 2)\n";
  cerr << " -p <col>       predicted values in column <col> (default 2)\n";
  cerr << " -z             strip leading zero's from identifiers\n";
  cerr << " -q <float>     cutoff(s) for experimental class divisions\n";
  cerr << " -n <float>     cutoff(s) for predicted    class divisions\n";
  cerr << " -n R:x0,x1,dx  range of predicted class division values\n";
  cerr << " -s             write the sum of SENS + PPV\n";
  cerr << " -w             suppress writing NPV values\n";
  cerr << " -D <fname>     write a file of OBS PRED DIFF softed by diff\n";
  cerr << " -X <fname>     write an output file for xmgrace\n";
  cerr << " -T <fname>     existing xmgrace template to use\n";
  cerr << " -M ...         various other options\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

typedef float activity_type_t;

static Activity_Data_From_File<activity_type_t> activity_data;

template class Activity_Data_From_File<activity_type_t>;

template int Activity_Data_From_File<activity_type_t>::construct_from_command_line(const Command_Line & cl, char, int);

static int predicted_column = 1;

class ID_and_Activities
{
  private:
    IWString _id;
    activity_type_t _obs;
    activity_type_t _pred;

  public:
    ID_and_Activities(const const_IWSubstring & s, activity_type_t o, activity_type_t p) :
      _id(s), _obs(o), _pred(p) {}

    const IWString & id() const { return _id;}
    activity_type_t obs()  const { return _obs;}
    activity_type_t pred() const { return _pred;}
};

template class resizable_array_p<ID_and_Activities>;
template class resizable_array_base<ID_and_Activities *>;

//typdef class ID_and_Activities_template<activity_type_t> ID_and_Activities;

class ID_and_Activities_Comparator
{
  private:
  public:
    int operator() (const ID_and_Activities *, const ID_and_Activities *) const;
};

int
ID_and_Activities_Comparator::operator() (const ID_and_Activities * idac1, 
                                          const ID_and_Activities * idac2) const
{
  activity_type_t d1 = fabs(idac1->obs() - idac1->pred());
  activity_type_t d2 = fabs(idac2->obs() - idac2->pred());

  if (d1 < d2)
    return 1;
  else if (d1 > d2)
    return -1;
  return 0;
}

template void resizable_array_base<ID_and_Activities*>::iwqsort<ID_and_Activities_Comparator>(ID_and_Activities_Comparator&);
template void iwqsort<ID_and_Activities*, ID_and_Activities_Comparator>(ID_and_Activities**, int, ID_and_Activities_Comparator&);
template void iwqsort<ID_and_Activities*, ID_and_Activities_Comparator>(ID_and_Activities**, int, ID_and_Activities_Comparator&, void*);
template void compare_two_items<ID_and_Activities*, ID_and_Activities_Comparator>(ID_and_Activities**, ID_and_Activities_Comparator&, void*);
template void move_in_from_left<ID_and_Activities*, ID_and_Activities_Comparator>(ID_and_Activities**, int&, int&, int, ID_and_Activities_Comparator&, void*);
template void move_in_from_right<ID_and_Activities*, ID_and_Activities_Comparator>(ID_and_Activities**, int&, int&, ID_and_Activities_Comparator&);
template void swap_elements<ID_and_Activities*>(ID_and_Activities*&, ID_and_Activities*&, void*);

static int
all_values_the_same (const resizable_array<activity_type_t> & c)
{
  int n = c.number_elements();

  activity_type_t c0 = c[0];

  for (int i = 1; i < n; i++)
  {
    if (c[i] != c0)
      return 0;
  }

  return 1;
}

/*
  We need to change some of the records in the template file
*/

static const IWString world("@    world 0, 0.2, 100, 1");
static const IWString title("@    title \"5HT2A\"");
static const IWString s0legend("@    s0 legend  ");
static const IWString s1legend("@    s1 legend  ");
static const IWString s2legend("@    s2 legend  ");
static const IWString tick_major("@    xaxis  tick major 20");
static const IWString subtitle("@    subtitle \"\"");

static int
echo_xmgrace_template (iwstring_data_source & input,
                       const Command_Line & cl,
                       const resizable_array<activity_type_t> & experimental_cutoff_value,
                       const resizable_array_p<ID_and_Activities> & zdata,
                       IWString_and_File_Descriptor & output)
{
  activity_type_t min_activity = zdata[0]->obs();
  activity_type_t max_activity = zdata[0]->obs();

  int n = zdata.number_elements();

  for (int i = 1; i < n; i++)
  {
    activity_type_t o = zdata[i]->obs();

    if (o < min_activity)
      min_activity = o;
    else if (o > max_activity)
      max_activity = o;
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (world == buffer)
    {
      output << "@    world " << min_activity << ", 0.2, " << max_activity << ", ";
      if (write_sum_sens_ppv_xmgrace)
        output << "1.5\n";
      else
        output << "1\n";
    }
    else if (title == buffer)
    {
      char dirname[300];
      if (! getcwd(dirname, sizeof(dirname))) {
        cerr << "Cannot get pwd\n";
      }
      const_IWSubstring tmp(dirname);
      int last_slash = tmp.rindex('/');
      if (last_slash > 0)
        tmp.set(dirname + last_slash + 1, tmp.length() - last_slash - 1);

      output << "@    title \"" << tmp << "\"\n";
    }
    else if (subtitle == buffer && all_values_the_same(experimental_cutoff_value))
    {
      output << "@    subtitle \"Cutoff " << experimental_cutoff_value[0] << "\"\n";
    }
    else if (tick_major == buffer)
    {
      float tm = (max_activity - min_activity) / 5.0;
      output << "@    xaxis  tick major " << tm << "\n";
    }
    else if (cl.number_elements() > 1 && buffer.starts_with(s0legend))
      ;
    else if (cl.number_elements() > 1 && buffer.starts_with(s1legend))
      ;
    else if (cl.number_elements() > 1 && buffer.starts_with(s2legend))
      ;
    else
      output << buffer << '\n';
  }

  if (cl.number_elements() > 1)   // write file names as legends
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      output << "@    s" << (ngraphs * i) << " legend \"" << cl[i] << "\"\n";
    }
  }

  output.write_whole_blocks_shift_unwritten();

  return 1;
}

/*static int
echo_xmgrace_template (const char * fname,
                       const resizable_array_p<ID_and_Activities> & zdata,
                       IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open xmgrace template file '" << fname << "'\n";
    return 0;
  }

  return echo_xmgrace_template(input, zdata, output);
}*/

static int
write_xmgrace_set(const resizable_array<float> & x,
                  const resizable_array<float> & y,
                  int set_number,
                  IWString_and_File_Descriptor & output)
{
  int n = x.number_elements();
  assert (n == y.number_elements());

  output << "@target G0.S" << set_number << "\n";
  output << "@type xy\n";

  for (int i = 0; i < n; i++)
  {
    output << x[i] << ' ' << y[i] << '\n';
  }
  output << "&\n";

  return output.size();
}

static int
do_xmgrace_output (const resizable_array<float> & cutoff,
                   int ndx,
                   IWString_and_File_Descriptor & output)
{
  int n = cutoff.number_elements();
  assert (n == all_sensitivity.number_elements());
  assert (n == all_ppv.number_elements());

  write_xmgrace_set(cutoff, all_sensitivity, ndx, output);
  ndx++;

  write_xmgrace_set(cutoff, all_ppv, ndx, output);
  ndx++;

  if (display_npv_curves_in_xmgrace_output)
  {
    assert (n == all_npv.number_elements());
    write_xmgrace_set(cutoff, all_npv, ndx, output);
    ndx++;
  }

  if (display_specificity_curves_in_xmgrace_output)
  {
    assert (n == all_specificity.number_elements());
    write_xmgrace_set(cutoff, all_specificity, ndx, output);
    ndx++;
  }

  if (write_sum_sens_ppv_xmgrace)
  {
    resizable_array<float> tmp;
    tmp.resize(n);

    for (int i = 0; i < n; i++)
    {
      tmp.add(all_sensitivity[i] + all_ppv[i]);
    }

    write_xmgrace_set(cutoff, tmp, ndx, output);
    ndx++;
  }

  output.flush();

  return 1;
}

static int
ppv (const resizable_array_p<ID_and_Activities> & zdata,
     activity_type_t experimental_cutoff_value,
     activity_type_t prediction_cutoff_value,
     float & score,
     IWString_and_File_Descriptor & output)
{
  int n = zdata.number_elements();

  int e1 = 0;    // items in experimental class 1
  int e2 = 0;    // items in experimental class 2

  int p1 = 0;    // items in predicted class 1
  int p2 = 0;    // items in predicted class 2

  int true_positives = 0;
  int true_negatives = 0;
  int false_positives = 0;
  int false_negatives = 0;

  for (int i = 0; i < n; i++)
  {
    const ID_and_Activities * idac = zdata[i];

    activity_type_t obs  = idac->obs();

    int expt_class;

    if (obs < experimental_cutoff_value)
    {
      e1++;
      expt_class = -1;
    }
    else
    {
      e2++;
      expt_class = 1;
    }
    
    int pred_class;

    activity_type_t pred = idac->pred();

    if (pred < prediction_cutoff_value)
    {
      p1++;
      pred_class = -1;
    }
    else
    {
      p2++;
      pred_class = 1;
    }

    if (pred_class == expt_class)
    {
      if (1 == pred_class)
        true_positives++;
      else
        true_negatives++;
    }
    else if (-1 == pred_class)
      false_negatives++;
    else
      false_positives++;
  }

  set_default_iwstring_float_concatenation_precision(3);

  int total_actives = true_positives + false_negatives;

  if (0 == total_actives)
  {
    cerr << "At expt cuttof " << experimental_cutoff_value << ", no actives\n";
    return 0;
  }

  if (0 == true_negatives && 0 == false_positives)
  {
    cerr << "At expt cuttof " << experimental_cutoff_value << ", no in-actives\n";
    return 0;
  }

  float sensitivity = static_cast<float>(true_positives) / static_cast<float>(true_positives + false_negatives);
  float ppv_value;
  if (0 == true_positives + false_positives)
    ppv_value = 0.0;
  else
    ppv_value = static_cast<float>(true_positives) / static_cast<float>(true_positives + false_positives);
  
  float specificity = static_cast<float>(true_negatives) / static_cast<float>(true_negatives + false_positives);
  float npv;
  if (0 == true_negatives && 0 == false_negatives)
    npv = 0.0;
  else
    npv = static_cast<float>(true_negatives) / static_cast<float>(true_negatives + false_negatives);

  output << "tot " << n << " active " << total_actives << ' ' << experimental_cutoff_value;
//if (experimental_cutoff_value != prediction_cutoff_value)
    output << ',' << prediction_cutoff_value;

  output << ": sensitivity " << sensitivity << " PPV " << ppv_value;

  if (write_sum_sens_ppv)
  {
    set_default_iwstring_float_concatenation_precision(4);
    output << ' ' << (ppv_value + sensitivity) << ' ' << (ppv_value + sensitivity + npv + specificity);
    set_default_iwstring_float_concatenation_precision(3);
  }

  output << '\n';

  if (write_npv)
  {
    output << "tot " << n << " active " << total_actives << ' ' << experimental_cutoff_value;
//  if (experimental_cutoff_value != prediction_cutoff_value)
      output << ',' << prediction_cutoff_value;

    output << ": specificity " << specificity << " NPV " << npv << '\n';
  }

  if (stream_for_xmgrace.is_open())
  {
    all_sensitivity.add(sensitivity);
    all_ppv.add(ppv_value);
    if (display_npv_curves_in_xmgrace_output)
      all_npv.add(npv);
    if (display_specificity_curves_in_xmgrace_output)
      all_specificity.add(specificity);
  }

  score = (ppv_value + sensitivity + npv + specificity);

  return 1;
}

static int
write_difference_file(resizable_array_p<ID_and_Activities> & zdata,
                      IWString_and_File_Descriptor & output)
{
  ID_and_Activities_Comparator idac;

  zdata.iwqsort(idac);

  output << "Id obs pred diff\n";

  int n = zdata.number_elements();

  for (int i = 0; i < n; i++)
  {
    const ID_and_Activities * idac = zdata[i];

    activity_type_t obs = idac->obs();
    activity_type_t pred = idac->pred();
    activity_type_t diff = fabs(obs - pred);

    output << idac->id() << ' ' << obs << ' ' << pred << ' ' << diff << '\n';

    if (output.size() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  output.flush();

  return 1;
}

static int
write_difference_file(const char * fname,
                      resizable_array_p<ID_and_Activities> & zdata)
{
  IWString_and_File_Descriptor output;

  if (! output.open(fname))
  {
    cerr << "Cannot open activity difference file '" << fname << "'\n";
    return 0;
  }

  return write_difference_file(zdata, output);
}

static int
create_dataitem (const const_IWSubstring & buffer,
                 resizable_array_p<ID_and_Activities> & zdata,
                 const Activity_Data_From_File<activity_type_t> & activity_data)
{
  const_IWSubstring id;
  int i = 0;

  if (! buffer.nextword(id, i))
  {
    cerr << "Cannot extract identifier\n";
    return 0;
  }

  float expt;

  if (! activity_data.get_activity(id, expt))
  {
    cerr << "No experimental data for '" << id << "'\n";
    return 0;
  }

  const_IWSubstring token;

  for (int col = 1; buffer.nextword(token, i); col++)
  {
    if (col != predicted_column)
      continue;

    float pred;
    if (token.numeric_value(pred))
    {
      ID_and_Activities * ida = new ID_and_Activities(id, expt, pred);
      zdata.add(ida);
      return 1;
    }
    else if (missing_value == token)
      return 1;
    else
    {
      cerr << "Invalid predicted data '" << token << "'\n";
      return 0;
    }
  }

  cerr << "Did not find predicted column " << (predicted_column + 1) << '\n';
  return 0;
}

static int
read_predicted_data (iwstring_data_source & input,
                     resizable_array_p<ID_and_Activities> & zdata)

{
  zdata.resize_keep_storage(0);

  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    if (! create_dataitem (buffer, zdata, activity_data))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  int n = zdata.size();

  if (verbose)
    cerr << "Read " << n << " predicted data values\n";

  if (0 == n)
  {
    cerr << "No predicted data read\n";
    return 0;
  }


  return 1;
}

static int
read_predicted_data (const char * fname,
                     resizable_array_p<ID_and_Activities> & zdata)
{
  zdata.resize_keep_storage(0);

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_predicted_data(input, zdata);
}

static int
report_data_extremeties (const resizable_array_p<ID_and_Activities> & zdata,
                         int ndx,
                         int number_files_with_predicted_data,
                         std::ostream & output)
{
  Accumulator<double> acc_obs, acc_pred;

  int n = zdata.number_elements();

  for (int i = 0; i < n; i++)
  {
    const ID_and_Activities * ida = zdata[i];
    acc_obs.extra(ida->obs());
    acc_pred.extra(ida->pred());
  }

  if (0 == ndx)
    output << "Experimental values between " << acc_obs.minval() << " and " << acc_obs.maxval() << " ave " << static_cast<float>(acc_obs.average()) << '\n';

  if (number_files_with_predicted_data > 1)
    output << ndx << ' ';

  output << "Predicted    values between " << acc_pred.minval() << " and " << acc_pred.maxval() << " ave " << static_cast<float>(acc_pred.average()) << '\n';

  return output.good();
}

/*
  we can only do xmgrace output if either
    all experimental cutoffs the same value
    expt == pred
*/

static int
cutoffs_consistent_with_xmgrace_output(const resizable_array<activity_type_t> & experimental_cutoff_value,
                                       const resizable_array<activity_type_t> & prediction_cutoff_value)
{
  int n = experimental_cutoff_value.number_elements();

  if (n != prediction_cutoff_value.number_elements())
    return 0;

  if (all_values_the_same(experimental_cutoff_value))
    return 1;

//cerr << "Checking both arrays with same contents\n";
  for (int i = 0; i < n; i++)
  {
//  cerr << experimental_cutoff_value[i] << " vs " << prediction_cutoff_value[i] << endl;
    if (experimental_cutoff_value[i] != prediction_cutoff_value[i])
      return 0;
  }

  return 1;
}

static int
fill_to_equal_sizes(const resizable_array<activity_type_t> & c1,
                    resizable_array<activity_type_t> & c2)
{
  for (int i = c2.number_elements(); i < c1.number_elements(); i++)
  {
    c2.add(c1[i]);
  }

  return 1;
}

/*
  Range looks like
  NNN
  XXX,YYY
  XXX,YYY,DDD
*/

static int
parse_range (const const_IWSubstring & r,
             resizable_array<activity_type_t> & cvalues)
{
  int nw = r.nwords(',');

  if (1 == nw)
  {
    activity_type_t a;
    if (! r.numeric_value(a))
    {
      cerr << "INvalid numeric specifier for range '" << r << "'\n";
      return 0;
    }

    cvalues.add(a);
    return 1;
  }

  const_IWSubstring token;
  int i = 0;

  r.nextword(token, i, ',');

  activity_type_t a0, a1;
  activity_type_t dx = 1.0;

  if (! token.numeric_value(a0))
  {
    cerr << "Invalid start of range '" << r << "'\n";
    return 0;
  }

  r.nextword (token, i, ',');

  if (! token.numeric_value(a1))
  {
    cerr << "Invalid end of range '" << r << "'\n";
    return 0;
  }

  if (a1 <= a0)
  {
    cerr << "Range is invalid " << a0 << " to " << a1 << endl;
    return 0;
  }

  if (r.nextword (token, i, ','))
  {
    if (! token.numeric_value(dx) || dx <= 0.0)
    {
      cerr << "Invalid step value '" << r << "'\n";
      return 0;
    }

    if (dx > (a1 - a0))
    {
      cerr << "IMpossible range '" << r << "', " << a0 << " to " << a1 << " by " << dx << endl;
      return 0;
    }
  }

  for (int i = 0; ; i++)
  {
    activity_type_t c = a0 + static_cast<activity_type_t>(i) * dx;

//  cerr << "Adding " << c << endl;
    cvalues.add(c);

    if (c >= a1)
      break;
  }

  return cvalues.number_elements();
}

/*
  Lots of possibilities about what a range looks like

  3-4
  -3-4
  -3--2
*/

/*static int
parse_range(const_IWSubstring r,    // local copy
            resizable_array<activity_type_t> & cvalues)
{
  float dx = 1.0;

  int icomma = r.index(',');

  if (icomma > 0)
  {
    const_IWSubstring b4comma, aftercomma;
    r.split (b4comma, ',', aftercomma);

    if (! aftercomma.numeric_value(dx) || dx <= 0.0)
    {
      cerr << "Invalid delta '" << r << "'\n";
      return 0;
    }

    r = b4comma;
  }

  int lhs_sign;

  if (r.starts_with('-'))
  {
    lhs_sign = -1;
    r++;
  }
  else
    lhs_sign = 1;

  int idash = r.index('-');

  const_IWSubstring rstart, rend;

  r.from_to(0, idash - 1, rstart);
  r.from_to(idash + 1, r.length() - 1, rend);

  if (verbose > 1)
    cerr << "Components '" << rstart << "' and '" << rend << "'\n";

  int r1, r2;

  if (! rstart.numeric_value(r1))
    ;
  else if (! rend.numeric_value(r2))
    ;
  else if ((lhs_sign * r1) >= r2)
    ;
  else
  {
    r1 = lhs_sign * r1;

    if (verbose > 1)
      cerr << "Range is " << r1 << " to " << r2 << ", dx " << dx << endl;

    if (1.0 == dx)
    {
      for (int i = r1; i <= r2; i++)
      {
        cvalues.add(static_cast<activity_type_t>(i));
      }
    }
    else
    {
      int istop = static_cast<int>((r2 - r1) / dx) + 1;

      for (int i = 0; i < istop; i++)
      {
        activity_type_t v = r1 + i * dx;
        cvalues.add(v);
      }
    }

    return 1;
  }

  return 0;
}*/

static int
get_cutoff_values (const Command_Line & cl,
                   char flag,
                   resizable_array<activity_type_t> & cvalues)
{
  int i = 0;
  const_IWSubstring n;
  while (cl.value(flag, n, i++))
  {
    if (n.starts_with("R:"))
    {
      n.remove_leading_chars(2);

      if (! parse_range(n, cvalues))
      {
        cerr << "Bad range '" << n << "'\n";
        return 0;
      }

      return 1;
    }

    int j = 0;
    const_IWSubstring token;
    while (n.nextword(token, j, ','))
    {
      activity_type_t nv;
      if (! token.numeric_value(nv))
      {
        cerr << "Invalid class specification '" << n << "'\n";
        return 0;
      }

      cvalues.add_if_not_already_present(nv);
    }
  }

  return cvalues.number_elements();
}

static void
display_dash_M_options (std::ostream & os)
{
  os << " -M sumx      display plot of Sens + PPV in the xmgrace output\n";
  os << " -M nonpv     omit display of NPV curves in xmgrace output\n";
  os << " -M spec      display specificity curves in xmgrace output\n";

  exit(0);
}

static int
ppv (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vE:e:p:D:n:q:X:T:wsM:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('E'))
  {
    cerr << "Must specify experimental data via the -E option\n";
    usage(3);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('w'))
  {
    write_npv = 0;
    if (verbose)
      cerr << "Will suppress writing NPV values\n";
  }

  if (cl.option_present('s'))
  {
    write_sum_sens_ppv = 1;
    if (verbose)
      cerr << "Will write the sum of sensitivity+PPV\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', predicted_column) || predicted_column < 2)
    {
      cerr << "Invalid predicted values column\n";
      usage (4);
    }

    if (verbose)
      cerr << "Predicted values in column " << predicted_column << endl;

    predicted_column--;
  }

  if (cl.option_present('e'))
  {
    int e;
    if (! cl.value('e', e) || e < 1)
    {
      cerr << "The experimental column must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Experimental data in column " << e << " of the activity file\n";

    activity_data.set_activity_column(e - 1);
  }

  if (cl.option_present('z'))
  {
    activity_data.set_strip_leading_zeros(1);
  }

  if (! activity_data.construct_from_command_line (cl, 'E', verbose))
  {
    cerr << "Cannot read experimental data\n";
    return 5;
  }

  if (cl.option_present('n'))
    ;
  else if (cl.option_present('q'))
    ;
  else
  {
    cerr << "Must specify one or both of the -n or -q options to establich a cutoff value\n";
    usage(4);
  }

  resizable_array<activity_type_t> prediction_cutoff_value;

  if (! cl.option_present('n'))
    ;
  else if (! get_cutoff_values (cl, 'n', prediction_cutoff_value))
  {
    cerr << "Cannot determine -n values\n";
    usage(5);
  }

  resizable_array<activity_type_t> experimental_cutoff_value;

  if (! cl.option_present('q'))
    ;
  else if (! get_cutoff_values (cl, 'q', experimental_cutoff_value))
  {
    cerr << "Cannot determine -q values\n";
    usage(5);
  }

  if (experimental_cutoff_value.number_elements() == prediction_cutoff_value.number_elements())
    ;
  else if (1 == experimental_cutoff_value.number_elements())
  {
    int nextra = prediction_cutoff_value.number_elements() - 1;

    activity_type_t p0 = experimental_cutoff_value[0];

    for (int i = 0; i < nextra; i++)
    {
      experimental_cutoff_value.add(p0);
    }
  }
  else if (0 == experimental_cutoff_value.number_elements())
    fill_to_equal_sizes(prediction_cutoff_value, experimental_cutoff_value);
  else if (0 == prediction_cutoff_value.number_elements())
    fill_to_equal_sizes(experimental_cutoff_value, prediction_cutoff_value);
  else
  {
    cerr << "Not sure what to do with " << experimental_cutoff_value.number_elements() << " experimental cutoffs and " << prediction_cutoff_value.number_elements() << " predicted ones?\n";
    usage(4);
  }

  if (verbose > 2)
  {
    cerr << experimental_cutoff_value.number_elements() << " experimental cutoffs and " << prediction_cutoff_value.number_elements() << " prediction cutoff values\n";
  }

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (m.starts_with("sumx"))
      {
        write_sum_sens_ppv_xmgrace = 1;
        ngraphs += 1;
        if (verbose)
          cerr << "Will write sens+ppv sums to the xmgrace output\n";
      }
      else if ("nonpv" == m)
      {
        display_npv_curves_in_xmgrace_output = 0;

        ngraphs--;

        if (verbose)
          cerr << "Will not display NPV curves in xmgrace eoutput\n";
      }
      else if ("spec" == m)
      {
        display_specificity_curves_in_xmgrace_output = 1;
        ngraphs++;
      }
      else if ("help" == m)
      {
        display_dash_M_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_M_options(cerr);
      }
    }
  }

// Initialise any stuff for xmgrace.

  iwstring_data_source xmgrace_template;

  if (cl.option_present('X'))
  {
    if (prediction_cutoff_value.number_elements() < 2)
    {
      cerr << "The -X option doesn't make sense with just 1 datapoint\n";
      return 4;
    }

    if (! cutoffs_consistent_with_xmgrace_output(experimental_cutoff_value, prediction_cutoff_value))
    {
      cerr << "Output for xmgrace doesn't make sense with both multiple -q values specified\n";
      usage(3);
    }

    const char * x = cl.option_value('X');

    if (! stream_for_xmgrace.open(x))
    {
      cerr << "Cannot initialise stream for xmgrace '" << x << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Output for xmgrace written to '" << x << "'\n";

    if (cl.option_present('T'))
    {
      const char * t = cl.option_value('T');

      if (! xmgrace_template.open(t))
      {
        cerr << "Cannot open xmgrace template file '" << t << "'\n";
        return 4;
      }
    }
  }

// We are ready to do work now

  float highest_score = static_cast<float>(0.0);
  float highest_score_cutoff = static_cast<float>(0.0);
  float middle_of_range = (prediction_cutoff_value[0] + prediction_cutoff_value.last_item()) * 0.5;

  IWString_and_File_Descriptor output(1);

  resizable_array_p<ID_and_Activities> zdata;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! read_predicted_data(cl[i], zdata))
      return i + 1;

    int n = zdata.number_elements();
    if (n < 3)
    {
      cerr << "Too few points " << n << ", cannot continue\n";
      return 8;
    }

    if (! cl.option_present('D'))   // only write out difference file on first dataset
      ;
    else if (0 == i)
    {
      const char * d = cl.option_value('D');

      write_difference_file(d, zdata);   // ignore any error messages
    }
    else
      cerr << "Only writes differences with first predicted set\n";

    if (verbose)
      report_data_extremeties(zdata, i, cl.number_elements(), cerr);

    for (int j = 0; j < experimental_cutoff_value.number_elements(); j++)
    {
      float score;
      if (! ppv(zdata, experimental_cutoff_value[j], prediction_cutoff_value[j], score, output))
      {
        cerr << "Fatal error processing '" << cl[i] << "'\n";
        return i + 1;
      }

      score = static_cast<int>(score * 100)/100.0;    // round

      if (score < highest_score)
        ;
      else if (score > highest_score)
      {
        highest_score = score;
        highest_score_cutoff = prediction_cutoff_value[j];
      }
      else if (fabs(middle_of_range - prediction_cutoff_value[j]) <
               fabs(middle_of_range - highest_score_cutoff))
      {
        highest_score = score;
        highest_score_cutoff = prediction_cutoff_value[j];
      }

      if (output.size() > 32768)
        output.write_whole_blocks_shift_unwritten();
    }

    cerr << "Highest score " << highest_score << " at " << highest_score_cutoff << '\n';

    if (stream_for_xmgrace.is_open())
    {
      if (0 == i && xmgrace_template.is_open())
      {
        echo_xmgrace_template (xmgrace_template, cl, experimental_cutoff_value, zdata, stream_for_xmgrace);
        xmgrace_template.do_close();
      }

      do_xmgrace_output(prediction_cutoff_value, ngraphs * i, stream_for_xmgrace);
    }

    all_sensitivity.resize_keep_storage(0);
    all_ppv.resize_keep_storage(0);
    all_npv.resize_keep_storage(0);
    all_specificity.resize_keep_storage(0);
    zdata.resize_keep_storage(0);
  }

  output.flush();

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ppv(argc, argv);

  return rc;
}
