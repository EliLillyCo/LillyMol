/*
  We want to correct bias in a prediction
  We read in a file that contains
  ID Activity Pred1 Pred2 ....
  We do a regression of difference vs activity.
  We write the results to a file

  There are much better ways of doing this today.
*/

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int linear = 1;
static int quadratic = 0;

#define  QUADRATIC_MODEL "quadratic"
#define  LINEAR_MODEL "linear"

//    subroutine lsq2_ (xd,yd,ndata,npoly,c)

// extern "C" void lsq2_(const float *, const float *, const int *, const int *, float *);
extern "C" void sgefa_(const float *,const int *,const int *,int *,const int *);
extern "C" void sgesl_(const float *,const int *,const int *,int *,float * c ,const int *);

static IWString_and_File_Descriptor stream_for_diff_vs_activity;

static int header_records_to_skip = 1;

static IWString rfile;

static IWString rfile_plotfile;

static int prediction_column = 1;

static int replace_existing_predicted_values = 1;

static char output_separator = ' ';

static IWString current_file_name;

const char * month [] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

static int take_stratified_sample = 0;

/*
  Observed some data sets where there are a lot of values at one end. These
  are often poorly predicted and we get a graph of obs vs pred that might look
  like

  | *
  | *                    *
  | *                  *
  | *                *
  | *              *
  | *            *
  | *          *
  | *        *
  | *      *
  | *    *
  | * *
  | *
  _____________________________________________________

  For this reason, we can trim the highest and/or lowest experimental
  values.
  For now, we make this a boolean and drop all points that have the min/max activity
*/

static int drop_lower_expt = 0;
static int drop_upper_expt = 0;

#define IGNORED_VALUE -7.2341e-26f

#define MIDPOINT_NOT_SET 2776704.0f
static float midpoint = MIDPOINT_NOT_SET;

static IWString ylab;

/*
  If we are doing an unscaling, and the user has asked for R output, we need
  to accumulate the raw and transformed values
*/

static resizable_array<float> raw, transformed;

class Pred
{
  private:
    std::vector<float> _pred;

  public:
    Pred (float);

    int extra (const const_IWSubstring &);

    void extra (float v) { _pred.push_back(v);}

    int fill_arrays (float a, resizable_array<float> & activity, resizable_array<float> & predicted) const;
};

Pred::Pred (float f)
{ 
  _pred.push_back(f);

  return;
}

int
Pred::extra (const const_IWSubstring & s)
{
  float f;
  if (! s.numeric_value(f))
  {
    cerr << "Pred::extra:invalid value '" << s << "'\n";
    return 0;
  }

  _pred.push_back(f);

  return 1;
}

int
Pred::fill_arrays (float a,
                   resizable_array<float> & activity,
                   resizable_array<float> & predicted) const
{
  auto n = _pred.size();

  for (decltype(n) i = 0; i < n; ++i)
  {
    activity.add(a);
    predicted.add(_pred[i]);
  }

  return n;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Finds a least squares fit of prediction diff vs activity which can be used to adjust predictions?\n";
  cerr << " -A <fname>     file with observed values\n";
  cerr << " -m <midpoint>  when building, take a stratified sample of values below <midpoint>\n";
  cerr << "                when correcting, do not correct predicted values below <midpoint>\n";
  cerr << " -s <nb>        number of buckets to use for stratified sampling (default 100)\n";
  cerr << " -F <fname>     file containing list of predicted files (one per line)\n";
//cerr << " -q             do a quadratic fit\n";
  cerr << " -D <fname>     write a file of 'activity diff' pairs\n";
  cerr << " -R <fname>     create a plot of the fit\n";
  cerr << " -d <...>       specify R plot file (foo.pdf) for example\n";
  cerr << " -Y <ylab>      y axis label for R plot (default from activity file)\n";
  cerr << " -h <n>         header records to skip (default 1)\n";
  cerr << " -U <fname>     adjust a prediction based on a previously generated result from this programme\n";
  cerr << " -p <col>       when correcting bias, predictions are in column <col>\n";
  cerr << " -i             insert the corrected results into the output (default is to replace the column)\n";
  cerr << " -c             drop all experimental values with the minimum value\n";
  cerr << " -C             drop all experimental values with the minimum value\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
replace_lower_upper_values (IW_STL_Hash_Map<IWString, float> & obs,
                            int drop_lower_expt,
                            int drop_upper_expt)
{
  iwminmax<float> m(std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

  for (const auto& o : obs)
  {
    m.extra(o.second);
  }

  cerr << "Observed values range between " << m.minval() << " and " << m.maxval() << endl;

  int rc = 0;

  for (auto & o : obs)
  {
    if (drop_lower_expt && o.second == m.minval())
    {
      o.second = IGNORED_VALUE;
      rc++;
    }
    if (drop_upper_expt && o.second == m.maxval())
    {
      o.second = IGNORED_VALUE;
      rc++;
    }
  }

  return rc;
}

static float
linear_unscaled (const float * fit,
                 float v)
{
  if (MIDPOINT_NOT_SET != midpoint && v < midpoint)   // do not change
    return v;

//float delta = (1.0 - fit[1]) * v - fit[0];
  const float newv = (v - fit[0]) / fit[1];

  if (rfile.length())
  {
    raw.add(v);
    transformed.add(newv);
  }

  return newv;
}

static float
quadratic_unscaled (const float * fit,
                   float v)
{
  float a = fit[2];
  float b = fit[1];
  float c = fit[0] - v;

  float s = b*b - 4.0f * a * c;

  assert (s >= 0.0f);

  float v1 = (-b + sqrt(s)) / (2.0f * a);
  float v2 = (-b - sqrt(s)) / (2.0f * a);

//cerr << "From " << v << " get " << v1 << " and " << v2 << endl;

  if (v1 < v2)
    return v1;

  return v2;
}

/*
  We are doing an unscaling and need to write the proper header
*/

static int
append_to_column_header (const const_IWSubstring & buffer,
                         int pcol,
                         const char * to_append,
                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring token;
  int i = 0;

  for (int col = 0; buffer.nextword(token, i); col++)
  {
    if (col > 0)
      output << output_separator;

    output << token;

    if (pcol == col)
    {
      output << to_append;

      if (! replace_existing_predicted_values)
        output << output_separator << token;
    }
  }

  output << "\n";

  return 1;
}

static int
do_correction_record (const const_IWSubstring & buffer,
                      const float * fit,
                      int degree,
                      IWString_and_File_Descriptor & output)
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); col++)
  {
    if (col > 0)
      output << output_separator;

    if (col != prediction_column)
      output << token;
    else
    {
      float v;
      if (! token.numeric_value(v))
      {
        cerr << "Invalid numeric '" << token << "'\n";
        return 0;
      }

      if (1 == degree)
        output << linear_unscaled(fit, v);
      else if (2 == degree)
        output << quadratic_unscaled (fit, v);

      if (! replace_existing_predicted_values)
        output << output_separator << token;
    }
  }

  output << "\n";

  return 1;
}

static int 
do_correction (iwstring_data_source & input,
               const float * fit,
               int degree,
               IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  append_to_column_header(buffer, prediction_column, (1 == degree ? ".L" : ".Q"), output);

  while (input.next_record(buffer))
  {
    if (! do_correction_record (buffer, fit, degree, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int 
do_correction (const char * fname,
               const float * fit,
               int degree,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return do_correction (input, fit, degree, output);
}

static int
parse_fit_data (const IWString & buffer,
                float * fit,
                int degree)
{
  if (1 == degree && (2+1) == buffer.nwords())
    ;
  else if (2 == degree && (3+1) == buffer.nwords())
    ;
  else
  {
    cerr << "Incorrect token count in input, got '" << buffer.nwords() << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  for (int ndx = 0; buffer.nextword(token, i); ndx++)
  {
    if (! token.numeric_value(fit[ndx]))
    {
      cerr << "Non numeric input '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
write_r_plotting_instructions (const IWString & rfile_plotfile,
                         IWString_and_File_Descriptor & output)
{
  if (rfile_plotfile.ends_with(".pdf"))
    output << "pdf('" << rfile_plotfile << "')\n";
  else if (rfile_plotfile.ends_with(".png"))
    output << "png('" << rfile_plotfile << "')\n";
  else
  {
    cerr << "Unrecognised graphics file type '" << rfile_plotfile << "'\n";
    return 0;
  }

  return 1;
}

static int
write_rfile_transformed (const resizable_array<float> & raw,
                         resizable_array<float> & transformed,
                         IWString_and_File_Descriptor & output)
{
  if (rfile_plotfile.length() > 0)
    write_r_plotting_instructions (rfile_plotfile, output);

  int n = raw.number_elements();
  assert (n == transformed.number_elements());

  float min_raw = raw[0];
  float max_raw = min_raw;
  output << "raw=c(" << raw[0];
  for (int i = 1; i < n; i++)
  {
    float r = raw[i];
    output << ',' << r;
    if (r < min_raw)
      min_raw = r;
    else if (r > max_raw)
      max_raw = r;
  }
  output << ")\n";

  float min_transformed = transformed[0];
  float max_transformed = min_transformed;
  output << "transformed=c(" << transformed[0];
  for (int i = 1; i < n; i++)
  {
    float t = transformed[i];
    output << ',' << t;
    if (t < min_transformed)
      min_transformed = t;
    else if (t > max_transformed)
      max_transformed = t;
  }
  output << ")\n";
  output << "plot(raw,transformed,xlab='raw',ylab='adjusted',col='red',pch=3,cex=0.8, main='Prediction Bias Correction')\n";

  output << "lines(c(" << min_raw << ',' << max_raw << "),c(" << min_raw << ',' << max_raw << "), col='blue', lw=2)\n";

  time_t tnow = ::time(NULL);
  struct tm * tm = localtime(&tnow);

  output << "mtext(\"" << (tm->tm_year+1900) << '-' << month[tm->tm_mon] << '-' << (tm->tm_mday < 10 ? "0" : "") << tm->tm_mday << ' ' << tm->tm_hour << ':' << (tm->tm_min < 10 ? "0" : "") << tm->tm_min << ':' << (tm->tm_sec < 10 ? "0" : "") << tm->tm_sec << "\", col='darkgreen', side=4, cex=0.8)\n";
  output << "mtext('blue line is y=1.0*x (unchanged)',side=3,col='blue')\n";


  return 1;
}

static int
write_rfile_transformed (const resizable_array<float> & raw,
                         resizable_array<float> & transformed,
                         IWString & rfile)
{
  IWString_and_File_Descriptor output;
  if (! output.open(rfile.null_terminated_chars()))
  {
    cerr << "Cannot open R file '" << rfile << "' for transformed values\n";
    return 0;
  }

  return write_rfile_transformed (raw, transformed, output);
}


static int
do_linear_correction (const Command_Line & cl,
                      const IWString & buffer,
                      IWString_and_File_Descriptor & output)
{
  float fit[2];

  if (! parse_fit_data(buffer, fit, 1))
    return 0;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! do_correction(cl[i], fit, 1, output))
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return 0;
    }
  }

  if (rfile.length())
    write_rfile_transformed(raw, transformed, rfile);

  return 1;
}

static int
do_quadratic_correction (const Command_Line & cl,
                         const IWString & buffer,
                         IWString_and_File_Descriptor & output)
{
  float fit[3];

  if (! parse_fit_data(buffer, fit, 2))
    return 0;


  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! do_correction(cl[i], fit, 2, output))
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
do_correction (const Command_Line & cl,
               iwstring_data_source & input,
               IWString_and_File_Descriptor & output)
{
  IWString buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read info from unscaling file (-U)\n";
    return 0;
  }
  
  if (buffer.starts_with(QUADRATIC_MODEL))
    return do_quadratic_correction (cl, buffer, output);
  else if (buffer.starts_with(LINEAR_MODEL))
    return do_linear_correction (cl, buffer, output);
  else
  {
    cerr << "Unrecognised fit type '" << buffer << "' (-U option)\n";
    return 0;
  }

  return 1;
}

static int
do_correction (const Command_Line & cl,
               const char * ufile,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(ufile);

  if (! input.good())
  {
    cerr << "Cannot open -U file '" << ufile << "'\n";
    return 0;
  }

  return do_correction (cl, input, output);
}
static int
create_r_file2 (const resizable_array<float> & activity,
                const resizable_array<float> & predicted,
                const float * fit,
                int degree,
                IWString_and_File_Descriptor & output)
{
  if (rfile_plotfile.length() > 0)
    write_r_plotting_instructions (rfile_plotfile, output);

  int n = activity.number_elements();

  float a = activity[0];
  float amin = a;
  float amax = a;
  output << "obs = c(" << a;
  for (int i = 1; i < n; i++)
  {
    a = activity[i];

    output << ',' << a;
    output.write_if_buffer_holds_more_than(8192);

    if (a > amax)
      amax = a;
    else if (a < amin)
      amin = a;
  }
  output << ")\n";

  float p = predicted[0];
  float pmin = p;
  float pmax = p;

  output << "pred = c(" << p;
  for (int i = 1; i < n; i++)
  {
    p = predicted[i];

    output << ',' << p;
    output.write_if_buffer_holds_more_than(8192);

    if (p > pmax)
      pmax = p;
    else if (p < pmin)
      pmin = p;
  }
  output << ")\n";

  float zmin = amin;
  if (pmin < zmin)
    zmin = pmin;

  float zmax = amax;
  if (pmax > zmax)
    zmax = pmax;

  output << "plot(obs, pred,xlim=c(" << zmin << ',' << zmax << "), ylim=c(" << zmin << ',' << zmax << "), xlab='Obs', ylab='" << ylab << "',pch=4,cex=0.4,col='blue', main='" << current_file_name << "')\n";

  output << "fit=c(";
  if (1 == degree)
  {
    for (int i = 0; i < n; i++)
    {
      float x = fit[0] + fit[1] * activity[i];
      if (i > 0)
        output << ',';

      output << x;
      output.write_if_buffer_holds_more_than(8192);
    }
    output << ")\n";
    output << "mtext('linear: intercept " << fit[0] << " slope " << fit[1] << "', side=3, cex=0.8,col='black')\n";
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      float a = activity[i];
      float x = fit[0] + fit[1] * a + fit[2] * a * a;
      if (i > 0)
        output << ',';

      output << x;
      output.write_if_buffer_holds_more_than(8192);
    }
    output << ")\n";
    output << "mtext('quadratic " << fit[0] << ' ' << fit[1] << ' ' << fit[2] << "', side=3, cex=0.8,col='black')\n";
  }

  output << "lines(obs,fit,col='red',lwd=3)\n";

  time_t tnow = ::time(NULL);
  struct tm * tm = localtime(&tnow);

  output << "mtext(\"" << (tm->tm_year+1900) << '-' << month[tm->tm_mon] << '-' << (tm->tm_mday < 10 ? "0" : "") << tm->tm_mday << ' ' << tm->tm_hour << ':' << (tm->tm_min < 10 ? "0" : "") << tm->tm_min << ':' << (tm->tm_sec < 10 ? "0" : "") << tm->tm_sec << "\", col='darkgreen', side=4, cex=0.8)\n";

  return 1;
}

static int
create_r_file (const resizable_array<float> & activity,
               const resizable_array<float> & predicted,
               const float * fit,
               const int degree,
               const char * rfile)
{
  IWString_and_File_Descriptor output;
  
  if (! output.open(rfile))
  {
    cerr << "Cannot create R file '" << rfile << "'\n";
    return 0;
  }

  return create_r_file2 (activity, predicted, fit, degree, output);
}

static int
write_diff_vs_activity (const resizable_array<float> & activity,
                        const resizable_array<float> & predicted,
                        IWString_and_File_Descriptor & output)
{
  int n = activity.number_elements();

  output << "Obs Obs-Pred\n";

  for (int i = 0; i < n; i++)
  {
    output << activity[i] << output_separator << (activity[i]-predicted[i]) << "\n";
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
do_linear (const resizable_array<float> & activity,
           const resizable_array<float> & predicted,
           float & slope,
           float & intercept)
{
  Accumulator<double> acc_a, acc_d;

  int n = activity.number_elements();

  double product = 0.0;

  for (int i = 0; i < n; i++)
  {
    double a = activity[i];
    double d = predicted[i];

    acc_a.extra(a);
    acc_d.extra(d);

    product += a*d;
  }

  slope = (product - acc_a.sum()*acc_d.sum()/n) / (acc_a.sum_of_squares() - acc_a.sum()*acc_a.sum()/n);

  intercept = acc_d.average() - slope * acc_a.average();

  if (0 == rfile.length())
    return 1;

  float c[2];
  c[0] = intercept;
  c[1] = slope;

  return create_r_file(activity, predicted, c, 1, rfile.null_terminated_chars());
}

static int
do_linear (const resizable_array<float> & activity,
           const resizable_array<float> & predicted,
           IWString_and_File_Descriptor & output)
{
  float slope, intercept;

  do_linear(activity, predicted, slope, intercept);

  output << LINEAR_MODEL << ' ' << intercept << ' ' << slope << "\n";

  return 1;
}

static int
do_quadratic (const resizable_array<float> & activity,
              const resizable_array<float> & predicted,
              IWString_and_File_Descriptor & output)
{
  int ndata = activity.number_elements();
  int npoly = 2;
  float * c = new float[npoly+1]; std::unique_ptr<float> free_c(c);

//lsq2_(activity.rawdata(), diff.rawdata(), &ndata, &npoly, c);

  float sumx = 0.0;
  float sumx2 = 0.0;
  float sumx3 = 0.0;
  float sumx4 = 0.0;

  float sumdi = 0.0;
  float sumditi = 0.0;
  float sumditi2 = 0.0;

  for (int i = 0; i < ndata; i++)
  {
    float a = activity[i];

    sumx  +=a;
    sumx2 += a*a;
    sumx3 += a*a*a;
    sumx4 += a*a*a*a;

//  float d = (activity[i]-predicted[i]);
    float d = predicted[i];

    sumdi += d;
    sumditi += d * a;
    sumditi2 += d * a * a;
  }

  float aa[9];
  aa[0] = ndata;
  aa[1] = sumx;
  aa[2] = sumx2;
  aa[3] = sumx;
  aa[4] = sumx2;
  aa[5] = sumx3;
  aa[6] = sumx2;
  aa[7] = sumx3;
  aa[8] = sumx4;


  c[0] = sumdi;
  c[1] = sumditi;
  c[2] = sumditi2;
  
  int neq = 3;
  int info;

#ifdef DEBUG_QUADRATIC_FIT
  for (int i = 0; i < 9; i++)
  {
    cerr << "Before sgefa aa_ i = " << i << " " << aa[i] << endl;
  }
#endif

  int ipvt[3];

  sgefa_(aa,&neq,&neq,ipvt,&info);

#ifdef DEBUG_QUADRATIC_FIT
  cerr << "From sgefa_ info " << info << "\n";
  for (int i = 0; i < 9; i++)
  {
    cerr << "After sgefa_ aa " << i << " " << aa[i] << endl;
  }

  for (int i = 0; i < 3; i++)
  {
    cerr << "C before sgesl_ i " << i << ' ' << c[i] << endl;
  }
#endif

  int job = 0;
  sgesl_(aa,&neq,&neq,ipvt,c ,&job);

#ifdef DEBUG_QUADRATIC_FIT
  for (int i = 0; i < 9; i++)
  {
    cerr << " After sgesl_ aa i = " << i << " aa " << aa[i] << endl;
  }
#endif

  output << QUADRATIC_MODEL << ' ' << c[0] << ' ' << c[1] << ' ' << c[2] << "\n";

  if (0 == rfile.length())
    return 1;

  return create_r_file(activity, predicted, c, 2, rfile.null_terminated_chars());
}

#ifdef OLD_STUFF_ASDASFDG

static int
parse_input_record (const const_IWSubstring & buffer,
                    resizable_array<float> & activity,
                    resizable_array<float> & predicted)
{
  const_IWSubstring token;
  int i = 0;

  if (buffer.nwords() < 3)
  {
    cerr << "Input file records must contain at least three tokens\n";
    return 0;
  }

  (void) buffer.nextword(token, i);
  (void) buffer.nextword(token, i);

  float a;
  if (! token.numeric_value(a))
  {
    cerr << "Invalid activity value '" << token << "'\n";
    return 0;
  }

  while (buffer.nextword(token, i))
  {
    float v;
    if (! token.numeric_value(v))
    {
      cerr << "Invalid numeric predicted value '" << token << "'\n";
      return 0;
    }

    activity.add(a);
    predicted.add(v);
  }

  return 1;
}
#endif

/*
  We try to preserve the most active molecules
*/

typedef std::pair<float, float> Activity_Predicted;

class Activity_Predicted_Compare
{
  private:
  public:
    bool operator() (const Activity_Predicted & ap1, const Activity_Predicted & ap2) const { return ap1.first < ap2.first;}
};

/*
  We are either given a cutoff value, or we just assume half the range.
  We put the upper half of the data into buckets.
  then sample the same number of points into the lower half buckets
*/

static int
do_take_stratified_sample (resizable_array<float> & activity,
                           resizable_array<float> & predicted)
{
  const int n = activity.number_elements();

  float * na = new float[n]; std::unique_ptr<float> free_na(na);
  float * np = new float[n]; std::unique_ptr<float> free_np(np);

  Activity_Predicted * ap = new Activity_Predicted[n]; std::unique_ptr<Activity_Predicted> free_ap(ap);

  for (int i = 0; i < n; i++)
  {
    ap[i].first  = activity[i];
    ap[i].second = predicted[i];
  }

  Activity_Predicted_Compare apc;

  std::sort(ap, ap + n, apc);

//#define CHECK_SORTINGJ
#ifdef CHECK_SORTINGJ
  for (int i = 0; i < 100; i++)
  {
    cerr << " " << ap[i].first << endl;
  }
#endif

  float minval = ap[0].first;

  float dx = (ap[n-1].first - minval) / static_cast<float>(take_stratified_sample);

  if (0.0f == dx)
  {
    cerr << "Cannot take stratified sample, activity data is constant!!\n";
    return 0;
  }

  if (MIDPOINT_NOT_SET == midpoint)
    midpoint = (minval + ap[n-1].first) * 0.5;

  int points_above_midpoint = -1;
  int midpoint_bucket = -1;

  for (int i = 0; i < n; i++)
  {
    if (ap[i].first < midpoint)
      continue;

    points_above_midpoint = n-i;
    midpoint_bucket = static_cast<int>((ap[i].first - minval) / dx + 0.4999f);

    break;
  }

  if (points_above_midpoint <= 1)
  {
    cerr << "None or too few activity values above " << midpoint << endl;
    return 0;
  }

  if (n == points_above_midpoint)
  {
    cerr << "All points above midpoint " << midpoint << ", stratified sample not done\n";
    cerr << "Activity range " << minval << " to " << ap[n-1].first << endl;
    return 0;
  }

// we now need to find an equal number of stratified sample of points below the midpoint
// would be better to select randomly from each bucket, but the bias introduced here should be small

  int * items_in_bucket = new_int(take_stratified_sample); std::unique_ptr<int> free_items_in_bucket(items_in_bucket);

// We want the same number of points below as above the midpoint

  int points_per_bucket = points_above_midpoint / (take_stratified_sample - midpoint_bucket) + 1;

  if (verbose)
    cerr << "From " << n << " points, found " << points_above_midpoint << " above midpoint. Expect " << points_per_bucket << " items in each bucket\n";

  activity.resize_keep_storage(0);
  predicted.resize_keep_storage(0);

  for (int i = 0; i < (n - points_above_midpoint); i++)
  {
    int b = static_cast<int>((ap[i].first - minval) / dx + 0.4999f);

    if (items_in_bucket[b] >= points_per_bucket)
      continue;

    activity.add(ap[i].first);
    predicted.add(ap[i].second);

    items_in_bucket[b]++;
  }

  for (int i = (n - points_above_midpoint); i < n; i++)
  {
    activity.add(ap[i].first);
    predicted.add(ap[i].second);
  }

  if (verbose)
    cerr << "Found " << points_above_midpoint << " activity values above midpoint " << midpoint << " from " << n << " values sample " << activity.number_elements() << " values\n";

  return 1;
}

static int
prediction_bias (const IW_STL_Hash_Map<IWString, float> & obs,
                 const IW_STL_Hash_Map<IWString, Pred *> & pred,
                 IWString_and_File_Descriptor & output)
{
  auto n = pred.size();

  resizable_array<float> activity(n);
  resizable_array<float> predicted(n);

  if (n < 2)
  {
    cerr << "Not enough data\n";
    return 0;
  }

  for (auto i = pred.begin(); i != pred.end(); ++i)
  {
    const IWString & id = (*i).first;

    auto f = obs.find(id);

    float a = (*f).second;

    if (IGNORED_VALUE == a)
      continue;

    (*i).second->fill_arrays(a, activity, predicted);
  }

  if (verbose)
    cerr << "Identified " << activity.number_elements() << " obs/pred pairs\n";

  IW_STL_Hash_Map<IWString, int> id_to_ndx;

  int ndx = 0;
  for (auto i = pred.begin(); i != pred.end(); ++i)
  {
    id_to_ndx[(*i).first] = ndx;
    ndx++;
  }

  if (stream_for_diff_vs_activity.is_open())
    write_diff_vs_activity(activity, predicted, stream_for_diff_vs_activity);

  if (take_stratified_sample)
    do_take_stratified_sample(activity, predicted);

  if (linear)
    return do_linear (activity, predicted, output);
  else if (quadratic)
    return do_quadratic (activity, predicted, output);
  else
    return 0;
}

static int
read_predicted_values_record (const const_IWSubstring & buffer,
                       IW_STL_Hash_Map<IWString, Pred *> & pred)
{
  IWString id;
  int i = 0;

  if (! buffer.nextword(id, i))
    return 0;

  const_IWSubstring token;

  float a;

  if (! buffer.nextword(token, i) || ! token.numeric_value(a))
  {
    cerr << "Invalid predicted value '" << buffer << "'\n";
    return 0;
  }

  auto f = pred.find(id);

  if (f != pred.end())
    (*f).second->extra(a);
  else
    pred[id] = new Pred(a);

  return 1;
}

static int
read_predicted_values (iwstring_data_source & input,
                       IW_STL_Hash_Map<IWString, Pred *> & pred)
{
  const_IWSubstring buffer;

  if (! input.next_record (buffer))
  {
    cerr << "Cannot read header record from predicted file\n";
    return 0;
  }

  while (input.next_record (buffer))
  {
    if (! read_predicted_values_record (buffer, pred))
    {
      cerr << "Cannot process predicted values '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_predicted_values (const char * fname,
                       IW_STL_Hash_Map<IWString, Pred *> & pred)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open predicted values file '" << fname << "'\n";
    return 0;
  }

  return read_predicted_values (input, pred);
}

static int
read_predicted_values_from_list_of_files (iwstring_data_source & input,
                                          IW_STL_Hash_Map<IWString, Pred *> & pred)
{
  IWString fname;

  while (input.next_record(fname))
  {
    if (! read_predicted_values(fname.null_terminated_chars(), pred))
    {
      cerr << "Cannot read predicted values from '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_predicted_values_from_list_of_files (const char * fname,
                                          IW_STL_Hash_Map<IWString, Pred *> & pred)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open list of files '" << fname << "'\n";
    return 0;
  }

  return read_predicted_values_from_list_of_files (input, pred);
}

static int
read_activity_data_record (const const_IWSubstring & buffer,
                           IW_STL_Hash_Map<IWString, float> & activity)
{
  IWString id;
  int i = 0;

  if (! buffer.nextword(id, i))
    return 0;

  const_IWSubstring token;

  float a;

  if (! buffer.nextword(token, i) || ! token.numeric_value(a))
  {
    cerr << "Missing or invalid activity value\n";
    return 0;
  }

  auto f = activity.find(id);

  if (f != activity.end())
  {
    cerr << "Duplicate activity value for '" << id << "', ignored\n";
  }
  else
    activity[id] = a;

  return 1;
}

static int
read_activity_data (iwstring_data_source & input,
                    IW_STL_Hash_Map<IWString, float> & activity)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record from activity file\n";
    return 0;
  }

  int i = 0;
  if (! buffer.nextword(ylab, i) || ! buffer.nextword(ylab, i))
  {
    cerr << "Experimental value files must contain at least two tokens '" << buffer << "' invalid\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    if (! read_activity_data_record (buffer, activity))
    {
      cerr << "Cannot process activity data record '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_activity_data (const char * fname,
                    IW_STL_Hash_Map<IWString, float> & activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open activity file '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_activity_data(input, activity);
}

static int
prediction_bias (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vqD:h:R:U:p:is:m:A:F:Y:d:cC");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('q'))
  {
    quadratic = 1;
    linear = 0;

    if (verbose)
      cerr << "Will do a quadratic fit\n";
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "The header records to skip option (-h) must be a non negative number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will skip " << header_records_to_skip << " header records\n";
  }

  if (cl.option_present('i'))
  {
    replace_existing_predicted_values = 0;

    if (verbose)
      cerr << "Will insert an extra column of corrected results\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');

    if (! stream_for_diff_vs_activity.open(d))
    {
      cerr << "Cannot open stream for activity/difference pairs '" << d << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Will write activity difference pairs to '" << d << "'\n";
  }

  if (cl.option_present('R'))
  {
    cl.value('R', rfile);

    if (verbose)
      cerr << "Will create plotting instructions in '" << rfile << "'\n";

    if (cl.option_present('d'))
    {
      cl.value('d', rfile_plotfile);
    }
  }

  if (cl.option_present('c'))
  {
    drop_lower_expt = 1;
    if (verbose)
      cerr << "All experimental values at lowest value will be dropped\n";
  }

  if (cl.option_present('C'))
  {
    drop_upper_expt = 1;
    if (verbose)
      cerr << "All experimental values at upper value will be dropped\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', take_stratified_sample) || take_stratified_sample < 2)
    {
      cerr << "The number of stratified samples (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will take a stratified sample of the input, " << take_stratified_sample << " samples\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', midpoint))
    {
      cerr << "The midpoint value must be a valid float\n";
      usage(2);
    }

    if (verbose)
      cerr << "Equal sized stratified sample taken around activity value " << midpoint << endl;

    if (0 == take_stratified_sample)
      take_stratified_sample = 100;
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('U'))
  {
    if (cl.option_present('p'))
    {
      if (! cl.value('p', prediction_column) || prediction_column < 1)
      {
        cerr << "The prediction column (-p) must be a valid column number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Predictions in column " << prediction_column << endl;

      prediction_column--;
    }

    const char * u = cl.option_value('U');

    if (!  do_correction(cl, u, output))
      return 1;

    return 0;
  }

// We are perceiving bias

  IW_STL_Hash_Map<IWString, float> obs;

  if (! cl.option_present('A'))
  {
    cerr << "Must specify activity file via the -A option\n";
    usage(2);
  }
  else
  {
    const char * fname = cl.option_value('A');

    if (! read_activity_data(fname, obs))
    {
      cerr << "Cannot read activity data from '" << fname << "'\n";
      return 0;
    }

    if (drop_lower_expt || drop_upper_expt)
    {
      int tmp = replace_lower_upper_values(obs, drop_lower_expt, drop_upper_expt);
      if (verbose)
        cerr << "Dropped " << tmp << " experimental values at either upper or lower ends of experimental range\n";
    }
  }

  if (cl.option_present('Y'))
  {
    cl.value('Y', ylab);

    if (verbose)
      cerr << "R plot will appear with Y label '" << ylab << "'\n";
  }

  IW_STL_Hash_Map<IWString, Pred *> pred;

  if (cl.option_present('F'))
  {
    const char * fname = cl.option_value('F');
    if (! read_predicted_values_from_list_of_files(fname, pred))
    {
      cerr << "Cannot process list of files in '" << fname << "'\n";
      return 0;
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! read_predicted_values(cl[i], pred))
      {
        cerr << "Cannot read predicted values from '" << cl[i] << "'\n";
        return i+1;
      }
    }
  }

// Now make sure we have activity data for all the predicted data

  for (auto i = pred.begin(); i != pred.end(); ++i)
  {
    if (! obs.contains((*i).first))
    {
      cerr << "Predicted data for '" << (*i).first << "' but no experimental data\n";
      return 2;
    }
  }

  if (verbose)
    cerr << "Read data for " << obs.size() << " measured values and " << pred.size() << " predicted values\n";

  if (! prediction_bias(obs, pred, output))
  {
    cerr << "Fatal error processing bias\n";
    return 2;
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

  int rc = prediction_bias(argc, argv);

  return rc;
}
