/*
  Imposes a correlation structure on columns of numbers.
  Works by computing the correlation structure within the data.
  Then, iteratively, re-orders, randonly selected, pairs of rows,
  and measures the new correlation structure. If this has moved
  us towards the objective, retain the move.
  The efficiency comes from how the new correlation coefficient
  is calculated - by a delta operation on the previous calculation.

  Based on a Julia package that I hope to release soon.
*/

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <time.h>
using std::cerr;
using std::endl;

#define IW_TABULAR_DATA_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwmisc/iw_tabular_data.h"

const char * prog_name = nullptr;

namespace correlate {

int verbose = 0;

char input_separator = ',';
char output_separator = ',';

int global_maxiter = 0;

int maxsec = 0;

Report_Progress report_progress;

IWString fmt(" %.4f");

int print_correlation_matrix_on_completion = 1;

IWString rfile;
IWString_and_File_Descriptor stream_for_plot;

void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Enforces a specified correlation structure on separately generated data\n";
  cerr << " -n <niter>     number of iterations to try (def 100 * nrows)\n";
  cerr << " -c <cor>       two columns in input, single correlation target\n";
  cerr << " -C <fname>     multiple columns in input, correlation matrix\n";
  cerr << " -t <tol>       scalar tolerance\n";
  cerr << " -T <fname>     tolerance matrix in a file\n";
  cerr << " -i <char>      input column separator (default comma)\n";
  cerr << " -o <char>      output column separator (default same as input)\n";
  cerr << " -H             input file has a header\n";
  cerr << " -r <n>         report progress every <n> iterations\n";
  cerr << " -x <seconds>   abandon calculation after <seconds> seconds\n";
  cerr << " -S <stem>      file name stem for output\n";
  cerr << " -R <fname>     file for R plots of changes with time (2D only)\n";
  cerr << " -f <fmt>       sprintf format for printing correlation matrix (def " << fmt << ")\n";
  cerr << " -b             do NOT print the correlation matrix with verbose upon completion\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template <typename T>
int
write_R_array(const char * vname, const T * v, const int n,
              IWString_and_File_Descriptor & output)
{
  output << vname << "=c(" << v[0];

  for (int i = 2; i < n; i += 2)
  {
    output << ',' << v[i];
  }

  output << ")\n";

  return 1;
}

template <typename T>
int
write_R_data(const T * xy, const int nrows, const T mycor, 
             const int ndx,
             IWString_and_File_Descriptor & output)
{
  IWString fname = rfile;
  fname << ndx << ".png";
  output << "png('" << fname << "',height=height,width=width)\n";
  write_R_array("x", xy, nrows,output);
  write_R_array("y", xy + 1, nrows,output);
  std::stringstream tmp;
  tmp << std::setprecision(2) << mycor;
  int col = static_cast<int>(mycor * 100.0 + 0.4999);
  if (col < 1)
    col = 1;

  output << "plot(x,y,main='Cor " << tmp.str() << "',col=" << "rainbow(100,start=0.6,end=1.0)[" << col << "],pch=pch,cex=cex)\n";
  output << "dev.off()\n";
  output.write_if_buffer_holds_more_than(32768);
  return 1;
}

template <typename T>
class Row_Changer
{
  private:
    const int _ncols;
    const T * _x;
    std::function<T(T,T)> _op;

  public:
    Row_Changer(const T * x, const int ncols, std::function<T(T, T)> o) : _ncols(ncols), _x(x), _op(o) {}

    int operator() (T * v) const
    {
      for (int i = 0; i < _ncols; ++i)
      {
        v[i] = _op(v[i], _x[i]);
      }
      return 1;
    }
};

template <typename T>
void
compute_cor(const T * data,
            const int nrows,
            const int ncols,
            T * cor)
{
  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      T xy = 0.0;
      T x2 = 0.0;
      T y2 = 0.0;
      for (int r = 0; r < nrows; ++r)
      {
        const T v1 = data[r*ncols+c1];
        const T v2 = data[r*ncols+c2];
        xy += v1*v2;
        x2 += v1*v1;
        y2 += v2*v2;
      }

      cor[c1 * ncols + c2] = xy / (sqrt(x2) * sqrt(y2));
    }
  }
}

template <typename T>
void
count_not_converged(const T * cdiff, const T * tol,
                    const int nrows, const int ncols,
                    int & nc,
                    T & maxdiff)
{
  nc = 0;
  maxdiff = 0.0;

  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      const auto d = cdiff[c1 * ncols + c2];
//    cerr << "count_not_converged:col " << c1 << " and " << c2 << " d = " << d << " tol " << tol[c1*ncols+c2] << endl;
      if (d > tol[c1 * ncols + c2])
      {
        if (d > maxdiff)
          maxdiff = d;
        nc++;
      }
    }
  }

//cerr << "count_not_converged nc " << nc << endl;

  return;
}

//for c1 in 1:ncols
//  for c2 in (c1+1):ncols
//    d = cdiff[c1,c2]
//    if d > tol[c1,c2]
//      if d > maxdiff
//        maxdiff = d
//      end
//      rc += 1
//    end
//  end
//end

//return (rc, maxdiff)


#ifdef JULIA_PRINT_WS
function print_with_space(output::Any, fmt::AbstractString, v::Float64)
  if v >= 0.0
    @printf(output, " ")
  end
  @printf(output, " %.4f", v)
end
#endif

template <typename T, typename O>
void
print_with_space(O & output, const char * fmt, const T v)
{
  if (v > 0.0)
    output << ' ';

  char buffer[32];
  sprintf(buffer, fmt, v);
  output << buffer;
  
  return;
}

void
choose_two_ordered(std::uniform_int_distribution<int> & u,
                   std::mt19937_64 & rng,
                   int & row1,
                   int & row2)
{
  row1 = u(rng);
  row2 = u(rng);

  if (row1 < row2)
    return;

  if (row1 > row2)
  {
    std::swap(row1, row2);
    return;
  }

  while (1)
  {
    row2 = u(rng);
    if (row2 == row1)
      continue;

    if (row1 < row2)
      return;

    std::swap(row1, row2);
    return;
  }
}

template <typename T>
int
improves_correlations(const T * data,
                      const int nrows,
                      const int ncols,
                      const int row1, 
                      const int row2,
                      const double * target_correlation,
                      const T * xx, 
                      const T * xy,
                      const T * cdiff,
                      const T * tol)
{
  double current_diff = 0.0;
 
  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      const auto d = cdiff[c1 * ncols + c2];
      if (d > tol[c1 * ncols + c2])
        current_diff += d;
    }
  }

#ifdef JULIA_IMPL
  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      d = cdiff[c1,c2]
      if d > tol[c1,c2]
        current_diff += d
      end
    end
  end
  @printf(STDERR, "Rows %d and %d, cirrent diff %f\n", row1, row2, current_diff)
#endif

//#define DEBUG_IMPROVES_CORRELATIONS
#ifdef DEBUG_IMPROVES_CORRELATIONS
  cerr << "Rows " << row1 << " and " << row2 << " current_diff " << current_diff << endl;
#endif

  T lowest = current_diff;
  int col_to_switch = -1;
 
  for (int c = 0; c < ncols; ++c)
  {
    double sum = 0.0;
    for (int x1 = 0; x1 < ncols; ++x1)
    {
      for (int x2 = x1 + 1; x2 < ncols; ++x2)
      {
        if (x1 == c || x2 == c)
        {
          double t = xy[x1*ncols+x2] - data[row1*ncols+x1]*data[row1*ncols+x2] - data[row2*ncols+x1]*data[row2*ncols+x2] +
                                       data[row1*ncols+x2]*data[row2*ncols+x1] + data[row1*ncols+x1]*data[row2*ncols+x2];
          double tc = t / (xx[x1] * xx[x2]);
          double d = fabs(tc - target_correlation[x1*ncols+x2]);
//        cerr << "c = " << c << " cols " << x1 << " and " << x2 << " d " << d << "  tol " << tol[x1*ncols+x2] << " gt " << (d > tol[x1*ncols+x2]) << endl;
          if (d > tol[x1*ncols+x2])
            sum += d;
        }
        else    // not involved in proposed switch, count based on cdiff
        {
           const double d = cdiff[x1*ncols+x2];
           if (d > tol[x1*ncols+x2])
             sum += d;
        }
      }
    }

#ifdef DEBUG_IMPROVES_CORRELATIONS
    cerr << "switching column " << c << "? sum " << sum << " cmp " << lowest << " " << (sum < lowest) << endl;
#endif

    if (sum < lowest)
    {
      lowest = sum;
      col_to_switch = c;
    }
  }

  return col_to_switch;
}

#ifdef JULIA_IMPROVES_XX
  for c in 1:ncols     # what happens if we swap column C
    sum = 0.0
    for x1 = 1:ncols      # check all correlations involving C
      for x2 = (x1+1):ncols
        if x1 == c || x2 == c
          t = xy[x1,x2] - data[row1,x1]*data[row1,x2] - data[row2,x1]*data[row2,x2] + 
                          data[row1,x2]*data[row2,x1] + data[row1,x1]*data[row2,x2]
          tc = t /(xx[x1] * xx[x2])
//        jj =  xy[x1,x2]/(xx[x1] * xx[x2])

//        @printf(STDERR, "switching col %d, testing %d and %d, c %f current %f target %f\n", c, x1, x2, tc, jj, target_correlation[x1,x2])
          d = fabs(tc - target_correlation[x1,x2])
          if d > tol[x1,x2]
            sum += d
          end
        else
          d = cdiff[x1,x2]
          if d > tol[x1,x2]
            sum += d
          end
        end
      end
    end
#   @printf(STDERR, "switching column %d, sum %f\n", c, sum)
    if sum < lowest
      lowest = sum
      col_to_switch = c
    end
  end
 
# @printf(STDERR, "Existing diff %f, new %f, column %d\n", current_diff, lowest, col_to_switch)

  return col_to_switch
#endif


template <typename T>
bool
correlate(IW_Tabular_Data<T> & zdata,
          const IW_Tabular_Data<double> & corm,
          const IW_Tabular_Data<double> & tolm)
{
  const int nrows = zdata.nrows();
  const int ncols = zdata.ncols();

//cerr << nrows << " rows and " << ncols << " cols\n";

  const T * target_correlation = corm.data();
  const T * tol = tolm.data();

  T * xx = new T[ncols]; std::unique_ptr<T[]> free_xx(xx);
  std::fill_n(xx, ncols, 0.0);

  T * xy = new T[ncols * ncols]; std::unique_ptr<T[]> free_xy(xy);
  std::fill_n(xy, ncols * ncols, 0.0);

  T * data = zdata.data();

  for (int c = 0; c < ncols; ++c)
  {
    for (int r = 0; r < nrows; ++r)
    {
      xx[c] += (data[r * ncols + c] * data[r * ncols + c]);
    }
  }

//for c in 1:ncols
//  write(STDERR, "mean of column $(c) $(mean(data[:,c]))\n")
//  for r in 1:nrows
//    xx[c] += (data[r,c] * data[r,c])
//  end
//end

  for (int c = 0; c < ncols; ++c)
  {
    xx[c] = sqrt(xx[c]);
  }

//for c in 1:ncols
//  xx[c] = sqrt(xx[c])
//end

  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = (c1+1); c2 < ncols; ++c2)
    {
      for (int r = 0; r < nrows; ++r)
      {
        xy[c1 * ncols + c2] += (data[r * ncols + c1] * data[r * ncols + c2]);
//      cerr << " c1 " << c1 << " c2 " << c2 << " r " << r << " value " << xy[c1*ncols+c2] << endl;
      }
    }
  }

//for c1 in 1:ncols
//  for c2 in (c1+1):ncols
//    for r in 1:nrows
//      xy[c1,c2] += (data[r,c1] * data[r,c2])
//    end
//  end
//end

  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      xy[c2 * ncols + c1] = xy[c1 * ncols + c2];
    }
  }

#ifdef DEBUG_CORRELATE
  for (int i = 0; i < (ncols*ncols); ++i)
  {
    cerr << "xy[" << i << "] = " << xy[i] << endl;
  }
#endif

//for c1 in 1:ncols
//  for c2 in (c1+1):ncols
//    xy[c2,c1] = xy[c1,c2]
//  end
//end

  T * mycor = new T[ncols * ncols]; std::unique_ptr<T[]> free_mycor(mycor);

  for (int c1 = 0; c1 < ncols; ++c1)
  {
    mycor[c1 * ncols + c1] = 1.0;
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      mycor[c1 * ncols + c2] = mycor[c2 * ncols + c1] = xy[c1 * ncols + c2] / (xx[c1] * xx[c2]);
    }
  }

//#define DEBUG_CORRELATE
#ifdef DEBUG_CORRELATE
  T * junk = new T[ncols*ncols];std::unique_ptr<T[]> free_junk(junk);
  compute_cor(data, nrows, ncols, junk);
  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1+1; c2 < ncols; ++c2)
    {
      cerr << " c1 " << c1 << " c2 " << c2 << " cor " << mycor[c1*ncols+c2] << " computed " << junk[c1*ncols+c2] << endl;
    }
  }
 
  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = 0; c2 < ncols; ++c2)
    {
      cerr << ' ' << mycor[c1*ncols+c2];
    }
    cerr << endl;
  }
#endif

//for c1 in 1:ncols
//  mycor[c1,c1] = 1.0
//  for c2 in (c1+1):ncols
//    mycor[c1,c2] = mycor[c2,c1] = xy[c1,c2] / (xx[c1] * xx[c2])
//  end
//end

//write(STDERR, "My correlation $(mycor)\n")
//write(STDERR, "   correlation $(cor(data,mean=0))\n")

  T * cdiff = new T[ncols * ncols]; std::unique_ptr<T[]> free_cdiff(cdiff);    // diff between where we are and target
  std::fill_n(cdiff, ncols * ncols, 0.0);            // not necessary

  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = c1 + 1; c2 < ncols; ++c2)
    {
      if (c1 != c2)
        cdiff[c1 * ncols + c2] = fabs(target_correlation[c1 * ncols + c2] - mycor[c1 * ncols + c2]);

//    if (c1 != c2)
//      cerr << c1 << ' ' << c2 << " cdiff " << cdiff[c1 * ncols + c2] << " target " << target_correlation[c1 * ncols + c2] << " mycor " << mycor[c1 * ncols + c2] << endl;
    }
  }

//cdiff = abs(target_correlation - mycor)

  int notconverged;
  T maxdiff;

  count_not_converged(cdiff, tol, nrows, ncols, notconverged, maxdiff);

//(nc,md) = count_not_converged(cdiff, tol)

  if (0 == notconverged)
  {
    cerr << "Already converged\n";
    return true;
  }

  T next_plot = 1.0;
  T plot_delta = 0.0;
  int next_plot_index = 0;

  if (stream_for_plot.is_open())
  {
    if (2 != ncols)
    {
      cerr << "Plotting only possible with two columns, not " << ncols << endl;
      return 0;
    }

    write_R_data(data, nrows, mycor[1], next_plot_index, stream_for_plot);
    plot_delta = 0.01;
    next_plot = mycor[1] + plot_delta;
    next_plot_index = 1;
  }

  int maxiter;
  if (0 == global_maxiter)
    maxiter = 200 * nrows;
  else
    maxiter = global_maxiter;

  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> uniform_rows(0, nrows - 1);

  const auto tzero = time(NULL);

  int switches_done = 0;
  int switches_last_report = 0;
  int last_switch = 0;

  int * column_switched = new int[ncols]; std::unique_ptr<int[]> free_column_switched(column_switched);
  std::fill_n(column_switched, ncols, 0);

  bool converged = false;

  int iter = 0;

  for ( ; iter < maxiter && ! converged; ++iter)
  {
    if (report_progress())
    {
      count_not_converged(cdiff, tol, nrows, ncols, notconverged, maxdiff);
      cerr << iter << " iterations, " << switches_done << " switches, " << (switches_done - switches_last_report) << " since last rpt. " << notconverged << " not converged, max diff " << maxdiff << endl;

#ifdef REPORT_CDIFF
      for (int i = 0; i < ncols; ++i)
      {
        for (int j = 0; j < ncols; ++j)
        {
          cerr << ' ' << cdiff[i*ncols+j];;
        }
        cerr << endl;
      }
#endif

      switches_last_report = switches_done;

      if (maxsec > 0)
      {
        const auto t = time(NULL);
        if (t - tzero >= maxsec)
        {
          cerr << "Timeout after " << (t - tzero) << " seconds\n";
          return false;
        }
      }
    }

    int row1, row2;
    choose_two_ordered(uniform_rows, rng, row1, row2);

//  (row1,row2) = choose_two_ordered(nrows)

    const auto col_to_switch = improves_correlations(data, nrows, ncols, row1, row2, target_correlation, xx, xy, cdiff, tol);

    if (col_to_switch < 0 )    
    {
      if (iter - last_switch > (maxiter / 10))
      {
        cerr << "No progress in " << (iter - last_switch) << " steps\n";
        break;
      }
      continue;
    }

#ifdef DEBUG_CORRELATE
    cerr << "Rows " << row1 << " and " << row2 << ", try switching column " << col_to_switch << endl;
#endif

    column_switched[col_to_switch] += 1;

    converged = true;      // until we find otherwise

    for (int c1 = 0; c1 < ncols; ++c1)
    {
      for (int c2 = c1+1; c2 < ncols; ++c2)
      {
        if (c1 == col_to_switch || c2 == col_to_switch)
        {
          xy[c1*ncols+c2] += -data[row1*ncols+c1] * data[row1*ncols+c2] - data[row2*ncols+c1] * data[row2*ncols+c2] +
                              data[row1*ncols+c1] * data[row2*ncols+c2] + data[row1*ncols+c2] * data[row2*ncols+c1];
          mycor[c1*ncols+c2] = xy[c1*ncols+c2] / (xx[c1] * xx[c2]);
          cdiff[c1*ncols+c2] = fabs(target_correlation[c1*ncols+c2] - mycor[c1*ncols+c2]);
        }

        if (converged && (cdiff[c1 * ncols + c2] > tol[c1 * ncols + c2]))
          converged = false;

#ifdef DEBUG_CORRELATE
        cerr << iter << ' ' << converged << " cdiff " << cdiff[c1*ncols+c2] << " and tol " << tol[c1*ncols+c2] << " switches " << switches_done << endl;
#endif
      }
    }

#ifdef DEBUG_CORRELATE_CHECK_COR
    compute_cor(data, nrows, ncols, junk);
    for (int c1 = 0; c1 < ncols; ++c1)
    {
      for (int c2 = c1+1; c2 < ncols; ++c2)
      {
        cerr << " c1 " << c1 << " c2 " << c2 << " cor " << mycor[c1*ncols+c2] << " computed " << junk[c1*ncols+c2] << " diff " << (mycor[c1*ncols+c2]-junk[c1*ncols+c2]) << endl;
      }
    }
#endif

    const T tmp = data[row1*ncols+col_to_switch];
    data[row1*ncols+col_to_switch] = data[row2*ncols+col_to_switch];
    data[row2*ncols+col_to_switch] = tmp;

    switches_done++;
    last_switch = iter;

    if (2 == ncols && mycor[1] >= next_plot)
    {
      write_R_data(data, nrows, mycor[1], next_plot_index, stream_for_plot);
      next_plot += plot_delta;
      next_plot_index++;
    }

    if (converged)
      break;

//  for c1 in 1:ncols
//    for c2 in 1:ncols
//      if c1 == c2 continue end
//      if c1 == col_to_switch || c2 == col_to_switch
//        xy[c1,c2] += - data[row1,c1]*data[row1,c2] - data[row2,c1]*data[row2,c2] + 
//                       data[row1,c1]*data[row2,c2] + data[row1,c2]*data[row2,c1]
//        mycor[c1,c2] = xy[c1,c2]/(xx[c1] * xx[c2])
//        cdiff[c1,c2] = abs(target_correlation[c1,c2] - mycor[c1,c2])
//        if converged && cdiff[c1,c2] > tol[c1,c2]
//          converged = false
//        end
//      elseif cdiff[c1,c2] > tol[c1,c2]
//        converged = false
//      end
//    end
//  end

//  data[row1,col_to_switch],data[row2,col_to_switch] = data[row2,col_to_switch],data[row1,col_to_switch]
//  switches_done += 1
  }

  if (! verbose)
  {
    if (report_progress.active())
    {
      if (converged)
        cerr << "Converged";
      else
        cerr << "Not converged";
      cerr << ' ' << iter << " iterations, " << switches_done << " switches\n";
    }

    return converged;
  }

  if (converged)
    cerr << "Converged. ";
  else
    cerr << "Not converged. \n";

  cerr << "Performed " << iter << " iterations, made " << switches_done << " switches\n";

#ifdef DEBUG_CORRELATE
  for (int i = 0; i <(ncols*ncols); ++i)
  {
    cerr << "mycor[" << i << "] = " << mycor[i] << endl;
  }
#endif

  double sum_errors = 0.0;
  int not_converged = 0;
  for (int c1 = 0; c1 < ncols; ++c1)
  {
    for (int c2 = 0; c2 < ncols; ++c2)
    {
      T v;
      if (c2 > c1)
        v = mycor[c1 * ncols + c2];
      else
        v = mycor[c2 * ncols + c1];

      if (print_correlation_matrix_on_completion)
        print_with_space(cerr, fmt.null_terminated_chars(), v);
      if (c2 > c1)
      {
        const auto d = fabs(target_correlation[c1 * ncols + c2] - v);
//      std::cout << "col " << c1 << " and " << c2 << " target_correlation " << target_correlation[c1*ncols+c2] << " mycor " << v << " tol " << tol[c1*ncols+c2] << " gt? " << (d > tol[c1*ncols+c2]) << endl;
        sum_errors += d;
        if (d > tol[c1 * ncols + c2])
          not_converged += 1;
      }
    }

    if (print_correlation_matrix_on_completion)
    {
      cerr << "        ";
      for (int c2 = 0; c2 < ncols; ++c2)
      {
        print_with_space(cerr, fmt.null_terminated_chars(), target_correlation[c1 * ncols + c2]);
      }
      cerr << endl;
    }
  }
  cerr << "Sum of errors " << sum_errors << ", " << not_converged << " items not converged\n";
  for (int c = 0; c < ncols; ++c)
  {
    cerr << "Switched column " << c << " " << column_switched[c] << " times\n";
  }

  const auto t = time(NULL);
  cerr << "Calculation took " << (t - tzero) << " seconds\n";

  return converged;
}

#ifdef THE_JULIA
    for c1 in 1:ncols
      for c2 in 1:ncols
        print_with_space(STDERR, " %3f", mycor[c1,c2])
        if c2 > c1
          sum_errors += fabs(target_correlation[c1,c2] - mycor[c1,c2])
          if fabs(target_correlation[c1,c2] - mycor[c1,c2]) > tol[c1,c2]
//          @printf(STDERR, "col %d vs col %d target %f mycor %f not converged\n", c1, c2, target_correlation[c1,c2], mycor[c1,c2])
            not_converged += 1
          end
        end
      end
      @printf(STDERR, "     ")
      for c2 in 1:ncols
        print_with_space(STDERR, " %3f", target_correlation[c1,c2])
      end
      @printf(STDERR, "\n")
    end
    @printf(STDERR, "Sum of errors %f, %d items not converged\n", sum_errors, not_converged)
    for c in 1:ncols
      @printf(STDERR, "Switched column %d %d times\n", c, column_switched[c])
    end
    @printf(STDERR, "Calculation took %.2f seconds\n", time()-tzero)
  end
#endif

template <typename T>
int
correlate(IW_Tabular_Data<T> & zdata,
          IW_Tabular_Data<double> & cor,
          IW_Tabular_Data<double> & tol,
          IWString_and_File_Descriptor & output)
{
  const int nc = zdata.ncols();
  const int nr = zdata.nrows();

  if (nc != cor.ncols())
  {
    cerr << "Column count mismatch, data has " << nc << " but corelation matrix has " << cor.ncols() << " impossible\n";
    return 0;
  }

// mean shift each column

  const T * data = zdata.data();

  T * mean = new T [nc]; std::unique_ptr<T[]> free_mean(mean);

  std::fill_n(mean, nc, 0.0);

  for (int i = 0; i < nr; ++i)
  {
    const T * row = data + i * nc;

    for (int j = 0; j < nc; ++j)
    {
      mean[j] += row[j];
#ifdef DEBUG_MEAN_DETERMINATION
      cerr << " processed " << row[j] << " tot now " << mean[j] << endl;
#endif
    }
  }

  for (int i = 0; i < nc; ++i)
  {
    mean[i] = mean[i] / static_cast<double>(nr);
  }

#ifdef DEBUG_MEAN_DETERMINATION
  for (int i = 0; i < nc; ++i)
  {
    cerr << "mean[" << i << "] = " << mean[i] << endl;
  }
#endif

#ifdef ECHO_DATA
  for (int i = 0; i < (nc*nr); ++i)
  {
    cerr << "data[" << i << "] = " << data[i] << endl;
  }
#endif

  Row_Changer<T> rch_minus(mean, nc, std::minus<T>());
  zdata.change_each_row(rch_minus);

#ifdef CHECK_MEAN
  double qq = 0.0;
  for (int i = 0; i < nr; ++i)
  {
    qq += data[i * nc + 1];
  }
  cerr << "First row mean " << qq << endl;
#endif

  correlate(zdata, cor, tol);

  Row_Changer<T> rch_plus(mean, nc, std::plus<T>());
  zdata.change_each_row(rch_plus);

  return zdata.do_write(output, output_separator);
}

int
append_likely_suffix(IWString & fname,
                     const char output_separator)
{
  if (! fname.ends_with('.'))
    fname += '.';

  if (' ' == output_separator || '\t' == output_separator)
    fname << "txt";
  else if (',' == output_separator)
    fname << "csv";
  else if ('|' == output_separator)
    fname << "vbar";
  else
    fname << "dat";

  return 1;
}

template <typename T>
void
fill_matrix(IW_Tabular_Data<T> & data,
            const int ncols,
            const T diagonal,
            const T off_diagonal)
{
  data.resize(ncols, ncols);
//cerr << "matrix resized to " << data.ncols() << " and " << data.nrows() << " rows\n";

  for (int i = 0; i < ncols; ++i)
  {
    for (int j = 0; j < ncols; ++j)
    {
      if (i == j)
        data.set(i, j, diagonal);
      else
        data.set(i, j, off_diagonal);
    }
  }

  return;
}

template <typename T>
int
fill_in_cor_and_tol(const int ncols,
                    const double scalar_correlation,
                    IW_Tabular_Data<T> & cor,
                    const double scalar_tolerance,
                    IW_Tabular_Data<double> & tol)
{
//cerr << scalar_correlation << " tol " << scalar_tolerance << endl;
  if (scalar_correlation >= -1.0)
  {
    fill_matrix(cor, ncols, 1.0, scalar_correlation);
  }

  if (scalar_tolerance >= 0.0)
  {
    fill_matrix(tol, ncols, 0.0, scalar_tolerance);
  }

  return 1;
}

int
correlate (int argc, char ** argv)
{
//Command_Line_v2 cl(argc, argv, "-v-t=f-T=sfile-i=s-o=s-H-c=f-C=sfile-S=s-n=ipos-r=ipos-x=ipos");
  Command_Line cl(argc, argv, "vt:T:i:o:Hc:C:S:n:r:x:f:bR:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', global_maxiter) || global_maxiter < 1)
    {
      cerr << "The number of iterations must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will perform " << global_maxiter << " iterations\n";
  }

  if (cl.option_present('b'))
  {
    print_correlation_matrix_on_completion = 0;

    if (verbose)
      cerr << "Will NOT print the correlation matrix upon exit\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress reporter\n";
      return 1;
    }
  }

  int header_record = 0;

  if (cl.option_present('H'))
  {
    header_record = 1;

    if (verbose)
      cerr << "Input assumed to have header record\n";
  }

  if (cl.option_present('x'))
  {
    if (! cl.option_present('r'))
    {
      cerr << "In order to have a max time, you also need to set a reporting interval with the -r option\n";
      usage(1);
    }

    if (! cl.value('x', maxsec) || maxsec < 0)
    {
      cerr << "The maximum seconds to devote (-x) must be a whole non negative number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will abandon calculation after " << maxsec << " seconds\n";
  }

  if (cl.option_present('f'))
  {
    cl.value('f', fmt);

    if (verbose)
      cerr << "output format '" << fmt << "'\n";
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "Invalid input separator '" << i << "'\n";
      return 1;
    }

    input_separator = i[0];

    if (verbose)
      cerr << "Input separator set to '" << input_separator << "'\n";

    output_separator = input_separator;
  }

  if (cl.option_present('o'))
  {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o))
    {
      cerr << "Invalid output separator '" << o << "'\n";
      return 1;
    }

    output_separator = o[0];

    if (verbose)
      cerr << "Output separator set to '" << output_separator << "'\n";
  }

  IWString output_stem;

  if (cl.option_present('S'))
  {
    cl.value('S', output_stem);

    if (verbose)
      cerr << "Output created with file name stem " << output_stem << endl;
  }

  const int nfiles = cl.number_elements();

  if (0 == nfiles)
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (nfiles > 1 && 0 == output_stem.length())
  {
    cerr << "If processing multiple files, must specify -S for output file name stem\n";
    usage(1);
  }

  IW_Tabular_Data<double> cor;
  double scalar_correlation = -5.0;

  if (cl.option_present('c'))
  {
    double c;
    if (! cl.value('c', c) || c > 1.0 || c < -1.0)
    {
      cerr << "The correlation value (-c) must be between -1 and 1\n";
      return 1;
    }

    scalar_correlation = c;

    if (verbose)
      cerr << "Scalar correlation target " << scalar_correlation << endl;
  }
  else if (cl.option_present('C'))
  {
    const char * fname = cl.option_value('C');

    if (! cor.build(fname, input_separator))
    {
      cerr << "Cannot read correlation matrix '" << fname << "'\n";
      return 1;
    }
  }
  else
  {
    cerr << "Must specify desired correlation via the -c or -C option\n";
    usage(1);
  }

  IW_Tabular_Data<double> tol;

  double scalar_tolerance = -1.0;

  if (cl.option_present('t'))
  {
    if (! cl.value('t', scalar_tolerance) || scalar_tolerance < 0.0 || scalar_tolerance > 2.0)
    {
      cerr << "The tolerance (-t) must be a valid tolerance, within [0, 2]\n";
      usage(1);
    }

    if (verbose)
      cerr << "Tolerance set to " << scalar_tolerance << endl;
  }
  else if (cl.option_present('T'))
  {
    const char * fname = cl.option_value('T');
    if (! tol.build(fname, input_separator))
    {
      cerr << "Cannot read tolerance matrix '" << fname << "'\n";
      return 1;
    }
  }
  else
  {
    cerr << "Must specify desired tolerance via the -t or -T option\n";
    usage(1);
  }

  if (scalar_tolerance >= 0)
  {
    const int ncols = cor.ncols();
    tol.resize(ncols, ncols);
    for (int i = 0; i < ncols; ++i)
    {
      for (int j = 0; j < ncols; ++j)
      {
        if (i == j)
          tol.set(i, j, 1.0);
        else
          tol.set(i, j, scalar_tolerance);
      }
    }
  }

  if (cl.option_present('R'))
  {
    rfile = cl.string_value('R');

    IWString tmp(rfile);
    if (tmp.ends_with(".r") || tmp.ends_with(".R"))
      ;
    else
      tmp << ".r";

    if (! stream_for_plot.open(tmp.null_terminated_chars()))
    {
      cerr << "Cannot open stream for R plots '" << tmp << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "R plotting functions written to '" << tmp << "'\n";

    stream_for_plot << "width=600\n";
    stream_for_plot << "height=600\n";
    stream_for_plot << "pch=4\n";
    stream_for_plot << "cex=0.5\n";
  }

  int converged = 0;
  int file_number = 0;

  int rc = 0;
  for (int i = 0; i < nfiles; i++)
  {
    IW_Tabular_Data<double> zdata;
    if (header_record)
      zdata.set_has_header(1);

    if (! zdata.build(cl[i], input_separator))
    {
      cerr << "Cannot read data from '" << cl[i] << "'\n";
      return 1;
    }

    fill_in_cor_and_tol(zdata.ncols(), scalar_correlation, cor, scalar_tolerance, tol);

    if (1 == nfiles && 0 == output_stem.length())
    {
      IWString_and_File_Descriptor output(1);
      if (correlate(zdata, cor, tol, output))
        converged++;
    }
    else
    {
      IWString fname;
      fname << output_stem << file_number;
      append_likely_suffix(fname, output_separator);

      file_number++;

      IWString_and_File_Descriptor output;
      if (! output.open(fname.null_terminated_chars()))
      {
        cerr << "Cannot open '" << fname << "'\n";
        return 1;
      }
      if (correlate(zdata, cor, tol, output))
        converged++;
    }

    if (! rc)
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Processed " << nfiles << " files, " << converged << " converged\n";
  }

  if (converged == nfiles)
    return 0;

  return converged + 1;
}

}  // namespace correlate

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = correlate::correlate(argc, argv);

  return rc;
}
