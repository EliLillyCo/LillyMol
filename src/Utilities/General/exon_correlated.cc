/*
  Scans a exon file from Rick Higgs looking for correlated exons
*/

#include <stdlib.h>
#include <random>
#include <algorithm>
#include <limits>
#include <memory>
#include <cmath>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"

using std::cerr;
using std::endl;
using std::numeric_limits;
using std::unique_ptr;

const char * prog_name = NULL;

static int verbose = 0;

static char input_separator = '\t';
static char output_separator = ' ';

static int header_records_to_skip = 1;

static IWString * column_name = nullptr;

static Report_Progress_Template<unsigned long int> report_progress;

static int nwrite = numeric_limits<int>::max();

static int compare_exons_from_same_gene = 1;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Looks for correlations in n exon file\n";
  cerr << " -r <nrow>      number of rows to process\n";
  cerr << " -c <ncol>      number of columns to process\n";
  cerr << " -i <space,tab,comma>  input  word separator\n";
  cerr << " -o <space,tab,comma>  output word separator\n";
  cerr << " -I <string>    exon identifier is in column with name <string>\n";
  cerr << " -G <string>    gene identifier is in column with name <string>\n";
  cerr << " -F <string>    exon data starts in column with header starting with <string>\n";
  cerr << " -t <r2>        discard all computed correlations below <r2>\n";
  cerr << " -p <n>         report progress every <n> computations done\n";
  cerr << " -w <n>         only write the <n> highest correlations found\n";
  cerr << " -x             do NOT compare exons that come from the same gene\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Column_Names_and_Numbers
{
  private:
    IWString _identifier_column_name;
    IWString _first_column_name;
    IWString _gene_identifier_column_name;

    int _identifier_column;
    int _first_column;
    int _gene_identifier_column;

  public:
    Column_Names_and_Numbers ();

    void set_identifier_column_name (const const_IWSubstring & s) { _identifier_column_name = s;}
    void set_first_column_name (const const_IWSubstring & s) { _first_column_name = s;}
    void set_gene_identifier_column_name (const const_IWSubstring & s) { _gene_identifier_column_name = s;}

    int identify_columns (const IWString * header, int n);

    int identifier_column () const { return _identifier_column;}
    int first_column () const { return _first_column;}
    int gene_identifier_column () const { return _gene_identifier_column;}
};

Column_Names_and_Numbers::Column_Names_and_Numbers ()
{
  _identifier_column = -1;
  _first_column = -1;
  _gene_identifier_column = -1;

  return;
}

static int
name_can_be_converted_to_column (const IWString & s,
                                 int & c)
{
  if (! s.numeric_value(c) || c < 1)
  {
    cerr << "Column_Names_and_Numbers::identify_columns:no match for '" << s << "'\n";
    return 0;
  }

  c--;    // convert to array index

  return 1;
}

int
Column_Names_and_Numbers::identify_columns (const IWString * header, int ncols)
{
  for (auto col = 0; col < ncols; ++col)
  {
    const IWString & s = header[col];

    if (_identifier_column < 0 && s == _identifier_column_name)
      _identifier_column = col;
    else if (_first_column < 0 && s == _first_column_name)
      _first_column = col;
    else if (_gene_identifier_column < 0 && s == _gene_identifier_column_name)
      _gene_identifier_column = col;
  }

  if (_identifier_column >= 0)
    ;
  else if (name_can_be_converted_to_column(_identifier_column_name, _identifier_column))
    ;
  else
    return 0;

  if (_first_column >= 0)
    ;
  else if (name_can_be_converted_to_column(_first_column_name, _first_column))
    ;
  else
    return 0;

  if (_gene_identifier_column >= 0)
    ;
  else if (name_can_be_converted_to_column(_gene_identifier_column_name, _gene_identifier_column))
    ;
  else
    return 0;

  return 1;
}

template <typename T>
class Association_Result
{
  private:
    int _id1;
    int _id2;
    T   _v;              // initially R2

  public:
    Association_Result ();
    Association_Result (int i1, int i2, T v) : _id1(i1), _id2(i2), _v(v) {};

    void set_id1 (int s) { _id1 = s;}
    void set_id2 (int s) { _id2 = s;}
    void set_value (T v) { _v = v;}

    void set (int id1, int id2, T v) { _id1 = id1; _id2 = id2; _v = v;}

    int id1 () const { return _id1;}
    int id2 () const { return _id2;}
    T    v () const { return _v;}
};

template <typename T>
Association_Result<T>::Association_Result ()
{
  _id1 = -1;
  _id2 = -1;
  _v = static_cast<T>(0.0);

  return;
}

static std::random_device rd;

static void
choose_rows (int * rows_to_test,
             const int nrows,
             const int points_to_test,
             const int d,
             const int p)
{
  std::fill(rows_to_test, rows_to_test + nrows, d);

  std::uniform_int_distribution<int> u(0, nrows);
  std::mt19937_64 rng(rd());

  int rows_selected = 0;

  while (rows_selected < points_to_test)
  {
    auto r = u(rng);

    if (p == rows_to_test[r])
      continue;

    rows_to_test[r] = p;
    rows_selected++;
  }

  return;

  int nsel = 0;
  for (auto i = 0; i < nrows; ++i)
  {
    if (rows_to_test[i])
    {
      cerr << points_to_test << " selected row " << i << endl;
      nsel++;
    }
  }

  return;
}

static void
choose_rows (int * rows_to_test,
             const int nrows,
             const int points_to_test)
{
  if (points_to_test < nrows / 2)
    choose_rows(rows_to_test, nrows, points_to_test, 0, 1);
  else
    choose_rows(rows_to_test, nrows, nrows - points_to_test, 1, 0);
}

template <typename T>
T
final_r2_computation (const int n,
                      const T r,
                      const T sumc1,
                      const T sumc2,
                      const T sum_squares_c1,
                      const T sum_squares_c2)
{
  const T c1bar = sumc1 / static_cast<T>(n);
  const T c2bar = sumc2 / static_cast<T>(n);

  T nx1bx2b = static_cast<T> (n) * c1bar * c2bar;

  T v1 = sum_squares_c1 - static_cast<T>(n) * c1bar * c1bar;
  T v2 = sum_squares_c2 - static_cast<T>(n) * c2bar * c2bar;

  if (static_cast<T>(0.0) == v1 || static_cast<T>(0.0) == v2)
    return static_cast<T>(0.0);

  T rho = (r - nx1bx2b) / sqrt (v1 * v2);

  return rho;
}

template <typename T>
float
do_test (const int nrows,
         const int ncols,
         const T * d,
         const int c1,
         const int c2,
         int points_to_test,   // the number of rows we will sample
         int * rows_to_test)
{

  choose_rows(rows_to_test, nrows, points_to_test);

  T r = 0.0f;
  T sumf1 = 0.0;
  T sumf2 = 0.0;
  T sum_squares_f1 = 0.0;
  T sum_squares_f2 = 0.0;

  for (auto i = 0; i < nrows; ++i)
  {
    if (0 == rows_to_test[i])
      continue;

    const T f1 = d[i * ncols + c1];
    const T f2 = d[i * ncols + c2];

    r += f1 * f2;
//  cerr << "Added " << f1 << " and " << f2 << endl;
    sumf1 += f1;
    sumf2 += f2;
    sum_squares_f1 += f1 * f1;
    sum_squares_f2 += f2 * f2;
  }

  return final_r2_computation(points_to_test, r, sumf1, sumf2, sum_squares_f1, sum_squares_f2);
}

#ifdef COLUMN_ORDER
template <typename T>
T
rsquared (const int nrows,
          const int ncols,
          const T * d,
          const int c1, 
          const int c2)
{
  T r = 0.0;
  T sumc1 = 0.0;
  T sumc2 = 0.0;
  T sum_squares_c1 = 0.0;
  T sum_squares_c2 = 0.0;

  for (auto i = 0; i < nrows; ++i)
  {
    const T f1 = d[i * ncols + c1];
    const T f2 = d[i * ncols + c2];

    r += f1 * f2;
//  cerr << "f1 " << f1 << " f2 " << f2 << " r " << r << endl;

    sumc1 += f1;
    sumc2 += f2;

    sum_squares_c1 += f1 * f1;
    sum_squares_c2 += f2 * f2;
  }

  return final_r2_computation (nrows, r, sumc1, sumc2, sum_squares_c1, sum_squares_c2);
}
#endif

template <typename T>
T
rsquared (const T * d,
          const int ncols,
          const int r1, 
          const int r2)
{
  T r = 0.0;
  T sumc1 = 0.0;
  T sumc2 = 0.0;
  T sum_squares_c1 = 0.0;
  T sum_squares_c2 = 0.0;

  const T * p1 = d + r1 * ncols;
  const T * p2 = d + r2 * ncols;

  for (auto i = 0; i < ncols; ++i)
  {
    const T f1 = *p1;
    const T f2 = *p2;

    p1++;
    p2++;

    r += f1 * f2;
//  cerr << "f1 " << f1 << " f2 " << f2 << " r " << r << endl;

    sumc1 += f1;
    sumc2 += f2;

    sum_squares_c1 += f1 * f1;
    sum_squares_c2 += f2 * f2;
  }

  return final_r2_computation (ncols, r, sumc1, sumc2, sum_squares_c1, sum_squares_c2);
}

template <typename T>
int
do_test (const int nrows,
         const int ncols,
         const T * d,
         const int pairs_to_test,
         const std::vector<int> & points_to_test)
{
  auto max_points = *(std::max_element(points_to_test.begin(), points_to_test.end()));

  int * working_array = new int[max_points + max_points]; unique_ptr<int> free_working_array(working_array);

  std::uniform_int_distribution<int> u(0, ncols - 1);

  std::mt19937_64 rng(rd());

  for (auto i = 0; i < pairs_to_test; ++i)
  {
    auto c1 = u(rng);
    auto c2 = u(rng);
    while (c1 == c2)
    {
      c2 = u(rng);
    }

    const T r2 = rsquared(nrows, ncols, d, c1, c2);

    cerr << "Between columns " << c1;
    if (nullptr != column_name)
      cerr << ' ' << column_name[c1];
    cerr << " and " << c2;
    if (nullptr != column_name)
      cerr << ' ' << column_name[c2];
    cerr << " rho " << r2 << endl;

    for (auto p = points_to_test.begin(); p != points_to_test.end(); ++p)
    {
      cerr << *p << " rows: " << r2;
      for (auto j = 0; j < 10; ++j)
      {
        const T rp = do_test(nrows, ncols, d, c1, c2, *p, working_array);
        cerr << ' ' << rp;
      }
      cerr << endl;
    }
  }

  return 1;
}

template <typename T>
int
correlated_columns (const int nrows,
                    const int ncols,
                    const T * d,
                    IWString_and_File_Descriptor & output)
{
  return 1;
}

template <typename T>
class  ASR_Sorter
{
  private:
  public:
    int operator () (const Association_Result<T> * a1, const Association_Result<T> * a2) const
    {
      if (a1->v() < a2->v())
        return -1;
      else if (a1->v() > a2->v())
        return 1;
      else
        return 0;
    }
};

//static ASR_Sorter<float> asr_sorter;

template <typename T>
int
do_output (const int nwrite,
           Association_Result<T> * r,
           int nstored,
           const IWString * id,
           const int * gene_id,
           std::ostream & output)
{
  auto n = nstored;
  if (n > nwrite)
    n = nwrite;

  for (auto i = 0; i < n; ++i)
  {
    const Association_Result<T> & ri = r[i];

    output << gene_id[ri.id1()] << output_separator <<
              id[ri.id1()] << output_separator <<
              gene_id[ri.id2()] << output_separator <<
              id[ri.id2()] << output_separator <<
              ri.v() << "\n";
  }

  return 1;
}

template <typename T>
void
discard_low_numbers (Association_Result<T> * r,
                     int n,
                     const int nkeep)
{
  std::random_shuffle(r, r + n);

  std::partial_sort (r, r + nkeep, r + n, [] (const Association_Result<T> & a1, const Association_Result<T> & a2) { return fabs(a1.v()) > fabs(a2.v());});

  return;
}

template <typename T>
int
exon_correlated (const T * d,
                 const int nrows,
                 const int ncols,
                 Association_Result<T> * r2,
                 T min_correlation,
                 const IWString * id,
                 const int * gene_id)
{
  //unsigned long int computations_done = 0;
  const unsigned long int todo = static_cast<unsigned long int> (nrows) * static_cast<unsigned long int>(nrows-1) / 2;

  Accumulator<T> acc;

  unsigned int nstored = 0;

  for (auto i = 0; i < nrows; ++i)
  {
    const auto gi = gene_id[i];

    for (auto j = i + 1; j < nrows; ++j)
    {
      if (report_progress())
        cerr << i << " completed " << report_progress.times_called() << " of " << todo << " computations " << static_cast<float>(report_progress.times_called()) / static_cast<float>(todo) << ", stored " << nstored << " results, max " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;

      if (! compare_exons_from_same_gene && gi == gene_id[j])
        continue;

      const T v = rsquared(d, ncols, i, j);

      if (fabs(v) < min_correlation)
        continue;

      if (verbose || report_progress.active())
        acc.extra(v);

      r2[nstored].set(i, j, v);
      nstored++;

      if (nstored == static_cast<unsigned int>(2 * nwrite))
      {
        discard_low_numbers(r2, nstored, nwrite);
        nstored = nwrite;
      }
    }
  }

  std::random_shuffle(r2, r2 + nstored);

  if (nstored > static_cast<unsigned int>(nwrite))
    std::partial_sort (r2, r2 + nwrite, r2 + nstored, [] (const Association_Result<T> & a1, const Association_Result<T> & a2) { return fabs(a1.v()) > fabs(a2.v());});
//r2.iwqsort(asr_sorter);

  if (verbose)
    cerr << "Found " << acc.n() << " associations above threshold\n";

  return do_output (nwrite, r2, nstored, id, gene_id, std::cout);
}

template <typename T>
void
echo_the_data (const int nrows,
               const int ncols,
               const T * d,
               std::ostream & output)
{
  for (auto r = 0; r < nrows; ++r)
  {
    const T * rowptr = d + r * ncols;

    for (auto c = 0; c < ncols; ++c)
    {
      if (c > 0)
        output << ' ';

      output << rowptr[c];
    }

    output << '\n';
  }

  return;
}

template <typename T>
int
read_row (const const_IWSubstring & buffer,
          T * rowptr,
          const int ncols,
          const Column_Names_and_Numbers & cnn,
          IWString & id,
          int & gene_id)
{
  const_IWSubstring token;
  int i = 0;

//cerr << "Processing row, first_column " << first_column << " identifier_column " << identifier_column << endl;

  for (auto col = 0; col < cnn.first_column(); ++col)
  {
    if (! buffer.nextword(token, i, input_separator))    // skip over identifier
      return 0;

    if (col == cnn.identifier_column())
      id = token;
    else if (col == cnn.gene_identifier_column())
    {
      if (! token.numeric_value(gene_id) || gene_id < 0)
      {
        cerr << "read_row::invalid numeric for gene id '" << token << "'\n";
        return 0;
      }
    }
  }

  for (auto col = 0; col < ncols; ++col)    // now 'col' is used for an index into the data array
  {
    if (! buffer.nextword(token, i, input_separator))
    {
      cerr << "Cannot fetch column " << col << endl;
      return 0;
    }

    if (! token.numeric_value(rowptr[col]))
    {
      cerr << "Invalid numeric detected '" << token << "' in column " << col << endl;
      return 0;
    }
  }

  return 1;
}

static int
process_header_record (const const_IWSubstring & buffer,
                       Column_Names_and_Numbers & cnn)
{
  int nw = buffer.nwords(input_separator);

  if (nw <5)
  {
    cerr << "Insufficient columns\n";
    return 0;
  }

  if (verbose)
    cerr << "Header contains " << nw << " tokens\n";

  column_name = new IWString[nw + 1];    // because the loop needs to fail

  int i = 0;

  for (auto col = 0; buffer.nextword(column_name[col], i, input_separator); col++)
  {
  }

  return cnn.identify_columns(column_name, nw);
}

template <typename T>
int
read_the_data (iwstring_data_source & input,
               int & nrows,
               int & ncols,
               T * & d,
               IWString * & ids,
               int * & gene_id,
               Column_Names_and_Numbers & cnn)
{
  const_IWSubstring buffer;

  for (auto i = 0; i < header_records_to_skip; ++i)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot read header\n";
      return 0;
    }

    if (0 == i)
    {
      if (! process_header_record (buffer, cnn))
      {
        cerr << "Cannot process header record '" << buffer << "'\n";
        return 0;
      }
    }
  }

  if (nrows < 0)
  {
    nrows = input.records_remaining();
    if (nrows < 3)
    {
      cerr << "Insufficient rows in file " << nrows << endl;
      return 0;
    }

    if (verbose)
      cerr << "File contains " << nrows << " rows\n";
  }

  ids = new IWString[nrows];
  gene_id = new int[nrows];

  for (auto r = 0; r < nrows; ++r)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Premature EOF\n";
      return 0;
    }

    if (ncols < 0)
    {
      ncols = buffer.nwords(input_separator) - cnn.first_column() + 1;
      if (ncols < 2)
      {
        cerr << "Insufficient columns in input\n";
        return 0;
      }

      ncols--;
    }

    if (nullptr == d)
    {
      if (verbose)
        cerr << "Allocating data for " << nrows << " rows x " << ncols << " columns\n";
      d = new T[nrows * ncols];
    }

    T * rowptr = d + r * ncols;

    if (! read_row(buffer, rowptr, ncols, cnn, ids[r], gene_id[r]))
    {
      cerr << "Fatal error reading line " << input.lines_read() << endl;
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
read_the_data (const char * fname, 
               int & nrows,
               int & ncols,
               T * & d,
               IWString * & ids,
               int * & gene_id,
               Column_Names_and_Numbers & cnn)
{
  iwstring_data_source input(fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_the_data (input, nrows, ncols, d, ids, gene_id, cnn);
}

static int
correlated_columns (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vc:r:i:o:t:T:I:F:p:w:G:x");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('w'))
  {
    cerr << "Must specify how many results to write via the -w option\n";
    usage(1);
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', nwrite) || nwrite < 1)
    {
      cerr << "The number of results to write (-w) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will write as many as " << nwrite << " of the best results\n";
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "Unrecognised input delimiter '" << i << "'\n";
      usage(1);
    }

    input_separator = i[0];
  }

  if (cl.option_present('o'))
  {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o))
    {
      cerr << "Unrecognised output delimiter '" << o << "'\n";
      usage(1);
    }

    output_separator = o[0];
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, don't know how to handle multiple input files\n";
    return 1;
  }

  if (cl.option_present('p'))
  {
    if (! report_progress.initialise (cl, 'p', verbose))
      return 1;
  }

  int nrows = -1;
  int ncols = -1;

  if (cl.option_present('r'))
  {
    if (! cl.value('r', nrows) || nrows < 3)
    {
      cerr << "The number of rows to process must be a whole +ve number > 2\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will process the first " << nrows << " rows of the file\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', ncols) || ncols < 2)
    {
      cerr << "The number of columns to process must be a whole +ve number > 1\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will process the first " << ncols << " columns of the file\n";
  }

  Column_Names_and_Numbers cnn;

  if (! cl.option_present('I'))
  {
    cerr << "Must specify identifier column name via the -I option\n";
    usage(1);
  }

  IWString s;
  cl.value('I', s);
  cnn.set_identifier_column_name(s);

  if (! cl.option_present('F'))
  {
    cerr << "Must specify first column that contains data via the -F option\n";
    usage(1);
  }

  cl.value('F', s);
  cnn.set_first_column_name(s);

  if (! cl.option_present('G'))
  {
    cerr << "Must specify gene identifier column via the -G option\n";
    usage(1);
  }

  cl.value('G', s);
  cnn.set_gene_identifier_column_name(s);

  float * d = nullptr;
  IWString * ids = nullptr;
  int * gene_id = nullptr;

  if (! read_the_data (cl[0], nrows, ncols, d, ids, gene_id, cnn))
  {
    cerr << "Cannot read '" << cl[0] << "'\n";
    return 2;
  }

  if (verbose)
    cerr << "Read " << nrows << " rows and " << ncols << " columns of data from '" << cl[0] << "'\n";

  if (cl.option_present('e'))
    echo_the_data (nrows, ncols, d, cerr);

  float min_correlation = static_cast<float>(0.0);

  if (cl.option_present('t'))
  {
    if (! cl.value('t', min_correlation) || min_correlation >= 1.0f)
    {
      cerr << "The minimum correlation coefficient to keep (-t) must be a valid correlation coefficient value\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will discard all correlation values below " << min_correlation << endl;
  }

  if (cl.option_present('x'))
  {
    compare_exons_from_same_gene = 0;

    if (verbose)
      cerr << "Will only compare exons that come from different genes\n";
  }

  if (cl.option_present('d'))
  {
//  resizable_array_p<Association_Result<double>> r2;
//  exon_correlated(d, nrows, ncols, r2, static_cast<double>(min_correlation));
  }
  else
  {
    Association_Result<float> * r2 = new Association_Result<float>[nwrite + nwrite];
    exon_correlated(d, nrows, ncols, r2, static_cast<float>(min_correlation), ids, gene_id);
  }

  if (verbose)
  {
  }

  if (nullptr != d)
    delete [] d;

  if (nullptr != column_name)
    delete [] column_name;

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = correlated_columns(argc, argv);

  return rc;
}
