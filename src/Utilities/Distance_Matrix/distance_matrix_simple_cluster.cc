/*
  Implementation of Murtagh's clustering software
*/

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <random>

#include "Foundational/cmdline/cmdline.h"

#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int iopt = 1;   // what kind of clustering to do

#include "iwhcl.h"

extern "C" void hc_(const int *,  // N
                const int *,  // LEN
                const int *,  // IOPT
                int *,        // IA
                int *,        // IB
                float *,      // CRIT
                int *,        // MEMBR
                int *,        // NN
                float *,      // DISNN
                int *,        // FLAG
                float *);     // DIS

static double resolution = 0.0;

static int nswap = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Implements Murtagh clustering functions\n";
  cerr << " -C <type>      what type of clustering to do\n";
  cerr << " -C wards       Wards\n";
  cerr << " -C slink       single linkage\n";
  cerr << " -C clink       complete linkage\n";
  cerr << " -C alink       average linkage\n";
  cerr << " -C mcq         Mcquitty's method\n";
  cerr << " -C gower       Gower's method\n";
  cerr << " -C centroid    centroid method\n";
  cerr << " -r <dist>      truncate all distances to <dist> resolution\n";
  cerr << " -s <nswap>     perform <nswap> arbitrary pair swaps to achieve randomisation\n";
  cerr << " -t ...         type of distance matrix (values 'float' and 'byte' only)\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

template <typename T>
int
do_write(IWDistanceMatrixBase<T> & dm,
         const char * fname)
{
  IWString_and_File_Descriptor output;
  if (! output.open(fname))
  {
    cerr << "do_write:cannot open '" << fname << "'\n";
    return 0;
  }

  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    const IWString & idi = dm.id(i);

    for (int j = 0; j < n; j++)
    {
      if (i == j)
        continue;

      const IWString & idj = dm.id(j);

      T v = dm.zvalue(i, j);

      output << idi << ' ' << idj << ' ' << v << '\n';
    }
  }

  return 1;
}

template int do_write<float>(IWDistanceMatrixBase<float>&, char const*);
template int do_write<unsigned char>(IWDistanceMatrixBase<unsigned char>&, char const*);

template <typename T>
int
do_swap_indices(IWDistanceMatrixBase<T> & dm,
                int nswap)
{
  int n = dm.number_molecules();

  std::random_device rd;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> u(0, n-1);

//do_write(dm, "before");
  for (int i = 0; i < nswap; i++)
  {
    int i1 = u(generator);
    int i2 = u(generator);

    if (i2 == i1)
      continue;

//  cerr << "Swapping " << i1 << " (" << dm.id(i1) << ") and " << i2 << " (" << dm.id(i2) << ")\n";
    dm.swap_items(i1, i2);
  }

//do_write(dm, "after");
  return 1;
}

template int do_swap_indices<float>(IWDistanceMatrixBase<float>&, int);
template int do_swap_indices<unsigned char>(IWDistanceMatrixBase<unsigned char>&, int);

template <typename T>
class Truncator
{
  private:
    T _resolution;

  public:
    Truncator(T r) : _resolution(r) {}

    T operator()(T) const;
};

template <typename T>
T
Truncator<T>::operator() (T v) const
{
  int i = static_cast<int> (v / _resolution + 0.4999);

  return i * _resolution;
}

template float Truncator<float>::operator()(float) const;
template unsigned char Truncator<unsigned char>::operator()(unsigned char) const;

template <typename T>
int
do_impose_resolution (IWDistanceMatrixBase<T> & dm,
                      T resolution)
{
  Truncator<T> t(resolution);

  dm.change_values(t);

  return 1;
}

template int do_impose_resolution<float>(IWDistanceMatrixBase<float>&, float);
template int do_impose_resolution<unsigned char>(IWDistanceMatrixBase<unsigned char>&, unsigned char);

template <typename T>
int
distance_matrix_simple_cluster (IWDistanceMatrixBase<T> & dm,
                                ostream & output)
{
  if (resolution > 0.0)
    do_impose_resolution(dm, static_cast<T>(resolution));

  if (nswap > 0)
    do_swap_indices(dm, nswap);

  typedef float hc_float_t;

  unsigned int n = dm.number_molecules();

  unsigned int len = n * (n - 1) / 2;

  hc_float_t * membr = new hc_float_t[n]; std::unique_ptr<hc_float_t[]> free_membr(membr);

  int * ia = new int[n]; std::unique_ptr<int[]> free_ia(ia);
  int * ib = new int[n]; std::unique_ptr<int[]> free_ib(ib);

  hc_float_t * crit = new hc_float_t[n]; std::unique_ptr<hc_float_t[]> free_crit(crit);

  unsigned int * nn = new unsigned int[n]; std::unique_ptr<unsigned int[]> free_nn(nn);

  T * disnn = new T[n]; std::unique_ptr<T[]> free_disnn(disnn);

  int * flag = new int[n]; std::unique_ptr<int[]> free_flag(flag);

  T * diss = const_cast<T *>(dm.rawdata());

  iwhcl(n, len, iopt, ia, ib, crit, membr, nn, disnn, flag, diss);

  for (unsigned int i = 0; i < n - 1; i++)
  {
    output << (ia[i] - 1) << ' ' << (ib[i] - 1) << ' ' << crit[i] << '\n';
  }

  if (verbose)
    cerr << "Clustering complete\n";

  return output.good();
}

template int distance_matrix_simple_cluster (IWDistanceMatrixBase<float> &, ostream &);
template int distance_matrix_simple_cluster (IWDistanceMatrixBase<unsigned char> &, ostream &);

template <typename T>
int
distance_matrix_simple_cluster (const char * fname,
     ostream & output)
{
  IWDistanceMatrixBase<T> dm;

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "distance_matrix_simple_cluster:cannot open '" << fname << "'\n";
    return 0;
  }

  if (1 == sizeof(T))
  {
  }

  if (! dm.do_read(fname))
  {
    cerr << "Cannot open distance matrix '" << fname << "'\n";
    return 0;
  }

  return distance_matrix_simple_cluster(dm, output);
}

template int distance_matrix_simple_cluster<float>(const char *, ostream &);
template int distance_matrix_simple_cluster<unsigned char>(const char *, ostream &);


static int
distance_matrix_simple_cluster (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vC:S:t:r:s:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present('r'))
  {
    if (! cl.value('r', resolution) || resolution <= 0.0)
    {
      cerr << "The resolution (-t) must be a positive floating point value\n";
      usage(4);
    }

    if (verbose)
      cerr << "Distances all truncated to " << resolution << " resolution\n";

    if (resolution >= 1.0)
      cerr << "Warning resolution may be out of range " << resolution << endl;
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', nswap) || nswap < 1)
    {
      cerr << "The nswap option (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will swap " << nswap << " items in the distance matrix\n";
  }

  if (cl.option_present ('C'))
  {
    const_IWSubstring c = cl.string_value ('C');

    if (c.numeric_value (iopt) && iopt >= 1 && iopt <= 7)
    {
      if (verbose)
        cerr << "Clustering method set to numeric code " << iopt << endl;
    }
    else if (c.starts_with ("ward"))
    {
      iopt = 1;
      if (verbose)
        cerr << "Wards clustering, iopt = 1\n";
    }
    else if ("slink" == c)
    {
      iopt = 2;
      if (verbose)
        cerr << "Single linkage, iopt = 2\n";
    }
    else if ("clink" == c)
    {
      iopt = 3;
      if (verbose)
        cerr << "Complete linkage, iopt = 3\n";
    }
    else if ("alink" == c)
    {
      iopt = 4;
      if (verbose)
        cerr << "Average linkage (group average), iopt = 4\n";
    }
    else if ("mcq" == c)
    {
      iopt = 5;
      if (verbose)
        cerr << "Mcquitty's method, iopt = 5\n";
    }
    else if ("gower" == c)
    {
      iopt = 6;
      if (verbose)
        cerr << "Gower's method, iopt = 6\n";
    }
    else if (c.starts_with ("cent"))
    {
      iopt = 7;
      if (verbose)
        cerr << "Centroid method, iopt = 7\n";
    }
    else
    {
      cerr << "Unrecognised clustering type '" << c << "'\n";
      usage(8);
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, only processes one distance matrix at a time\n";
    usage(4);
  }

  int rc = 0;
  
  if (cl.option_present ('t'))
  {
    const_IWSubstring t = cl.string_value ('t');
    if ("float" == t)
      rc = distance_matrix_simple_cluster<float> (cl[0], cout);
    else if ("byte" == t)
      rc = distance_matrix_simple_cluster<unsigned char> (cl[0], cout);
    else
    {
      cerr << "Unrecognised -t qualifier '" << t << "'\n";
      usage(5);
    }
  }
  else
  {
    rc = distance_matrix_simple_cluster<float> (cl[0], cout);
  }

  if (0 == rc)
    return 4;

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_simple_cluster (argc, argv);

  return rc;
}
