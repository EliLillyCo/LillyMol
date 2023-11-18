/*
  Scores an svm model
*/

#include <stdlib.h>
#include <pthread.h>
#include <iostream>
#include <memory>

#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/normalisation.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "gfp_bit_subset.h"

using std::cerr;
using std::endl;

using std::unordered_map;

static int tdts_read = 0;

static double threshold_b = 0.0;

static double kernel_multiplier = 1.0;

static similarity_type_t lower_similarity_cutoff = static_cast<similarity_type_t>(0.0);

static similarity_type_t lower_singleton_similarity_cutoff = static_cast<similarity_type_t>(0.0); 

static GFP_Bit_Subset bit_subset;

static Accumulator<double> items_adjusted_by_singletons;

typedef Accumulator<double> double_accumulator;

/*
  In order to identify the most important bits, we maintain statistics
  about each of the bits. We later gather this info and sort it
*/

// static hash_map<unsigned int, double_accumulator *> global_bit_contributions;
static unordered_map<unsigned int, double_accumulator *> global_bit_contributions;

class Fingerprint_and_Weight : public IW_General_Fingerprint
{
  private:
    double _weight;

  public:
    Fingerprint_and_Weight();

    int construct_from_tdt (IW_TDT &, int &);

    double weighted_similarity (IW_General_Fingerprint & fp, float &);

    int extract_numeric_value_from_id (int w);

    double weight() const { return _weight;}
};

static IWString svmfp_weight_tag("SVMFPW<");

static IWString pcn("PCN<");

static IWString smiles_tag ("$SMI<");
static int output_smiles = 0;

static int pool_size = 0;
static Fingerprint_and_Weight * pool = nullptr;

/*
  The singletons don't really have a weight, they have the 
  activity of the singleton. We just use that class as a
  fingerprint and a number
*/

static Fingerprint_and_Weight * singleton = nullptr;
static int number_singletons = 0;

static similarity_type_t minimum_similarity_for_singleton_consideration = static_cast<similarity_type_t>(0.0);

static Accumulator<double> acc_scores;

static NColumn normalisation;

static int normalisation_active = 0;

static int write_shortest_distance = 0;

static int gather_influence_data = 0;

static IWString_and_File_Descriptor stream_for_vector_contributions;
static IWString_and_File_Descriptor stream_for_bit_contributions;

Fingerprint_and_Weight::Fingerprint_and_Weight()
{
  _weight = 0.0;

  return;
}

double
Fingerprint_and_Weight::weighted_similarity (IW_General_Fingerprint & fp,
                                             float & max_similarity)
{
//similarity_type_t sim = this->tanimoto(fp);
  similarity_type_t sim = this->equal_weight_tanimoto(fp);

  if (sim < lower_similarity_cutoff)
    return 0.0;

  if (sim > max_similarity)
    max_similarity = sim;

  return _weight * static_cast<double>(sim);
}

int
Fingerprint_and_Weight::construct_from_tdt (IW_TDT & tdt,
                                            int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal))
    return 0;

  if (0 == svmfp_weight_tag.length())    // reading singletons
    return 1;

  if (! tdt.dataitem_value(svmfp_weight_tag, _weight))
  {
    cerr << "Fingerprint_and_Weight::construct_from_tdt:cannot extract weight, tag '" << svmfp_weight_tag << "'\n";
    cerr << tdt;
    fatal = 1;
    return 0;
  }

// When models are built with the -j option, we can get weights outside this range

/*if (_weight < -1.0 || _weight > 1.0)
  {
    cerr << "Fingerprint_and_Weight::construct_from_tdt:invalid weight " <<_weight << '\n';
    fatal = 1;
    return 0;
  }*/

//cerr << id() << " weight " << _weight << endl;

  return 1;
}


int
Fingerprint_and_Weight::extract_numeric_value_from_id (int w)
{
  if (_id.nwords() < w + 1)
  {
    cerr << "Fingerprint_and_Weight::extract_numeric_value_from_id:not enough tokens in name '" << _id << "' need " << (w + 1) << endl;
    return 0;
  }

  const_IWSubstring token = _id.word(w);

  if (! token.numeric_value(_weight))
  {
    cerr << "Fingerprint_and_Weight::extract_numeric_value_from_id:invalid numeric, word " << (w + 1) << " '" << _id << "'\n";
    return 1;
  }

  if (normalisation_active)
    normalisation.scale(_weight, _weight);

  return 1;
}

const char * prog_name = nullptr;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Evaluates svmfp models\n";
  cerr << " -U <fname>     support vectors\n";
  cerr << " -s <n>         number of support vectors present\n";
  cerr << " -b <float>     threshold b parameter\n";
  cerr << " -p <precision> default output precision\n";
  cerr << " -S <tag>       support vector weight tag (default '" << svmfp_weight_tag << "')\n";
  cerr << " -N <fname>     activity normalisation data file\n";
  cerr << " -d             include shortest distance to support vector in output\n";
  cerr << " -i             include largest negative and positive influence values in output \n";
  cerr << " -I <fname>     file of singletons, from ./nn_identify_active_outliers\n";
  cerr << " -g <dist>      minimum distance for singleton consideration\n";
  cerr << " -o ...         options for normalisation, enter '-o help' for info\n";
  cerr << " -T <dist>      upper distance threshold - ignore molecule > <dist>\n";
  cerr << " -X <fname>     bit cross reference file created by gfp_to_svm_lite\n";
  cerr << " -D <fname>     write individual vector contributions to <fname>\n";
  cerr << " -B <fname>     write individual bit contributions to <fname>\n";
  cerr << " -m             write smiles rather than descriptors\n";
  cerr << " -k <double>    kernel multiplier\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

int
allocate_pool (int n,
               Fingerprint_and_Weight * & fp)
{
  if (NULL != fp)
  {
    cerr << "Warning, fingerprint array already allocated!!\n";
    delete [] fp;
  }

  fp = new Fingerprint_and_Weight[n];

  if (NULL == fp)
  {
    cerr << "Memory failure, could not allocate " << n << " fingerprints\n";
    return 0;
  }

  return n;
}

static int
read_singletons (iwstring_data_source & input)
{
  std::unique_ptr<re2::RE2> rx = std::make_unique<RE2>("^PCN<");

  int s = input.grep(*rx);

  if (0 == s)
  {
    cerr << "No data in singletons file\n";
    return 0;
  }

  if (! allocate_pool(s, singleton))
    return 0;

  number_singletons = s;

  IW_TDT tdt;

  int items_in_pool = 0;

  int fatal;

  while (tdt.next(input))
  {
    if (! singleton[items_in_pool].construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build pool item\n";
        cerr << tdt;
        return 0;
      }
      continue;
    }

    if (! singleton[items_in_pool].extract_numeric_value_from_id(1))
    {
      cerr << "Invalid singleton activity '" << singleton[items_in_pool].id() << "'\n";
      return 0;
    }

    if (bit_subset.active())
      bit_subset.reduce_to_subset(pool[items_in_pool]);

    items_in_pool++;
  }

  return 1;
}

static int
read_singletons(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open singletons file '" << fname << "'\n";
    return 0;
  }

  return read_singletons (input);
}

static int
read_pool (iwstring_data_source & input)
{
  if (0 == svmfp_weight_tag.length())
  {
    cerr << "Cannot read pool, no weight tag\n";
    return 0;
  }

  int items_in_pool = 0;

  IW_TDT tdt;

  int fatal;   // scope here for efficiency

  while (tdt.next(input))
  {
    if (! pool[items_in_pool].construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build pool item\n";
        cerr << tdt;
        return 0;
      }
      continue;
    }

    if (bit_subset.active())
      bit_subset.reduce_to_subset(pool[items_in_pool]);

    items_in_pool++;

    if (items_in_pool == pool_size)
      break;
  }

  if (tdt.next(input))
  {
    cerr << "Cannot read support vectors, -s option WRONG\n";
    return 0;
  }

  pool_size = items_in_pool;

  return items_in_pool;
}

static int
read_pool (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    pool_size = count_tdts_in_file(input, pcn);
    if (0 == pool_size)
      return 0;

    if (! allocate_pool(pool_size, pool))
      return 0;
  }

  return read_pool(input);
}

/*
  We only allow one singleton to influence the score, so first task
  is to identify the most similar singleton
*/

static int
update_score_for_proximity_to_singletons (IW_General_Fingerprint & fp,
                                          double & rc)
{
  iwmaxid<double, int> maxsim(0.0, -1);

  for (int i = 0; i < number_singletons; i++)
  {
    if (rc >= singleton[i].weight())   // already more active, cannot be incremented
      continue;

    similarity_type_t tmp = singleton[i].tanimoto(fp);

    if (tmp < lower_singleton_similarity_cutoff)
      continue;

    maxsim.try_this(tmp, i);
  }

  int highest_sim = maxsim.which_is_max();

  if (highest_sim < 0)
    return 0;

// A linear interpolation

  double zextra = (singleton[highest_sim].weight() - rc) * (maxsim.maxval() - lower_singleton_similarity_cutoff) / (1.0 - lower_singleton_similarity_cutoff);

  items_adjusted_by_singletons.extra (zextra);

  rc += zextra;

  return 1;
}

class Bit_and_Contribution
{
  private:
    unsigned int _bit;
    double       _contribution;
  public:
    Bit_and_Contribution();

    Bit_and_Contribution(unsigned int b, double c) : _bit(b), _contribution(c) {};

    void set_bit (unsigned int b) { _bit = b;}
    unsigned int bit_number() const { return _bit;}
    void extra (double e) { _contribution += e;}
    double contribution () const { return _contribution;}
};

Bit_and_Contribution::Bit_and_Contribution()
{
  _bit = 0;
  _contribution = 0.0;

  return;
}

class BC_Comparator
{
  private:
  public:
    int operator() (const Bit_and_Contribution &, const Bit_and_Contribution &) const;
};

int
BC_Comparator::operator() (const Bit_and_Contribution & bc1, const Bit_and_Contribution & bc2) const
{
  double c1 = bc1.contribution();
  double c2 = bc2.contribution();

  if (c1 < c2)
    return -1;

  if (c1 > c2)
    return 1;

  return 0;
}

static BC_Comparator bc_comparator;

static int
analyse_bit_contributions (const IWString & smiles,
                           const IWString & id, 
                           const Sparse_Fingerprint & sfp,
                           IWString_and_File_Descriptor & output)
{
  const int n = sfp.nbits();

  Bit_and_Contribution * bc = new Bit_and_Contribution[n]; std::unique_ptr<Bit_and_Contribution[]> free_bc(bc);

  int i = 0;
  unsigned int b;
  int c;

  int ndx = 0;   // index into bc array

  while (sfp.next_bit_set(i, b, c))
  {
    double double_c = static_cast<double>(c);
    bc[ndx].set_bit(b);

    for (int j = 0; j < pool_size; j++)
    {
      const Sparse_Fingerprint & sfpj = pool[j].sparse_fingerprint(0);

      int jc = sfpj.count_for_bit(b);

      if (0 == jc)
        ;
      else if (jc > c)
        bc[ndx].extra(pool[j].weight() * double_c / static_cast<double>(jc));
      else if (jc == c)
        bc[ndx].extra(pool[j].weight());
      else
        bc[ndx].extra(pool[j].weight() * static_cast<double>(jc) / double_c);
    }

    ndx++;
  }

  assert (ndx == n);

  iwqsort(bc, ndx, bc_comparator);

  double sum_contributions = 0.0;

  output << id;
  for (int i = 0; i < n; i++)
  {
    unsigned int b = bc[i].bit_number();

    output << ' ' << b << ' ' << static_cast<float>(bc[i].contribution());

    sum_contributions += bc[i].contribution();

    // hash_map<unsigned int, double_accumulator *>::iterator f = global_bit_contributions.find(b);
    unordered_map <unsigned int, double_accumulator *>::iterator f = global_bit_contributions.find(b);

    if (f == global_bit_contributions.end())
    {
      global_bit_contributions[b] = new double_accumulator;
      global_bit_contributions[b]->extra(bc[i].contribution());
    }
    else
      (*f).second->extra(bc[i].contribution());
  }

  output << ' ' << static_cast<float>(sum_contributions);


  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

template void iwqsort<Bit_and_Contribution, BC_Comparator>(Bit_and_Contribution*, int, BC_Comparator&);
template void iwqsort<Bit_and_Contribution, BC_Comparator>(Bit_and_Contribution*, int, BC_Comparator&, void*);
template void compare_two_items<Bit_and_Contribution, BC_Comparator>(Bit_and_Contribution*, BC_Comparator&, void*);
template void move_in_from_left<Bit_and_Contribution, BC_Comparator>(Bit_and_Contribution*, int&, int&, int, BC_Comparator&, void*);
template void move_in_from_right<Bit_and_Contribution, BC_Comparator>(Bit_and_Contribution*, int&, int&, BC_Comparator&);
template void swap_elements<Bit_and_Contribution>(Bit_and_Contribution&, Bit_and_Contribution&, void*);

static int
analyse_bit_contributions (IW_General_Fingerprint & fp,
                           const IWString & smiles,
                           IWString_and_File_Descriptor & output)
{
  if (bit_subset.active())
    bit_subset.reduce_to_subset(fp);

  int n = number_fingerprints() + fp.number_sparse_fingerprints();

  if (1 != n)
  {
    cerr << "Sorry, bit_contributions cannot process " << n << " fingerprints\n";
    return 0;
  }

  if (1 != fp.number_sparse_fingerprints())
  {
    cerr << "Sorry, bit_contributions restricted to sparse fingerprints only, see Ian\n";
    return 0;
  }

  const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(0);

  return analyse_bit_contributions (smiles, fp.id(), sfp, output);
}

class Double_Accumulator_Compartor
{
  private:
  public:
    int operator() (const Accumulator<double> *, const Accumulator<double> *) const;
};

int
Double_Accumulator_Compartor::operator() (const Accumulator<double> * a1,
                                          const Accumulator<double> * a2) const
{
  double ave1 = a1->average();
  double ave2 = a2->average();

  if (ave1 < ave2)
    return -1;

  if (ave1 > ave2)
    return 1;

  return 0;
}

class Global_Bit_Summary
{
  private:
    unsigned int _bit;
    int _number_molecules;
    double _average_contribution;

  public:
    Global_Bit_Summary();

    void set_bit_number (unsigned int b) { _bit = b;}
    void set_number_molecules (int n) { _number_molecules = n;}
    void set_average_contribution (double s) { _average_contribution = s;}

    unsigned int bit_number() const { return _bit;}
    int number_molecules() const { return _number_molecules;}
    double average_contribution() const { return _average_contribution;}
};

Global_Bit_Summary::Global_Bit_Summary()
{
  _bit = 0;
  _number_molecules = 0;
  _average_contribution = 0.0;

  return;
}

class Global_Bit_Summary_Comparator
{
  private:
  public:
    int operator() (const Global_Bit_Summary &, const Global_Bit_Summary &) const;
};

int
Global_Bit_Summary_Comparator::operator() (const Global_Bit_Summary & gb1,
                                           const Global_Bit_Summary & gb2) const
{
  double a1 = gb1.average_contribution();
  double a2 = gb2.average_contribution();

  if (a1 < a2)
    return -1;

  if (a1 > a2)
    return 1;

  return 0;
}

static int
// write_global_bit_contributions (const hash_map<unsigned int, double_accumulator *> & global_bit_contributions,
write_global_bit_contributions (const unordered_map<unsigned int, double_accumulator *> & global_bit_contributions,
                                IWString_and_File_Descriptor & output)
{
  int threshold = pool_size / 100;
  if (0 == threshold)
    threshold = 2;

  int n = 0;

  // typedef hash_map<unsigned int, double_accumulator *> H;
  typedef unordered_map<unsigned int, double_accumulator *> H;

  for (H::const_iterator i = global_bit_contributions.begin(); i != global_bit_contributions.end(); i++)
  {
    const double_accumulator * a = (*i).second;

    if (a->n() >= static_cast<unsigned int>(threshold))
      n++;
  }

  if (0 == n)
  {
    cerr << "write_global_bit_contributions:no bits meet threshold " << threshold << endl;
    return 0;
  }

  Global_Bit_Summary * tmp = new Global_Bit_Summary[n]; std::unique_ptr<Global_Bit_Summary[]> free_tmp(tmp);

  n = 0;

  for (H::const_iterator i = global_bit_contributions.begin(); i != global_bit_contributions.end(); i++)
  {
    const double_accumulator * a = (*i).second;

    if (a->n() < static_cast<unsigned int>(threshold))
      continue;

    tmp[n].set_bit_number((*i).first);
    tmp[n].set_average_contribution((*i).second->average());
    tmp[n].set_number_molecules((*i).second->n());

    n++;
  }

  Global_Bit_Summary_Comparator gbsc;

  iwqsort (tmp, n, gbsc);

  for (int i = 0; i < n; i++)
  {
    output << tmp[i].bit_number() << ' ' << tmp[i].number_molecules() << ' ' << static_cast<float>(tmp[i].average_contribution()) << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

struct Do_Part_of_Pool
{
  int istart;
  int istop;
  double result;
};

static Do_Part_of_Pool * dpp = nullptr;

static int nthreads = 0;
static pthread_t * child_thread = nullptr;

static pthread_barrier_t start_barrier;
static pthread_barrier_t end_barrier;

static IW_General_Fingerprint * needle = nullptr;

static double 
thread_svmfp_score(int istart, int istop, 
                   IW_General_Fingerprint & fp)
{
  double rc = 0.0;

  float max_similarity = 0.0;

  if (1.0 == kernel_multiplier)
  {
    for (int i = istart; i != istop; i++)
    {
      rc += pool[i].weighted_similarity(fp, max_similarity);
    }
  }
  else
  {
    for (int i = istart; i != istop; i++)
    {
      rc += kernel_multiplier * pool[i].weighted_similarity(fp, max_similarity);
    }
  }

  return rc;
}

static void *
do_child_thread (void * v)
{
  Do_Part_of_Pool * mydpp = reinterpret_cast<Do_Part_of_Pool *>(v);

//cerr << "Child waiting for start\n";
  pthread_barrier_wait(&start_barrier);

//cerr << "Child got past start_barrier\n";

  while (NULL != needle)
  {
    mydpp->result = thread_svmfp_score(mydpp->istart, mydpp->istop, *needle);
//  cerr << "Thread result " << mydpp->result << endl;

    pthread_barrier_wait(&end_barrier);
    pthread_barrier_wait(&start_barrier);
  }

  pthread_barrier_wait(&end_barrier);
  pthread_exit(NULL);
  return NULL;
}

static int
initialise_threads(int nthreads)
{
  assert (nthreads > 0);

  if (nthreads > pool_size)
  {
    cerr << pool_size << " fingerprints, cannot use " << nthreads << " threads\n";
    return 0;
  }

  int items_per_thread = pool_size / nthreads;
  if (items_per_thread < 1)
  {
    cerr << "Pool too small, " << pool_size << " for " << nthreads << " threads\n";
    items_per_thread = 1;
  }

  if (verbose)
    cerr << pool_size << " fingerprints and " << nthreads << " threads means " << items_per_thread << " items per thread\n";

  dpp = new Do_Part_of_Pool[nthreads];

  for (int i = 0; i < nthreads; i++)
  {
    dpp[i].istart = i * items_per_thread;
    if (nthreads - 1 == i)
      dpp[i].istop = pool_size;
    else
      dpp[i].istop = i * items_per_thread + items_per_thread;
  }

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init (&attr);
  pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    Do_Part_of_Pool & dppi = dpp[i];

    if (0 != pthread_create (&(child_thread[i]), &attr, do_child_thread, &dppi))
    {
      perror ("cannot create child process");
      return 0;
    }
  }

  if (pthread_barrier_init(&start_barrier, NULL, nthreads + 1) < 0)
    return 0;

  if (pthread_barrier_init(&end_barrier, NULL, nthreads + 1) < 0)
    return 0;

  return nthreads;
}

static double
svmfp_score_threaded (IW_General_Fingerprint & fp)
{
  needle = &fp;

//cerr << "Post to start barrier\n";
  pthread_barrier_wait(&start_barrier);
//cerr << "Post to end barrier\n";
  pthread_barrier_wait(&end_barrier);

  double rc = 0;
  for (int i = 0; i <nthreads; i++)
  {
    rc += dpp[i].result;
  }

//cerr << "Result is " << rc << endl;

  return - threshold_b + rc;
}

static int
svmfp_score (IW_General_Fingerprint & fp,
             const IWString & smiles,
             IWString_and_File_Descriptor & output)
{
  if (bit_subset.active())
  {
    bit_subset.reduce_to_subset(fp);
//  int t = bit_subset.reduce_to_subset(fp);
//  cerr << "Reduced " << t << " bits from '" << fp.id() << "'\n";
  }

  double rc = - threshold_b;

  float max_similarity = 0.0;

// Use pointers to these objects to avoid costly constructors.

  iwminid<double, int> * ineg = nullptr;
  iwmaxid<double, int> * ipos = nullptr;

  if (stream_for_vector_contributions.is_open())
  {
    append_first_token_of_name (fp.id(), stream_for_vector_contributions);

    double tmp;
    for (int i = 0; i < pool_size; i++)
    {
      Fingerprint_and_Weight & pi = pool[i];

      tmp = pi.weighted_similarity(fp, max_similarity);
      stream_for_vector_contributions << ' ' << static_cast<float>(tmp);
      rc += tmp;
    }

    stream_for_vector_contributions << ' ' << static_cast<float>(rc);

    if (! normalisation_active)
      stream_for_vector_contributions << '\n';

    if (stream_for_vector_contributions.size() > 32768)
      stream_for_vector_contributions.write_whole_blocks_shift_unwritten();
  }
  else if (gather_influence_data)
  {
    ineg = new iwminid<double, int>(9.9e+09, -1);
    ipos = new iwmaxid<double, int>(-9.9e+09, -1);

    double tmp;   // scope here for efficiency

    for (int i = 0; i < pool_size; i++)
    {
      tmp = pool[i].weighted_similarity(fp, max_similarity);
      ineg->try_this(tmp, i);
      ipos->try_this(tmp, i);
      rc += tmp;
    }
  }
  else if (nthreads > 0)
    rc = svmfp_score_threaded(fp);
  else if (1.0 == kernel_multiplier)
  {
    for (int i = 0; i < pool_size; i++)
    {
      rc += pool[i].weighted_similarity(fp, max_similarity);
    }
  }
  else
  {
    for (int i = 0; i < pool_size; i++)
    {
      rc += kernel_multiplier * pool[i].weighted_similarity(fp, max_similarity);
    }
  }

  if (number_singletons)
    update_score_for_proximity_to_singletons (fp, rc);

  if (output_smiles)
    output << smiles << ' ' << fp.id();
  else
    append_first_token_of_name(fp.id(), output);

  output << ' ';

  if (normalisation_active)
  {
    double unscaled;
    normalisation.unscale(rc, unscaled);
//  cerr << "Unscale '" << fp.id() << "' " << rc << " to " << unscaled << endl;
    output << static_cast<float>(unscaled);

    if (stream_for_vector_contributions.active())
      stream_for_vector_contributions << ' ' << static_cast<float>(unscaled) << '\n';

    if (verbose)
      acc_scores.extra(unscaled);
  }
  else
  {
    output << static_cast<float>(rc);

    if (verbose)
      acc_scores.extra(rc);
  }

  if (write_shortest_distance)
    output << ' ' << (static_cast<float>(1.0) - max_similarity);

  if (gather_influence_data)
  {
    int j = ineg->which_is_min();
    output << ' ' << pool[j].id() << ' ' << ineg->minval();

    j = ipos->which_is_max();
    output << ' ' << pool[j].id() << ' ' << ipos->maxval();

    delete ineg;
    delete ipos;
  }

  output << '\n';

  return 1;
}

static int
svmfp_score (IW_TDT & tdt,
             IWString_and_File_Descriptor & output)
{
  IW_General_Fingerprint fp;

  int fatal;
  if (! fp.construct_from_tdt (tdt, fatal))
  {
    if (fatal)
      return 0;

    return 1;
  }

  IWString smiles;

  if (! output_smiles)
    ;
  else if (tdt.dataitem_value(smiles_tag, smiles))
    ;
  else
  {
    cerr << "Cannot extract smiles '" << smiles_tag << "' from tdt\n";
    return 0;
  }

  if (stream_for_bit_contributions.is_open())
    return analyse_bit_contributions(fp, smiles, stream_for_bit_contributions);

  return svmfp_score (fp, smiles, output);
}

static int
svmfp_score (iwstring_data_source & input,
             IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    tdts_read++;

    if (! svmfp_score (tdt, output))
      return 0;

    if (output.size() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  return 1;
}

static int
svmfp_score (const char * fname,
             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return svmfp_score(input, output);
}

static int
read_normalisation_data(iwstring_data_source & input,
                        NColumn & normalisation)
{
  IWString buffer;

  while (input.next_record(buffer))
  {
    int fatal;

    if (parse_normalisation_file_record (buffer, fatal))
      continue;

    if (fatal)
      return 0;

    if (buffer.starts_with('#'))
      continue;

    if (0 == buffer.length())
      continue;

    buffer.truncate_at_first ('#');

    buffer.remove_leading_words(1);

    return normalisation.establish_range_from_pre_existing_data(buffer);
  }

  cerr << "No normalisation data in file\n";
  return 0;
}

static int
read_normalisation_data(const IWString & n,
                        NColumn & normalisation)
{
  iwstring_data_source input(n);

  if (! input.good())
  {
    cerr << "Cannot open normalisation data file '" << n << "'\n";
    return 0;
  }

  return read_normalisation_data (input, normalisation);
}

/*
  The way svm_learn is run, all fingerprints, of all types get mixed together
  with equal weight. But all the gfp_* programmes run by assigning each fingerprint
  equal weight by default. In order to reproduce the results from svm_learn,
  we need to adjust the weights

  WRONG! this never worked properly. Maybe there's a way to figure out how
  to do this, but it isn't obvious
*/

#ifdef DOES_NOT_WORK4
static int
adjust_fingerprint_weights_according_to_number_of_bits (const GFP_Bit_Subset & bit_subset)
{
  int total_bits = bit_subset.number_properties ();

  for (int i = 0; i < bit_subset.number_fixed_width_fingerprints(); i++)
  {
    total_bits += bit_subset.bits_used_in_fixed_width_fingerprint(i);
  }

  for (int i = 0; i < bit_subset.number_sparse_fingerprints(); i++)
  {
    total_bits += bit_subset.bits_in_sparse_fingerprint(i);
  }

  double w = static_cast<double>(bit_subset.number_properties()) / static_cast<double>(total_bits);

  set_property_weight_integer(w);

  if (verbose > 1)
    cerr << "Property weight " << w << endl;

  for (int i = 0; i < bit_subset.number_fixed_width_fingerprints(); i++)
  {
    int nb = bit_subset.bits_used_in_fixed_width_fingerprint(i);

    w = static_cast<double>(nb) / static_cast<double>(total_bits);

    set_fixed_width_fingerprint_weight(i, w);

    if (verbose)
      cerr << "Fixed width fingerprint " << i << " weight " << w << endl;
  }

  for (int i = 0; i < bit_subset.number_sparse_fingerprints(); i++)
  {
    int nb = bit_subset.bits_in_sparse_fingerprint(i);

    w = static_cast<double>(nb) / static_cast<double>(total_bits);

    set_sparse_fingerprint_weight(i, w);

    if (verbose)
      cerr << "Sparse fingerprint " << i << " weight " << w << endl;
  }

  return 1;
}
#endif

static int
svmfp_score (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vP:p:S:n:P:F:b:s:H:N:dio:U:T:I:g:D:X:mB:k:h:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('U'))
  {
    cerr << "Must specify the support vectors via the -U option\n";
    usage(2);
  }

  if (! cl.option_present('s'))
  {
    cerr << "Must specify the number of support vectors via the -s option\n";
    usage(2);
  }

  if (! cl.option_present('b'))
  {
    cerr << "Must specify the threshold b parameter via the -b option\n";
    usage(2);
  }

  if (! cl.value('b', threshold_b))
  {
    cerr << "Invalid threshold b parameter (-b)\n";
    return 4;
  }

  if (verbose)
    cerr << "Threshold set to " << threshold_b << '\n';

  if (cl.option_present('k'))
  {
    if (! cl.value('k', kernel_multiplier) || kernel_multiplier <= 0.0)
    {
      cerr << "The kernel multiplier (-k) must be a +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "kernel multiplier " << kernel_multiplier << endl;
  }

  if (! verbose)
    NColumn_set_report_out_of_range_values(0);

  if (cl.option_present('T'))
  {
    if (! cl.value('T', lower_similarity_cutoff) || lower_similarity_cutoff < 0.0 ||
          lower_similarity_cutoff > 1.0)
    {
      cerr << "The upper distance cutoff (-T) must be a valid distance\n";
      usage(3);
    }

    lower_similarity_cutoff = static_cast<float>(1.0) - lower_similarity_cutoff;
  }

  if (cl.option_present('d'))
  {
    write_shortest_distance = 1;

    if (verbose)
      cerr << "Will include the shortest distance to a support vector in the output\n";
  }

  if (cl.option_present('i'))
  {
    gather_influence_data = 1;

    if (verbose)
      cerr << "Will identify the support vector with greatest positive and greatest negative influence\n";
  }

  if (cl.option_present('p'))
  {
    int p;
    if (! cl.value('p', p) || p < 1)
    {
      cerr << "Invalid default output precision (-p)\n";
      usage(3);
    }

    if (verbose)
      cerr << "Default floating point output precision " << p << '\n';

    set_default_iwstring_float_concatenation_precision(p);
  }

  set_sparsefp_warn_empty_data(0);

  if (cl.option_present('o'))
  {
    if (! parse_normalisation_options (cl, 'o', verbose))
    {
      cerr << "Cannot discern normalisation options (-o)\n";
      return 7;
    }
  }

  if (cl.option_present('m'))
  {
    output_smiles = 1;

    if (verbose)
      cerr << "Will write smiles output\n";
  }

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    if (! bit_subset.do_read(x))
    {
      cerr << "Cannot read bit subset information '" << x << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Bit subset information read from '" << x << "'\n";

    if (verbose > 1)
      bit_subset.debug_print(cerr);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }

//if (bit_subset.active())
//  adjust_fingerprint_weights_according_to_number_of_bits (bit_subset);

  if (cl.option_present('N'))
  {
    IWString n = cl.string_value('N');

    if (! dash_s(n.null_terminated_chars()))
    {
      cerr << "Missing or empty normalisation file '" << n << "'\n";
      return 3;
    }
    if (! read_normalisation_data(n, normalisation))
    {
      cerr << "Invalid normalisation data '" << n << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Read normalisation data from '" << n << "'\n";

    normalisation_active = 1;
  }

  if (cl.option_present('S'))
  {
    cl.value('S', svmfp_weight_tag);

    if (verbose)
      cerr << "Weight tag in '" << svmfp_weight_tag << "'\n";

    if (! svmfp_weight_tag.ends_with('<'))
      svmfp_weight_tag.add('<');
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', nthreads) || nthreads < 1)
    {
      cerr << "The number of threads option (-h) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will run with " << nthreads << " threads\n";
  }

  if (cl.option_present('D'))
  {
    if (cl.option_present('i'))
    {
      cerr << "Sorry, the -D option is mutually incompatible with the -i option\n";
      usage(3);
    }

    const char * d = cl.option_value('D');

    if (! stream_for_vector_contributions.open(d))
    {
      cerr << "Cannot open stream for individual vector contributions '" << d << "'\n";
      return 9;
    }

    if (verbose)
      cerr << "INdividual vector contributions written to '" << d << "'\n";
  }

  if (cl.option_present('B'))
  {
    if (cl.option_present('i'))
    {
      cerr << "Sorry, the -D option is mutually incompatible with the -i option\n";
      usage(3);
    }

    const char * b = cl.option_value('B');

    if (! stream_for_bit_contributions.open(b))
    {
      cerr << "Cannot open stream for individual bit contributions '" << b << "'\n";
      return 9;
    }

    if (verbose)
      cerr << "INdividual bit contributions written to '" << b << "'\n";

    set_default_iwstring_float_concatenation_precision(4);
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The pool size value (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (! allocate_pool(pool_size, pool))
    {
      cerr << "Cannot allocate pool " << pool_size << '\n';
      return 3;
    }
  }

  if (cl.option_present('U'))
  {
    IWString p = cl.string_value('U');
    if (! read_pool(p))
    {
      cerr << "Cannot read pool '" << p << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read pool\n";
  }

  if (nthreads > 0)
  {
    if (! initialise_threads(nthreads))
    {
      cerr << "Cannot initialise threads\n";
      return 4;
    }
  }

  if (stream_for_vector_contributions.is_open())
  {
    stream_for_vector_contributions << "ID";
    for (int i = 0; i < pool_size; i++)
    {
      stream_for_vector_contributions << ' ';
      append_first_token_of_name (pool[i].id(), stream_for_vector_contributions);
      stream_for_vector_contributions.write_if_buffer_holds_more_than(32768);
    }
    if (normalisation_active)
      stream_for_vector_contributions << " unscaled_pred pred\n";
    else
      stream_for_vector_contributions << " pred\n";
  }

  if (cl.option_present('I'))
  {
    const char * g = cl.option_value('I');

    IWString save_weight_tag(svmfp_weight_tag);

    svmfp_weight_tag.resize_keep_storage(0);

    if (! read_singletons(g))
    {
      cerr << "Cannot read active singletons from '" << g << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read " << number_singletons << " active singletons from '" << g << "'\n";

    svmfp_weight_tag = save_weight_tag;

    if (cl.option_present('g'))
    {
      if (! cl.value('g', minimum_similarity_for_singleton_consideration) ||
            minimum_similarity_for_singleton_consideration < 0.0 || minimum_similarity_for_singleton_consideration > 1.0)
      {
        cerr << "The minimum similarity for singleton comparison (-g) must be a valid similarity\n";
        usage(4);
      }

      if (verbose)
        cerr << "Will only consider singletons closer than " << minimum_similarity_for_singleton_consideration << '\n';

//    convert to similarity

      minimum_similarity_for_singleton_consideration = 1.0 - minimum_similarity_for_singleton_consideration;
    }
  }

  IWString_and_File_Descriptor output(1);

  if (! output_smiles)
  {
    if (cl.option_present('H'))
    {
      output << "Name SVMFP:";
      output << cl.string_value('H');
    }
    else
      output << "Name gfp_SVM_Predicted";

    if (write_shortest_distance)
      output << " mindist";

    if (gather_influence_data)
      output << " Neginf Infneg Posinf Infpos";

    output << '\n';
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! svmfp_score(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << tdts_read << " fingerprints for scoring\n";
    if (acc_scores.n() > 1)
      cerr << "Predicted values between " << static_cast<float>(acc_scores.minval()) << " and " <<
              static_cast<float>(acc_scores.maxval()) << " ave " << static_cast<float>(acc_scores.average()) << '\n';
    if (items_adjusted_by_singletons.n() > 0)
    {
      cerr << items_adjusted_by_singletons.n() << " items adjusted for proximimity to any of " << number_singletons << " singletons\n";
      cerr << "Increments between " << static_cast<float>(items_adjusted_by_singletons.minval()) << " and " << static_cast<float>(items_adjusted_by_singletons.maxval()) << " ave " << static_cast<float>(items_adjusted_by_singletons.average_if_available_minval_if_not()) << '\n';
    }
  }

  if (nthreads)
  {
    needle = nullptr;
    pthread_barrier_wait(&start_barrier);
    pthread_barrier_wait(&end_barrier);

    pthread_barrier_destroy(&start_barrier);
    pthread_barrier_destroy(&end_barrier);
  }

  if (global_bit_contributions.size())
    write_global_bit_contributions (global_bit_contributions, output);

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = svmfp_score(argc, argv);

  return rc;
}
