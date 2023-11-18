/*
  Implementation of Isolation Forests

  http://cs.nju.edu.cn/zhouzh/zhouzh.files/publication/icdm08b.pdf

  Likely better versions of this should now be available.
*/

#define USE_OMP
#ifdef USE_OMP
#include <omp.h>
#endif

#include <random>
#include <type_traits>
#include <cmath>
#include <memory>
#include <algorithm>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/report_progress.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int has_header = 1;
static int has_rownames = 1;

static char input_separator = ' ';
static char output_separator = ' ';

static int print_tree = 0;

static std::random_device rd;

static int run_consistency_check = 0;

static int out_of_range = 0;

static int produce_sorted_output = 0;

static int include_out_of_range_values_in_regular_output = 1;

static int include_raw_data_with_output = 0;

static int include_original_data_in_output = 0;

static int extra_columns_for_out_of_range_info = 0;

static int include_header_in_output = 1;

static IWString_and_File_Descriptor stream_for_out_of_range;

static Report_Progress report_progress;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Implementation of Isolation Forests\n";
  cerr << " -n <ntrees>    number of trees to build (default 100)\n";
  cerr << " -p <psi>       psi parameter from paper (default value 256)\n";
  cerr << " -T <type>      type of numbers to use (double, int, float)\n";
  cerr << " -C <fname>     reference collection - input determines outliers wrt <fname>\n";
//cerr << " -s             produce sorted output\n";
  cerr << " -x             exclude out of range values from output\n";
  cerr << " -X <fname>     write out of range values to <fname>\n";
  cerr << " -O             add an extra column with out of range info\n";
  cerr << " -P             print the tree after it is formed\n";
  cerr << " -e             echo input data with the score\n";
  cerr << " -E             echo comparator data also - score will be written 0.5\n";
  cerr << " -c             run a consistency check on the tree once built\n";
  cerr << " -y             do NOT produce a header record in the output\n";
  cerr << " -u             input files to NOT contain header records\n";
  cerr << " -w             input file do NOT contain row names\n";
  cerr << " -i <sep>       input  separator (default space, tab, comma,... recognised)\n";
  cerr << " -o <sep>       output separator (default space, tab, comma,... recognised)\n";
  cerr << " -z <prec>      default floating point precision for output\n";
#if defined(OMP)
  cerr << " -h <nthreads>  number of OMP threads to use for scoring\n";
#else
  cerr << " -r <n>         report progress every <n> items scored\n";
#endif
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Members_as_Bits
{
  private:
    IW_Bits_Base _b;

  public:
    Members_as_Bits();

    int set_size(const int s);

    void set_all_items_in_set();
    void clear() { _b.clear();}

    int first_member() const;

    int is_member(const int s) const { return _b.is_set(s);}
    void add_member(const int s) { _b.set(s);}
    int remove_member(const int s) { _b.set(s, 0); return 1;}    // should check to see if set

    int next_member(int & s) const { return _b.next_on_bit(s);}

    int number_members() const { return _b.nset();}
};

Members_as_Bits::Members_as_Bits()
{
}

int
Members_as_Bits::set_size(const int s)
{
  int b = s;
  if (0 != b % 8)
    b = (b / 8 + 1) * 8;

  _b.allocate_space_for_bits(b);
  _b.clear();

  return 1;
}

void
Members_as_Bits::set_all_items_in_set()
{
  _b.set_all();
}

int
Members_as_Bits::first_member() const
{
  return _b.first_bit();
}

template <typename T>
class Tabular_Data
{
  private:
    T * _v;
    int _nrows;
    int _ncols;

    T * _colmin;
    T * _colmax;

    IWString * _colname;
    IWString * _rowname;

//  private functions

    int _initialise_column_names(const const_IWSubstring & buffer);
    int _parse_input_record(const const_IWSubstring & buffer, const int zrow, T * v);

  public:
    Tabular_Data();
    ~Tabular_Data();

    int nrows() const { return _nrows;}
    int ncols() const { return _ncols;}

    T colmin(const int c) const { return _colmin[c];}
    T colmax(const int c) const { return _colmax[c];}

    int build(const char * fname);
    int build(iwstring_data_source &);

    int same_descriptor_names(const Tabular_Data<T> &) const;

    const IWString & column_name(const int s) const { return _colname[s];}
    const IWString & row_name(const int s) const { return _rowname[s];}

    T operator[](const int s) { return _v[s];}

    T operator() (const int r, const int c) const { return _v[r * _ncols + c];}

    const T * row(const int s) const { return _v + s * _ncols;}

    int first_direction_out_of_range(const T * v) const;
    int number_attributes_out_of_range(const T * x) const;
    int number_attributes_out_of_range(const T * x, T & max_out_of_range) const;

    int create_DataFrame(const char * vname, IWString_and_File_Descriptor & output) const;
};

template <typename T>
Tabular_Data<T>::Tabular_Data()
{
  _v = nullptr;
  _nrows = 0;
  _ncols = 0;

  _colmin = nullptr;
  _colmax = nullptr;

  _colname = nullptr;
  _rowname = nullptr;

  return;
}

template <typename T>
Tabular_Data<T>::~Tabular_Data()
{
  if (nullptr != _v)
    delete [] _v;

  if (nullptr != _rowname)
    delete [] _rowname;

  if (nullptr != _colname)
    delete [] _colname;

  if (nullptr != _colmin)
    delete [] _colmin;

  if (nullptr != _colmax)
    delete [] _colmax;

  return;
}

template <typename T>
int
Tabular_Data<T>::same_descriptor_names(const Tabular_Data<T> & rhs) const
{
  if (! has_header)   // assume neither one has a header
    return 1;

  if (_ncols != rhs._ncols)
    return 0;

  for (int i = 0; i < _ncols; ++i)
  {
    if (_colname[i] != rhs._colname[i])
      return 0;
  }

  return 1;
}

/*
  We need something to keep track of the tree as it is formed.
*/

#ifdef NEED_TSF_QWE
template <typename T>
class Tree_so_Far
{
  private:
    T _split;
    int _col;
    std::binary_function<T, T, bool> _f;

  public:
    Tree_so_Far();

    void set_column(const int s) { _col = s;}
    void set_split_point(const T s) { _split = s;}
    void set_function(std::binary_function<T, T, bool> & s) { _f = s;}

    int col() const { return _col;}
    T split() const { return _split;}

    bool operator() (const T s) const { return _f(s, _split);}
};
#endif

template <typename T>
class IF_Result
{
  public:
    int    _ndx;
    double _score;
    int    _out_of_bounds;
    T      _max_out_of_bounds;
  public:
    IF_Result();
};

template <typename T>
IF_Result<T>::IF_Result()
{
  _ndx = 0;
  _score = 0.0;
  _out_of_bounds = 0;
  _max_out_of_bounds = 0.0;
}

template <typename T, typename M>
class IFTree
{
  private:
    int _split_column;
    T _split_value;
    int _mydepth;

    int _n;    // if we are a terminal node, we will have a number of items

    std::mt19937_64 _rng;

    const IFTree * _parent;

    M _members;

    IFTree * _left;
    IFTree * _right;

//  private functions

    int _determine_range(const IFTree<T, M> * requestor, const int col, T & zmin, T & zmax) const;
    template <typename C> int _cutoffs_consistent(C cmp, const int parent_col, const T from_parent) const;

  public:
    IFTree();
    ~IFTree();

    void set_parent(const IFTree<T, M> * s) { _parent = s;}

    void set_size(const int s);

    void set_depth(const int s) { _mydepth = s;}

    int print_tree(const char prefix, std::ostream & output) const;

    int build(const Tabular_Data<T> & zdata, const int current_tree_height, const int height_limit);

    int number_members() const { return _members.number_members();}

    int split_column() const { return _split_column;}

    void add_member(const int s) { _members.add_member(s);}

    double score(const T *) const;

    int cutoffs_consistent() const;

    void gather_tree_depths(extending_resizable_array<int> & acc) const;
};

template <typename T, typename M>
IFTree<T, M>::IFTree() : _rng(rd())
{
  _split_column = -1;
  _mydepth = 0;

  _n = 0;

  _parent = nullptr;
  _left = nullptr;
  _right = nullptr;

  return;
}

template <typename T, typename M>
IFTree<T, M>::~IFTree()
{
  if (nullptr != _left)
    delete _left;

  if (nullptr != _right)
    delete _right;

  return;
}

template <typename T, typename M>
void
IFTree<T, M>::set_size(const int s)
{
  _members.set_size(s);

  return;
}

template <typename T, typename M>
void
IFTree<T, M>::gather_tree_depths(extending_resizable_array<int> & acc) const
{
  if (_members.number_members() > 0)
  {
    acc[_mydepth]++;
    return;
  }

  if (nullptr != _left)
    _left->gather_tree_depths(acc);

  if (nullptr != _right)
    _right->gather_tree_depths(acc);

  return;
}

template <typename T, typename M>
int
IFTree<T, M>::cutoffs_consistent() const
{
  if (nullptr != _left)
  {
    if (! _left->_cutoffs_consistent(std::greater<T>(), _split_column, _split_value))
    {
      cerr << "IFTree::cutoffs_consistent:left child inconsistent. Levl " << _mydepth << " col " << _split_column << " value " <<  _split_value << endl;
      return 0;
    }
  }

  if (nullptr != _right)
  {
    if (! _right->_cutoffs_consistent(std::less<T>(), _split_column, _split_value))
    {
      cerr << "IFTree::cutoffs_consistent:right child inconsistent. Levl " << _mydepth << " col " << _split_column << " value " <<  _split_value << endl;
      return 0;
    }
  }

  return 1;
}

template <typename T, typename M> template <typename C>
int
IFTree<T, M>::_cutoffs_consistent(C cmp,
                                 const int parent_column,
                                 const T parent_split) const
{
  if (_split_column != parent_column)    // do not need to check anything
    ;
  else if (! cmp(parent_split, _split_value))
  {
    cerr << "IFTree::cutoffs_consistent:inconsistent. Level " << _mydepth << " col " << _split_column << " from parent " << parent_split << " my split " << _split_value << endl;
    return 0;
  }

  if (nullptr != _left)
  {
    if (! _left->_cutoffs_consistent(std::greater<T>(), _split_column, _split_value))
    {
      cerr << "IFTree::current_tree_height:left inconsistent. Level " << _mydepth << " col " << _split_column << " value " << _split_value << endl;
      return 0;
    }
  }

  if (nullptr != _right)
  {
    if (! _right->_cutoffs_consistent(std::less<T>(), _split_column, _split_value))
    {
      cerr << "IFTree::current_tree_height:right inconsistent. Level " << _mydepth << " col " << _split_column << " value " << _split_value << endl;
      return 0;
    }
  }

  return 1;
}

/*
  Equations from the paper
*/

static double
harmonic_number(const double x)
{
  if (x <= 0.0)   // not really sure what to do...
    return 1.0;

  return log(x) + 0.5772156649;
}

static double
c_function_1(const int n)
{
  const double x = static_cast<double>(n);

  return 2.0 * harmonic_number(x-1.0) - (2.0 * (x-1.0)/x);
}

//#define DEBUG_IFT_SCORE

template <typename T, typename M>
double
IFTree<T, M>::score(const T * d) const
{
  const int n = _members.number_members();

#ifdef DEBUG_IFT_SCORE
  if (n > 0)
    cerr << "IFTree::score:terminal node, depth " << _mydepth << " contains " << n << " values\n";
#endif

  if (n > 0)
    return static_cast<double>(_mydepth) + c_function_1(n);

  if (nullptr == _left)    // tree never got any points added, and did not grow
    return 0.0;

#ifdef DEBUG_IFT_SCORE
  cerr << "IFTree::score:level " << _mydepth << " checking column " << _split_column << " input value " << d[_split_column] << " my split " << _split_value << endl;
#endif

  if (d[_split_column] < _split_value)
    return _left->score(d);
  else
    return _right->score(d);
}

#ifdef BAD_IDEA
template <typename T>
void
determine_range(const Tabular_Data<T> & zdata,
                const Tree_so_Far<T> * tree_so_far,
                const int current_tree_height,
                const int col,
                T & zmin,
                T & zmax)
{
  zmin = zdata.colmin(col);
  zmax = zdata.colmax(col);

  if (0 == current_tree_height)
    return;

  for (int i = 0; i < current_tree_height; ++i)
  {
    if (tree_so_far[i].col() != col)    // split did not involve our column, no impact
      continue;

    if (tree_so_far[i](zmin))
      zmin = tree_so_far[i].split();
    else if (tree_so_far[i](zmax))
      zmax = tree_so_far[i].split();
  }
}
#endif

/*
  A descendant wants to split, and they need to know the range of the variable they are
  considering. If we split on that variable, we can adjust zmin or zmax
*/

template <typename T, typename M>
int
IFTree<T, M>::_determine_range(const IFTree<T, M> * requestor,   // one of our children
                               const int col,
                               T & zmin,
                               T & zmax) const

{
  if (col != _split_column)   // we cannot adjust the range
  {
    if (nullptr == _parent)    // done
      return 1;

    return _parent->_determine_range(this, col, zmin, zmax);
  }

// We are splitting on the column.

  if (requestor == _left)
    zmax = _split_value;
  else if (requestor == _right)
    zmin = _split_value;
  else
  {
    cerr << "IFTree::determine_range:not a child\n";
    return 0;
  }

  return 1;            // the first ancestor who split this column will be the limit for a descendant
}

template <typename T, typename M>
int
IFTree<T, M>::print_tree(const char prefix, std::ostream & output) const
{
  for (int i = 0; i < _mydepth; ++i)
  {
    output << ' ';
  }

  output << prefix << ' ';

  const int n = _members.number_members();

  if (n)
  {
    output << "terminal node with " << n << " items, depth " << _mydepth << '\n';
    return 1;
  }

  if (_split_column < 0)
  {
    output << "empty\n";
    return 1;
  }

  output << "split col " << _split_column << " value " << _split_value << endl;

  if (nullptr != _left)
    _left->print_tree('L', output);

  if (nullptr != _right)
    _right->print_tree('R', output);

  return 1;
}

//#define DEBUG_IFT_BUILD

template <typename T, typename M>
int
IFTree<T, M>::build(const Tabular_Data<T> & zdata,
                    const int current_tree_height,
                    const int height_limit)
{
  _mydepth = current_tree_height;

  if (0 == _mydepth)
  {
    _members.set_size(zdata.nrows());
    _members.set_all_items_in_set();
  }

  std::uniform_int_distribution<int> uc(0, zdata.ncols() - 1);

  _split_column = uc(_rng);

  T zmin = zdata.colmin(_split_column);
  T zmax = zdata.colmax(_split_column);

  if (nullptr != _parent)
    _parent->_determine_range(this, _split_column, zmin, zmax);

  if (zmin == zmax)  // no splitting possible
  {
#ifdef DEBUG_IFT_BUILD
    cerr << "IFTree::build:no variability in column " << _split_column << ' ' << zdata.column_name(_split_column) << " at level " << _mydepth << " constant " << zmin << ", have " << _members.number_members() << " members\n";
#endif
    return 1;
  }

  std::uniform_real_distribution<double> uv(zmin, zmax);
  _split_value = uv(_rng);

#ifdef DEBUG_IFT_BUILD
  cerr << "Tree at dpeth " << _mydepth << " with " << _members.number_members() << " members, split on column " << _split_column << " at value " << _split_value << " min " << zmin << " max " << zmax << endl;
#endif

  if (std::is_same<T, double>::value || std::is_same<T, float>::value)
    ;
  else if (_split_value > 0.0)
    _split_value = static_cast<T>(_split_value + 0.4999);
  else
    _split_value = static_cast<T>(_split_value - 0.4999);

  _left  = new IFTree<T, M>;
  _right = new IFTree<T, M>;

  _left->set_parent(this);
  _right->set_parent(this);

  _left->set_size(zdata.nrows());
  _right->set_size(zdata.nrows());

  _left->set_depth(_mydepth+1);   // set the depth in case they do not get any members
  _right->set_depth(_mydepth+1);

  int i = 0;
  int j;
  while ((j = _members.next_member(i)) >= 0)
  {
    if (j >= zdata.nrows())
      break;

    _members.remove_member(j);

#ifdef DEBUG_IFT_BUILD
    if (0 == _mydepth)
      cerr << "Where to put " << zdata(j, _split_column) << " compared to " << _split_value << endl;
#endif

    if (zdata(j, _split_column) < _split_value)
      _left->add_member(j);
    else
      _right->add_member(j);
  }

  _members.clear();

#ifdef DEBUG_IFT_BUILD
  cerr << "Tree at depth " << _mydepth << " LHS " << _left->number_members() << " RHS " << _right->number_members() << endl;
#endif

  if (current_tree_height >= height_limit)
    return 1;

  if (_left->number_members() > 1)
    _left->build(zdata, current_tree_height+1, height_limit);

  if (_right->number_members() > 1)
    _right->build(zdata, current_tree_height+1, height_limit);

  return 1;
}

template <typename T, typename M>
class Isolation_Forest
{
  typedef IFTree<T, M> IFTree_Type;

  private:
    IFTree<T, M> * _tree;
    int _ntrees;
    double _psi;

    int _nrows;     // we need to keep track of how many rows were in the data set used to build us

//  private functions

    int _build_tree(const Tabular_Data<T> & zdata, IFTree_Type & t, const int height_limit);

  public:
    Isolation_Forest();
    ~Isolation_Forest();

    int set_ntrees(const int n);
    void set_psi(const double s) { _psi = s;}

    int print_tree(std::ostream &) const;


    int run_consistency_check() const;

    int gather_tree_depths(extending_resizable_array<int> &) const;

    int build_forest(const Tabular_Data<T> & v);

    double score(const T *) const;
};

template <typename T, typename M>
Isolation_Forest<T, M>::Isolation_Forest()
{
  _tree = nullptr;
  _ntrees = 0;
  _psi = 256.0;

  _nrows = 0;

  return;
}

template <typename T, typename M>
Isolation_Forest<T, M>::~Isolation_Forest()
{
  if (nullptr != _tree)
    delete [] _tree;

  return;
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::set_ntrees(const int n)
{
  if (nullptr != _tree)
    delete [] _tree;

  _tree = new IFTree_Type[n];
  _ntrees = n;

  return 1;
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::build_forest(const Tabular_Data<T> & zdata)
{
  if (0 == _ntrees)
  {
    cerr << "Isolation_Forest::build:ntrees not specified\n";
    return 0;
  }

  _nrows = zdata.nrows();

  const int height_limit = static_cast<int>(log2(_psi) + 0.5);

  for (int i = 0; i < _ntrees; ++i)
  {
    if (verbose > 1)
      cerr << "Building tree " << i << endl;

    _build_tree(zdata,_tree[i], height_limit);
    if (_tree[i].split_column() < 0)
      cerr << "NEgative split column on tree " << i << endl;
  }

  return 1;
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::_build_tree(const Tabular_Data<T> & zdata,
                                 IFTree_Type & t,
                                 const int height_limit)
{
  return t.build(zdata, 0, height_limit);
}

template <typename T, typename M>
double
Isolation_Forest<T, M>::score(const T * d) const
{
  double rc = 0.0;

  for (int i = 0; i < _ntrees; ++i)
  {
    if (_tree[i].split_column() < 0)
      cerr << "Isolation_Forest::score:bad split column, tree " << i << " col " << _tree[i].split_column() << endl;
  }

  for (int i = 0; i < _ntrees; ++i)
  {
    rc += _tree[i].score(d);
  }

#ifdef DEBUG_SCORE_FOREST
  cerr << "Score from " << _ntrees << " trees " << rc << endl;
#endif

  const double Eh = rc / static_cast<double>(_ntrees);

  const double s = pow(2.0, - Eh/ c_function_1(_nrows));

#ifdef DEBUG_SCORE_FOREST
  cerr << "eh " << Eh << " nrows " << _nrows << " c_function_1" << c_function_1(_nrows) << " score " << s << endl;
#endif

  return s;
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::print_tree(std::ostream & output) const
{
  output << "Isolation_Forest with " << _ntrees << " trees\n";
  for (int i = 0; i < _ntrees; ++i)
  {
    output << i << '\n';
    _tree[i].print_tree(' ', output);
  }

  return output.good();
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::gather_tree_depths(extending_resizable_array<int> & acc) const
{
  cerr << "gather_tree_depths::across " << _ntrees << " trees\n";

  for (int i = 0; i < _ntrees; ++i)
  {
    _tree[i].gather_tree_depths(acc);
  }

  return 1;
}

template <typename T, typename M>
int
Isolation_Forest<T, M>::run_consistency_check() const
{
  int failures = 0;

  for (int i = 0; i < _ntrees; ++i)
  {
    if (! _tree[i].cutoffs_consistent())
    {
      cerr << "Tree " << i << " inconsistent\n";
      _tree[i].print_tree(' ', cerr);
      failures++;
    }
  }

  return failures;
}

template <typename T>
int
Tabular_Data<T>::build(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

template <typename T>
int 
get_next_token(const const_IWSubstring & buffer,
               T & s,
               int & i)
{
  if (' ' == input_separator)
    return buffer.nextword(s, i);
  else
    return buffer.nextword_single_delimiter(s, i, input_separator);
}

template <typename T>
int
Tabular_Data<T>::_parse_input_record(const const_IWSubstring & buffer,
                                     const int zrow,
                                     T * v)
{
  const_IWSubstring token;
  int i = 0;

  if (has_rownames)
  {
    if (! get_next_token(buffer, _rowname[zrow], i))
    {
      cerr << "Tabular_Data::_parse_input_record:cannot extract row name " << zrow << endl;
      return 0;
    }
  }

  int c = 0;
  for (; get_next_token(buffer, token, i); ++c)
  {
//  cerr << "Examining " << token << " for column " << c << endl;
    if (! token.numeric_value(v[c]))
    {
      cerr << "Tabular_Data::_parse_input_record:invalid numeric '" << token << " col " << c << " row " << zrow << endl;
      return 0;
    }
  }

  if (c == _ncols)
    return 1;

  if (0 == _ncols)
  {
    if (0 == c)
    {
      cerr << "Tabular_Data::_parse_input_record:no columns in input\n";
      return 0;
    }

    _ncols = c;

    return 1;
  }

  cerr << "Tabular_Data::_parse_input_record:column mismatch, expected " << _ncols << " got " << c << endl;
  return 0;
}


template <typename T>
int
Tabular_Data<T>::build(iwstring_data_source & input)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    cerr << "Tabular_Data::build:cannot read header record\n";
    return 0;
  }

  if (has_header)
  {
    if (! _initialise_column_names(buffer))
      return 0;
  }
  else
  {
    _ncols = buffer.nwords_single_delimiter(input_separator);
    if (! input.seekg(0))
    {
      cerr << "Tabular_Data::build:cannot seek back to start of file\n";
      return 0;
    }

    if (has_rownames)
      _ncols--;
  }

  _nrows = input.records_remaining();

  if (0 == _nrows)
  {
    cerr << "Tabular_Data::build:no data\n";
    return 0;
  }

  if (nullptr != _v)
    delete [] _v;

  _v = new T[_nrows*_ncols];

  if (verbose)
    cerr << "Tabular data contains " << _nrows << " rows by " << _ncols << " columns\n";

  if (0 == _ncols)
  {
    cerr << "Tabular_Data::build:no columns in input\n";
    return 0;
  }

  if (has_rownames)
    _rowname = new IWString[_nrows];

  for (int r = 0; input.next_record(buffer); ++r)
  {
    if (! _parse_input_record(buffer, r, _v + _ncols * r))
    {
      cerr << "Tabular_Data::build:fatal error on line " << r << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  _colmin = new T[_ncols];
  _colmax = new T[_ncols];

  for (int c = 0; c < _ncols; ++c)
  {
    _colmin[c] = _v[c];
    _colmax[c] = _v[c];
  }

  for (int r = 1; r < _nrows; ++r)
  {
    for (int c = 0; c < _ncols; ++c)
    {
      const auto x = _v[r * _ncols + c];

      if (x < _colmin[c])
        _colmin[c] = x;
      else if (x > _colmax[c])
        _colmax[c] = x;
    }
  }

#ifdef ECHO_DATA_QWEQQWE
  for (int i = 0; i < _ncols; ++i)
  {
    cerr << "Col " << i << " " << _colname[i] << " btw " << _colmin[i] << " and " << _colmax[i] << endl;
  }
#endif

  return 1;
}

template <typename T>
int
Tabular_Data<T>::_initialise_column_names(const const_IWSubstring & buffer)
{
  if (' ' == input_separator)
    _ncols = buffer.nwords();
  else
    _ncols = buffer.nwords_single_delimiter(input_separator);

  if (0 == _ncols)
  {
    cerr << "Tabular_Data::_initialise_column_names:no columns\n";
    return 0;
  }

  if (1 == _ncols && 0 == has_rownames)
  {
    cerr << "Tabular_Data::_initialise_column_names:just one column, no data\n";
    return 0;
  }

  if (has_rownames)
    _ncols--;

  if (nullptr != _colname)
    delete [] _colname;

  _colname = new IWString[_ncols];

  int i = 0;
  const_IWSubstring s;
  if (has_rownames)
    get_next_token(buffer, s, i);

  for (int c = 0; get_next_token(buffer, s, i); ++c)
  {
    _colname[c] = s;
  }

  return 1;
}

template <typename T>
int
Tabular_Data<T>::first_direction_out_of_range(const T * x) const
{
  for (int i = 0; i < _ncols; ++i)
  {
//  cerr << "Checking " << x[i] << " vs min " << _colmin[i] << " max " << _colmax[i] << endl;

    if (x[i] < _colmin[i])
      return -1;
    if (x[i] > _colmax[i])
      return 1;
  }

  return 0;
}
template <typename T>
int
Tabular_Data<T>::number_attributes_out_of_range(const T * x,
                                                T & max_out_of_range) const
{
  int rc = 0;

  for (int i = 0; i < _ncols; ++i)
  {
    if (x[i] >= _colmin[i] && x[i] <= _colmax[i])
      continue;

    T oor;
    if (x[i] < _colmin[i])
      oor = _colmin[i] - x[i];
    else
      oor = x[i] - _colmax[i];

    if (_colmax[i] > _colmin[i])
      oor = oor / (_colmax[i] - _colmin[i]);

    if (0 == rc)
      max_out_of_range = oor;
    else if (oor > max_out_of_range)
      max_out_of_range = oor;

    rc++;
  }

  return rc;
}

template <typename T>
int
write_r_variable(const char * vname,
                 const T * v,
                 const int n,
                 IWString_and_File_Descriptor & output)
                
{
  return 1;
}

template <typename T>
int
Tabular_Data<T>::create_DataFrame(const char * vname,
                                  IWString_and_File_Descriptor & output) const 
{
  T * tmp = new T[_nrows]; std::unique_ptr<T[]> free_tmp(tmp);

  for (int c = 0; c < _ncols; ++c)
  {
    output << _colname[c] << "=c(" << _v[c];
    for (int r = 1; r < _nrows; ++r)
    {
      output << ',' << _v[r * _ncols + c];
    }
    output << ")\n";
  }

  output << vname << "=data.frame(";
  for (int c = 0; c < _ncols; ++c)
  {
    if (c > 0)
      output << ',';
    output << _colname[c];
  }
  output << ")\n";

  return 1;
}

static void
report_results(const Accumulator<double> & acc,
               std::ostream & output)
{
  output << out_of_range << " values out of range of reference set\n";

  if (0 == acc.n())
    output << "No values in range of comparator set\n";
  else
    output << "Anomaly values btw " << static_cast<float>(acc.minval()) << " and " << static_cast<float>(acc.maxval()) << " ave " << static_cast<float>(acc.average()) << endl;

  return;
}

class R_Plot
{
  private:
    resizable_array<double> _cutoff;
    resizable_array_p<IWString> _colour;

    int _include_training_set;

    int _pch_train;
    int _pch_test;

    float _cex_train;
    float _cex_test;

    int _plot_area_is_test_data;

    IWString_and_File_Descriptor _output;

//  private functions

    int _display_r_plot_directives(std::ostream &) const;
    int _parse_cutoff_colours(const IWString & s_cutoff, const IWString & s_colour);
    int _write_colour_for_value(const double d);

  public:
    R_Plot();

    int initialise(const Command_Line & cl, const char flag, const int verbose);

    template <typename T> int make_plot(const Tabular_Data<T> & training_set,
                                        const Tabular_Data<T> & zdata,
                                        const std::pair<int, double> * scored);

    template <typename T> int make_plot(const Tabular_Data<T> & training_set,
                                        const Tabular_Data<T> & zdata,
                                        const IF_Result<T> * scored);

    int active() const { return _output.is_open();}
};

R_Plot::R_Plot()
{
  _include_training_set = 1;

  _pch_train = 20;
  _pch_test = 4;

  _cex_train = 0.5;
  _cex_test = 1.0;

  _cutoff.add(0.50);
  _colour.add(new IWString("blue"));
  _colour.add(new IWString("red"));

  _plot_area_is_test_data = 0;

  return;
}

static R_Plot r_plot;    // too many arguments if we pass it down that way...

int
R_Plot::_display_r_plot_directives(std::ostream & output) const
{
  output << "Creates R plot file for Isolation Trees\n";
  output << " cutoff=<x>      cutoff(s) between outliers and non-outliers\n";
  output << " col=<col>       colour(s) for regions defined by cutoff values\n";
  output << " pch_train=<n>   pch value for training set points\n";
  output << " pch_test=<n>    pch value for training set points\n";
  output << " acol=<col>      colour for\n";
  output << " bcol=<col>      colour for\n";
  output << " <fname>         otherwise unrecognised is the file name\n";

  return 1;
}

int
R_Plot::initialise(const Command_Line & cl,
                   const char flag,
                   const int verbose)
{
  IWString fname;
  IWString s_cutoff;
  IWString s_colour;

  const_IWSubstring s;
  for (int i = 0; cl.value(flag, s, i); ++i)
  {
    if (s.starts_with("cutoff="))
    {
      s.remove_leading_chars(7);
      s_cutoff = s;
    }
    else if (s.starts_with("col="))
    {
      s.remove_leading_chars(4);
      s_colour = s;
    }
    else if (s.starts_with("pch_train="))
    {
      s.remove_leading_chars(10);
      if (! s.numeric_value(_pch_train) || _pch_train < 1)
      {
        cerr << "The pch_train value must be a whole +ve number\n";
        return 0;
      }
    }
    else if (s.starts_with("pch_test="))
    {
      s.remove_leading_chars(9);
      if (! s.numeric_value(_pch_test) || _pch_test < 1)
      {
        cerr << "The pch_test value must be a whole +ve number\n";
        return 0;
      }
    }
    else if (s.starts_with("acol="))
    {
      s.remove_leading_chars(5);
      _colour[1] = new IWString(s);
    }
    else if (s.starts_with("bcol="))
    {
      s.remove_leading_chars(5);
      _colour[0] = new IWString(s);
    }
    else if ("ptest" == s)
    {
      _plot_area_is_test_data = 1;
    }
    else if (0 == fname.length())
    {
      fname = s;
    }
    else if ("help" == s)
    {
      _display_r_plot_directives(cerr);
      exit(1);
    }
    else
    {
      cerr << "R_Plot::initialise:unrecognised input " << s << "\n";
      return 1;
    }
  }

  if (0 == fname.length())
  {
    cerr << "R_Plot::initialise:no output file specified\n";
    return 0;
  }

  if (0 == s_cutoff.length() && 0 == s_colour.length())
    ;
  else if (0 == s_cutoff.length() || 0 == s_colour.length())
  {
    cerr << "Must specify both the cutoff values and corresponding colours\n";
    return 0;
  }
  else
  {
    if (! _parse_cutoff_colours(s_cutoff, s_colour))
    {
      cerr << "R_Plot::initialise:cannot parse colour/cutoff specifications\n";
      return 0;
    }
  }

  if (! fname.ends_with(".r"))
    fname << ".r";

  if (! _output.open(fname.null_terminated_chars()))
  {
    cerr << "R_Plot::initialise:cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "R output written to '" << fname << "'\n";

  return 1;
}

int
R_Plot::_parse_cutoff_colours(const IWString & s_cutoff,
                              const IWString & s_colour)
{
  const int ncu = s_cutoff.nwords(',');
  const int nco = s_colour.nwords(',');
  if (ncu + 1 != nco)
  {
    cerr << "Must specify <n> cutoff values and <n+1> colours\n";
    return 0;
  }

  _cutoff.resize_keep_storage(0);
  _colour.resize_keep_storage(0);

  const_IWSubstring token;
  int i = 0;
  while (s_cutoff.nextword(token, i, ','))
  {
    double v;
    if (! token.numeric_value(v) || v <= 0.0 || v >= 1.0)
    {
      cerr << "R_Plot::_parse_cutoff_colours:invalid cutoff '" << token << "'\n";
      return 0;
    }

    if (_cutoff.number_elements() && v <= _cutoff.last_item())
    {
      cerr << "R_Plot::_parse_cutoff_colours:cutoffs out of order " << v << endl;
      return 0;
    }

    _cutoff.add(v);
  }

  i = 0;
  while (s_colour.nextword(token, i, ','))
  {
    _colour.add(new IWString(token));
  }

  return 1;
}

template <typename T>
int
R_Plot::make_plot(const Tabular_Data<T> & training_set,
                  const Tabular_Data<T> & zdata,
                  const std::pair<int, double> * scored)
{
  if (2 != training_set.ncols())
  {
    cerr << "do_produce_R_plot:sorry, only works with 2 columns of data\n";
    return 0;
  }

  const int nrows = zdata.nrows();

  if (_include_training_set)
  {
    training_set.create_DataFrame("t", _output);
    _output << "plot(t,xlab='" << training_set.column_name(0) << "',ylab='" << training_set.column_name(1) << "',cex=" << _cex_train << ",pch=" << _pch_train << ",las=1";
    if (_plot_area_is_test_data)
    {
      _output << ",xlim=c(" << zdata.colmin(0) << ',' << zdata.colmax(0) << "),ylim=c(" << zdata.colmin(1) << ',' << zdata.colmax(1) << ")";
    }
    _output << ")\n";
  }

  zdata.create_DataFrame("e", _output);

  _output << "col=c(";
  for (int i = 0; i < nrows; ++i)
  {
    if (i > 0)
      _output << ',';

    _write_colour_for_value(scored[i].second);
  }
  _output << ")\n";

  _output << "points(e,col=col,pch=" << _pch_test << ",cex=" << _cex_test << ",lwd=2)\n";

  int outliers = 0;
  for (int i = 0; i < nrows; ++i)
  {
    if (scored[i].second > _cutoff[0])
      outliers++;
  }
  _output << "title('Isolation Forest\\n" << outliers << " outliers " << _cutoff[0] << "')\n";

  return 1;
}

template <typename T>
int
R_Plot::make_plot(const Tabular_Data<T> & training_set,
                  const Tabular_Data<T> & zdata,
                  const IF_Result<T> * scored)
{
  if (2 != training_set.ncols())
  {
    cerr << "do_produce_R_plot:sorry, only works with 2 columns of data\n";
    return 0;
  }

  const int nrows = zdata.nrows();

  if (_include_training_set)
  {
    training_set.create_DataFrame("t", _output);
    _output << "plot(t,xlab='" << training_set.column_name(0) << "',ylab='" << training_set.column_name(1) << "',cex=" << _cex_train << ",pch=" << _pch_train << ",las=1";
    if (_plot_area_is_test_data)
    {
      _output << ",xlim=c(" << zdata.colmin(0) << ',' << zdata.colmax(0) << "),ylim=c(" << zdata.colmin(1) << ',' << zdata.colmax(1) << ")";
    }
    _output << ")\n";
  }

  zdata.create_DataFrame("e", _output);

  _output << "col=c(";
  for (int i = 0; i < nrows; ++i)
  {
    if (i > 0)
      _output << ',';

    _write_colour_for_value(scored[i]._score);
  }
  _output << ")\n";

  _output << "points(e,col=col,pch=" << _pch_test << ",cex=" << _cex_test << ",lwd=2)\n";

  int outliers = 0;
  for (int i = 0; i < nrows; ++i)
  {
    if (scored[i]._score > _cutoff[0])
      outliers++;
  }
  _output << "title('Isolation Forest\\n" << outliers << " outliers " << _cutoff[0] << "')\n";

  return 1;
}

int
R_Plot::_write_colour_for_value(const double d)
{
  if (d <= -1.0 || d >= 1.0)   // out of range, gets the fullest outlier colour
  {
    _output << '\'' << *(_colour.last_item()) << '\'';
    return 1;
  }

  if (d <= _cutoff[0])
  {
    _output << '\'' << *(_colour[0]) << '\'';
    return 1;
  }

  for (int i = 1; i < _cutoff.number_elements(); ++i)
  {
    if (d <= _cutoff[i])
    {
      _output << '\'' << (*_colour[i]) << '\'';
      return 1;
    }
  }

  _output << '\'' << *(_colour.last_item()) << '\'';

  return 1;
}

template <typename T>
int
handle_out_of_range(const Tabular_Data<T> & training_set,
                    const Tabular_Data<T> & zdata,
                    const int ndx,
                    const int number_out_of_bounds,
                    const int max_out_of_bounds)
{
  if (stream_for_out_of_range.is_open())
  {
    stream_for_out_of_range << zdata.row_name(ndx) << output_separator << number_out_of_bounds << output_separator << max_out_of_bounds << '\n';
    stream_for_out_of_range.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

#ifdef OBSOLETE_NOT_USED_____123
template <typename T, typename M>
int
do_isolation_forest_sorted(const Tabular_Data<T> & training_set,
                           Isolation_Forest<T, M> & ifrst,
                           const Tabular_Data<T> & zdata,
                           IWString_and_File_Descriptor & output)
{
  Accumulator<double> acc;
  const int nrows = zdata.nrows();

  typedef std::pair<IWString, double> Scored;

  auto * scored = new Scored[nrows]; std::unique_ptr<Scored[]> free_scored(scored);

  int ndx = 0;

  for (int i = 0; i < nrows; ++i)
  {
    const T * d = zdata.row(i);

    const auto x = training_set.first_direction_out_of_range(d);
    if (0 != x)
    {
      handle_out_of_range(training_set, zdata, i, x);
      out_of_range++;
//    continue;
    }

    const double s = ifrst.score(d);

    acc.extra(s);

    if (s > 0.5)
      above5++;

    scored[ndx].first = zdata.row_name(i);
    scored[ndx].second = s;
    ndx++;
  }

  if (0 == ndx)
  {
    report_results(acc, cerr);
    return 1;
  }

  std::sort(scored, scored + ndx, [](const Scored & s1, const Scored & s2) { return s1.second > s2.second;});

  for (int i = 0; i < ndx; ++i)
  {
    output << scored[i].first << output_separator << scored[i].second << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  report_results(acc, cerr);

  return 1;
}
#endif

template <typename T>
int
write_row(const T * v,
          const int n,
          IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < n; ++i)
  {
    output << output_separator << v[i];
  }

  return 1;
}

template <typename T>
int
do_write_header(const Tabular_Data<T> & zdata,
                IWString_and_File_Descriptor & output)
{
  output << "ID";
  if (include_raw_data_with_output)
  {
    if (has_header)
    {
      for (int i = 0; i < zdata.ncols(); ++i)
      {
        output << output_separator << zdata.column_name(i);
      }
    }
    else
    {
      for (int i = 0; i < zdata.ncols(); ++i)
      {
        output << output_separator << 'Y' << i;
      }
    }
  }
  output << output_separator << "OutlierScore";
  if (extra_columns_for_out_of_range_info)
    output << output_separator << "N_Out_of_Range" << output_separator << "Max_Out_of_Range";
  output << '\n';

  return 1;
}

template <typename T>
int
do_write_original_data(const Tabular_Data<T> & training_set,
                       IWString_and_File_Descriptor & output)
{
  for (int r = 0; r < training_set.nrows(); ++r)
  {
    if (has_rownames)
      output << training_set.row_name(r);
    write_row(training_set.row(r), training_set.ncols(), output);
    output << output_separator << "0.5\n";
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

template <typename T, typename M>
int
do_isolation_forest(const Tabular_Data<T> & training_set,
                    Isolation_Forest<T, M> & ifrst,
                    const Tabular_Data<T> & zdata,
                    IWString_and_File_Descriptor & output)
{
  if (include_header_in_output)
    do_write_header(zdata, output);

  if (include_original_data_in_output)
    do_write_original_data(training_set, output);

//if (produce_sorted_output)
//  return do_isolation_forest_sorted(training_set, ifrst, zdata, output);

  const int nrows = zdata.nrows();

  IF_Result<T> * scores = new IF_Result<T>[nrows]; std::unique_ptr<IF_Result<T>[]> free_scores(scores);

#pragma omp parallel for schedule(dynamic,256)
  for (int i = 0; i < nrows; ++i)
  {
    IF_Result<T> & r = scores[i];
    
    r._ndx = i;

    const T * d = zdata.row(i);

    r._out_of_bounds = training_set.number_attributes_out_of_range(d, r._max_out_of_bounds);

    r._score = ifrst.score(d);

#if ! defined(OMP)
    if (report_progress())
      cerr << "Scored " << i << " of " << nrows << " points\n";
#endif
  }

  if (produce_sorted_output)
    std::sort(scores, scores + nrows, [] (const IF_Result<T> & s1, const IF_Result<T> & s2) { return s1._score > s2._score;});

  Accumulator<double> acc;

  int above5 = 0;     // the number of scores above 0.5

  extending_resizable_array<int> all_scores;     // converted to int

  for (int i = 0; i < nrows; ++i)
  {
    const IF_Result<T> & r = scores[i];

    const int ndx = r._ndx;

    if (r._out_of_bounds > 0)
    {
      handle_out_of_range(training_set, zdata, ndx, r._out_of_bounds, r._max_out_of_bounds);
      out_of_range++;
      if (! include_out_of_range_values_in_regular_output)
        continue;
    }

    const double s = r._score;

    acc.extra(s);
    all_scores[static_cast<int>(100.0 * s + 0.4999)]++;
    if (s > 0.5)
      above5++;

    if (has_rownames)
      output << zdata.row_name(ndx);
    else
      output << ndx;

    if (include_raw_data_with_output)
      write_row(zdata.row(ndx), zdata.ncols(), output);

    output << output_separator << s;
    if (extra_columns_for_out_of_range_info)
      output << output_separator << r._out_of_bounds << output_separator << r._max_out_of_bounds;
    output << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  if (verbose)
  {
    report_results(acc, cerr);
    cerr << above5 << " of " << acc.n() << " scores above 0.5 " << static_cast<float>(above5) / static_cast<float>(acc.n()) << endl;

    for (int i = 0; i < all_scores.number_elements(); ++i)
    {
      if (0 == all_scores[i])
        continue;

      float x = static_cast<float>(i) / 100.0f;

      cerr << all_scores[i] << " values near " << x << endl;
    }
  }

  if (r_plot.active())
    r_plot.make_plot(training_set, zdata, scores);

  return 1;
}


template <typename T, typename M>
int
do_isolation_forest(const char * reference_set,
                    const int ntrees,
                    const double psi,
                    const char * comparison_set,
                    IWString_and_File_Descriptor & output)
{
  Tabular_Data<T> zdata;
  if (! zdata.build(reference_set))
  {
    cerr << "do_isolation_forest:cannot read data '" << reference_set << "'\n";
    return 0;
  }

  if (verbose)
  {
    cerr << "From " << reference_set << endl;
    for (int c = 0; c < zdata.ncols(); ++c)
    {
      cerr << "Column " << c << " ";
      if (has_header)
        cerr << zdata.column_name(c);
      cerr << " btw " << zdata.colmin(c) << " and " << zdata.colmax(c) << endl;
    }
  }

  Isolation_Forest<T, M> ifrst;

  ifrst.set_ntrees(ntrees);
  ifrst.set_psi(psi);

  if (! ifrst.build_forest(zdata))
  {
    cerr << "do_isolation_forest:cannot build forest\n";
    return 0;
  }

  if (print_tree)
  {
    ifrst.print_tree(cerr);
    extending_resizable_array<int> depths;

    ifrst.gather_tree_depths(depths);

    Accumulator_Int<int> acc;

    for (int i = 0; i < depths.number_elements(); ++i)
    {
      cerr << depths[i] << " terminal nodes at depth " << i << '\n';
      if (depths[i] > 0)
        acc.extra(i, depths[i]);
    }

    cerr << "Depths btw " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
  }

  if (run_consistency_check)
  {
    const auto failures = ifrst.run_consistency_check();

    if (failures)
    {
      cerr << "Isolation forest inconsistent\n";
      return 0;
    }
  }

  if (reference_set == comparison_set)
    cerr << "Ref and CMP the same\n";

  if (reference_set == comparison_set)
    return do_isolation_forest(zdata, ifrst, zdata, output);

  Tabular_Data<T> to_score;
  if (! to_score.build(comparison_set))
  {
    cerr << "Cannot build comparator set '" << comparison_set << "'\n";
    return 0;
  }

  if (! to_score.same_descriptor_names(zdata))
  {
    cerr << "Descriptor name mismatch between reference set and comparator\n";
    return 0;
  }

  if (verbose)
    cerr << "Comparison set contains " << to_score.nrows() << " rows of data\n";

  return do_isolation_forest(zdata, ifrst, to_score, output);
}

static int
parse_separator(const_IWSubstring & s,
                char & sep)
{
  if (1 == s.length())
    sep = s[0];
  else if ("space" == s)
    sep = ' ';
  else if ("tab" == s)
    sep = '\t';
  else if ("comma" == s)
    sep = ',';
  else
  {
    return 0;
  }

  return 1;
}

template int do_isolation_forest<double, Members_as_Bits>(char const*, int, double, char const*, IWString_and_File_Descriptor&);
template class Isolation_Forest<double, Members_as_Bits>;
template class IFTree<double, Members_as_Bits>;

static int
isolation_forest (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vn:C:i:o:p:PcsyxX:OR:eEuwz:r:h:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('i'))
  {
    const_IWSubstring s = cl.string_value('i');
    if (! parse_separator(s, input_separator))
    {
      cerr << "Unrecognised input separator specification '" << s << "'\n";
      return 1;
    }
  }

  if (cl.option_present('o'))
  {
    const_IWSubstring s = cl.string_value('o');
    if (! parse_separator(s, output_separator))
    {
      cerr << "Unrecognised output separator specification '" << s << "'\n";
      return 1;
    }
  }

  if (cl.option_present('P'))
  {
    print_tree = 1;

    if (verbose)
      cerr << "Will print the tree\n";
  }

  int ntrees = 100;
  if (cl.option_present('n'))
  {
    if (! cl.value('n', ntrees) || ntrees < 1)
    {
      cerr << "The number of trees mst be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will generate " << ntrees << " trees\n";
  }

  double psi = 256.0;
  if (cl.option_present('p'))
  {
    if (! cl.value('p', psi) || psi <= 1.0)
    {
      cerr << "The psi value (-p) must be a value above 1.0\n";
      usage(1);
    }

    if (verbose)
      cerr << "psi set to " << psi << endl;
  }

#ifdef OMP
  if (cl.option_present('h'))
  {
    int h;
    if (! cl.value('h', h) || h < 0)
    {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve number\n";
      usage(2);
    }

    omp_set_num_threads(h);
  }
  else
  {
    omp_set_num_threads(1);
  }
#endif

  if (cl.option_present('c'))
  {
    run_consistency_check = 1;

    if (verbose)
      cerr << "Will run a consitency check on the tree\n";
  }

  if (cl.option_present('s'))
  { 
    if (cl.option_present('w'))
    {
      cerr << "If the input lacks row names (-w), sorted output (-s) not possible\n";
      return 1;
    }

    produce_sorted_output = 1;

    if (verbose)
      cerr << "Will produce sorted output\n";
  }

  if (cl.option_present('x'))
  {
    include_out_of_range_values_in_regular_output = 0;

    if (verbose)
      cerr << "Will exclude out of range values from the output\n";
  }

  if (cl.option_present('O'))
  {
    extra_columns_for_out_of_range_info = 1;

    if (verbose)
      cerr << "Will include extra columns with out of range information\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 1;
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    if (cl.option_present('w'))
    {
      cerr << "If data does not have row names (-w) does not make sense to write outliers (-X)\n";
      return 1;
    }

    if (cl.option_present('O'))
    {
      cerr << "Cannot use both the -X and -O options\n";
      return 1;
    }

    if (! stream_for_out_of_range.open(x))
    {
      cerr << "Cannot open stream for out of range values '" << x << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Will write out of range values to " << x << "'\n";

    stream_for_out_of_range << "ID FirstOutOfRange NumberOutOfRange\n";

    include_out_of_range_values_in_regular_output = 0;
  }

  if (cl.option_present('R'))
  {
    if (! r_plot.initialise(cl, 'R', verbose))
    {
      cerr << "Cannot initialise R data\n";
      return 1;
    }
  }

  if (cl.option_present('e'))
  {
    if (cl.option_present('s') || cl.option_present('x') || cl.option_present('X'))
    {
      cerr << "Sorry, the -e/-x and -s options are not compatible\n";
      usage(1);
    }

    include_raw_data_with_output = 1;

    if (verbose)
      cerr << "Raw data included with output\n";
  }

  if (cl.option_present('E'))
  {
    include_original_data_in_output = 1;

    if (verbose)
      cerr << "Will include training data with output\n";
  }

  const char * reference_set;
  const char * comparison_set;
  if (cl.option_present('C'))
  {
    reference_set = cl.option_value('C');
    comparison_set = cl[0];

    if (verbose)
      cerr << "Assessment will be with respect to " << reference_set << endl;
  }
  else
  {
    reference_set = cl[0];
    comparison_set = cl[0];
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('y'))
  {
    include_header_in_output = 0;
  }

  if (cl.option_present('u'))
  {
    has_header = 0;

    if (verbose)
      cerr << "Input assumed to lack header record\n";
  }

  if (cl.option_present('z'))
  {
    int z;
    if (! cl.value('z', z) || z < 2)
    {
      cerr << "The default floating point precision option (-z) must be a whole +ve number\n";
      usage(1);
    }

    set_default_iwstring_float_concatenation_precision(z);
    set_default_iwstring_double_concatenation_precision(z);
  }

  if (cl.option_present('w'))
  {
    has_rownames = 0;

    if (verbose)
      cerr << "Input assumed to lack row names\n";
  }

  int rc = 0;

  if (cl.option_present('T'))
  {
    const_IWSubstring t = cl.string_value('T');
    if ("int" == t)
    {
      rc = do_isolation_forest<int, Members_as_Bits>(reference_set, ntrees, psi, comparison_set, output);
    }
    else if ("float" == t)
    {
      rc = do_isolation_forest<float, Members_as_Bits>(reference_set, ntrees, psi, comparison_set, output);
    }
    else
    {
      cerr << "Unrecognised -T qualifier '" << t << "'\n";
      return 1;
    }
  }
  else
  {
    rc = do_isolation_forest<double, Members_as_Bits>(reference_set, ntrees, psi, comparison_set, output);
  }

  if (verbose)
  {
  }

  if (rc)
    return 0;

  return 1;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = isolation_forest(argc, argv);

  return rc;
}
