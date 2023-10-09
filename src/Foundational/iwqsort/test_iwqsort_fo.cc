/*
  Test function object stuff
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwrandom.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "iwqsort.h"

using std::cerr;
using std::ostream;

static int verbose = 0;

static int range = 100;

static int n = 10;

static void
usage (int rc)
{
  exit (rc);
}

// Found no advantage to using friend declarations

//#define USE_FRIEND_DECLARATIONS

#ifdef USE_FRIEND_DECLARATIONS
class Foo_Comparitor_Pointer;
class Foo_Comparitor;
#endif

class Foo
{
#ifdef USE_FRIEND_DECLARATIONS
  friend Foo_Comparitor_Pointer;
  friend Foo_Comparitor;
#endif

  public:
    static int _nFoo;

  private:

    int _unique_id;

    IWString _string_rep;

    int _sort_key;

  public:
    Foo ();
    ~Foo ();

    int debug_print (ostream &) const;

    int unique_id () const { return _unique_id;}

    int sort_key () const { return _sort_key;}

    const IWString & string_rep () const { return _string_rep;}
};

Foo::Foo ()
{
  _unique_id = _nFoo;
  _nFoo++;

  _string_rep << _unique_id;

  (void) _string_rep.null_terminated_chars ();     // force a terminator

  _sort_key = intbtwij (0, range);

  return;
}

Foo::~Foo ()
{
  if (-77373 == _sort_key)
    cerr << "Deleting already deleted Foo " << _unique_id << endl;

  _sort_key = -77373;

  return;
}

int
Foo::debug_print (ostream & os) const
{
  os << _string_rep << " id " << _unique_id << " sort key " << _sort_key << endl;

  return os.good ();
}

ostream &
operator << (ostream & os, const Foo & foo)
{
  os << foo.unique_id () << " sort key " << foo.sort_key ();

  return os;
}

ostream & 
operator << (ostream & os, const Foo * foo)
{
  os << foo->unique_id () << " sort key " << foo->sort_key ();

  return os;
}

int Foo::_nFoo = 0;

class Foo_Comparitor_Pointer
{
  private:
  public:
    int operator () (const Foo *, const Foo *);
};

class Foo_Comparitor
{
  private:
  public:
    int operator () (const Foo &, const Foo &);
};

#ifdef USE_FRIEND_DECLARATIONS
int
Foo_Comparitor_Pointer::operator () (const Foo * f1, const Foo * f2)
{
  if (f1->_sort_key < f2->_sort_key)
    return -1;

  if (f1->_sort_key > f2->_sort_key)
    return 1;

  return 0;
}

int
Foo_Comparitor::operator () (const Foo & f1, const Foo & f2)
{
  if (f1._sort_key < f2._sort_key)
    return -1;

  if (f1._sort_key > f2._sort_key)
    return 1;

  return 0;
}
#else
int
Foo_Comparitor_Pointer::operator () (const Foo * f1, const Foo * f2)
{
  if (f1->sort_key () < f2->sort_key ())
    return -1;

  if (f1->sort_key () > f2->sort_key ())
    return 1;

  return 0;
}

int
Foo_Comparitor::operator () (const Foo & f1, const Foo & f2)
{
  if (f1.sort_key () < f2.sort_key ())
    return -1;

  if (f1.sort_key () > f2.sort_key ())
    return 1;

  return 0;
}

#endif

template void iwqsort<Foo *, Foo_Comparitor_Pointer>(Foo * *, int, Foo_Comparitor_Pointer &);
template void iwqsort<Foo, Foo_Comparitor>(Foo *, int, Foo_Comparitor &);

static int
test_array_of_pointers (int n)
{
  Foo ** f = new Foo *[n];

  if (nullptr == f)
  {
    cerr << "Memory failure trying to allocate " << n << " Foo pointers\n";
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    f[i] = new Foo;
  }

  Foo_Comparitor_Pointer foo_comparitor_pointer;
  iwqsort (f, n, foo_comparitor_pointer);

  if (verbose > 1)
    f[0]->debug_print (cerr);

  int failures = 0;

  for (int i = 1; i < n; i++)
  {
    if (verbose > 1)
    {
      cerr << " i = " << i << ' ';
      f[i]->debug_print (cerr);
    }

    if (f[i]->sort_key () < f[i - 1]->sort_key ())
    {
      cerr << "Out of order, i = " << i << " i - 1 is " << f[i - 1]->sort_key () << " i is " << f[i]->sort_key () << endl;
      failures++;
    }

    delete f[i - 1];
  }

  delete f[n - 1];
  
  delete [] f;

  if (failures)
    return 0;

  if (verbose)
    cerr << "Array of pointers test successful\n";

  return 1;
}

static int
test_array (int n)
{
  Foo * f = new Foo[n];

  if (nullptr == f)
  {
    cerr << "Bad news, cannot allocate array of " << n << " Foo objects\n";
    return 0;
  }

  Foo_Comparitor foo_comparitor;
  iwqsort (f, n, foo_comparitor);

  if (verbose > 1)
    f[0].debug_print (cerr);

  int failures = 0;
  for (int i = 1; i < n; i++)
  {
    if (verbose > 1)
    {
      cerr << " i = " << i << ' ';
      f[i].debug_print (cerr);
    }

    if (f[i].sort_key () < f[i - 1].sort_key ())
    {
      cerr << "Out of order, i = " << i << " i - 1 is " << f[i - 1].sort_key () << " i is " << f[i].sort_key () << endl;
      failures++;
    }
  }

  delete [] f;

  if (failures)
    return 0;

  if (verbose)
    cerr << "Array of objects test successful\n";

  return 1;
}

static int
test_iwqsort_fo (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', n) || n < 2)
    {
      cerr << "Invalid array size (-n option)\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will process arrays of size " << n << endl;
  }

  iw_random_seed ();

  int failures = 0;
  if (! test_array_of_pointers (n))
    failures++;

  if (! test_array (n))
    failures++;

  cerr << Foo::_nFoo << " objects created\n";
  
  if (failures)
  {
    cerr << failures << " of two tests failed\n";
    return 5;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_iwqsort_fo (argc, argv);

  return rc;
}
