#include <stdlib.h>
#include <utility>
#include <unistd.h>
#include <assert.h>
#include <random>
#include <memory>

#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWARAY_EACH_IMPLEMENTATION
#include "iwaray.h"

int verbose = 0;

static int test_failures = 0;

static int nfoos = 0;    // global counter of the number of Foo objects created

#ifdef QQQQQ
static int
speed_test_index(const int ntest,
                 const int n)
{
  resizable_array<int> v(n);

  for (int i = 0; i < n; ++i)
  {
    v.add(i);
  }

  int rc = 1;

  for (int i = 0; i < ntest; ++i)
  {
    int sum = 0;
    for (int j = 0; j < n; ++j)
    {
      sum += v[j];
    }

    if (sum < n + 10)
      rc = 0;
  }

  return rc;
}

static int
speed_test_iterator(const int ntest,
                    const int n)
{
  resizable_array<int> v(n);

  for (int i = 0; i < n; ++i)
  {
    v.add(i);
  }

  int rc = 1;

  for (int i = 0; i < ntest; ++i)
  {
    int sum = 0;
    for (auto j : v)
    {
      sum += j;
    }

    if (sum < n + 10)
      rc = 0;
  }

  return rc;
}
#endif

static int
test_move ()
{
  const int n = 10;

  resizable_array<int> * x = new resizable_array<int>[n]; std::unique_ptr<resizable_array<int>[]> free_x(x);

  int * initial_size = new int[n]; std::unique_ptr<int[]> free_initial_size(initial_size);

  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<unsigned int> s(1,10);

  for (int i = 0; i < n; ++i)
  {
    initial_size[i] = s(rng);

    x[i].resize(initial_size[i]);
    for (int j = 0; j < initial_size[i]; ++j)
    {
      x[i].add(i + j);
    }
  }

// shift everything down one 

  for (int i = 0; i < (n-1); ++i)
  {
    x[i] = std::move(x[i+1]);
  }

  int ok = 1;

  for (int i = 0; i < (n-1); ++i)
  {
    if (initial_size[i+1] != x[i].number_elements())
    {
      cerr << "Move failed, i = " << i << " expected " << initial_size[i+1] << " items, got " << x[i].number_elements() << endl;
      ok = 0;
    }
  }

  if (0 != x[n-1].number_elements())
  {
    cerr << "Move failed, last array item not empty, has " << x[n-1].number_elements() << " items\n";
    ok = 0;
  }

  if (! ok)
  {
    test_failures++;
    return 0;
  }

  return 1;
}

static int
test_remove_two_items ()
{
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> s(10,1000);

  for (int i = 0; i < 100; ++i)
  {
    int asize = s(rng);
    resizable_array<int> a(asize);

    std::uniform_int_distribution<int> js(0, asize * 100);

    while (a.number_elements() < asize)
    {
      a.add_if_not_already_present(js(rng));
    }

    std::uniform_int_distribution<int> us(0, asize - 1);

    const int ndx1 = us(rng);
    int ndx2 = us(rng);
    while (ndx1 == ndx2)
    {
      ndx2 = us(rng);
    }

    const int v1 = a[ndx1];
    const int v2 = a[ndx2];

    resizable_array<int> acopy(a);

    a.remove_two_items(v1, v2);

    if (a.number_elements() + 2 != acopy.number_elements())
    {
      cerr << "remove_two_items failed to remove two elements, array size " << asize << ", values " << v2 << " and " << v2 << endl;
      for (int j = 0; j < asize; ++j)
      {
        cerr << " i = " << j << " value " << a[j];
        if (j == ndx1 || j == ndx2)
          cerr << " *";
        cerr << endl;
      }
      test_failures++;
      return 0;
    }

    int offset = 0;
    for (int j = 0; j < a.number_elements(); ++j)
    {
      if (j == ndx1 || j == ndx2)
        offset++;

      if (a[j] == acopy[j + offset])
        continue;

      cerr << "remove_two_items mismatch at j = " << j << " ndx1 " << ndx1 << " ndx2 " << ndx2 << endl;
      test_failures++;
      return 0;
    }

  }

  return 1;
}

static int
test_remove()
{
  resizable_array<int> a;
  int asize = 101;

// fill the array with disinct elements (very important)

  for (int i = 0; i < asize; i++)
  {
    a.add(i);
  }

  const int nremove = 11;
  resizable_array<int> to_remove;
  to_remove.resize(nremove);

  for (int i = 0; i <= 100; i += 10)
  {
    to_remove.add(i);
  }

// Make an array of what should be removed - this is why the values must be distinct

  resizable_array<int> should_be_removed;
  should_be_removed.resize(nremove);
  for (int i = 0; i < nremove; i++)
  {
    int j = to_remove[i];
    should_be_removed.add(a[j]);
//  cerr << "Will remove element " << j << " which is " << a[j] << endl;
  }

// Unfortunately, I could not get remove_items to compile as a resizable array,
// so make an array

  int * tmp = new int[nremove]; std::unique_ptr<int[]> free_tmp(tmp);
  for (int i = 0; i < nremove; i++)
  {
    tmp[i] = to_remove[i];
  }

  a.remove_items(tmp, nremove);

  if (a.number_elements() != (asize - nremove))
  {
    cerr << "Remove items, incorrect count. Array started with " << asize << " members\n";
    cerr << "Removed " << nremove << " but now has " << a.number_elements() << " members\n";
    test_failures++;
    return 99;
  }

  int ok = 1;
  for (int i = 0; i < a.number_elements(); i++)
  {
    if (should_be_removed.contains(a[i]))
    {
      cerr << "OOps, element " << i << " (" << a[i] << ") still in array\n";
      ok = 0;
    }
  }

  if (! ok)
  {
    test_failures++;
    return 0;
  }

  return 1;
}

class Foo
{
  private:
    int _id;

  public:
    Foo();
    ~Foo();

    int id() const { return _id;}

    void voidfun     ();
    void constvoidfun() const;
    int  intfun      ();
    int  constintfun () const;
};

static void
voiddofoo (Foo & f)
{
  cerr << "Void dofoo accessing " << f.id() << '\n';
}

Foo::Foo()
{
  _id = nfoos++;
}

Foo::~Foo()
{
  if (-5 == _id)
    cerr << "Deleting an already deleted foo object\n";
//else
//  cerr << "Deleting foo " << _id << endl;

  _id = -5;
}

void
Foo::constvoidfun() const
{
  cerr << "Const voidfun for object " << _id << endl;
}

void
Foo::voidfun()
{
  cerr << "voidfun for object " << _id << endl;
}

int
Foo::intfun()
{
  cerr << "intfun for object " << _id << endl;

  return _id;
}

int
Foo::constintfun() const
{
  cerr << "const intfun for object " << _id << endl;

  return _id;
}

template void resizable_array_p<Foo>::each(void (&)(Foo&));

//template class resizable_array<int>;
//template class resizable_array_base<int>;

static int
test_int()
{
  resizable_array<int> mtest;

  for (int i = 0; i < 100; i++)
  {
    mtest.add(i);
  }

  for (int i = 0; i < mtest.number_elements(); i++)
  {
    int j = mtest[i];
    if (verbose)
      std::cout << "item " << i << " is " << j << "\n";
    assert (i == j);
  }

  const int items = mtest.number_elements();

  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> u(0, items - 2);

// choose a random element to remove

  const int ii = u(rng);

  mtest.remove_item(ii);

  assert (items - 1 == mtest.number_elements());
  assert (ii + 1 == mtest.item(ii));

  mtest.remove_item(0);
  assert (items - 2 == mtest.number_elements());

//int j = intbtwij(0, mtest.number_elements());
//assert (j + 2 == mtest[j]);

  assert (0 == mtest.remove_all(-1));

//assert (1 == mtest.remove_all(mtest[j]));
  
//assert (items - 3 == mtest.number_elements());

  mtest.resize(0);

  assert (0 == mtest.number_elements());

  const int tsize = 10;

  for (int i = 0; i <= tsize; i++)
    if (0 == i % 2)
      mtest.add(i);

  assert (tsize / 2 + 1 == mtest.number_elements());

  for (int i = 0; i < tsize; i++)
  {
    if (0 != i % 2)
      mtest.insert_before(i, i);
  }

  assert (tsize + 1 == mtest.number_elements());

  int ok = 1;

  for (int i = 0; i < tsize; i ++)
  {
    if (i != mtest[i])
    {
      cerr << "Test failed, item " << i << " is " << mtest.item(i) << "\n";
      ok = 0;
    }
  }

  if (! ok)
  {
    test_failures++;
    return 0;
  }

  return 1;
}

//#ifdef __GNUG__
template class resizable_array_p<Foo>;
template class resizable_array_base<Foo *>;
//#endif

static int
test_ptr()
{
  resizable_array_p<Foo> myfoos;
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> u(2, 20);

  const int nf = u(rng);
  cerr << "Array sized to " << nf << '\n';

  for (int i = 0; i < nf; i++)
  {
    Foo * tmp = new Foo;
    myfoos.add(tmp);
  }

  assert (nf == myfoos.number_elements());

  myfoos.each(voiddofoo);

  return 0;
}

static int 
test_iwaray (int argc, char **argv)
{
#ifndef _WIN32
  int o;
  while ((o = getopt(argc, argv, "v")) != EOF)
    switch(o)
    {
      case 'v':
        verbose = 1;
        break;

      default:
        cerr << "Unrecognised option '" << o << "'\n";
        return 1;
    }
#endif

  test_int();

  test_ptr();

  test_remove();

  test_remove_two_items();

  test_move();

//if (! speed_test_iterator(10000000, 514))
//  return 1;

#ifdef _WIN32
  int i;
  cin >> i;
#endif

  if (test_failures > 0)
    cerr << test_failures << " test failures\n";
  else
    cerr << "all tests successful\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_iwaray(argc, argv);

  return rc;
}
