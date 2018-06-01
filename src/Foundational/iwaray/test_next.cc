#include <stdlib.h>

/*
  Tester for the next function
*/

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_next_.h"

class Foo
{
  private:
    int _f;
  public:
    Foo (int);

    int f () const { return _f;}
};

Foo::Foo (int f) : _f(f)
{
}

template resizable_array_p<Foo>;
template resizable_array_base<Foo *>;

class Foo_Comparitor
{
  private:
    int _n;
  public:
    Foo_Comparitor (int x) : _n (x) {};

    int operator () (const Foo &) const;

    int n () const { return _n;}
};

int
Foo_Comparitor::operator () (const Foo & f) const
{
  return _n == f.f ();
}

class Multiples
{
  protected:
    int _n;
  public:
    Multiples (int n) : _n (n) {}

    int n () const { return _n;}
};

class Detect_Multiples : public Multiples
{
  private:
  public:
    Detect_Multiples (int);

    int operator () (int) const;
};

Detect_Multiples::Detect_Multiples (int x) : Multiples (x)
{
}

int
Detect_Multiples::operator () (int i) const
{
  return i % _n;
}

class Count_Multiples : public Multiples
{
  private:
    int _sum;

  public:
    Count_Multiples (int);

    int operator () (int);

    int sum () const { return _sum;}
};

Count_Multiples::Count_Multiples (int x) : Multiples (x)
{
  _sum = 0;

  return;
}

int
Count_Multiples::operator () (int i) 
{
  int rc = i % _n;

  if (0 == rc)
    _sum++;

  return rc;
}

class Increment
{
  private:
    int _n;
  public:
    Increment (int x) : _n (x) {};

    int n () const { return _n;}

    void operator () (int &) const;
};

void
Increment::operator () (int & i) const
{
  i += _n;
}


static int verbose = 1;

int
main ()
{
  int n = 100;

  resizable_array<int> x (n);

  for (int i = 0; i < n; i++)
  {
    x.add (i);
  }

  int matches = 0;

  Detect_Multiples detect_multiples (5);

  resizable_array_base<int>::next_iterator i = x.next_begin (); 
  while (x.next (detect_multiples, i))
  {
    if (verbose)
      cerr << "Matched " << i << endl;
    matches++;
  }

  cerr << "Found " << matches << " multiples of " << detect_multiples.n () << endl;

  Count_Multiples count_multiples (10);

  i = x.next_begin ();

  while (x.next (count_multiples, i))
  {
    cerr << "Match for i = " << i << " value " << x[i] << endl;
  }

  cerr << "Found " << count_multiples.sum () << " multiples of " << count_multiples.n () << endl;

  Increment increment (5);

  x.each (increment);

  for (int i = 0; i < x.number_elements (); i++)
  {
    if (x[i] != i + 5)
    {
      cerr << "Each failure, i = " << i << " got " << x[i] << " expected " << (i + 5) << endl;
    }
  }

  resizable_array_p<Foo> s;

  for (int i = 0; i < n; i++)
  {
    Foo * tmp = new Foo (i);

    s.add (tmp);
  }

  Foo_Comparitor foo_comparitor (5);

  resizable_array_p<Foo>::next_iterator j = s.next_begin ();
  while (s.next (foo_comparitor, j))
  {
    cerr << "Got Foo match for " << foo_comparitor.n () << " at j = " << j << endl;
  }

  return 0;
}
