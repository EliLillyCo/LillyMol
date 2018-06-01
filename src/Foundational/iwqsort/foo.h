#ifndef IWFOO_H
#define IWFOO_H

/*
  Sometimes we want to experiment with varying the size of things
*/

#define FOO_ARRAY_SIZE 3

template <typename T>
class Foo
{
  private:
    T _value;
    T _notused;

    static T _static_variable;

    T _foo_array[FOO_ARRAY_SIZE];

  public:
    Foo ();
    
    void set_value (T v) { _value = v;};

    T zvalue () const { return _value;}

    int iwqsortcompare (const Foo &) const;

    static int iwqsort_mfn (Foo<T> & f1, Foo<T> & f2) { return f1.iwqsortcompare (f2);}
//  static int (Foo<T>::*iwqsort_mfn) (const Foo<T> &) const;

//  void set_compare_ascending  () { iwqsort_mfn = &(Foo::compare_ascending);}
//  void set_compare_descending () { iwqsort_mfn = &(Foo::compare_descending);}

    static int compare_ascending  (const Foo<T> &);
    static int compare_descending (const Foo<T> &);
};

extern int foo_comparitor_int   (const Foo<int> &, const Foo<int> &);
extern int foo_comparitor_float (const Foo<float> &, const Foo<float> &);

#endif
