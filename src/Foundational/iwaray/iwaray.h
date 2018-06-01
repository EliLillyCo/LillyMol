#ifndef IWARAY_H
#define IWARAY_H

//#include <tbb/scalable_allocator.h>

#include <iostream>
#include "iwconfig.h"

#if (__GNUC__ == 3) && (__GNUC_MINOR__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

#if (__GNUC__ >= 4)
// Avoid redefinition warning
#ifndef IW_TWO_PHASE_TEMPLATES
#define IW_TWO_PHASE_TEMPLATES
#endif
#endif

/*
  Frequently we have arrays of objects. They always have
    number of array elements allocated.
    number of array elements used.
    pointer to allocated storage for array elements.

  When the number of elements is reduced to 0, we free the array -
  it is the responsibility of the caller to free the upper structure.
*/

#define IWARAY_MAGIC_NUMBER 4444444

template <class T> class resizable_array_iterator;
template <class T> class resizable_array_base;

/*
  We want the iterator on next to be the index of the most recently processed
  item, so we need to start it with a negative value. Hide that from the user
  by allowing a call to next_begin () which returns an arbitrary negative number
*/

#define IWARAY_NEXT_ITERATOR_START -76704

template <typename T>
class resizable_array_base 
{
  friend
    class resizable_array_iterator<T>;

  protected:
    T * _things;
    int _number_elements;
    int _elements_allocated;
    int _magic;
  public:
    resizable_array_base        ();
    ~resizable_array_base       ();

    int ok() const;
    int debug_print            (std::ostream &) const;

    inline int number_elements  () const { return _number_elements; }
    inline unsigned int size    () const { return _number_elements; }
    inline int elements_allocated () const { return _elements_allocated;}

    int ok_index (int _i) const { return _i >= 0 && _i < _number_elements;}

    inline int array_is_full    () const 
                            { return _number_elements == _elements_allocated; }

//  If we know we are going to be adding a certain number of items, we can
//  tell the resizable_array_base to resize itself to allow that number of extra items

    int  make_room_for_extra_items (int);

    int  resize                 (int);
    int  add                    (T);
    void  add_no_check_space     (T extra) { _things[_number_elements] = extra; _number_elements++;}
    void operator +=            (T);
    void operator +=            (const resizable_array_base<T> &);
    int  add_if_space           (T);
    int  add                    (const T *, int);
    int  remove_item            (int);
    int  erase                  (int i) { (void) remove_item (i); return i;}
    int  erase                  (int, int);
    int  remove_first           (const T);
    int  remove_all             (const T);
    int  remove_items           (const int *, int);     // dense array of items to be deleted
    int  remove_items           (const int *);          // non zero items removed
    int  remove_two_items       (const T i1, const T i2);     // two specific things to be excised. Must be only one occurrence of each
    int  index                  (const T) const;
    int  rindex                 (const T) const;
    int  contains               (const T) const;
#ifdef NDEBUG
    T &  operator []            (int ndx) const { return _things[ndx];}
    T    item                   (int ndx) const { return _things[ndx];}
#else
    T    item                   (int) const;
    T &  operator []            (int) const;
#endif

#ifdef NDEBUG
    T    first                  () const { return _things[0];}
    T    last_item              () const { return _things[_number_elements - 1];}
#else
    T    first                  () const;
    T    last_item              () const;
#endif
    void sort                   ( int  (*) (const T *, const T *));
    template <typename C>
    void iwqsort                (C &);
    template <typename C>
    void iwqsort_lambda         (C);
    T    pop                    ();
    void seti                   (int, const T);
    int  insert_at_beginning    (const T item_to_insert);
    int  insert_before          (int, T);
    int  insert_after           (int, T);
    int  insert_in_order        (T, int (*) (const T *, const T *));
    int  swap_elements          (int, int);
    T    next_after_wrap        (int &, int = 1) const;   // returns next item
    int  next_index_after_wrap  (int, int = 1) const;     // returns next index
    void each                   (void (*) (T));
    T * rawdata                 () const { return _things;}    // dangerous, but good for performance

    typedef int next_iterator;

    int next_begin              () const { return IWARAY_NEXT_ITERATOR_START;}

          T * begin()  const { return _things;}
    const T * cbegin() const { return _things;}
          T * end ()   const { return _things + _number_elements;}
    const T * cend ()  const { return _things + _number_elements;}
};

template <typename T>
class resizable_array_p : public resizable_array_base<T *>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array_base<T *>::_number_elements;
    using resizable_array_base<T *>::_elements_allocated;
    using resizable_array_base<T *>::_things;
#endif

  public:
    resizable_array_p            ();
    resizable_array_p            (int);
    resizable_array_p            (T *);
    ~resizable_array_p           ();
 
    int  resize                   (int);
    int  resize_no_delete         (int);
    int  resize_keep_storage      (int);
    int  resize_keep_storage_no_delete (int);

    void seti                     (int, T *);

    resizable_array_p<T> & operator = (resizable_array_p<T> &&);

    int  remove_item              (int);
    T *  remove_no_delete         (int);
    int  remove_no_delete         (T *);
    int  remove                   (T *);
    int  remove_items             (const int *, int);     // dense array of items to be deleted
    int  remove_items             (const int *);          // non zero items removed
    template <typename O>
    int  remove_items_fn          (O op);          // remove all items for which op evaluates as true

    void chop                     (int = 1);    // remove from end

    int  transfer_in              (resizable_array_p<T> &);
    int  transfer_in              (resizable_array_p<T> &, int);

    template <typename F> void each (F &) const;
    template <typename F> void each (F &);
    template <typename F> void each (const F &) const;
    template <typename F> void each (const F &);

    template <typename F> void each_lambda (F) const;
};

/*
  No provision for objects with constructors.
*/

template <typename T>
class resizable_array: public resizable_array_base <T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array_base<T>::_number_elements;
    using resizable_array_base<T>::_elements_allocated;
    using resizable_array_base<T>::_things;
#endif

  public:
    resizable_array          ();
    resizable_array          (int);
    resizable_array          (int, const T);
    resizable_array          (const resizable_array<T> &);
    ~resizable_array         ();

    int  resize_keep_storage    (int);
    int  add_if_not_already_present (const T);
    void chop                (int = 1);     // remove from the end
    int  set_all             (const T);
    int  insert_in_order      (const T, int = 1);
    int  insert_in_order_if_not_already_present (const T, int = 1);
    int  operator ==          (const resizable_array<T> &) const;
    int  operator ==          (int) const;
    resizable_array<T> & operator = (const resizable_array<T> &);
    resizable_array<T> & operator = (resizable_array<T> &&);
    int  copy                 (const resizable_array<T> &, int = -1);
    void reverse             ();
    int extend               (int, T = T (0));
    T   max_val              () const;
    T   min_val              () const;
    template <typename O>
    int  remove_items_fn     (O op);          // remove all items for which op evaluates as true

    int add_non_duplicated_elements (const resizable_array<T> &);

    template <typename F> void each (F &);
    template <typename F> void each (const F &);
    template <typename F> void each (F &) const;
    template <typename F> void each (const F &) const;

    template <typename F> void each_lambda (F) const;

//  mimic standard library

    unsigned int count (const T &) const;
    void fill (const T & v) {set_all(v);}
};

/*
  Whenever the [] operator of an extending resizable array is accessed,
  the array is extended if necessary
*/

template <typename T>
class extending_resizable_array : public resizable_array<T>
{
  protected:
#ifdef IW_TWO_PHASE_TEMPLATES
    using resizable_array_base<T>::_number_elements;
    using resizable_array_base<T>::_elements_allocated;
    using resizable_array_base<T>::_things;
#endif

    T _initialiser;

  public:
    extending_resizable_array (T = T (0));

// Very dangerous function. What if the array has already been
// extended using the previous initialiser!!!

    void set_initialiser (T s) { _initialiser = s;}

    int extend (int);

    T & operator[] (int);
    const T & operator [] (int i) const { return _things[i];}
};

/*
  This derived type supports next () and previous () operators, and
  has a notion of the current item.
*/

template <typename T>
class resizable_array_iterator
{
  private:
    int _index;
    const resizable_array_base<T> * _array;

  public:
    resizable_array_iterator (const resizable_array_base<T> *);
    ~resizable_array_iterator ();

    void reset_iterator () { _index = 0;}
    void set_iterator (int);
    int  iterator_index () { return _index;}

//  int next     (T &);
//  int previous (T &);

    int push ();
};

template <typename T>
class iwaray
{
  protected:
    int _number_elements;
    T * _things;

  public:
    iwaray ();
    iwaray (int);
    ~iwaray ();

    int resize (int);

    int number_elements () const { return _number_elements;}

    T & operator [] (int);
    const T & operator [] (int) const;

    iwaray<T> & operator = (T);
};

#define OK_INDEX(i) ((i) >= 0 && (i) < _number_elements)

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(RESIZABLE_ARRAY_IMPLEMENTATION)

#include <stdlib.h>
#include <string.h>

using std::endl;
using std::cerr;

template <typename T>
resizable_array_base<T>::resizable_array_base ()
{
  _elements_allocated = _number_elements = 0;
  _magic = IWARAY_MAGIC_NUMBER;
  _things = NULL;
}

template <typename T>
resizable_array_p<T>::resizable_array_p ()
{
  return;
}

template <typename T>
resizable_array<T>::resizable_array ()
{
  return;
}

#include <assert.h>

template <typename T>
resizable_array_p<T>::resizable_array_p (int n)
{
  assert (n >= 0);

  this->resize(n);

  return;
}

template <typename T>
resizable_array_p<T>::resizable_array_p (T * item)
{
  this->add(item);

  return;
}

template <typename T>
resizable_array<T>::resizable_array (int n)
{
  assert (n >= 0);

  this->resize(n);

  return;
}

template <typename T>
resizable_array<T>::resizable_array (int n, const T initialiser)
{
  assert (n >= 0);
    
  this->resize(n);

  for (int i = 0; i < _elements_allocated; i++)
  {
    _things[i] = initialiser;
  }
                                
  _number_elements = _elements_allocated;

  return;
}

template <typename T>
resizable_array<T>::resizable_array (const resizable_array<T> & rhs)
{
  _number_elements = rhs._number_elements;
  _elements_allocated = rhs._number_elements;
  if (_number_elements > 0)
  {
    _things = new T[rhs._number_elements];
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i] = rhs._things[i];
    }
  }
  else
    _things = NULL;

  return;
}

template <typename T>
resizable_array_base<T>::~resizable_array_base ()
{
  assert (ok());
   
  _number_elements = _elements_allocated = -1;   //paranoia !!
  _magic = -1;

  if (_things)
  {
    delete [] _things;
    _things = NULL;
  }

  return;
}

template <typename T>
resizable_array_p<T>::~resizable_array_p ()
{
  assert (this->ok());

  this->resize(0);
}

template <typename T>
resizable_array<T>::~resizable_array ()    // should not exist at all
{
}
 
template <typename T>
int
resizable_array_base<T>::index (const T x) const
{
  assert (IWARAY_MAGIC_NUMBER == _magic);
   
  for (int i = 0; i < _number_elements; i++)
  {
    if (x == _things[i])
      return i;
  }

  return -1;
} 

template <typename T>
int
resizable_array_base<T>::rindex (const T x) const
{
  assert (IWARAY_MAGIC_NUMBER == _magic);
   
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    if (x == _things[i])
      return i;
  }

  return -1;
}

template <typename T>
int
resizable_array_base<T>::contains (const T x) const
{
  assert (IWARAY_MAGIC_NUMBER == _magic);
   
  for (int i = 0; i < _number_elements; i++)
  {
    if (x == _things[i])
      return i + 1;
  }

  return 0;
}

#if ! defined(NDEBUG)
template <typename T>
T &
resizable_array_base<T>::operator [] (int ii) const
{
//assert (IWARAY_MAGIC_NUMBER == _magic);

  assert (ii >= 0 && ii < _number_elements);

  return _things[ii];
}
#endif

/*
  The difference between add and add_if_space is that add ()
  always works, by growing the list if needed. 
*/

#define INITIAL_NON_EMPTY_SIZE 3

template <typename T>
int
resizable_array_base<T>::add (T extra)
{
  assert (IWARAY_MAGIC_NUMBER == _magic);

  if (_number_elements == _elements_allocated)
  {
    if (0 == _number_elements)
      resize(INITIAL_NON_EMPTY_SIZE);
    else
      resize(_number_elements + _number_elements / 2 + 1);
  }

  _things[_number_elements] = extra;
  _number_elements++;

  return _number_elements;
}

template <typename T>
int
resizable_array_base<T>::add_if_space (T extra)
{
  assert (IWARAY_MAGIC_NUMBER == _magic);

  if (_number_elements == _elements_allocated)
    return -1;

  _things[_number_elements] = extra;
  _number_elements++;

  return _number_elements;
}

template <typename T>
int
resizable_array_base<T>::add (const T * extra, int nextra)
{
  assert (IWARAY_MAGIC_NUMBER == _magic);
  assert (NULL != extra);
  if (0 == nextra)
    return 1;

  if (_number_elements + nextra > _elements_allocated)
    resize(_number_elements + nextra);

  memcpy(_things + _number_elements, extra, nextra * sizeof(T));

//for (int i = 0; i < nextra; i++)
//{
//  _things[_number_elements + i] = extra[i];
//}

  _number_elements += nextra;

  return 1;
}

/*
  audit function for one of these objects
*/

template <typename T>
int
resizable_array_base<T>::ok() const
{
  if (IWARAY_MAGIC_NUMBER != _magic)
    return 0;

  if (_elements_allocated < _number_elements)
    return 0;

  if (_number_elements > 0)
    return NULL != _things;

  if (_number_elements < 0)
    return 0;

// the cases of 0 == _number_elements

  if (0 == _elements_allocated && NULL == _things)
    return 1;

  if (_elements_allocated && _things)
    return 1;

  return 0;
}

/*
  For efficiency, remove_item does NOT do a resize.
*/

template <typename T>
int
resizable_array_base<T>::remove_item (int item_to_remove)
{
  assert (ok());
  assert (item_to_remove >= 0 && item_to_remove < _number_elements);

  for (int i = item_to_remove; i < _number_elements - 1; i++)
  {
    _things[i] = _things[i + 1];
  }

  _number_elements--;

  if (0 == _number_elements)
  {
    if (NULL != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _things = NULL;
  }

  return _number_elements;
}

template <typename T>
int
resizable_array_base<T>::erase (int istart, int istop)
{
  assert (istart >= 0 && istart <= istop && istop < _number_elements);

  int j = istart;
  for (int i = istop + 1; i < _number_elements; i++)
  {
    _things[j++] = _things[i];
  }

  _number_elements = j;

  return istart;
}

template <typename T>
T *
resizable_array_p<T>::remove_no_delete (int item_to_remove)
{
  T * rc = _things[item_to_remove];

  resizable_array_base<T *>::remove_item(item_to_remove);

  return rc;
}

template <typename T>
int
resizable_array_p<T>::remove_item (int item_to_remove)
{
  assert (this->ok());
  assert (item_to_remove >= 0 && item_to_remove < _number_elements);

  delete _things[item_to_remove];

  return resizable_array_base<T *>::remove_item(item_to_remove);
}

template <typename T>
int
resizable_array_base<T>::remove_first (const T thing_to_remove)
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_item(i);
      return 1;
    }
  }

  return 0;
}

template <typename T>
int
resizable_array_p<T>::remove (T * thing_to_remove)
{
  assert (NULL != thing_to_remove);
  assert (this->ok());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_item(i);
      return 1;
    }
  }

  return 0;
}

template <typename T>
int
resizable_array_p<T>::remove_no_delete (T * thing_to_remove)
{
  assert (NULL != thing_to_remove);
  assert (this->ok());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_no_delete(i);
      return 1;
    }
  }

  return 0;
}


template <typename T>
int
resizable_array_base<T>::remove_all (const T item_to_remove)
{
  assert (ok());

  int j = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (item_to_remove != _things[i])
      _things[j++] = _things[i];
  }

  int items_removed = _number_elements - j;
  _number_elements = j;

  if (0 == _number_elements)
  {
    if (NULL != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _things = NULL;
  }

  return items_removed;
}

/*
  Chop is a little unusual because it does not adjust the storage
  to match. This is for efficiency.
*/

template <typename T>
void
resizable_array<T>::chop (int nchop)
{
  assert (this->ok());
  assert (nchop > 0);

  assert (nchop <= _number_elements);

  _number_elements -= nchop;

  return;
}

template <typename T>
void
resizable_array_p<T>::chop (int nchop)
{
  assert (this->ok());
  assert (nchop > 0);

  assert (nchop <= _number_elements);

  for (int i = _number_elements - nchop; i < _number_elements; i++)
  {
    delete _things[i];
  }

  _number_elements -= nchop;

  return;
}

template <typename T>
int
resizable_array_base<T>::resize (int new_size)
{
  assert (ok());
  assert (new_size >= 0);

  if (0 == new_size)
  {
    if (NULL != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _number_elements = 0;
    _things = NULL;
    return 1;
  }

  if (new_size == _elements_allocated)
    return 1;
  
  if (0 == _elements_allocated)
  {
    assert (NULL == _things);
    _things = new T[new_size];
    if (NULL == _things)
    {
      cerr << "resizable_array_base<T>::resize: malloc failure, size " << new_size << endl;
      return 1;
    }

    _elements_allocated = new_size;
    return 1;
  }

// Create a new array, and copy existing data.

  T * new_things = new T[new_size];

  if (NULL == new_things)
  {
    cerr << "resizable_array_base<T>::resize: malloc failure, size " << new_size << endl;
    return 1;
  }

  int ncopy = _number_elements;
  if (ncopy > new_size)
    ncopy = new_size;

  for (int i = 0; i < ncopy; i++)
  {
    new_things[i] = _things[i];
  }

  delete [] _things;
  _things = new_things;

  _elements_allocated = new_size;

  if (_number_elements > new_size)
    _number_elements = new_size;

  return 1;
}

template <typename T>
int
resizable_array_p<T>::resize_no_delete (int new_size)
{
  return resizable_array_base<T *>::resize(new_size);
}

template <typename T>
int
resizable_array_p<T>::resize_keep_storage (int new_size)
{
  if (new_size < _number_elements)    // expected to be the most common case
  {
    for (int i = new_size; i < _number_elements; i++)
    {
      delete _things[i];
    }

    _number_elements = new_size;
  }
  else if (new_size > _elements_allocated)
  {
    resizable_array_base<T *>::resize(new_size);
  }

  return 1;
}

template <typename T>
int
resizable_array_p<T>::resize_keep_storage_no_delete (int new_size)
{
  assert (new_size >= 0);

  if (new_size <= _number_elements)
  {
    _number_elements = new_size;
    return 1;
  }

  cerr << "resize_keep_storage_no_delete:new size larger!, currently " << _number_elements << " request " << new_size << ", ignored\n";

  return 1;
}

template <typename T>
int
resizable_array_p<T>::resize (int new_size)
{
  assert (this->ok());
  assert (new_size >= 0);

//  First delete any elements which will be lost. Depending on the value
//  of  new_size, this loop may not execute at all.

  for (int i = new_size; i < _number_elements; i++)
  {
    delete _things[i];
  }

  return resize_no_delete(new_size);
}

/*
  When repeatedly re-using these objects, we often want to get rid of all
  the current contents, without giving up our storage. This can do that.
*/

template <typename T>
int
resizable_array<T>::resize_keep_storage (int new_size)
{
  assert (this->ok());
  assert (new_size >= 0);

  if (new_size <= _number_elements)  // expected to be the most common case.
  {
    _number_elements = new_size;
    return 1;
  }

  return this->resize(new_size);
}

template <typename T>
int
resizable_array_base<T>::debug_print (std::ostream & os) const
{
  os << "Resizable array has " << _elements_allocated << " allocated, " <<
        _number_elements << " used\n";
  
  if (! ok())
    cerr << "Warning, OK fails, this = " << (void *) this << " magic = " << _magic << endl;

  return ok();
}

#if ! defined(NDEBUG)
template <typename T>
T
resizable_array_base<T>::item (int i) const
{
  assert (ok());

  assert (i >= 0 && i < _number_elements);

  return _things[i];
}
#endif

template <typename T>
void
resizable_array_base<T>::seti (int i, const T item)
{
  assert (ok());

  assert (i >= 0 && i < _number_elements);

  _things[i] = item;

  return;
}

template <typename T>
void
resizable_array_p<T>::seti (int i, T * item)
{
  assert (this->ok());

  assert (i >= 0 && i < _number_elements);

  delete _things[i];

  _things[i] = item;

  return;
}

template <typename T>
int
resizable_array<T>::set_all (const T fill_value)
{
  assert (this->ok());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = fill_value;
  }

  return _number_elements;
}

#define CAST_FOR_QSORT ( int (*)(const void *, const void *) )

template <typename T>
void
resizable_array_base<T>::sort (int (* comparitor)(const T *, const T *))
{
  assert (ok());
  assert (NULL != comparitor);

  if (_number_elements <= 1)
    return;

// This generates an error, but dammit, I just cannot figure out how to
// get this right

  qsort((void *) _things, _number_elements, sizeof(T), CAST_FOR_QSORT comparitor);

  return;
}

template <typename T>
int
resizable_array_base<T>::insert_before (int where_to_insert,
                                        T item_to_insert)
{
  assert (ok());

  assert (where_to_insert >= 0);
  if (_number_elements > 0)
  {
    assert (where_to_insert < _number_elements);
  }

  if (_number_elements == _elements_allocated)
  {
    if (0 == _number_elements)
      resize(INITIAL_NON_EMPTY_SIZE);
    else
      resize(_number_elements + _number_elements / 2 + 1);
  }

  for (int i = _number_elements; i > where_to_insert; i--)
  {
    _things[i] = _things[i - 1];
  }

  _things[where_to_insert] = item_to_insert;

  _number_elements++;

  return _number_elements;
}

template <typename T>
int
resizable_array_base<T>::insert_after (int where_to_insert,
                                       T item_to_insert)
{     
  assert (ok());
  if (_number_elements > 0)
  {
    assert (where_to_insert < _number_elements);
  }

  if (_number_elements == _elements_allocated)
  {
    if (0 == _number_elements)
      resize (INITIAL_NON_EMPTY_SIZE);
    else
      resize (_number_elements + _number_elements / 2 + 1);
  }

  for (int i = _number_elements; i > where_to_insert + 1; i--)
  {
    _things[i] = _things[i - 1];
  }

  _things[where_to_insert + 1] = item_to_insert;

  _number_elements++;

  return _number_elements;
}

template <typename T>
int
resizable_array_base<T>::insert_at_beginning (const T item_to_insert)
{
  assert (ok());

  if (0 == _number_elements )
    return this->add(item_to_insert);
  else
    return this->insert_before(0, item_to_insert);
}

template <typename T>
int
resizable_array_base<T>::insert_in_order (T item_to_insert, 
                       int (* comparitor)(const T *, const T *))
{
  assert (ok());

  if (0 == _number_elements)
    return add(item_to_insert);

  int i;
  for (i = 0; i < _number_elements; i++)
  {
    if (comparitor(&item_to_insert, &_things[i]) <= 0)
      break;
  }

  if (i == _number_elements)
    return this->add(item_to_insert);

  return this->insert_before(i, item_to_insert);
}

template <typename T>
int
resizable_array<T>::insert_in_order (const T item_to_insert, int increasing_order)
{
  if (0 == _number_elements)
    return this->add(item_to_insert);

  if (increasing_order)
  {
    if (item_to_insert < _things[0])
      return this->insert_at_beginning(item_to_insert);

    if (item_to_insert >= _things[_number_elements - 1])
      return this->add(item_to_insert);

    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert < _things[i])
        return this->insert_before(i, item_to_insert);
    }

    return this->add(item_to_insert);
  }
  else     // decreasing order, largest first.
  {
    if (item_to_insert > _things[0])
      return this->insert_at_beginning(item_to_insert);

    if (item_to_insert <= _things[_number_elements - 1])
      return this->add(item_to_insert);

    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert > _things[i])
        return this->insert_before(i, item_to_insert);
    }

    return this->add(item_to_insert);
  }
}

template <typename T>
int
resizable_array<T>::insert_in_order_if_not_already_present (const T item_to_insert, int increasing_order)
{
  if (0 == _number_elements)
    return this->add(item_to_insert);

  if (increasing_order)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert < _things[i])
        return this->insert_before(i, item_to_insert);

      if (item_to_insert == _things[i])
        return 0;
    }

    return this->add(item_to_insert);
  }
  else     // decreasing order, largest first.
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert > _things[i])
        return this->insert_before(i, item_to_insert);

      if (item_to_insert == _things[i])
        return 0;
    }

    return this->add(item_to_insert);
  }

  return 0;     // no path to here
}

template <typename T>
T
resizable_array_base<T>::pop ()
{
  assert (ok());
  assert (_number_elements > 0);

  return _things[--_number_elements];
}

#if ! defined(NDEBUG)
template <typename T>
T
resizable_array_base<T>::last_item () const
{
  assert (ok());
  assert (_number_elements > 0);

  return _things[_number_elements - 1];
}

template <typename T>
T
resizable_array_base<T>::first () const
{
  assert (ok());
  assert (_number_elements > 0);

  return _things[0];
}
#endif

template <typename T>
void
resizable_array_base<T>::operator += (const resizable_array_base<T> & oo)
{
  assert (ok());

  int noo = oo.number_elements();
  resize (_number_elements + noo);

  for (int i = 0; i < noo; i++)
  {
    this->add(oo.item(i));
  }

  return;
}

template <typename T>
void
resizable_array_base<T>::operator += (T extra)
{
  (void) this->add(extra);

  return;
}

/*
  Transfer points from one resizable_array to another
*/

template <typename T>
int
resizable_array_p<T>::transfer_in (resizable_array_p<T> & oo)
{
  assert (this->ok());

  int noo = oo.number_elements();
  this->resize(_number_elements + noo);

  int j = _number_elements;
  for (int i = 0; i < noo; i++)
  {
    _things[j++] = oo._things[i];
  }

  _number_elements += noo;

  oo.resize_no_delete(0);

  return _number_elements;
}

/*
  Transfer a single item from one array to this one
  Make sure that we remove it from the other one
*/

template <typename T>
int
resizable_array_p<T>::transfer_in (resizable_array_p<T> & oo,
                                   int item_to_transfer)
{
  assert (this->ok());
  assert (oo.ok_index(item_to_transfer));

  T * tmp = oo._things[item_to_transfer];
  this->add(tmp);
  oo.remove_no_delete(item_to_transfer);

  return 1;
}

template <typename T>
T
resizable_array_base<T>::next_after_wrap (int & i, int direction) const
{
  assert (direction);

  i += direction;
  if (i >= _number_elements)
    i -= _number_elements;
  else if (i < 0)
    i = _number_elements + i;

  return _things[i];
}

template <typename T>
int
resizable_array_base<T>::next_index_after_wrap (int i, int direction) const
{
  assert (direction);

  i += direction;
  if (i >= _number_elements)
    i -= _number_elements;
  else if (i < 0)
    i = _number_elements + i;

  return i;
}

/*
  These are equal if all elements are equal
*/

template <typename T>
int
resizable_array<T>::operator == (const resizable_array<T> & other) const
{
  if (other.number_elements() != _number_elements)
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] != other._things[i])
      return 0;
  }

  return 1;
}

template <typename T>
resizable_array<T> & 
resizable_array<T>::operator = (const resizable_array<T> & other)
{   
  int nother = other._number_elements;

  if (_elements_allocated < nother)
    this->resize(nother);

  for (int i = 0; i < nother; i++)
  {
    _things[i] = other._things[i];
  }

  _number_elements = nother;

  return *this;
}

template <typename T>
resizable_array_p<T> & 
resizable_array_p<T>::operator = (resizable_array_p<T> && other)
{   
  resize(0);

  _things = other._things;
  _number_elements = other._number_elements;
  _elements_allocated = other._elements_allocated;

  other._things = NULL;
  other._number_elements = 0;
  other._elements_allocated = 0;

//cerr << "resizable_array_p::operator = move\n";

  return *this;
}

template <typename T>
resizable_array<T> & 
resizable_array<T>::operator = (resizable_array<T> && other)
{   
//cerr << "resizable_array::operator = rvalue\n";
  if (NULL != _things)
    delete [] _things;

//cerr << "operator = && called, my size " << _number_elements << " rhs " << other._number_elements << endl;

  _things = other._things;
  _elements_allocated = other._elements_allocated;
  _number_elements = other._number_elements;

  other._elements_allocated = 0;
  other._number_elements = 0;
  other._things = NULL;

  return *this;
}

template <typename T>
int
resizable_array<T>::copy (const resizable_array<T> & other, int ncopy)
{
  if (-1 == ncopy)
    ncopy = other._number_elements;
  else
  {
    assert (other.ok_index(ncopy));
  }

  _number_elements = 0;
  if (_elements_allocated < ncopy)
    this->resize(ncopy);

  for (int i = 0; i < ncopy; i++)
  {
    _things[i] = other._things[i];
  }

  _number_elements = ncopy;

  return 1;
}


template <typename T>
void
resizable_array<T>::reverse ()
{
  int n = _number_elements / 2;
  for (int i = 0; i < n; i++)
  {
    T tmp = _things[i];
    int j = _number_elements - 1 - i;
    _things[i] = _things[j];
    _things[j] = tmp;
    j--;
  }

  return;
}

template <typename T>
int
resizable_array<T>::extend (int new_size, T fill_value)
{
  if (new_size <= _number_elements)
    return 0;

  if (new_size > _elements_allocated)
    this->resize(new_size);
  
  for (int i = _number_elements; i < new_size; i++)
  {
    _things[i] = fill_value;
  }

  return _number_elements = new_size;
}

template <typename T>
T
resizable_array<T>::max_val () const
{
  if (0 == _number_elements)
    return T(0);

  T qmax = _things[0];
  for (int i = 1; i < _number_elements; i++)
  {
    if (_things[i] > qmax)
      qmax = _things[i];
  }

  return qmax;
}

template <typename T>
T
resizable_array<T>::min_val () const
{
  if (0 == _number_elements)
    return T(0);

  T qmin = _things[0];
  for (int i = 1; i < _number_elements; i++)
  {
    if (_things[i] < qmin)
      qmin = _things[i];
  }

  return qmin;
}

template <typename T>
int
resizable_array_base<T>::swap_elements (int i1, int i2)
{
  assert (ok_index(i1));
  assert (ok_index(i2));
  
  T tmp = _things[i1];
  _things[i1] = _things[i2];
  _things[i2] = tmp;

  return 1;
}

template <typename T>
int
resizable_array<T>::add_non_duplicated_elements (const resizable_array<T> & qq)
{
  int rc = 0;

  int nq = qq._number_elements;
  for (int i = 0; i < nq; i++)
  {
    T qi = qq._things[i];
    int foundqi = 0;
    for (int j = 0; j < _number_elements && 0 == foundqi; j++)    // do we already have it?
    {
      if (_things[j] == qi)
        foundqi++;
    }
    if (0 == foundqi)
    {
      this->add(qi);
      rc++;
    }
  }

  return rc;
}

template <typename T>
int
resizable_array<T>::add_if_not_already_present (const T qq)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (qq == _things[i])
      return 0;
  }

// No match found, add qq

  return this->add(qq);
}

/*template <typename T>
void
resizable_array_p<T>::each (void (T:: *xx) ())
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->xx();
  }

  return;
}*/

template <typename T>
void
resizable_array_base<T>::each (void (* something)(T))
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    something(_things[i]);
  }

  return;
}

template <typename T>
int
operator < (int i, const resizable_array<T> & qq)
{
  return i < qq.number_elements();
}

template <typename T>
int
operator <= (int i, const resizable_array<T> & qq)
{
  return i <= qq.number_elements();
}

template <typename T>
int
operator > (int i, const resizable_array<T> & qq)
{
  return i > qq.number_elements();
}

template <typename T>
int
operator >= (int i, const resizable_array<T> & qq)
{
  return i >= qq.number_elements();
}

template <typename T>
int
operator < (int i, const resizable_array_p<T> & qq)
{
  return i < qq.number_elements();
}

template <typename T>
int
operator <= (int i, const resizable_array_p<T> & qq)
{
  return i <= qq.number_elements();
}

template <typename T>
int
operator > (int i, const resizable_array_p<T> & qq)
{
  return i > qq.number_elements();
}

template <typename T>
int
operator >= (int i, const resizable_array_p<T> & qq)
{
  return i >= qq.number_elements();
}

template <typename T>
int
operator == (int i, const resizable_array_p<T> & qq)
{
  return i == qq.number_elements();
}

template <typename T>
int
resizable_array_base<T>::remove_items (const int * items_to_remove, int nremove)
{
  int jptr = 0;

// From and to pointers within our own array

  int fptr = 0;
  int tptr = 0;

  int items_removed = 0;
  while (fptr < _number_elements)
  {
    if (fptr == items_to_remove[jptr])
    {
      items_removed++;
      jptr++;
      fptr++;
      if (jptr == nremove)
        break;
    }
    else
    {
      _things[tptr] = _things[fptr];
      tptr++;
      fptr++;
    }
  }


  while (fptr < _number_elements)
  {
    _things[tptr] = _things[fptr];
    tptr++;
    fptr++;
  }

  _number_elements = tptr;
  
  if (jptr == nremove)
    return items_removed;

// Looks like a serious problem

  cerr << "resizable_array_base<T>::remove_items: array not exhausted, WAS IT SORTED?\n";
  for (int i = 0; i < nremove; i++)
  {
    cerr << "  " << i << ' ' << items_to_remove[i] << endl;
  }

  abort();

  return 0;
}

template <typename T>
int
resizable_array_base<T>::remove_items (const int * items_to_remove)
{
  int rc = 0;

  int tptr = 0;      // where we put items in re-filling the array

  for (int i = 0; i < _number_elements; ++i)
  {
    if (! items_to_remove[i])
    {
      _things[tptr] = _things[i];
      tptr++;
    }
    else
      rc++;
  }

  _number_elements = tptr;

  return rc;
}

template <typename T>
int
resizable_array_p<T>::remove_items (const int * items_to_remove, int nremove)
{
  for (int i = 0; i < nremove; i++)
  {
    int j = items_to_remove[i];
    if (! this->ok_index(j))
    {
      cerr << "resizable_array_p<T>::remove_items: cannot remove item " << j << endl;
      abort();
      return 0;
    }

    delete _things[j];
  }

  return resizable_array_base<T *>::remove_items (items_to_remove, nremove);
}

template <typename T>
int
resizable_array_p<T>::remove_items (const int * items_to_remove)
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (items_to_remove[i])
      delete _things[i];
  }

  return resizable_array_base<T *>::remove_items(items_to_remove);
}

template <typename T> template <typename O>
int
resizable_array_p<T>::remove_items_fn (O opq)
{
  int tptr = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if (opq(_things[i]))    // being removed
      delete _things[i];
    else
    {
      if (i != tptr)
        _things[tptr] = _things[i];
      tptr++;
    }
  }

  if (tptr == _number_elements)
    return 0;

  const int rc = _number_elements - tptr;

  _number_elements = tptr;

  return rc;
}

template <typename T> template <typename O>
int
resizable_array<T>::remove_items_fn (O opq)
{
  int tptr = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if (opq(_things[i]))    // being removed
      ;
    else
    {
      if (i > tptr)
        _things[tptr] = _things[i];
      tptr++;
    }
  }

  if (tptr == _number_elements)
    return 0;

  const int rc = _number_elements - tptr;

  _number_elements = tptr;

  return rc;
}

/*
  this makes the assumption that there is only one occurrence of each in the array.
*/

template <typename T>
int
resizable_array_base<T>::remove_two_items (const T i1, const T i2)
{
  int ndx1 = -1;
  int ndx2 = -1;
  int nfound = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if (i1 == _things[i])
      ndx1 = i;
    else if (i2 == _things[i])
      ndx2 = i;
    else
      continue;

    nfound++;
    if (2 == nfound)
      break;
  }

  if (2 != nfound)
    return 0;

  if (ndx1 > ndx2)
  {
    int tmp = ndx1;
    ndx1 = ndx2;
    ndx2 = tmp;
  }

//cerr << " found at " << ndx1 << " and " << ndx2 << endl;
  for (int i = ndx1; i < ndx2; ++i)
  {
    _things[i] = _things[i+1];
  }

  _number_elements -= 2;

  for (int i = ndx2; i < _number_elements; ++i)
  {
    _things[i] = _things[i+2];
  }

  return 1;
}

template <typename T>
int
resizable_array<T>::operator== (int rhs) const
{
  return _number_elements == rhs;
}

template <typename T>
int
operator == (int lhs, const resizable_array<T> & rhs)
{
  return lhs == rhs.number_elements();
}

template <typename T>
int
resizable_array_base<T>::make_room_for_extra_items (int e)
{
  assert (ok());
  assert (e >= 0);

  if (_number_elements + e <= _elements_allocated)    // already enough room
    return 1;

  int save_number_elements = _number_elements;

  int rc = resize(_number_elements + e);

  _number_elements = save_number_elements;

  return rc;
}

template <typename T>
unsigned int
resizable_array<T>::count (const T & needle) const
{
  unsigned int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (needle == _things[i])
      rc++;
  }

  return rc;
}

#endif

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(EXTENDING_RESIZABLE_ARRAY_IMPLEMENTATION)

#include <assert.h>

template <typename T>
extending_resizable_array<T>::extending_resizable_array (T ini)
{
  _initialiser = ini;

  return;
}

template <typename T>
T &
extending_resizable_array<T>::operator [] (int i)
{
  assert (i >= 0);

  if (i < _number_elements)
    return _things[i];

  if (i < 10)
    resizable_array<T>::extend(10, _initialiser);
  else if (i < 100)
    resizable_array<T>::extend(100, _initialiser);
  else
    resizable_array<T>::extend(i + i, _initialiser);

  return _things[i];
}

template <typename T>
int
extending_resizable_array<T>::extend (int new_size)
{
  return resizable_array<T>::extend(new_size, _initialiser);
}

#endif

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARAY_IMPLEMENTATION)

template <typename T>
iwaray<T>::iwaray ()
{
  _number_elements = 0;
  _things = NULL;

  return;
}

template <typename T>
iwaray<T>::iwaray (int qsize)
{
  assert (qsize > 0);

  _things = new T [qsize];

  _number_elements = qsize;
}

template <typename T>
iwaray<T>::~iwaray ()
{
  if (NULL != _things)
    delete [] _things;

  _number_elements = -1;

  return;
}

template <typename T>
int
iwaray<T>::resize (int qsize)
{
  if (_number_elements)
    delete [] _things;

  _things = new T[qsize];
  _number_elements = qsize;

  return 1;
}

template <typename T>
const T &
iwaray<T>::operator [] (int i) const
{
  assert (i >= 0 && i < _number_elements);

  return _things[i];
}

template <typename T>
T &
iwaray<T>::operator [] (int i)
{
  assert (i >= 0 && i < _number_elements);

  return _things[i];
}

#endif

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARAYQ_ASSIGN)

template <typename T>
iwaray<T> &
iwaray<T>::operator = (T rhs)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = rhs;
  }

  return *this;
}

#endif

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARAY_EACH_IMPLEMENTATION)

template <typename T> template <typename F>
void
resizable_array_p<T>::each (F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each (const F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(*(_things[i]));
  }

  return;
}
template <typename T> template <typename F>
void
resizable_array_p<T>::each (F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each (const F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (const F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (const F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f(_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each_lambda (F f) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    f(_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each_lambda (F f) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    f(_things[i]);
  }

  return;
}

#endif

#endif
