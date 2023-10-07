#ifndef IWARAY_IMPLEMENTATION_H
#define IWARAY_IMPLEMENTATION_H

//#define IWARAY_USE_IWMALLOC
#ifdef IWARAY_USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iwaray.h"

template <typename T>
resizable_array_base<T>::resizable_array_base ()
{
  _elements_allocated = _number_elements = 0;
  _magic = IWARAY_MAGIC_NUMBER;
  _things = nullptr;
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

  this->resize (n);

  return;
}

template <typename T>
resizable_array_p<T>::resizable_array_p (T * item)
{
  add (item);

  return;
}

template <typename T>
resizable_array<T>::resizable_array (int n)
{
  assert (n >= 0);

  this->resize (n);

  return;
}

template <typename T>
resizable_array<T>::resizable_array (int n, const T initialiser)
{
  assert (n >= 0);
    
  this->resize (n);

  for (int i = 0; i < _elements_allocated; i++)
  {
    _things[i] = initialiser;
  }
                                
  _number_elements = _elements_allocated;

  return;
}

template <typename T>
resizable_array_base<T>::~resizable_array_base ()
{
  assert (ok ());
   
  _number_elements = _elements_allocated = -1;   //paranoia !!
  _magic = -1;

  if (_things)
  {
    delete [] _things;
    _things = nullptr;
  }

  return;
}

template <typename T>
resizable_array_p<T>::~resizable_array_p ()
{
  assert (this->ok ());

  this->resize (0);
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

template <typename T>
T &
resizable_array_base<T>::operator [] (int ii) const
{
//assert (IWARAY_MAGIC_NUMBER == _magic);

  assert (ii >= 0 && ii < _number_elements);

  return _things[ii];
}

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
      resize (INITIAL_NON_EMPTY_SIZE);
    else
      resize (_number_elements + _number_elements);
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
  assert (nullptr != extra && nextra > 0);

  if (_number_elements + nextra > _elements_allocated)
    resize (_number_elements + nextra);

  memcpy (_things + _number_elements, extra, nextra * sizeof (T));

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
resizable_array_base<T>::ok () const
{
  if (IWARAY_MAGIC_NUMBER != _magic)
    return 0;

  if (_elements_allocated < _number_elements)
    return 0;

  if (_number_elements > 0)
    return nullptr != _things;

  if (_number_elements < 0)
    return 0;

// the cases of 0 == _number_elements

  if (0 == _elements_allocated && nullptr == _things)
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
  assert (ok ());
  assert (item_to_remove >= 0 && item_to_remove < _number_elements);

  for (int i = item_to_remove; i < _number_elements - 1; i++)
  {
    _things[i] = _things[i + 1];
  }

  _number_elements--;

  if (0 == _number_elements)
  {
    if (nullptr != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _things = nullptr;
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
int
resizable_array_p<T>::remove_no_delete (int item_to_remove)
{
  return resizable_array_base<T *>::remove_item (item_to_remove);
}

template <typename T>
int
resizable_array_p<T>::remove_item (int item_to_remove)
{
  assert (this->ok ());
  assert (item_to_remove >= 0 && item_to_remove < _number_elements);

  delete _things[item_to_remove];

  return resizable_array_base<T *>::remove_item (item_to_remove);
}

template <typename T>
int
resizable_array_base<T>::remove_first (const T thing_to_remove)
{
  assert (ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_item (i);
      return 1;
    }
  }

  return 0;
}

template <typename T>
int
resizable_array_p<T>::remove (T * thing_to_remove)
{
  assert (nullptr != thing_to_remove);
  assert (this->ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_item (i);
      return 1;
    }
  }

  return 0;
}

template <typename T>
int
resizable_array_p<T>::remove_no_delete (T * thing_to_remove)
{
  assert (nullptr != thing_to_remove);
  assert (this->ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    if (thing_to_remove == _things[i])
    {
      (void) remove_no_delete (i);
      return 1;
    }
  }

  return 0;
}


template <typename T>
int
resizable_array_base<T>::remove_all (const T item_to_remove)
{
  assert (ok ());

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
    if (nullptr != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _things = nullptr;
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
  assert (this->ok ());
  assert (nchop > 0);

  assert (nchop <= _number_elements);

  _number_elements -= nchop;

  return;
}

template <typename T>
void
resizable_array_p<T>::chop (int nchop)
{
  assert (this->ok ());
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
  assert (ok ());
  assert (new_size >= 0);

  if (0 == new_size)
  {
    if (nullptr != _things)      // free anything already there
      delete [] _things;
    _elements_allocated = 0;
    _number_elements = 0;
    _things = nullptr;
    return 1;
  }

  if (new_size == _elements_allocated)
    return 1;
  
  if (0 == _elements_allocated)
  {
    assert (nullptr == _things);
    _things = new T[new_size];
    if (nullptr == _things)
    {
      cerr << "resizable_array_base<T>::resize: malloc failure, size " << new_size << endl;
      return 1;
    }

    _elements_allocated = new_size;
    return 1;
  }

// Create a new array, and copy existing data.

  T * new_things = new T[new_size];

  if (nullptr == new_things)
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
  return resizable_array_base<T *>::resize (new_size);
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
    resizable_array_base<T *>::resize (new_size);
  }

  return 1;
}

template <typename T>
int
resizable_array_p<T>::resize (int new_size)
{
  assert (this->ok ());
  assert (new_size >= 0);

//  First delete any elements which will be lost. Depending on the value
//  of  new_size, this loop may not execute at all.

  for (int i = new_size; i < _number_elements; i++)
  {
    delete _things[i];
  }

  return resize_no_delete (new_size);
}

/*
  When repeatedly re-using these objects, we often want to get rid of all
  the current contents, without giving up our storage. This can do that.
*/

template <typename T>
int
resizable_array<T>::resize_keep_storage (int new_size)
{
  assert (this->ok ());
  assert (new_size >= 0);

  if (new_size <= _number_elements)  // expected to be the most common case.
  {
    _number_elements = new_size;
    return 1;
  }

  return this->resize (new_size);
}

#include <iostream>

template <typename T>
int
resizable_array_base<T>::debug_print (ostream & os) const
{
  os << "Resizable array has " << _elements_allocated << " allocated, " <<
        _number_elements << " used\n";
  
  if (! ok ())
    cerr << "Warning, OK fails, this = " << (void *) this << " magic = " << _magic << endl;

  return ok ();
}

template <typename T>
T
resizable_array_base<T>::item (int i) const
{
  assert (ok ());

  assert (i >= 0 && i < _number_elements);

  return _things[i];
}

template <typename T>
void
resizable_array_base<T>::seti (int i, const T item)
{
  assert (ok ());

  assert (i >= 0 && i < _number_elements);

  _things[i] = item;

  return;
}

template <typename T>
void
resizable_array_p<T>::seti (int i, T * item)
{
  assert (this->ok ());

  assert (i >= 0 && i < _number_elements);

  delete _things[i];

  _things[i] = item;

  return;
}

template <typename T>
int
resizable_array<T>::set_all (const T fill_value)
{
  assert (this->ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = fill_value;
  }

  return _number_elements;
}

#define CAST_FOR_QSORT ( int (*) (const void *, const void *) )

template <typename T>
void
resizable_array_base<T>::sort (int (* comparitor) (const T *, const T *))
{
  assert (ok ());
  assert (nullptr != comparitor);

  if (_number_elements <= 1)
    return;

// This generates an error, but dammit, I just cannot figure out how to
// get this right

  qsort ((void *) _things, _number_elements, sizeof (T), CAST_FOR_QSORT comparitor);

  return;
}

template <typename T>
int
resizable_array_base<T>::insert_before (int where_to_insert,
                                        T item_to_insert)
{
  assert (ok ());

  assert (where_to_insert >= 0);
  if (_number_elements > 0)
  {
    assert (where_to_insert < _number_elements);
  }

  if (_number_elements == _elements_allocated)
  {
    if (0 == _number_elements)
      resize (INITIAL_NON_EMPTY_SIZE);
    else
      resize (_number_elements + _number_elements);
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
  assert (ok ());
  if (_number_elements > 0)
  {
    assert (where_to_insert < _number_elements);
  }

  if (_number_elements == _elements_allocated)
  {
    if (0 == _number_elements)
      resize (INITIAL_NON_EMPTY_SIZE);
    else
      resize (_number_elements + _number_elements);
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
  assert (ok ());

  if (0 == _number_elements )
    return add (item_to_insert);
  else
    return insert_before (0, item_to_insert);
}

template <typename T>
int
resizable_array_base<T>::insert_in_order (T item_to_insert, 
                       int (* comparitor) (const T *, const T *))
{
  assert (ok ());

  if (0 == _number_elements)
    return add (item_to_insert);

  int i;
  for (i = 0; i < _number_elements; i++)
  {
    if (comparitor (&item_to_insert, &_things[i]) <= 0)
      break;
  }

  if (i == _number_elements)
    return add (item_to_insert);

  return insert_before (i, item_to_insert);
}

template <typename T>
int
resizable_array<T>::insert_in_order (const T item_to_insert, int increasing_order)
{
  if (0 == _number_elements)
    return add (item_to_insert);

  if (increasing_order)
  {
    if (item_to_insert < _things[0])
      return insert_at_beginning (item_to_insert);

    if (item_to_insert >= _things[_number_elements - 1])
      return add (item_to_insert);

    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert < _things[i])
        return insert_before (i, item_to_insert);
    }

    return add (item_to_insert);
  }
  else     // decreasing order, largest first.
  {
    if (item_to_insert > _things[0])
      return insert_at_beginning (item_to_insert);

    if (item_to_insert <= _things[_number_elements - 1])
      return add (item_to_insert);

    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert > _things[i])
        return insert_before (i, item_to_insert);
    }

    return add (item_to_insert);
  }
}

template <typename T>
int
resizable_array<T>::insert_in_order_if_not_already_present (const T item_to_insert, int increasing_order)
{
  if (0 == _number_elements)
    return add (item_to_insert);

  if (increasing_order)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert < _things[i])
        return insert_before (i, item_to_insert);

      if (item_to_insert == _things[i])
        return 0;
    }

    return add (item_to_insert);
  }
  else     // decreasing order, largest first.
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (item_to_insert > _things[i])
        return insert_before (i, item_to_insert);

      if (item_to_insert == _things[i])
        return 0;
    }

    return add (item_to_insert);
  }

  return 0;     // no path to here
}

template <typename T>
T
resizable_array_base<T>::pop ()
{
  assert (ok ());
  assert (_number_elements > 0);

  return _things[--_number_elements];
}

template <typename T>
T
resizable_array_base<T>::last_item () const
{
  assert (ok ());
  assert (_number_elements > 0);

  return _things[_number_elements - 1];
}

template <typename T>
void
resizable_array_base<T>::operator += (const resizable_array_base<T> & oo)
{
  assert (ok ());

  int noo = oo.number_elements ();
  resize (_number_elements + noo);

  for (int i = 0; i < noo; i++)
  {
    add (oo.item (i));
  }

  return;
}

template <typename T>
void
resizable_array_base<T>::operator += (T extra)
{
  (void) add (extra);

  return;
}

/*
  Transfer points from one resizable_array to another
*/

template <typename T>
int
resizable_array_p<T>::transfer_in (resizable_array_p<T> & oo)
{
  assert (this->ok ());

  int noo = oo.number_elements ();
  this->resize (_number_elements + noo);

  int j = _number_elements;
  for (int i = 0; i < noo; i++)
  {
    _things[j++] = oo._things[i];
  }

  _number_elements += noo;

  oo.resize_no_delete (0);

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
  assert (this->ok ());
  assert (oo.ok_index (item_to_transfer));

  T * tmp = oo._things[item_to_transfer];
  add (tmp);
  oo.remove_no_delete (item_to_transfer);

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

/*
  These are equal if all elements are equal
*/

template <typename T>
int
resizable_array<T>::operator == (const resizable_array<T> & other) const
{
  if (other.number_elements () != _number_elements)
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
    this->resize (nother);

  for (int i = 0; i < nother; i++)
  {
    _things[i] = other._things[i];
  }

  _number_elements = nother;

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
    assert (other.ok_index (ncopy));
  }

  _number_elements = 0;
  if (_elements_allocated < ncopy)
    this->resize (ncopy);

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
    this->resize (new_size);
  
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
    return T (0);

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
    return T (0);

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
  assert (ok_index (i1));
  assert (ok_index (i2));
  
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
      add (qi);
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

  return add (qq);
}

/*template <typename T>
void
resizable_array_p<T>::each (void (T:: *xx) ())
{
  assert (ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->xx ();
  }

  return;
}*/

template <typename T>
void
resizable_array_base<T>::each (void (* something) (T))
{
  assert (ok ());

  for (int i = 0; i < _number_elements; i++)
  {
    something (_things[i]);
  }

  return;
}

template <typename T>
int
operator < (int i, const resizable_array<T> & qq)
{
  return i < qq.number_elements ();
}

template <typename T>
int
operator <= (int i, const resizable_array<T> & qq)
{
  return i <= qq.number_elements ();
}

template <typename T>
int
operator > (int i, const resizable_array<T> & qq)
{
  return i > qq.number_elements ();
}

template <typename T>
int
operator >= (int i, const resizable_array<T> & qq)
{
  return i >= qq.number_elements ();
}

template <typename T>
int
operator < (int i, const resizable_array_p<T> & qq)
{
  return i < qq.number_elements ();
}

template <typename T>
int
operator <= (int i, const resizable_array_p<T> & qq)
{
  return i <= qq.number_elements ();
}

template <typename T>
int
operator > (int i, const resizable_array_p<T> & qq)
{
  return i > qq.number_elements ();
}

template <typename T>
int
operator >= (int i, const resizable_array_p<T> & qq)
{
  return i >= qq.number_elements ();
}

template <typename T>
int
operator == (int i, const resizable_array_p<T> & qq)
{
  return i == qq.number_elements ();
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

  abort ();

  return 0;
}

template <typename T>
int
resizable_array_p<T>::remove_items (const int * items_to_remove, int nremove)
{
  for (int i = 0; i < nremove; i++)
  {
    int j = items_to_remove[i];
    if (! this->ok_index (j))
    {
      cerr << "resizable_array_p<T>::remove_items: cannot remove item " << j << endl;
      abort ();
      return 0;
    }

    delete _things[j];
  }

  return resizable_array_base<T *>::remove_items (items_to_remove, nremove);
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
  return lhs == rhs.number_elements ();
}

template <typename T>
int
resizable_array_base<T>::make_room_for_extra_items (int e)
{
  assert (ok ());
  assert (e >= 0);

  if (_number_elements + e <= _elements_allocated)    // already enough room
    return 1;

  int save_number_elements = _number_elements;

  int rc = resize (_number_elements + e);

  _number_elements = save_number_elements;

  return rc;
}

#endif
