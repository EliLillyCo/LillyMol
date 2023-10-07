#ifndef IWARAY_NP_IMPLEMENTATION
#define IWARAY_NP_IMPLEMENTATION

#include "assert.h"
#include <iomanip>
#include "iwaray.h"

template <typename T>
resizable_array_iterator<T>::resizable_array_iterator (const resizable_array_base<T> * aa)
{
  assert (aa->ok ());

  _array = aa;

  _index = 0;

  return;
}

template <typename T>
resizable_array_iterator<T>::~resizable_array_iterator ()
{
  _index = -1;
  _array = nullptr;
}

template <typename T>
void
resizable_array_iterator<T>::set_iterator (int new_index)
{
  int ne = _array->number_elements ();    // will invoke OK function

  assert (new_index >= 0 && new_index <= ne);

  _index = new_index;

  return;
}

template <typename T>
int
resizable_array_iterator<T>::push ()
{
  (void) _array->number_elements ();

  if (_index > 0)
  {
    _index--;
    return 1;
  }

  return 0;
}

template <typename T>
int
resizable_array_iterator<T>::next (T & result)
{
  int ne = _array->number_elements ();

  if (_index >= 0 && _index < ne)
  {
    result = _array->_things[_index];
    _index++;
    return 1;
  }

  return 0;
}

template <typename T>
int
resizable_array_iterator<T>::previous (T & result)
{
  int ne = _array->number_elements ();

  if (_index > 0 && _index <= ne)
  {
    _index--;
    result = _array->_things[_index];
    return 1;
  }
  
  return 0;
}

#ifdef NOT_IMPLEMENTED_JJ
template <typename T>
resizable_array_p_np<T>::resizable_array_p_np () : resizable_array_iterator<T *> (this)
{
}

template <typename T>
resizable_array_p_np<T>::resizable_array_p_np (int initial_size) :
            resizable_array_p<T> (initial_size) , resizable_array_iterator<T *> (this)
{
}

template <typename T>
resizable_array_np<T>::resizable_array_np () : resizable_array_iterator<T> (this)
{
}

template <typename T>
resizable_array_np<T>::resizable_array_np (int initial_size) : resizable_array<T> (initial_size),
            resizable_array_iterator<T> (this)
{
}
#endif

#endif
