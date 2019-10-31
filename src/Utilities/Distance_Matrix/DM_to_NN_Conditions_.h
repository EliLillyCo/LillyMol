#ifndef IW_DM_TO_NN_CONDITIONS_IMPLEMENTATION_H
#define IW_DM_TO_NN_CONDITIONS_IMPLEMENTATION_H

#include "IWDistanceMatrixBase.h"
#include "iwdmsupport.h"

template <typename T>
DM_to_NN_Conditions<T>::DM_to_NN_Conditions ()
{
  _min_neighbours = 0;
  _max_neighbours = 0;

  _min_distance = static_cast<T> (0);
  _max_distance = static_cast<T> (0);

  _dash_ho = 0;

  return;
}

template <typename T>
int
DM_to_NN_Conditions<T>::ok () const
{
  if (_min_neighbours < 0 || _max_neighbours < 0)
    return 0;

// Convoluted logic to avoid compiler warnings about testing negatives for unsigned variables

  if (_min_distance > static_cast<T> (0))
    ;
  else if (static_cast<T> (0) == _min_distance)
    ;
  else
    return 0;

  if (_max_distance > static_cast<T> (0))
    ;
  else if (static_cast<T> (0) == _max_distance)
    ;
  else
    return 0;

  if (_min_neighbours > _max_neighbours)
    return 0;

  if (_min_distance > _max_distance)
    return 0;

  return 1;
}

template <typename T>
int
DM_to_NN_Conditions<T>::debug_print (std::ostream & os) const
{
  os << "DM_to_NN_Conditions::debug_print\n";
  if (_min_neighbours > 0)
    os << " min neighbours " << _min_neighbours << endl;

  if (_max_neighbours > 0)
    os << " max neighbours " << _max_neighbours << endl;

  if (_min_distance > 0)
    os << " min distance " << _min_distance << endl;

  if (_max_distance > 0)
    os << " max distance " << _max_distance << endl;

  if (_dash_ho)
    os << " -h -o\n";

  return os.good ();
}

template <typename T>
int
DM_to_NN_Conditions<T>::recognised (const IWString & n,
                                    int error_occurred)
{
  if (n.starts_with ("maxd="))
  {
    if (! parse_directive (n, _max_distance))
      error_occurred = 1;

    return 1;
  }

  if (n.starts_with ("mind="))
  {
    if (! parse_directive (n, _min_distance))
      error_occurred = 1;

    return 1;
  }
  else if (n.starts_with ("maxn"))
  {
    if (! parse_directive (n, _max_neighbours))
      error_occurred = 1;

    return 1;
  }
  else if (n.starts_with ("minn"))
  {
    if (! parse_directive (n, _min_neighbours))
      error_occurred = 1;

    return 1;
  }
  else if ("ho" == n)
  {
    _dash_ho = 1;

    return 1;
  }

  return 0;
}

#endif
