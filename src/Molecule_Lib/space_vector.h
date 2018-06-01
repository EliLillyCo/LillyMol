#ifndef IW_3DVECTOR_H
#define IW_3DVECTOR_H

#include <assert.h>
#include <cmath>

#include "iwmtypes.h"

class const_IWSubstring;

template <typename T>
class Space_Vector
{
  protected:
    T _x;
    T _y;
    T _z;
    void _default_values ();

  public:
    Space_Vector ();
    Space_Vector (T, T, T);
    Space_Vector (const Space_Vector<T> &);
    ~Space_Vector ();

    T x () const { return _x; }
    T y () const { return _y; }
    T z () const { return _z; }

    T & x () { return _x; }
    T & y () { return _y; }
    T & z () { return _z; }

    void setxyz (T, T, T);
    void setxyz (const Space_Vector<T> &);

//  function to help speed up whim

    void getxyz (double * c) const { c[0] = static_cast<double>(_x);
                                     c[1] = static_cast<double>(_y);
                                     c[2] = static_cast<double>(_z);}

    int read (const const_IWSubstring &, char=' ');

    void add (T, T, T);
    void translate (T tx, T ty, T tz) { Space_Vector<T>::add (tx, ty, tz);}
    void translate (const Space_Vector<T> &);
    void  normalise ();
    T norm () const;
    T length () const { return norm ();}
    T normsquared () const { return _x * _x + _y * _y + _z * _z;}

    Space_Vector<T> & operator =  (const Space_Vector<T> &);
    Space_Vector<T> operator +  (T) const;
    Space_Vector<T> operator +  (const Space_Vector<T> &) const;
    Space_Vector<T> operator -  (T) const;
    Space_Vector<T> operator -  (const Space_Vector<T> &) const;
    Space_Vector<T> operator *  (T) const;
    Space_Vector<T> operator /  (T) const;
    Space_Vector<T> & operator *= (const Space_Vector<T> &);          // cross product
    Space_Vector<T> operator *  (const Space_Vector<T> &) const;    // cross product
    Space_Vector<T> operator -  ()const ;

    void        operator += (T);
    void        operator += (const Space_Vector<T> &);
    void        operator -= (T);
    void        operator -= (const Space_Vector<T> &);
    void        operator *= (T);
    void        operator /= (T);
    T           operator ^ (const Space_Vector<T> &) const;       // dot product
    int         operator == (const Space_Vector<T> &)const ;
    int         operator != (const Space_Vector<T> &)const ;
    void        cross_product (const Space_Vector<T> &);

    angle_t     angle_between              (const Space_Vector<T> &) const;
    angle_t     angle_between_unit_vectors (const Space_Vector<T> &) const;

    T           distance (const Space_Vector<T> &) const;
    T           distance_squared (const Space_Vector<T> &) const;

    T           dot_product (const Space_Vector<T> &) const;

    Space_Vector<T> form_unit_vector (const Space_Vector<T> &) const;

    angle_t     angle_between (const Space_Vector<T> & a1, const Space_Vector<T> & a2) const;   // bond angle, we assume  a1 - this - a2, we return the angle a1-this-a2

    bool        closer_than (const Space_Vector<T> &, T) const;   // tries to be efficient
};

#include <iostream>

template <typename T> std::ostream & operator << (std::ostream &, const Space_Vector<T> &);

#if defined(SPACE_VECTOR_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

#include <stdlib.h>
#include <iostream>
#include <algorithm>

using std::cerr;
using std::endl;

#include "iwstring.h"

template <typename T>
void
Space_Vector<T>::_default_values ()
{
  _x = _y = _z = static_cast<T>(0.0);

  return;
}

template <typename T>
Space_Vector<T>::Space_Vector ()
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "In default constructor " << this << endl;
#endif

  _default_values();
}

template <typename T>
Space_Vector<T>::Space_Vector (T x, T y, T z)
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "Individual value constructor " << this << endl;
#endif

  _x = x;
  _y = y;
  _z = z;

  return;
}

template <typename T>
Space_Vector<T>::Space_Vector (const Space_Vector<T> & rhs)
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "In copy constructor " << this << " rhs is " << rhs << endl;
#endif

  _x = rhs._x;
  _y = rhs._y;
  _z = rhs._z;

  return;
}

template <typename T>
Space_Vector<T>::~Space_Vector ()
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "descructor called for " << this << endl;
#endif

//_default_values();   no reason to reset these
}

template <typename T>
Space_Vector<T> &
Space_Vector<T>::operator = (const Space_Vector<T> & rhs)
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "Assignment operator between " << this << " and " << &rhs << endl;
#endif

  _x = rhs._x;
  _y = rhs._y;
  _z = rhs._z;

  return *this;
}

template <typename T>
void
Space_Vector<T>::normalise ()
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "Normalising " << this << endl;
#endif

  double mynorm = sqrt(static_cast<double>(_x) * static_cast<double>(_x) + static_cast<double>(_y) * static_cast<double>(_y) + static_cast<double>(_z) * static_cast<double>(_z));

  if (mynorm > 0.0)
  {
    _x = static_cast<T>(_x / mynorm);
    _y = static_cast<T>(_y / mynorm);
    _z = static_cast<T>(_z / mynorm);
  }

  return;
}

template <typename T>
T
Space_Vector<T>::norm () const
{
  T norm = static_cast<T>(sqrt(static_cast<double>(_x) * static_cast<double>(_x) + static_cast<double>(_y) * static_cast<double>(_y) + static_cast<double>(_z) * static_cast<double>(_z)));

  return norm;
}

template <typename T>
void
Space_Vector<T>::setxyz (T x, T y, T z)
{
  _x = x;
  _y = y;
  _z = z;

  return;
}
template <typename T>
void
Space_Vector<T>::setxyz (const Space_Vector<T> & rhs)
{
  _x = rhs._x;
  _y = rhs._y;
  _z = rhs._z;

  return;
}

template <typename T>
T
Space_Vector<T>::operator ^ (const Space_Vector<T> & v2) const
{
  return _x * v2._x + _y * v2._y + _z * v2._z;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator * (const Space_Vector<T> & v2) const
{
  Space_Vector<T> result(_y * v2._z - _z * v2._y,
                         _z * v2._x - _x * v2._z,
                         _x * v2._y - _y * v2._x);

  return result;
}

template <typename T>
Space_Vector<T> &
Space_Vector<T>::operator *= (const Space_Vector<T> & v2)
{
  cross_product(v2);

  return *this;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator + (T extra) const
{
  Space_Vector<T> result(_x + extra, _y + extra, _z + extra);

  return result;
}

template <typename T>
void
Space_Vector<T>::operator += (T extra)
{
  _x += extra;
  _y += extra;
  _z += extra;
}

template <typename T>
void
Space_Vector<T>::operator += (const Space_Vector<T> & v2)
{
  _x += v2._x;
  _y += v2._y;
  _z += v2._z;

  return;
}

template <typename T>
void
Space_Vector<T>::add (T xx, T yy, T zz)
{
  _x += xx;
  _y += yy;
  _z += zz;
}

template <typename T>
void
Space_Vector<T>::translate (const Space_Vector<T> & whereto)
{
  _x += whereto._x;
  _y += whereto._y;
  _z += whereto._z;

  return;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator * (T factor) const
{
  Space_Vector<T> result( _x * factor, _y * factor, _z * factor);

  return result;
}

template <typename T>
void
Space_Vector<T>::operator *= (T factor)
{
  _x *= factor;
  _y *= factor;
  _z *= factor;
}

template <typename T>
void
Space_Vector<T>::operator /= (T factor)
{
  _x /= factor;
  _y /= factor;
  _z /= factor;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator / (T factor) const
{
  assert (static_cast<T>(0.0) != factor);

  Space_Vector<T> result(_x / factor, _y / factor, _z / factor);
  
  return result;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator + (const Space_Vector<T> & v2) const
{
  Space_Vector<T> result(_x + v2._x, _y + v2._y, _z + v2._z);

  return result;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator - (T extra) const
{
  Space_Vector<T> result(_x - extra, _y - extra, _z - extra);

  return result;
}

template <typename T>
void
Space_Vector<T>::operator -= (T extra)
{
  _x -= extra;
  _y -= extra;
  _z -= extra;

  return;
}

template <typename T>
void
Space_Vector<T>::operator -= (const Space_Vector<T> & v2)
{
  _x -= v2._x;
  _y -= v2._y;
  _z -= v2._z;

  return;
}

template <typename T>
Space_Vector<T> 
Space_Vector<T>::operator - (const Space_Vector<T> & v2) const
{
  Space_Vector<T> result(_x - v2._x, _y - v2._y, _z - v2._z);

  return result;
}

template <typename T>
int
Space_Vector<T>::operator == (const Space_Vector<T> & v2) const
{
  return (_x == v2._x && _y == v2._y && _z == v2._z);
}

// not equal if any components differ

template <typename T>
int
Space_Vector<T>::operator != (const Space_Vector<T> & v2) const
{
  return (_x != v2._x || _y != v2._y || _z != v2._z);
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::operator - () const
{
  Space_Vector<T> result(-_x, -_y, -_z);

  return result;
}

template <typename T>
angle_t
Space_Vector<T>::angle_between_unit_vectors (const Space_Vector<T> & v1) const
{
  double tmp = static_cast<double>(_x) * static_cast<double>(v1._x) + static_cast<double>(_y) * static_cast<double>(v1._y) + static_cast<double>(_z) * static_cast<double>(v1._z);    // the dot product

  if (fabs(tmp) <= 1.0)    // that's good
    ;
  else if (fabs(tmp) < 1.00001)
  {
//  cerr << "Space_Vector::angle_between_unit_vectors:numerical roundoff discarded " << (fabs(tmp) - 1.0) << endl;
    tmp = 1.0;
  }
  else
  {
    cerr << "Space_Vector::angle_between_unit_vectors:vector might not be a unit vector, tmp = " << tmp << endl;
    abort();
  }

  return static_cast<angle_t>(acos(tmp));
}

template <typename T>
angle_t
Space_Vector<T>::angle_between (const Space_Vector<T> & v1) const
{
  Space_Vector<T> lhs(_x, _y, _z);
  lhs.normalise();

  Space_Vector<T> rhs(v1._x, v1._y, v1._z);
  rhs.normalise();

  return lhs.angle_between_unit_vectors(rhs);
}

template <typename T>
void
Space_Vector<T>::cross_product (const Space_Vector<T> & v2)
{
  double xorig = static_cast<double>(_x);
  double yorig = static_cast<double>(_y);
  double zorig = static_cast<double>(_z);

  _x = static_cast<T>(yorig * static_cast<double>(v2._z) -  zorig * static_cast<double>(v2._y));
  _y = static_cast<T>(zorig * static_cast<double>(v2._x) -  xorig * static_cast<double>(v2._z));
  _z = static_cast<T>(xorig * static_cast<double>(v2._y) -  yorig * static_cast<double>(v2._x));

  return;
}

template <typename T>
int
Space_Vector<T>::read (const const_IWSubstring & buffer, char separator)
{
  if (buffer.nwords(separator) < 3)
  {
    cerr << "vector::read: must have at least 3 tokens '" << buffer << "' has " << buffer.nwords(separator) << endl;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword(token, i, separator);
  if (! token.numeric_value(_x))
  {
    cerr << "vector::read: cannot parse first token as numeric\n";
    return 0;
  }

  (void) buffer.nextword(token, i, separator);
  if (! token.numeric_value(_y))
  {
    cerr << "vector::read: cannot parse second token as numeric\n";
    return 0;
  }

  (void) buffer.nextword(token, i, separator);
  if (! token.numeric_value(_z))
  {
    cerr << "vector::read: cannot parse third token as numeric\n";
    return 0;
  }

  return 1;
}

template <typename T>
T
Space_Vector<T>::distance (const Space_Vector<T> & rhs) const
{
  return static_cast<T>(sqrt(static_cast<double>(_x - rhs._x) * static_cast<double>(_x - rhs._x) + 
                             static_cast<double>(_y - rhs._y) * static_cast<double>(_y - rhs._y) +
                             static_cast<double>(_z - rhs._z) * static_cast<double>(_z - rhs._z)));
}

template <typename T>
T
Space_Vector<T>::distance_squared (const Space_Vector<T> & rhs) const
{
  return static_cast<T>(static_cast<double>(_x - rhs._x) * static_cast<double>(_x - rhs._x) + 
                        static_cast<double>(_y - rhs._y) * static_cast<double>(_y - rhs._y) +
                        static_cast<double>(_z - rhs._z) * static_cast<double>(_z - rhs._z));
}

template <typename T>
T
Space_Vector<T>::dot_product (const Space_Vector<T> & rhs) const
{
  return _x * rhs._x + _y * rhs._y + _z * rhs._z;
}

template <typename T>
Space_Vector<T>
Space_Vector<T>::form_unit_vector (const Space_Vector<T> & rhs) const
{
#ifdef DEBUG_SPACE_VECTOR
  cerr << "In form_unit_vector from " << this << " to " << &rhs << endl;
#endif

  Space_Vector<T> rc(_x - rhs._x, _y - rhs._y, _z - rhs._z);

  rc.normalise();

  return rc;
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const Space_Vector<T> & qq)
{
  os << "(" << qq.x() << "," << qq.y() << "," << qq.z() << ")";

  return os;
}

#ifdef NO_LONGER_USED_ASDASDASD
static angle_t
internal_angle_between (double a,
                        double b,
                        double c)
{
  double x = (a * a - b * b - c * c) / (2.0 * b);

  double tmp = (b + x) / a;

  if (tmp > 1.0)
    return static_cast<angle_t>(0.0);
  else if (tmp < -1.0)
    return static_cast<angle_t>(M_PI * 0.5);

  return static_cast<angle_t>(acos(tmp));
}
#endif

template <typename T>
angle_t
Space_Vector<T>::angle_between (const Space_Vector<T> & a1,
                                const Space_Vector<T> & a2) const
{
  Space_Vector<double> ab(a1.x() - _x, a1.y() - _y, a1.z() - _z);
  Space_Vector<double> ac(a2.x() - _x, a2.y() - _y, a2.z() - _z);

  ab.normalise();
  ac.normalise();

  return ab.angle_between_unit_vectors(ac);
}

template <typename T>
bool
Space_Vector<T>::closer_than (const Space_Vector<T> & rhs, T d) const
{
  T sum = static_cast<T>(0.0);

  T q = fabs(_x - rhs._x);
  if (q > d)
    return false;

  sum += q * q;

  q = fabs(_y - rhs._y);
  if (q > d)
    return false;

  sum += q * q;

  if (sqrt(sum) > d)
    return false;

  q = fabs(_z - rhs._z);

  if (q > d)
    return false;

  sum += q * q;

  return sqrt(sum) <= d;
}

#endif

#endif
