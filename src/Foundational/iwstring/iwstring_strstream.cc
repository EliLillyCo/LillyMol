#include <stdlib.h>

#include "iwstring.h"

 
#if (__GNUC_MINOR__ == 95)
#include <strstream>

IWString &
IWString::operator = (ostrstream & zstream)
{
  int nchars = zstream.pcount ();
  char * c = zstream.str ();
  IWString::strncpy (c, nchars);
 
#ifdef __GNUG__
  zstream.freeze (0);   // allow zstream to delete the buffer, we have a copy
#else
  delete c;
#endif
 
  return * this;
}

/*
  Append the text in ZEXTRA
*/
 
IWString &  
operator << (IWString & s, ostrstream & zextra)
{
  assert (s.ok ());
 
  int nchars = zextra.pcount ();
  char * c = zextra.str ();
 
  s.strncat (c, nchars);
 
#ifdef __GNUG__
  zextra.freeze (0);
#else
  delete c;
#endif
 
  return s; 
}

#else
#include <sstream>

/*IWString &
IWString::operator = (ostringstream & zstream)
{
  std::string c = zstream.str ();

  IWString::operator = (c);

  return *this;
}

IWString &
operator << (IWString & s, ostringstream & zextra)
{
  assert (s.ok ());

  string c = zextra.str ();

  s.strncat (c.data (), c.length ());

  return s;
}*/

#endif

