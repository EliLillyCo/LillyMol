#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "iwstring.h"

/*
  Function to substitute one string for another in a string
  We abort on out of bounds.....
*/

int
string_substitute (char *buffer, int len, 
                   const char *old_string, const char *new_string)
{
  assert (NULL != buffer);
  assert (NULL != old_string);
  assert (NULL != new_string);

  char *c;
  if (NULL == (c = strstr (buffer, old_string)))    /* not present at all */
    return 0;

  char *new_buffer;
  if (NULL == (new_buffer = new char[len]))
  {
    abort ();
  }

  char *out_of_bounds = new_buffer + len;

  int number_changes = 0;
  char *t1 = buffer;
  char *t2 = new_buffer;
  while (NULL != c)
  {
    while (t1 < c)
      *t2++ = *t1++;
    if (t2 >= out_of_bounds)
      abort ();
    *t2 = '\0';
    (void) strcat (t2, new_string);
    t1 += strlen (old_string);
    t2 += strlen (new_string);

    number_changes++;
    c = strstr (t1, old_string);
  }

  strcat (new_buffer, t1);

/*
  At this stage, the string is transformed. Copy back into buffer
*/

  strcpy (buffer, new_buffer);
  delete new_buffer;

  return number_changes;
}
