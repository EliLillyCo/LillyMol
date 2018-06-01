#include <stdlib.h>
#include <stdio.h>

// Convert a string to hex

static int
hex_to_int (const char * s,
            int nchars,
            unsigned int * result)
{
  int i;
  int bshift;

  if (s[0] != '0' || s[1] != 'x')
    return 0;

  s++;
  nchars -=2;

  if (nchars > 8)
    return 0;

  result = 0;
  for (i = nchars - 1, bshift = 0; i >= 0; i--, bshift += 4)
  {
    register char c = s[i];
    unsigned int j;
    if (c >= '0' && c <= '9')
    {
      j = c - '0';
    }
    else if (c >= 'a' && c <= 'f')
    {
      j = c - 'a' + 0xa;
    }
    else
      return 0;

    j = j << bshift;
    result |= j;

  }

  return 1;
}

int
main ()
{
  char buffer[32];

  strcpy (buffer, "hello");

  while (strlen (buffer) > 1)
  {
    char * nl;
    int rc;
    unsigned int result;

    fscanf (stdin, "%s", buffer);
    nl = strchr (buffer, '\n');
    if (NULL == nl)
      continue;

    *nl = '\0';

    rc = hex_to_int (buffer, strlen (buffer), result);

    if (0 == rc)
      fprintf (stderr, "Conversion failed\n";
    else
      fprintf (stderr, "Conversion from hex successful, result is %x\n", result);
  }

  return 0;
}
