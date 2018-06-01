/*
  This is like Perl's split function.
  We return an array of the strings obtained by strtok'ing buffer.
*/

#include <stdlib.h>
#include <string.h>
#include <String.h>
#include <iwaray.h>
#include <iwstring.h>

/*
resizable_array_p<String> *
tokenise (const String &buffer, const char *pattern, int ntokens)
{
  assert (NULL != pattern);

  if (0 == buffer.length ())
    return NULL;

  if (NTOKENS_UNSPECIFIED != ntokens)
  {
    assert (ntokens > 0);
  }

  int len_pattern = strlen (pattern);
  assert (len_pattern > 0);

  resizable_array_p<String> * result = new resizable_array_p<String>;

  int i = 0;
  int j = buffer.index (pattern, i);

  while (j > 0)
  {
    String *tmp = new String (&buffer[i], j);
    result->add (tmp);

    j = i;
    j = buffer.index (buffer, i);
  }

  if (0 == result->number_elements ())
  {
    delete result;
    return NULL;
  }

  return result;
}*/

template class resizable_array_p<String>;
template class resizable_array_base<String *>;

resizable_array_p<String> *
tokenise (const char *buffer, const char *pattern, int ntokens)
{
  assert (NULL != buffer);
  int buflen = strlen (buffer);

  if (0 == buflen)
    return NULL;

  assert (NULL != pattern);
  int len_pattern = strlen (pattern);
  assert (len_pattern > 0);

  if (NTOKENS_UNSPECIFIED != ntokens)
  {
    assert (ntokens > 0);
  }

  resizable_array_p<String> * result = new resizable_array_p<String>;

  const char *c = strpbrk (buffer, pattern);
//fprintf (stderr, "c set to %x\n", c);

  while (c)
  {
    int nchars = c - buffer;
    if (nchars > 0)
    {
      String * new_token = new String (buffer, nchars);
      result->add (new_token);
//    cerr << "Added string '" << *new_token << "'\n";
    }

    buffer += nchars;
    buffer += strspn (buffer, pattern);

//  (void) fprintf (stderr, "Buffer incremented to %x\n", buffer);
    c = strpbrk (buffer, pattern);
//  fprintf (stderr, "c set to %x, points to '%s'\n", c, c);
  }

  buflen = strlen (buffer);
  if (buflen > 0)
  {
    String * last_token = new String (buffer, buflen);
    result->add (last_token);
  }

  if (0 == result->number_elements ())
  {
    delete result;
    result = NULL;
  }

  return result;
}

#include <data_source.h>

/*int
main ()
{
  data_source *input = new_data_source (stdin);

  input->set_strip_newlines ();

  const char *buffer;
  while (NULL != (buffer = input->next_record ()))
  {
    cout << "Buffer: " << buffer;

    resizable_array_p<String> *tk = tokenise (buffer, " ");
    if (NULL == tk)
    {
      (void) fprintf (stderr, "No tokens\n");
    }
    else
    {
      (void) fprintf (stdout, "%d tokens\n", tk->number_elements ());
      for (int i = 0; i < tk->number_elements (); i++)
        cout << i << "'" << tk->item (i) << "'\n";
      delete tk;
    }
  }
  return 0;
} */
