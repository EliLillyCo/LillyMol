#include <stdlib.h>

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <libgen.h>

#include <iostream>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

static int
nextdir_reading_file_of_directories (const IWString & pwd,
                                     int argc, char ** argv)
{
//cerr << "File of directories checking '" << pwd << "'\n";
  iwstring_data_source input(argv[1]);

  if (! input.good())
  {
    cerr << "Cannot open '" << argv[0] << "'\n";
    return 2;
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.ends_with('/'))
      buffer.chop();

    if (buffer == pwd)
      ;
    else if (buffer.starts_with(pwd) && isspace(buffer[pwd.length()]))
      ;
    else
      continue;

    if (! input.next_record(buffer))
    {
      fprintf(stderr, "End of file of directories\n");
      return 2;
    }

    IWString dirname(buffer);

    dirname.truncate_at_first(' ');

    fprintf(stdout, "../%s", dirname.null_terminated_chars());
    return 0;
  }

  fprintf(stderr, "No next directory\n");

  return 1;
}

class String_Comparitor
{
  private:
  public:
    int operator()(const IWString *, const IWString *) const;
};

int
String_Comparitor::operator()(const IWString * s1, const IWString * s2) const
{
  return s1->strcasecmp(*s2);
}

#define NEED_TEMPLATE_INSTANTIATION
#ifdef NEED_TEMPLATE_INSTANTIATION
template void resizable_array_base<IWString*>::iwqsort<String_Comparitor>(String_Comparitor&);
template void iwqsort<IWString*, String_Comparitor>(IWString**, int, String_Comparitor&);
template void iwqsort<IWString*, String_Comparitor>(IWString**, int, String_Comparitor&, void*);
template void compare_two_items<IWString*, String_Comparitor>(IWString**, String_Comparitor&, void*);
template void move_in_from_left<IWString*, String_Comparitor>(IWString**, int&, int&, int, String_Comparitor&, void*);
template void move_in_from_right<IWString*, String_Comparitor>(IWString**, int&, int&, String_Comparitor&);
template void swap_elements<IWString*>(IWString*&, IWString*&, void*);
#endif

//#define DEBUG_NEXTDIR

/*
  Returns the next sub-directory from where we are
*/

static int
nextdir(int argc, char ** argv)
{
  char tmp[512];    // first determine PWD

  if (NULL == getcwd(tmp, sizeof(tmp)))
  {
    fprintf(stderr, "Cannot get PWD!\n");
    return 2;
  }

  IWString pwd;

  if (argc > 1)
  {
    pwd = basename(tmp);
    return nextdir_reading_file_of_directories (pwd, argc, argv);
  }

  pwd << "../" << basename(tmp);

  DIR * parent_directory = opendir("..");

  if (NULL == parent_directory)
  {
    fprintf(stderr, "Cannot open parent directory\n");
    return 1;
  }

  resizable_array_p<IWString> directories;

  struct dirent * d;

  while (NULL != (d = readdir(parent_directory)))
  {
    if ('.' == (d->d_name)[0])    // we don't process these
      continue;

    IWString * tmp = new IWString(32);

    (*tmp) << "../" << d->d_name;

#ifdef DEBUG_NEXTDIR
    fprintf (stderr, "Examining directory '%s', dir? %d\n", d->d_name, dash_d(tmp->null_terminated_chars()));
#endif

    if (! dash_d(tmp->null_terminated_chars()))
      continue;

    directories.add(tmp);
  }

  int n = directories.number_elements();

#ifdef DEBUG_NEXTDIR
  fprintf(stderr, "Found %d directories\n", n);
#endif

  if (1 == n)
  {
    fprintf(stderr, "No next directory found, staying put\n");
    fprintf (stdout, ".");
  }

  String_Comparitor s;
  directories.iwqsort(s);

  int return_next_directory = 0;
  for (int i = 0; i < n; i++)
  {
    IWString & s = *(directories[i]);

#ifdef DEBUG_NEXTDIR
    fprintf (stderr, "Examining sorted '%s'\n", s.null_terminated_chars());
#endif

    if (return_next_directory)
    {
      fprintf(stdout, "%s", s.null_terminated_chars());
      return 0;
    }

    if (s == pwd)
    {
      return_next_directory = 1;
      continue;
    }
  }

// If we get to here, no next directory found

  fprintf(stderr, "No next directory found, staying put\n");

  fprintf(stdout, ".");

  return 0;
}

int
main(int argc, char ** argv)
{
  int rc = nextdir(argc, argv);

  return rc;
}
