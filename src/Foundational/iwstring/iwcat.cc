#include <stdlib.h>
#include "iostream.h"
#include "unistd.h"
#include "stdio.h"

#include "cmdline.h"
#include "iwstring_data_source.h"

int 
main (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    return 1;
  }

  iwstring_data_source input (cl[0]);


  if (! input.ok ())
  {
    cerr << "Cannot open '" << cl[0] << "'\n";
    return 3;
  }


#define BUFFERED_OUTPUT
#ifdef BUFFERED_OUTPUT
  IWString output_buffer;
  output_buffer.resize (8192 * 2);
#else
  char * newline = "\n";
#endif

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
#ifdef BUFFERED_OUTPUT
    output_buffer += buffer;
    output_buffer += '\n';
    if (output_buffer.length () > 8192)
    {
//    write (1, output_buffer.rawchars (), output_buffer.length ());
      cout.write (output_buffer.rawchars (), output_buffer.length ());
      output_buffer.resize_keep_storage (0);
    }
#else
    fwrite (buffer.rawchars (), 1, buffer.length (), stdout);
    fprintf (stdout, "\n");
//  write (1, buffer.rawchars (), buffer.length ());
//  write (1, newline, 1);
#endif
  }

#ifdef BUFFERED_OUTPUT
  if (output_buffer.nchars ())
    cout.write (output_buffer.rawchars (), output_buffer.length ());
#endif

  return 0;
}
