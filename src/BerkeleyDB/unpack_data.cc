#include <iostream>

#include "unpack_data.h"

using std::cerr;
using std::endl;

static int
get_pack_directive (const IWString & upformat,
                    int & ndx,
                    char & directive,
                    int & count)
{
  if (ndx >= upformat.length())
    return 0;

  count = 0;
  while (isdigit(upformat[ndx]))
  {
    int tmp = upformat[ndx] - '0';
    count = count * 10 + tmp;

    ndx++;

    if (ndx == upformat.length())    // must have a letter following
      return 0;
  }

  if (0 == count)
    count = 1;

  if (ndx == upformat.length())
    return 0;

  directive = upformat[ndx];

  ndx++;

  return 1;
}

/*
  This function is made complex by the possibility that p may not be
  aligned properly
*/

template <typename T>
int
write_as(const char * p, IWString_and_File_Descriptor & output)
{
  T tmp;

  char * ptmp = reinterpret_cast<char *>(&tmp);

  memcpy(ptmp, p, sizeof(T));

  output << tmp;

  return output.good();
}

template int write_as<unsigned int>(const char * p, IWString_and_File_Descriptor & output);
template int write_as<int>(const char * p, IWString_and_File_Descriptor & output);
template int write_as<unsigned char>(const char * p, IWString_and_File_Descriptor&);

template <typename T>
int
write_as(const char * & p, int count, IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < count; i++)
  {
    write_as<T>(p, output);
    p += sizeof(T);
  }

  return 1;
}

template int write_as<unsigned int>(const char * &, int, IWString_and_File_Descriptor &);
template int write_as<int>(const char * &, int, IWString_and_File_Descriptor &);
template int write_as<unsigned char>(const char * &, int, IWString_and_File_Descriptor&);

#ifdef NOW_DONE_WITH_TEMPLATE

static int
write_as_int(const char * p, IWString_and_File_Descriptor & output)
{
  int tmp;

  char * ptmp = reinterpret_cast<char *>(&tmp);

  memcpy(ptmp, p, 4);

  output << tmp;

  return output.good();
}

static int
write_as_int(const char * & p, int count, IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < count; i++)
  {
    write_as_int(p, output);
    p += sizeof(int);
  }

  return 1;
}

#endif

static int
write_as_hex(char p, IWString_and_File_Descriptor & output)
{
//unsigned char c = p;

//output << hex << int(c) << dec;

  return output.good();
}

static int
write_as_hex(const char * & p, int count, IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < count; i++)
  {
    write_as_hex(*p, output);

    p++;
  }

  return output.good();
}

static int
write_some_data (const char * & p,
                 char directive,
                 int count,
                 IWString_and_File_Descriptor & output)
{
  switch(directive)
  {
    case('i'):
    {
      write_as<int>(p, count, output);
      break;
    }
    case ('u'):
    {
      write_as<unsigned int>(p, count, output);
      break;
    }
    case ('h'):
    {
      write_as_hex(p, count, output);
      break;
    }
    case ('c'):
    {
      write_as<unsigned char>(p, count, output);
      break;
    }
    default:
    {
      cerr << "Unrecognised unpack directive '" << directive << "'\n";
      return 0;
    }
  }

  return output.good();
}

template <typename T>
int
Unpack_Binary_Data::write_unpacked_data(const char * p,
                                        const int len,
                                        T & output) const
{
  if (0 == _unpack_format.length())
  {
    cerr << "Unpack_Binary_Data::write_unpacked_data:no pattern\n";
    return 0;
  }

  int i = 0;     // index into pack string
  const char * pend = p + len; 

  char previous_directive = ' ';
  char directive;
  int count;
  while (get_pack_directive(_unpack_format, i, directive, count))
  {
    if (i > 0)
      output << ' ';

    write_some_data(p, directive, count, output);

    if (p >= pend)
      return 1;

    previous_directive = directive;
  }

// If our unpack directive didn't consume the data, repeat...

  while (p < pend)
  {
    output << ' ';

    if (! write_some_data(p, previous_directive, 1, output))
      return 0;
  }

  return 1;
}

template int Unpack_Binary_Data::write_unpacked_data(const char *, const int, IWString_and_File_Descriptor &) const;
