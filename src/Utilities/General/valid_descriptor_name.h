#ifndef IWVALID_CHARS_IN_DNAME
#define IWVALID_CHARS_IN_DNAME

#include "cmdline.h"

class Valid_Descriptor_Name
{
  private:
    int _max_name_length;
    char _valid_character[256];

  public:
    Valid_Descriptor_Name ();

    int valid (const IWString & s) const;
    int valid (const IWString & s, IWString &) const;

    int construct_from_command_line (Command_Line &, char, char, int);
};

extern int display_standard_valid_descriptor_name_options (std::ostream & os, char mflag, char vflag);

#endif
