#include <stdlib.h>

#define ACTIVITY_DATA_IMPLEMENATION_H

#include "activity_data_from_file.h"
#include "cmdline.h"

template class Activity_Data_From_File<float>;

template int Activity_Data_From_File<float>::construct_from_command_line(const Command_Line & ,
                                                        char , int);
