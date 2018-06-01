#include <stdlib.h>

#define REPORT_PROGRESS_IMPLEMENTATION

#include "report_progress.h"
#include "cmdline.h"

template class Report_Progress_Template<unsigned int>;
template int Report_Progress_Template<unsigned int>::initialise (Command_Line & cl, char flag, int verbose);
