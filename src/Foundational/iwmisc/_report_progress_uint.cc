
#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "report_progress.h"

template class Report_Progress_Template<unsigned int>;
template int Report_Progress_Template<unsigned int>::initialise(Command_Line & cl, char flag, int verbose);
