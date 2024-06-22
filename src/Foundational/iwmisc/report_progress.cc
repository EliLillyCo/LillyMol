#include <stdlib.h>
#include <limits>

#include "Foundational/cmdline/cmdline.h"

#define REPORT_PROGRESS_IMPLEMENTATION

#include "report_progress.h"

template class Report_Progress_Template<uint64_t>;
template int Report_Progress_Template<uint64_t>::initialise<Command_Line>(Command_Line&, char, int);
