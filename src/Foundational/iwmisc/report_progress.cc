#include <stdlib.h>
#include <limits>

#include "Foundational/cmdline/cmdline.h"

#define REPORT_PROGRESS_IMPLEMENTATION

#include "report_progress.h"

template class Report_Progress_Template<uint64_t>;
template int Report_Progress_Template<uint64_t>::initialise<Command_Line>(Command_Line&, char, int);

namespace report_progress_internal {
void
DisplayHelpMessage(std::ostream& output) {
  output << R"(ReportProgress triggers every <n> items processed.
Usually activated via the -r option, but not always.
The following options are recognised.
 -r time        at each occurrence also report time spent.
 -r human       write numbers as human readable, 1,000,000
)";

}

}  // namespace report_progress_internal
