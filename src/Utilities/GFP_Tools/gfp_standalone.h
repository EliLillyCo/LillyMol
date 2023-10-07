#ifndef GFP_STANDALONE_FILE_H
#define GFP_STANDALONE_FILE_H

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <memory>
#include <limits>
#include <fstream>
#include <queue>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#define IWQSORT_IMPLEMENTATION

#include "cmdline.h"
#include "iw_tdt.h"
#include "iwqsort.h"
#include "accumulator.h"
#include "iwdigits.h"

#include "gfp_standard.h"


class Extern_Gfp_Standalone
{
public:
	void extern_usage(int rc);

	int extern_gfp_standalone(int argc, char ** argv);

};

#endif
