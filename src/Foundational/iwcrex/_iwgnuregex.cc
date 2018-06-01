#include <stdlib.h>
#include <iostream>
using namespace std;

#define IWCREX_IMPLEMENTATION
#include "iwcrex.h"
#include "iwgnuregex.h"

template class IW_Regular_Expression_Template<IW_gnu_regex>;
