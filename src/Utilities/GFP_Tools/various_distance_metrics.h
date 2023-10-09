#ifndef IWVDM_H
#define IWVDM_H

#include "Foundational/iwbits/iwbits.h"

extern similarity_type_t fligner_verducci_blower (int nb,
                         int nset1,
                         int nset2,
                         int n00,
                         int n11);

/*
  The relative weights given to T1 and T0
*/

extern int set_fvb_ratios (similarity_type_t, similarity_type_t);

// Comma separated form

extern int set_fvb_ratios (const const_IWSubstring & fvb);

#endif
