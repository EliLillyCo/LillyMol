#ifndef BSQUARED_H
#define BSQUARED_H

extern int compute_b_squared (const float * predicted,
                   float * ybest,
                   int n,
                   double & Bsquared);

extern int compute_b_squared (const float * predicted,
                   float * ybest,
                   int n,
                   double * Bsquared);
#endif
