#ifndef POPCOUNT_H
#define POPCOUNT_H

#if defined(__x86_64__)
#include <nmmintrin.h>
#define POPCOUNT _mm_popcnt_u64
#elif defined(__i386__)
#include <nmmintrin.h>
#define POPCOUNT _mm_popcnt_u32
#elif __APPLE__
#define POPCOUNT __builtin_popcountll
#else
#include <bit>
#define POPCOUNT std::popcount
#endif

#endif /* POPCOUNT_H */
