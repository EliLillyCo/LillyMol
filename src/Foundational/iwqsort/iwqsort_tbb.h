#ifndef IW_TBB_QSORT_H
#define IW_TBB_QSORT_H

#include <stdlib.h>

#include "tbb/parallel_for.h"

template <typename RandomAccessIterator, typename Compare>
struct quick_sort_range
{
  static const size_t grainsize = 500;
  const Compare & comp;
  RandomAccessIterator begin;
  size_t size;

  quick_sort_range (RandomAccessIterator begin_,
                    size_t size_,
                    const Compare & comp_) :
                      comp(comp_), begin(begin_), size(size_) {}
  bool empty () const { return 0 == size;}
  bool is_divisible() const { return size >= grainsize;}

  quick_sort_range (quick_sort_range & range, tbb::split) : comp(range.comp) 
  {
    RandomAccessIterator array = range.begin;
    RandomAccessIterator key0 = range.begin;
    size_t m = range.size/2u;

    std::swap(array[0], array[m]);

    size_t i = 0;
    size_t j = range.size;

    // Partition interval [i+1,j-1] with key *key0

    for (;;)
    {
      __TBB_ASSERT(i < j, NULL);
      // loop must terminate since array[l]==*key0
      do
      {
        --j;
        __TBB_ASSERT(i <= j, "bad ordering relation");
      } while (comp(*key0, array[j]));
      do
      {
        __TBB_ASSERT(i<=j, NULL);
        if (i == j)
          goto partition;
        ++i;
      } while ( comp(array[i], *key0));
      if (i == j)
        goto partition;
        std::swap(array[i], array[j]);
    }

    partition:
    // Put the partition key where it belongs
    std::swap(array[j], *key0);
    // array[1..j) is less than or equal to key
    // array[j..r) is greater than or equal to key
    i = j + 1;
    begin = array + i;
    size = range.size - i;
    range.size = j;
  }
};

template<typename RandomAccessIterator, typename Compare>
struct
quick_sort_body
{
  void operator() (const quick_sort_range<RandomAccessIterator, Compare> & range) const
  {
    // SerialQuicksort(range.begin, range.size, range.comp)
    std::sort(range.begin, range.begin + range.size, range.comp);
  };
};

template <typename RandomAccessIterator, typename Compare>
void
parallel_quick_sort (RandomAccessIterator begin,
                     RandomAccessIterator end,
                     const Compare & comp)
{
  tbb::parallel_for(quick_sort_range<RandomAccessIterator, Compare>(begin, end-begin, comp),
               quick_sort_body<RandomAccessIterator, Compare>() );
}

template <typename RandomAccessIterator, typename Compare>
void parallel_sort (RandomAccessIterator begin,
                    RandomAccessIterator end,
                    const Compare & comp)
{
  const int min_parallel_size = 1000;
  if (end > begin)
  {
    if (end - begin < min_parallel_size)
      std::sort(begin, end, comp);
    else
      parallel_quick_sort(begin, end, comp);
  }
}


template <typename RandomAccessIterator, typename Compare>
void parallel_sort (RandomAccessIterator begin,
                    RandomAccessIterator end)
{
  parallel_sort(begin, end, std::less< typename std::iterator_traits<RandomAccessIterator>::value_type >());
}

template <typename T>
inline void
parallel_sort (T * begin, T * end)
{
  parallel_sort(begin, end, std::less<T>());
}

#endif
