#include <ruby.h>
#include <math.h>
#include <stdint.h>

inline static uint64_t pow2(int r)
{
  return (uint64_t)1<<r;
}

#define MASK(bits) (~(~0<<(bits)))
#define BIT(nth) (1<<(nth))

