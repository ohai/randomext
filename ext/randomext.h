#include <ruby.h>
#include <math.h>
#include <stdint.h>

double randomext_random_standard_normal(VALUE random);

inline static uint64_t pow2(int r)
{
  return (uint64_t)1<<r;
}

#define MASK(bits) (~(~0<<(bits)))
#define BIT(nth) (1<<(nth))

#define MAX2(n, m) (((n) < (m)) ? (m) : (n))
#define MIN2(n, m) (((n) < (m)) ? (n) : (m))

