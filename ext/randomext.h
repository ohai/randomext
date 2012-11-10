#include <ruby.h>
#include <math.h>
#include <stdint.h>

double randomext_random_standard_normal(VALUE random);
double randomext_sumlog(int from, int to);
double randomext_logcombination(int n, int m);

inline static uint64_t pow2(int r)
{
  return (uint64_t)1<<r;
}

inline static double random_open_interval(VALUE random)
{
  for (;;) {
    double u = rb_random_real(random);
    if (u != 0.0) return u;
  }
}

#define MASK(bits) (~(~0<<(bits)))
#define BIT(nth) (1<<(nth))

#define MAX2(n, m) (((n) < (m)) ? (m) : (n))
#define MIN2(n, m) (((n) < (m)) ? (n) : (m))

