#include "randomext.h"

#define R 3.442619855899
#define V 9.91256303526217e-3
#define K 7
#define M 64
#define N (1<<K)

#include "standard_normal_table.h"

inline static double sn(double x)
{
  return exp(-x*x/2);
}

static double sample_from_tail(VALUE random)
{
  for (;;) {
    double x = sqrt(R*R-2*log(1-rb_random_real(random)));
    if (x*rb_random_real(random) <= R)
      return x;
  }
}

double randomext_random_standard_normal(VALUE random)
{
  int i;
  uint64_t u;
  int sign;
  double ux;

  for (;;) {
    unsigned int u0 = rb_random_int32(random);
    i = u0 & MASK(K);
    sign = (u0 & BIT(K)) ? 1 : -1;
    u = ((uint64_t)(u0 >> (K+1)) << 32) | rb_random_int32(random);

    if (u < k[i])
      return sign*(u*w[i]);
    if (i == N-1)
      return sign*sample_from_tail(random);
    ux = u * w[i];
    if ( rb_random_real(random)*(f[i]-f[i+1]) <= sn(ux)-f[i+1])
      return sign*ux;
  }

}

/*
 * Draws a random sample from the standard normal distribution.
 *
 * Ziggurat method is used for random sampling.
 * 
 * @return [Float] a random sample 
 */
static VALUE random_standard_normal(VALUE self)
{
  return DBL2NUM(randomext_random_standard_normal(self));
}

void randomext_standard_normal_init(VALUE cRandom)
{
  rb_define_method(cRandom, "standard_normal", random_standard_normal, 0);
}
