#include "randomext.h"

#define R 3.442619855899
#define V 9.91256303526217e-3
#define K 7
#define M 64
#define N (1<<K)

static double w[N];
static int64_t k[N];
static double f[N];
static double x[N+1];

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
  int64_t u;
  int sign;
  double ux;
  
  for (;;) {
    unsigned int u0 = rb_random_int32(random);
    i = u0 & MASK(K);
    sign = (u0 & BIT(K)) ? 1 : -1;
    u = ((uint64_t)(u0 >> (K+1)) << 32) | rb_random_int32(random);

    if (u < k[i])
      return sign*u*w[i];
    if (i == N-1)
      return sign*sample_from_tail(random);
    ux = u * w[i];
    if ( rb_random_real(random)*(f[i]-f[i+1]) <= sn(ux)-f[i+1])
      return sign*ux;
  }

}

/*
 * call-seq: prng.standard_normal() -> float
 *
 * Draws a random sample from the standard normal distribution.
 *
 * Ziggurat method is used for random sampling.
 */
static VALUE random_standard_normal(VALUE self)
{
  return DBL2NUM(randomext_random_standard_normal(self));
}

static void init_table(void)
{
  int i;
  
  w[N-1] = V*exp(R*R/2)/pow2(M-K-1);
  w[N-2] = R/pow2(M-K-1);
  k[N-1] = ceil(R/w[N-1]);
  f[N-1] = sn(R);
  x[N] = V*sn(R);
  x[N-1] = R;

  for (i=N-2; i>=1; --i) {
    x[i] = sqrt(-2*log(sn(x[i+1])+V/x[i+1]));
    w[i-1] = x[i]/pow2(M-K-1);
    k[i] = ceil(x[i]/w[i]);
    f[i] = sn(x[i]);
  }
  k[0] = 0;
  f[0] = 1;
}

void randomext_standard_normal_init(VALUE cRandom)
{
  rb_define_method(cRandom, "standard_normal", random_standard_normal, 0);
  init_table();
}
