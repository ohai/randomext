#include "randomext.h"

#define R 3.442619855899
#define V 9.91256303526217e-3
#define K 7
#define M 64
#define N (1<<K)

static double* w = NULL;
static uint64_t* k;
static double* f;

inline static double sn(double x)
{
  return exp(-x*x/2);
}

static void init_snormal_table(void)
{
  int i;
  double xi;
  
  w = ALLOC_N(double, N);
  k = ALLOC_N(uint64_t, N);
  f = ALLOC_N(double, N);
  
  w[N-1] = V*exp(R*R/2)/pow2(M-K-1);
  w[N-2] = R/pow2(M-K-1);
  k[N-1] = floor(R/w[N-1]);
  f[N-1] = sn(R);
  xi = R;
  
  for (i=N-2; i>=1; --i) {
    xi = sqrt(-2*log(sn(xi)+V/xi));
    w[i-1] = xi/pow2(M-K-1);
    k[i] = floor(xi/w[i]);
    f[i] = sn(xi);
  }
  k[0] = 0;
  f[0] = 1;
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

  if (w == NULL)
    init_snormal_table();
  
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
