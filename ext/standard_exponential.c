#include "randomext.h"

#define K 8
#define n (1<<K)
#define v 0.00394965982258
#define r 7.697117470131
#define m 64

static double* w = NULL;
static double* f;
static uint64_t* k;

void init_exponentail_table(void)
{
  double xi;
  int i;
  
  w = ALLOC_N(double, n);
  f = ALLOC_N(double, n);
  k = ALLOC_N(uint64_t, n);

  w[n-1] = v*exp(r)/pow2(m-K);
  w[n-2] = r/pow2(m-K);
  k[n-1] = floor(r/w[n-1]);
  f[n-1] = exp(-r);
  xi = r;
  
  for (i=n-2; i >= 1; --i) {
    xi = -log(exp(-xi) + v/xi);
    w[i-1] = xi/pow2(m-K);
    k[i] = floor(xi/w[i]);
    f[i] = exp(-xi);
  }
  k[0] = 0;
  f[0] = 1;
}

static double standard_exponential(VALUE rng)
{
  uint32_t u1; 
  uint32_t i; 
  uint64_t u2; 

  if (w == NULL)
    init_exponentail_table();
  
retry:
  u1 = rb_random_int32(rng);
  i = MASK(K) & u1;
  u2 = (uint64_t)rb_random_int32(rng) | (((uint64_t)u1 >> K) << 32);
  
  if (u2 < k[i])
    return (double)u2 * w[i];
  if (i == n-1) {
    double u = rb_random_real(rng);
    return r - log(1-u);
  } else {
    double ux = u2*w[i];
    double u = rb_random_real(rng);
    if (u*(f[i]-f[i+1]) <= exp(-ux) - f[i+1])
      return ux;
    goto retry;
  }
}

static VALUE random_standard_exp(VALUE self)
{
  return DBL2NUM(standard_exponential(self));
}

void randomext_standard_exponential_init(VALUE cRandom)
{
  rb_define_method(cRandom, "standard_exponential", random_standard_exp, 0);
}
