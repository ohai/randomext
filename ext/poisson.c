#include "randomext.h"

static double poisson_distribution(int m, double lambda)
{
  int i;
  double sumlog = 0;
  
  for (i=1; i<=m; i++)
    sumlog += log(i);

  return exp(-lambda + m*log(lambda) - sumlog);
}

/* Returns Poi(x+1|lambda)/Poi(x|lambda) */
static inline double forward_ratio(int x, double lambda)
{
  return lambda/(x+1);
}

/* Returns Poi(x-1|lambda)/Poi(x|lambda) */
static inline double backward_ratio(int x, double lambda)
{
  return x/lambda;
}

/*
 * call-seq: prng.poisson(lambda) -> int
 *
 * Draws a random sample from a Poisson distribution.
 *
 * Inverse function method is used.
 */
static VALUE random_poisson_inv(VALUE self, VALUE l)
{
  double lambda = NUM2DBL(l);
  int mode = floor(lambda);
  int xu = mode;
  int xl = mode - 1;
  double pu = poisson_distribution(mode, lambda);
  double pl = pu * backward_ratio(xu, lambda);
  double u = rb_random_real(self);
  
  for (;;) {
    if (u <= pu)
      return INT2NUM(xu);
    u = u - pu;
    pu *= forward_ratio(xu, lambda);
    ++xu;

    if (xl >= 0) {
      if (u <= pl)
        return INT2NUM(xl);
      u = u - pl;
      pl *= backward_ratio(xl, lambda);
      --xl;
    }
  }
}

void randomext_poisson_init(VALUE cRandom)
{
  rb_define_method(cRandom, "poisson", random_poisson_inv, 1);
}
