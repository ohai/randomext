#include "randomext.h"

static double poisson_distribution(int m, double lambda)
{
  return exp(-lambda + m*log(lambda) - randomext_sumlog(1, m));
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
 * Draws a random sample from a Poisson distribution.
 *
 * Inverse function method is used.
 *
 * @overload poisson(lambda)
 * @param [Float] lambda mean
 * @return [Integer] a random sample in [0, INFINITY)
 */
static VALUE random_poisson_inv(VALUE self, VALUE l)
{
  double lambda = NUM2DBL(l);
  int mode, xu, xl;
  double pu, pl, u;
  
  if (lambda <= 0.0)
    rb_raise(rb_eArgError, "Random#poisson: lambda must be positive");

  mode = floor(lambda);
  xu = mode;
  xl = mode - 1;
  pu = poisson_distribution(mode, lambda);
  pl = pu * backward_ratio(xu, lambda);
  u = rb_random_real(self);
  
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
