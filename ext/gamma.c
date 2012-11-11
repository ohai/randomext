#include "randomext.h"

/*
 * @private
 */
static VALUE random_gamma(VALUE self, VALUE shape)
{
  double c, d;
  
  if (NUM2DBL(shape) < 1.0)
    rb_raise(rb_eArgError, "Random#_gamma: shape parameter must be >= 1.0");
  
  d = NUM2DBL(shape) - 1.0/3.0;
  c = 1/sqrt(9*d);

  for (;;) {
    double z, v, y, w, u;
    z = randomext_random_standard_normal(self);
    v = 1 + c*z;
    if (v <= 0) continue;
    w = v*v*v; y = d*w;
    u = random_open_interval(self);
    if (u > 1 - 0.0331*(z*z*z*z) && z*z/2 + d*log(w) - y + d < log(u))
      continue;
    return DBL2NUM(y);
  }
}

void randomext_gamma_init(VALUE cRandom)
{
  rb_define_private_method(cRandom, "_gamma", random_gamma, 1);
}
