#include "randomext.h"

static double poisson_distribution(double lambda, int m)
{
  int i;
  double sumlog = 0;
  
  for (i=1; i<=m; i++)
    sumlog += log(i);

  return exp(-lambda + m*log(lambda) - sumlog);
}

static VALUE random_poisson_inv(VALUE self, VALUE l)
{
  double lambda = NUM2DBL(l);
  int mode = floor(lambda);
  double d = poisson_distribution(lambda, mode);
  double pu, pl;
  int Xu, Xl;
  double U = rb_random_real(self);
  
  pu = pl = d; Xu = Xl = mode;
  for (;;) {
    double V = U - pu;
    if (V <= 0)
      return INT2NUM(Xu);

    U = V;
    if (Xl > 0) {
      pl *= Xl/lambda;
      --Xl;
      V = U - pl;
      if (V <= 0)
        return INT2NUM(Xl);
      U = V;
    }
    ++Xu; pu *= lambda/Xu;
  }
}

void randomext_poisson_init(VALUE cRandom)
{
  rb_define_method(cRandom, "poisson", random_poisson_inv, 1);
}
