#include "randomext.h"

static VALUE random_zipf(int argc, VALUE *argv, VALUE self)
{
  VALUE vN, vs, vq;
  int N;
  double s, q;
  double sum;
  int i;
  double u;
  
  rb_scan_args(argc, argv, "12", &vN, &vq, &vs);
  N = NUM2INT(vN);
  s = NIL_P(vs) ? 1.0 : NUM2DBL(vs);
  q = NIL_P(vq) ? 0.0 : NUM2DBL(vq);

  if (N <= 0 || s <= 0 || q < 0)
    rb_raise(rb_eArgError, "Random#zipf_mandelbrot: parameters must be N >0, s > 0, q >= 0");
  
  for (i=1, sum=0; i<=N; ++i)
    sum += 1.0/pow(i+q, s);

  u = rb_random_real(self);
  
  for (i=1; i<=N; ++i) {
    double p = 1.0/pow(i+q, s)/sum;
    if (u <= p)
      return INT2NUM(i);
    u -= p;
  }
  
  return INT2NUM(N);
}

static VALUE random_zeta(VALUE self, VALUE vs)
{
  double s = NUM2DBL(vs);
  double q = s - 1.0;
  double r = -1.0/q;
  double t = pow(2.0, q);

  for (;;) {
    double u = 1.0 - rb_random_real(self);
    double v = rb_random_real(self);
    int x = floor(pow(u, r));
    double w = pow(1.0 + 1.0/x, q);
    if (v*x*(w-1)/(t-1) <= w/t)
      return INT2NUM(x);
  }
}

void randomext_other_init(VALUE cRandom)
{
  rb_define_method(cRandom, "zipf_mandelbrot", random_zipf, -1);
  rb_define_method(cRandom, "zeta", random_zeta, 1);
}
