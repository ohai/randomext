#include <ruby.h>
#include <math.h>
#include <stdint.h>

#define R 3.442619855899
#define V 9.91256303526217e-3
#define K 7
#define M 64
#define N (1<<K)

#define MASK(bits) (~(~0<<(bits)))
#define BIT(nth) (1<<(nth))

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

static double rb_random_standard_normal(VALUE random)
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
  return DBL2NUM(rb_random_standard_normal(self));
}

inline static uint64_t pow2(int r)
{
  return (uint64_t)1<<r;
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

inline static double random_open_interval(VALUE random)
{
  for (;;) {
    double u = rb_random_real(random);
    if (u != 0.0) return u;
  }
}

static VALUE random_gamma(VALUE self, VALUE shape)
{
  double d = NUM2DBL(shape) - 1.0/3.0;
  double c = 1/sqrt(9*d);

  for (;;) {
    double z, v, y, w, u;
    z = rb_random_standard_normal(self);
    v = 1 + c*z;
    if (v <= 0) continue;
    w = v*v*v; y = d*w;
    u = random_open_interval(self);
    if (u > 1 - 0.0331*(z*z*z*z) && z*z/2 + d*log(w) - y + d < log(u))
      continue;
    return DBL2NUM(y);
  }
}

static VALUE random_binomial_inv(VALUE self, VALUE num, VALUE prob)
{
  int n = NUM2INT(num);
  double theta = NUM2DBL(prob);
  double s = theta/(1-theta);
  double a = (n+1)*s;
  double t = 1.0/s;
  double b = (n+1)*t;
  int m = ceil(theta*(n-1));
  double d;
  int i;
  double pu, pl; 
  int xu, xl; 
  double u = rb_random_real(self);
  d = pow((1-theta), n);
  for (i=1; i<=m; ++i)
    d *= (n-i+1)*s/i;
  pu = pl = d;
  xu = xl = m;
  
  for (;;) {
    double v = u - pu;
    //printf("u:%f pu:%e\n", u, pu);
    if (v <= 0)
      return INT2NUM(xu);
    u = v;
    for (;;) {
      if (xl > 0) {
        xl--;
        pl = (b/(n-xl) - t)*pl;
        v = u - pl;
        if (v <= 0)
          return INT2NUM(xl);
        u = v;
      } else if (xl == 0 && xu == n) {
        return INT2NUM(n);
      }
      if (xu < n) {
        ++xu;
        pu *= a/xu - s;
        break;
      } else if (xu == n && xl == 0) {
        return INT2FIX(0);
      }
    }
  }
}

static double binomial_distribution(int n, int k, double theta)
{
  double ret = 1.0;
  int i;
  for (i=1; i <= k; ++i) {
    ret *= n-k+i;
    ret /= i;
    ret *= theta;
    if (i <= n - k)
      ret *= 1 - theta;
  }
  for (i=0; i<n-2*k; ++i)
    ret *= 1 - theta;

  return ret;
}

static void fill_binomial_table(double p[], int n, double theta)
{
  int mode = ceil(theta*(n+1));
  double s = theta/(1-theta);
  int k;

  p[mode] = binomial_distribution(n, mode, theta);
  for (k = mode+1; k <= n; ++k) {
    p[k] = s*((n+1)/(double)k - 1)*p[k-1];
  }
  for (k = mode-1; k >= 0; --k) {
    p[k] = p[k+1]/(s*((n+1)/(double)(k+1)-1));
  }
}

static VALUE random_binomial_table(VALUE self, VALUE num, VALUE prob)
{
  int n = NUM2UINT(num);
  double theta = NUM2DBL(prob);
  double *p = ALLOCA_N(double, n+1);
  int *t = ALLOCA_N(int, n+1);
  double nt;
  int i, k, b;
  
  fill_binomial_table(p, n, theta);

  for (k=8,, nt=0; k < 16 && (double)nt/pow2(k) < 0.9; ++k) {
    int x;
    nt = 0;
    for (x=0; x<=n; ++x) {
      t[x] = ceil(b*p[x]);
      nt += t[x];
    }
  }
  
  return Qnil;
}

void Init_randomext_native()
{
  VALUE cRandom = rb_const_get(rb_cObject, rb_intern("Random"));

  rb_define_method(cRandom, "standard_normal", random_standard_normal, 0);
  rb_define_private_method(cRandom, "_gamma", random_gamma, 1);
  rb_define_method(cRandom, "binomial1", random_binomial_inv, 2);
  rb_define_method(cRandom, "binomial2", random_binomial_table, 2);
  init_table();
}
