#include "randomext.h"

/* Returns HGeo(x+1| N, M, n)/HGeo(x| N, M, n) */
inline static double forward_ratio(int x, int N, int M, int n)
{
  return (double)((M-x)*(n-x))/((x+1)*(N-M-n+x+1));
}

/* Returns HGeo(x-1| N, M, n)/HGeo(x| N, M, n) */
inline static double backward_ratio(int x, int N, int M, int n)
{
  return 1.0/forward_ratio(x-1, N, M, n);
}

/* Returns HGeo(x| N, M, n) */
static inline double hypergeometric_distribution(int x, int N, int M, int n)
{
  return exp(randomext_logcombination(M, x)
             + randomext_logcombination(N-M, n-x)
             - randomext_logcombination(N, n));
}

/*
 * Draws a random sample from a hypergeometric distribution.
 *
 * Inverse method is used.
 *
 * @overload hypergeometric(N, M, n)
 * @param [Integer] N a population (N >= 0)
 * @param [Integer] M the number of successes (0 <= M <= N)
 * @param [Integer] n the number of samples (0 <= n <= N)
 * @return [Integer] a random sample in [max(0, n-(N-M)), min(n, M)]
 */
static VALUE random_hypergoemtric_inv(VALUE self, VALUE vN, VALUE vM, VALUE vn)
{
  int N = NUM2INT(vN);
  int M = NUM2INT(vM);
  int n = NUM2INT(vn);
  int ok = (N >= 0) && (M >= 0) && (n >= 0) && (M <= N) && (n <= N);
  int mode, x_min, x_max, xu, xl;
  double pl, pu;
  double u;
  
  if (!ok)
    rb_raise(rb_eArgError,
             "Random#hypergeometric: paramters must be:"
             "(N >= 0) && (M >= 0) && (n >= 0) && (M <= N) && (n <= N)");

  mode = (M+1)*(n+1) / (N+2);
  x_min = MAX2(0, n - (N-M));
  x_max = MIN2(n, M);
  xu = mode;
  pu = hypergeometric_distribution(mode, N, M, n);
  xl = mode-1;
  pl = pu * backward_ratio(mode, N, M, n);
  u = rb_random_real(self);

  for (;x_min <= xl || xu <= x_max;) {
    if (xu <= x_max) {
      if (u <= pu)
        return INT2NUM(xu);
      u -= pu;
      pu *= forward_ratio(xu, N, M, n);
      ++xu;
    }
    if (xl >= x_min) {
      if (u <= pl)
        return INT2NUM(xl);
      u -= pl;
      pl *= backward_ratio(xl, N, M, n);
      --xl;
    }
  }
  
  return INT2NUM(x_min);
}

void randomext_hypergeometric_init(VALUE cRandom)
{
  rb_define_method(cRandom, "hypergeometric", random_hypergoemtric_inv, 3);
}
