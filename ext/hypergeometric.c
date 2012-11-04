#include "randomext.h"

#define SUMLOG_TABLE_MAX 1000

/* sumlog_table[0] = 0 */
/* sumlog_table[i] = log(1) + log(2) + ... + log(i) */
static double sumlog_table[SUMLOG_TABLE_MAX] = { -1.0 };

static void fill_sumlog_table(void)
{
  int i;
  
  sumlog_table[0] = 0.0;
  for (i=1; i<SUMLOG_TABLE_MAX; ++i) {
    sumlog_table[i] = sumlog_table[i-1] + log(i);
  }
}

static double sumlog(int from, int to)
{
  int i;
  double ret = 0.0;

  if (sumlog_table[0] < 0.0)
    fill_sumlog_table();
  
  if (1 <= from && to < SUMLOG_TABLE_MAX)
    return sumlog_table[to] - sumlog_table[from-1];

  for (i=from; i<=to; ++i)
    ret += log(i);
  return ret;
}

static double logcombination(int n, int m)
{
  return sumlog(n-m+1, n) - sumlog(1, m);
}

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
  return exp(logcombination(M, x)
             + logcombination(N-M, n-x)
             - logcombination(N, n));
}

static VALUE random_hypergoemtric_inv(VALUE self, VALUE vN, VALUE vM, VALUE vn)
{
  int N = NUM2INT(vN);
  int M = NUM2INT(vM);
  int n = NUM2INT(vn);
  int ok = (N >= 0) && (M >= 0) && (n >= 0) && (M <= N) && (n <= N);
  int mode = (M+1)*(n+1) / (N+2);
  int x_min = MAX2(0, n - (N-M));
  int x_max = MIN2(n, M);
  int xu = mode;
  double pu = hypergeometric_distribution(mode, N, M, n);
  int xl = mode-1;
  double pl = pu * backward_ratio(mode, N, M, n);
  double u = rb_random_real(self);

  if (!ok)
    rb_raise(rb_eArgError,
             "Paramters of hypergeometric distribution must be:"
             "(N >= 0) && (M >= 0) && (n >= 0) && (M <= N) && (n <= N)");
  
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

  fill_sumlog_table();
}
