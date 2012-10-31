#include <ruby.h>
#include <math.h>
#include <stdint.h>
#include "randomext.h"

#define BINOMIAL_M 64

typedef struct {
  int n;
  double theta;
  int N;
  int k;
  double *p;
  int *T;
  int *K;
  double *V;
} binomial_t;

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

static VALUE random_binomial_inv(VALUE self, VALUE num, VALUE prob)
{
  int n = NUM2INT(num);
  double theta = NUM2DBL(prob);
  double s = theta/(1-theta);
  double a = (n+1)*s;
  double t = 1.0/s;
  double b = (n+1)*t;
  int m = ceil(theta*(n-1));
  double d = binomial_distribution(n, m, theta);
  double pu, pl; 
  int xu, xl; 
  double u = rb_random_real(self);
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

static void fill_binomial_table(binomial_t *bin)
{
  double theta = bin->theta;
  double n = bin->n;
  int mode = floor(theta*(n+1));
  double s = theta/(1-theta);
  int k;

  bin->p[mode] = binomial_distribution(n, mode, theta);
  for (k = mode+1; k <= n; ++k) {
    bin->p[k] = s*((n+1)/(double)k - 1)*bin->p[k-1];
  }
  for (k = mode-1; k >= 0; --k) {
    bin->p[k] = bin->p[k+1]/(s*((n+1)/(double)(k+1)-1));
  }
}

static void fill_binomial_VK_table(binomial_t *bin, double ntheta[])
{
  double c = 1.0/(bin->n + 1);
  int i, j, n, x;
  
  bin->K = ALLOC_N(int, bin->n + 1);
  bin->V = ALLOC_N(double, bin->n + 1);

  for (i=0; i<=bin->n; ++i) {
    bin->K[i] = i;
    bin->V[i] = (i+1)*c;
  }

  for (n=0; n<bin->n; ++n) {
    for (x=0, i=j=0; x<=bin->n; ++x) {
      if (ntheta[x] < ntheta[i])
        i = x;
      if (ntheta[x] > ntheta[j])
        j = x;
    }
    bin->K[i] = j;
    bin->V[i] = i*c + ntheta[i];
    ntheta[j] = ntheta[j] - (c - ntheta[i]);
    ntheta[i] = c;
  }
}

static void fill_binomial_T_VT_table(binomial_t *bin)
{
  int k, i, x;
  int nt;
  double w;
  int *qt = ALLOC_N(int, bin->n + 2);
  double *theta = ALLOC_N(double, bin->n + 1);
  double *ntheta = ALLOC_N(double, bin->n + 1);
  double sum_theta = 0;
  
  for (k=7; ;k++) {
    int b = pow2(k);
    nt = 0;
    qt[0] = 0;
    for (x=0; x<=bin->n; x++) {
      qt[x+1] = floor(b*bin->p[x]) + qt[x];
    }
    if (k > 16 || qt[bin->n + 1] > 0.9*b)
      break;
  }

  bin->k = k;
  bin->N = pow2(k);
  bin->T = ALLOC_N(int, bin->N);
  w = pow(2, k - BINOMIAL_M);
  for (x=0; x<=bin->n; ++x) {
    for (i=qt[x]; i<qt[x+1]; ++i)
      bin->T[i] = x;
  }
  for (i=qt[bin->n + 1]; i<bin->N; ++i)
    bin->T[i] = -1;

  for (x=0; x<=bin->n; ++x) {
    theta[x] = pow2(k) * bin->p[x] - (qt[x+1]-qt[x]);
    sum_theta += theta[x];
  }
  
  for (x=0; x<=bin->n; ++x)
    ntheta[x] = theta[x]/sum_theta;

  fill_binomial_VK_table(bin, ntheta);

  xfree(qt);
  xfree(theta);
  xfree(ntheta);
}

static void binomial_free(binomial_t *bin)
{
  xfree(bin->p);
  xfree(bin->T);
  xfree(bin->K);
  xfree(bin->V);
  xfree(bin);
}
      
static VALUE binomial_alloc(VALUE klass)
{
  binomial_t *bin;
  VALUE obj = Data_Make_Struct(klass, binomial_t, 0, binomial_free, bin);
  bin->n = -1; 
  bin->p = NULL; bin->T = NULL; bin->K = NULL; bin->V = NULL;
  return obj;
}

static VALUE binomial_initialize(VALUE self, VALUE rng, VALUE num, VALUE prob)
{
  binomial_t *bin;
  
  rb_iv_set(self, "rng", rng);
  Data_Get_Struct(self, binomial_t, bin);
  bin->n = NUM2UINT(num);
  bin->theta = NUM2DBL(prob);
  bin->p = ALLOC_N(double, bin->n + 1);

  fill_binomial_table(bin);
  fill_binomial_T_VT_table(bin);
  
  return Qnil;
}

static VALUE binomial_rand(VALUE self)
{
  /* Assume BINOMIAL_M == 64 */
  VALUE rng = rb_iv_get(self, "rng");
  binomial_t *bin;
  uint32_t I0 = rb_random_int32(rng);
  uint32_t ILk;
  double U;
  int J;

  Data_Get_Struct(self, binomial_t, bin);
  
  ILk = I0 & MASK(bin->k);
  if (bin->T[ILk] >= 0)
    return INT2NUM(bin->T[ILk]);
  
  U = rb_random_real(rng);
  J = floor((bin->n + 1)*U);
  if (U < bin->V[J])
    return INT2NUM(J);
  else
    return INT2NUM(bin->K[J]);
}

#if 0
static VALUE binomial_debug_info(VALUE self)
{
  binomial_t *bin;
  int i;
  
  Data_Get_Struct(self, binomial_t, bin);
  

  printf("N=%d\n", bin->N);
  for (i=0; i<bin->N; ++i) {
    printf("%d ", bin->T[i]);
  }
  puts("");
  for (i=0; i<=bin->n; ++i) {
    printf("%f %d\n", bin->V[i], bin->K[i]);
  }
  return Qnil;
}
#endif

void randomext_binomial_init(VALUE cRandom)
{
  VALUE cBinomial = rb_define_class_under(cRandom, "Binomial", rb_cObject);
  
  rb_define_method(cRandom, "binomial1", random_binomial_inv, 2);
  rb_define_alloc_func(cBinomial, binomial_alloc);
  rb_define_method(cBinomial, "initialize", binomial_initialize, 3);
  rb_define_method(cBinomial, "rand", binomial_rand, 0);
  //rb_define_method(cBinomial, "debug_info", binomial_debug_info, 0);
}
