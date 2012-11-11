#include "randomext.h"

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

static double binomial_distribution(int k, int n, double theta)
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

/* Returns Bin(x+1| n, theta)/Bin(x| n, theta) */
inline static double forward_ratio(double x, double n, double theta)
{
  return ((n+1)/(x+1) - 1)*(theta/(1-theta));
}

/* Returns Bin(x-1| n, theta)/Bin(x| n, theta) */
inline static double backward_ratio(double x, double n, double theta)
{
  return ((n+1)/(n+1-x)-1)*((1-theta)/theta);
}

static void check_binomial_params(int n, double theta, const char* method_name)
{
  if (n < 1)
    rb_raise(rb_eArgError, "%s: n must be >= 1", method_name);
  if (theta <= 0.0 || 1.0 <= theta)
    rb_raise(rb_eArgError, "%s: n must be in (0, 1)", method_name);
  
}
/*
 * Draws a random sample from a binomial distribution
 *
 * Inverse function method is used.
 *
 * @overload binomial(n, theta)
 * @param [Integer] n the number of trials (n > 0)
 * @param [Float] theta success probability (0 < theta < 1)
 * @return [Integer] a random sample in 0..n
 */
static VALUE random_binomial_inv(VALUE self, VALUE num, VALUE prob)
{
  int n = NUM2INT(num);
  double theta = NUM2DBL(prob);
  int mode = floor(theta*(n+1));
  int xl = mode;
  int xu = mode+1;
  double pl = binomial_distribution(xl, n, theta);
  double pu = pl*forward_ratio(xl, n, theta);
  double u = rb_random_real(self);

  check_binomial_params(n, theta, "Random#binomial");
  
  for (;xl >=0 || xu <= n;) {
    if (xl >= 0) {
      if (u <= pl)
        return INT2NUM(xl);
      u = u - pl;
      pl *= backward_ratio(xl, n, theta);
      --xl;
    }
    if (xu <= n) {
      if (u <= pu)
        return INT2NUM(xu);
      u = u - pu;
      pu *= forward_ratio(xu, n, theta);
      ++xu;
    }
  }
  
  return INT2FIX(0);
}

static void fill_binomial_table(binomial_t *bin)
{
  double theta = bin->theta;
  double n = bin->n;
  int mode = floor(theta*(n+1));
  int k;

  bin->p[mode] = binomial_distribution(mode, n, theta);
  for (k = mode+1; k <= n; ++k) {
    bin->p[k] = bin->p[k-1]*forward_ratio(k-1, n, theta);
  }
  for (k = mode-1; k >= 0; --k) {
    bin->p[k] = bin->p[k+1]*backward_ratio(k+1, n, theta);
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

/*
 * Returns a random sampler from a binomial distribution.
 *
 * This sampler uses table plus square histgram method with
 * Robin Hoot method. This method constructs a table for each
 * distribution. If you once construct the table, you can
 * draw a random sample fast for large n (n >= 40), but the
 * cost of the table construction is expensive. Therefore,
 * if you need to draw many samples from the same binomial distribution,
 * you had better to use this class. Otherwise, you should use
 * Random#binomial.
 *
 * @overload initialize(rng, n, theta)
 * @param [Random] rng a random number generator
 * @param [Integer] n the number of trials (n > 0)
 * @param [Float] theta success probability (0 < theta < 1)
 * @return [Random::Binomial] a random number generator from the specified binomial distribution
 */
static VALUE binomial_initialize(VALUE self, VALUE rng, VALUE num, VALUE prob)
{
  binomial_t *bin;
  int n = NUM2INT(num);
  double theta = NUM2DBL(prob);

  check_binomial_params(n, theta, "Random::Binomial.new");
  
  rb_iv_set(self, "rng", rng);
  Data_Get_Struct(self, binomial_t, bin);
  bin->n = n;
  bin->theta = theta;
  bin->p = ALLOC_N(double, bin->n + 1);

  fill_binomial_table(bin);
  fill_binomial_T_VT_table(bin);
  
  return Qnil;
}

/*
 * Draws a sample from the binomimial distribution whose parameters
 * are specified in Random::Binomial.new.
 *
 * @return [Integer] a random sample
 */
static VALUE binomial_rand(VALUE self)
{
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

/*
 * @!attribute [r] n
 * @return [Integer] the parameter n
 */
static VALUE binomial_n(VALUE self)
{
  binomial_t *bin;
  Data_Get_Struct(self, binomial_t, bin);

  return INT2NUM(bin->n);
}

/*
 * @!attribute [r] theta
 * @return [Float] the parameter theta
 */
static VALUE binomial_theta(VALUE self)
{
  binomial_t *bin;
  Data_Get_Struct(self, binomial_t, bin);

  return DBL2NUM(bin->theta);
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
  
  rb_define_method(cRandom, "binomial", random_binomial_inv, 2);
  rb_define_alloc_func(cBinomial, binomial_alloc);
  rb_define_method(cBinomial, "initialize", binomial_initialize, 3);
  rb_define_method(cBinomial, "rand", binomial_rand, 0);
  rb_define_method(cBinomial, "n", binomial_n, 0);
  rb_define_method(cBinomial, "theta", binomial_theta, 0);
  //rb_define_method(cBinomial, "debug_info", binomial_debug_info, 0);
}
