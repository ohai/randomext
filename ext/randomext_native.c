#include "randomext.h"

#define SUMLOG_TABLE_SIZE_INIT 1024
#define SUMLOG_TABLE_SIZE_MAX (1024*32)

/* sumlog_table[0] = 0 */
/* sumlog_table[i] = log(1) + log(2) + ... + log(i) */
static double* sumlog_table = NULL;
static int sumlog_table_size = -1;

static inline void setup_sumlog_table(int need)
{
  int old_sumlog_table_size;
  int i;
  
  if (sumlog_table_size > need || need >= SUMLOG_TABLE_SIZE_MAX)
    return;
  
  if (sumlog_table == NULL) {
    sumlog_table_size = SUMLOG_TABLE_SIZE_INIT;
    sumlog_table = ALLOC_N(double, sumlog_table_size);
    sumlog_table[0] = 0;
    old_sumlog_table_size = 1;
  } else {
    old_sumlog_table_size = sumlog_table_size;
  }
  
  for (;sumlog_table_size < need; sumlog_table_size <<= 1)
    ;

  REALLOC_N(sumlog_table, double, sumlog_table_size);

  for (i = old_sumlog_table_size; i < sumlog_table_size; ++i) {
    sumlog_table[i] = sumlog_table[i-1] + log(i);
  }
}

double randomext_sumlog(int from, int to)
{
  int i;
  double ret = 0.0;

  setup_sumlog_table(to);
  if (to < SUMLOG_TABLE_SIZE_MAX)
    return sumlog_table[to] - sumlog_table[from-1];
  
  for (i=from; i<=to; ++i)
    ret += log(i);
  return ret;
}

double randomext_logcombination(int n, int m)
{
  return randomext_sumlog(n-m+1, n) - randomext_sumlog(1, m);
}

extern void randomext_standard_normal_init(VALUE cRandom);
extern void randomext_standard_exponential_init(VALUE cRandom);
extern void randomext_gamma_init(VALUE cRandom);
extern void randomext_binomial_init(VALUE cRandom);
extern void randomext_poisson_init(VALUE cRandom);
extern void randomext_hypergeometric_init(VALUE cRandom);
extern void randomext_other_init(VALUE cRandom);

void Init_randomext_native()
{
  VALUE cRandom = rb_const_get(rb_cObject, rb_intern("Random"));
  
  randomext_standard_normal_init(cRandom);
  randomext_standard_exponential_init(cRandom);
  randomext_gamma_init(cRandom);
  randomext_binomial_init(cRandom);
  randomext_poisson_init(cRandom);
  randomext_hypergeometric_init(cRandom);
  randomext_other_init(cRandom);
}
